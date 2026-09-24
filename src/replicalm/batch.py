r"""Process many tiles, resumably, across worker processes.

WHY THIS EXISTS
---------------
A tool that handles one file at a time is not an operational tool. A survey is
hundreds of tiles, a run over them takes days, and anything that takes days will
be interrupted -- a reboot, a full disk, one malformed header. So the unit of
work is the tile, the record of what happened is written after every tile, and a
second invocation picks up where the first stopped.

WHAT ISOLATION MEANS HERE
-------------------------
One tile failing must not end the run or contaminate its neighbours. Each tile
is processed in a worker process that does nothing else, so a crash inside PDAL
or GDAL costs that tile and no other, and the traceback is recorded rather than
printed and lost. `maxtasksperchild=1` retires each worker afterwards, which
also bounds the native-library memory growth that shows up over hundreds of
tiles.

WHY PROCESSES AND NOT THREADS
-----------------------------
The work is numpy, PDAL and GDAL, and the part that dominates -- the kriging
solve -- holds the GIL. Threads would serialise on exactly the expensive step.

ORDER
-----
Largest tiles first. With N workers and unequal tiles, starting the big ones
early keeps every worker busy to the end, instead of leaving one straggler
running alone after the others have drained.
"""
import argparse
import glob
import json
import multiprocessing as mp
import os
import time
import traceback


def _tile_key(path):
    return os.path.splitext(os.path.basename(str(path)))[0]


def _one(job):
    """Process a single tile in a worker. Returns a record, never raises."""
    path, out_dir, preset, cell_m, make_g1, keep_ground, rvt_script = job
    # Set before the first numpy import in this process, which is why it is here
    # and not at module level: LAPACK and the neighbour query each thread across
    # every core by default, so N workers would each claim all of them and spend
    # their time contending. Parallelism belongs at the tile level here.
    for var in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                "NUMEXPR_NUM_THREADS", "REPLICALM_QUERY_WORKERS"):
        os.environ.setdefault(var, "1")
    t0 = time.time()
    rec = {"tile": _tile_key(path), "las": str(path), "status": "failed"}
    try:
        from .config import PRESETS
        from . import pipeline
        cfg = PRESETS[preset]
        out = pipeline.process(path, out_dir, cfg=cfg, cell_m=cell_m,
                               make_g1=make_g1, keep_ground=keep_ground,
                               rvt_script=rvt_script, progress=None)
        rec.update(status="done", seconds=round(time.time() - t0, 1),
                   dem=out["dem"], ground_points=out["ground_points"],
                   density=round(out["density"], 3),
                   cell_m=round(out["cell_m"], 4),
                   search_radius_m=round(out["search_radius_m"], 3),
                   cells_after_trim=out["cells_after_trim"],
                   filled_cells=out.get("filled_cells", 0),
                   # How much of this tile was actually kriged. Without it a
                   # survey-wide run cannot say what fraction of its cells came
                   # from the covariance model and what fraction from inverse
                   # distance, which is a methodological fact rather than a
                   # diagnostic.
                   fallback_fraction=round(out.get("fallback_fraction", 0.0), 5),
                   variogram_range_m=(round(out["variogram_range_m"], 4)
                                      if out.get("variogram_range_m") else None),
                   max_points=cfg.max_points, min_points=cfg.min_points,
                   locked_baseline=out.get("locked_baseline"))
    except BaseException as exc:            # a worker must not die silently
        rec.update(seconds=round(time.time() - t0, 1),
                   error="%s: %s" % (type(exc).__name__, exc),
                   traceback=traceback.format_exc(limit=6))
    return rec


def _load(manifest_path):
    if os.path.exists(manifest_path):
        with open(manifest_path, encoding="utf-8") as fh:
            return json.load(fh)
    return {}


def _save(manifest_path, records):
    # Written through a temporary file: an interrupted write must never leave a
    # truncated manifest, because the manifest is what makes the run resumable.
    tmp = manifest_path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(records, fh, indent=1)
    os.replace(tmp, manifest_path)


def run(las_paths, out_dir, preset="clear", workers=None, cell_m=None,
        make_g1=False, keep_ground=False, rvt_script=None, redo=False,
        manifest=None, verbose=True):
    """Process every tile in `las_paths`. Returns the manifest dict."""
    os.makedirs(out_dir, exist_ok=True)
    manifest = manifest or os.path.join(out_dir, "batch_manifest.json")
    records = {} if redo else _load(manifest)

    todo = []
    for p in las_paths:
        done = records.get(_tile_key(p), {})
        # A tile counts as done only if it said so and its DEM is still there.
        # Trusting the manifest alone would skip tiles whose output was deleted.
        if (not redo and done.get("status") == "done"
                and done.get("dem") and os.path.exists(done["dem"])):
            continue
        todo.append(p)
    todo.sort(key=lambda p: os.path.getsize(p), reverse=True)

    workers = workers or max(1, min(8, (os.cpu_count() or 2) // 2))
    if verbose:
        print("%d tiles, %d already done, %d to process, %d workers"
              % (len(las_paths), len(las_paths) - len(todo), len(todo), workers),
              flush=True)
    if not todo:
        return records

    jobs = [(p, out_dir, preset, cell_m, make_g1, keep_ground, rvt_script)
            for p in todo]
    t0, n_ok, n_bad = time.time(), 0, 0
    with mp.Pool(processes=workers, maxtasksperchild=1) as pool:
        for i, rec in enumerate(pool.imap_unordered(_one, jobs), 1):
            records[rec["tile"]] = rec
            _save(manifest, records)
            if rec["status"] == "done":
                n_ok += 1
                msg = ("%d pts, %.2f/m2, cell %.2f m, %d cells, %.0f%% idw, "
                       "%.0f s"
                       % (rec["ground_points"], rec["density"], rec["cell_m"],
                          rec["cells_after_trim"],
                          100 * rec.get("fallback_fraction", 0.0),
                          rec["seconds"]))
            else:
                n_bad += 1
                msg = "FAILED %s" % rec["error"][:80]
            if verbose:
                elapsed = time.time() - t0
                eta = elapsed / i * (len(jobs) - i)
                print("[%4d/%d] %-34s %s  (elapsed %.1f h, eta %.1f h)"
                      % (i, len(jobs), rec["tile"], msg,
                         elapsed / 3600, eta / 3600), flush=True)
    if verbose:
        print("\n%d done, %d failed, %.2f h total"
              % (n_ok, n_bad, (time.time() - t0) / 3600))
        for k, r in sorted(records.items()):
            if r.get("status") != "done":
                print("  failed: %-34s %s" % (k, r.get("error", "")[:90]))
    return records


def main(argv=None):
    from .config import PRESETS
    p = argparse.ArgumentParser(
        prog="replicalm-batch",
        description="Process many tiles, resumably, across worker processes.")
    p.add_argument("inputs", nargs="+",
                   help="LAS/LAZ files, directories, or glob patterns")
    p.add_argument("-o", "--out", required=True, help="output directory")
    p.add_argument("--profile", "--preset", dest="preset", default="clear",
                   choices=sorted(PRESETS))
    p.add_argument("--workers", type=int, default=None,
                   help="worker processes (default: half the cores, max 8)")
    p.add_argument("--cell", type=float, default=None,
                   help="output cell size in metres (default: from density)")
    p.add_argument("--g1", action="store_true", help="also build the G1 image")
    p.add_argument("--keep-ground", action="store_true",
                   help="keep the classified point cloud for each tile")
    p.add_argument("--rvt-script", default=None)
    p.add_argument("--redo", action="store_true",
                   help="reprocess tiles already recorded as done")
    a = p.parse_args(argv)

    paths = []
    for item in a.inputs:
        if os.path.isdir(item):
            paths += sorted(glob.glob(os.path.join(item, "**", "*.la[sz]"),
                                      recursive=True))
        else:
            hits = sorted(glob.glob(item))
            paths += hits if hits else [item]
    paths = [q for q in dict.fromkeys(paths) if os.path.isfile(q)]
    if not paths:
        p.error("no LAS or LAZ files matched")

    run(paths, a.out, preset=a.preset, workers=a.workers, cell_m=a.cell,
        make_g1=a.g1, keep_ground=a.keep_ground, rvt_script=a.rvt_script,
        redo=a.redo)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
