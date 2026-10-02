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

BUT EACH PROCESS STILL WANTS SEVERAL THREADS
--------------------------------------------
Pinning every worker to one thread was right when workers outnumbered cores. It
is wrong here. Memory, not cores, decides how many tiles run at once -- three,
on a 128 GB machine -- so on 72 cores that left 69 of them idle. Measured on
this machine, going from one thread to sixteen takes the neighbour query from
3.05 s to 0.22 s and the batched solve from 0.48 s to 0.21 s; both saturate by
sixteen. So each worker is given cores divided by the number of tiles the memory
budget will actually admit, capped at sixteen because past that it buys
nothing and only invites contention.

ORDER
-----
Largest tiles first. With N workers and unequal tiles, starting the big ones
early keeps every worker busy to the end, instead of leaving one straggler
running alone after the others have drained.
"""
import argparse
import concurrent.futures as cf
import glob
import json
import multiprocessing as mp
import os
import time
import traceback

import numpy as np

# Peak resident memory per tile, as a multiple of its file size. Measured, not
# guessed: a 2.0 GB LAS held 11.6 GB while classifying, a ratio of 5.8. Eight
# was that with margin, and it was still too low -- on 2026-09-27 a run at an
# 88 GB budget admitted about six 2 GB tiles at once and the machine ran out
# before the first of them finished, while the same run at 64 GB held three
# tiles at 42-45 GB charged and survived two and a half hours. Ten is what that
# implies. Overestimating costs a slower run; underestimating has now killed
# three.
MEMORY_FACTOR = 10.0


def _tile_key(path):
    return os.path.splitext(os.path.basename(str(path)))[0]


def _one(job):
    """Process a single tile in a worker. Returns a record, never raises."""
    (path, out_dir, preset, cell_m, make_g1, keep_ground, rvt_script,
     threads) = job
    # Set before the first numpy import in this process, which is why it is here
    # and not at module level: LAPACK and the neighbour query read these once.
    # Left at one, N workers on a big machine use N cores of it -- see the note
    # at the top of this file.
    for var in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                "NUMEXPR_NUM_THREADS", "REPLICALM_QUERY_WORKERS"):
        os.environ[var] = str(threads)
    t0 = time.time()
    rec = {"tile": _tile_key(path), "las": str(path), "status": "failed"}
    # Peak resident memory for this worker, recorded rather than assumed.
    # MEMORY_FACTOR has been wrong twice and each time the run died hours in
    # with nothing to show why; a tile that records what it actually took lets
    # the next run be sized from measurement.
    try:
        import psutil
        _proc = psutil.Process()
    except Exception:
        _proc = None
    try:
        from .config import PRESETS
        from . import pipeline
        cfg = PRESETS[preset]
        out = pipeline.process(path, out_dir, cfg=cfg, cell_m=cell_m,
                               make_g1=make_g1, keep_ground=keep_ground,
                               rvt_script=rvt_script, progress=None)
        peak = None
        if _proc is not None:
            try:
                mi = _proc.memory_info()
                peak = getattr(mi, "peak_wset", None) or mi.rss
            except Exception:
                peak = None
        rec.update(status="done", seconds=round(time.time() - t0, 1),
                   peak_gb=(round(peak / (1 << 30), 2) if peak else None),
                   peak_over_las=(round(peak / max(os.path.getsize(path), 1), 2)
                                  if peak else None),
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


def _budget_bytes(explicit_gb=None):
    """How much memory the run may commit to tiles at once."""
    if explicit_gb:
        return int(explicit_gb * (1 << 30))
    try:
        import ctypes

        class MS(ctypes.Structure):
            _fields_ = [("dwLength", ctypes.c_ulong),
                        ("dwMemoryLoad", ctypes.c_ulong),
                        ("ullTotalPhys", ctypes.c_ulonglong),
                        ("ullAvailPhys", ctypes.c_ulonglong),
                        ("ullTotalPageFile", ctypes.c_ulonglong),
                        ("ullAvailPageFile", ctypes.c_ulonglong),
                        ("ullTotalVirtual", ctypes.c_ulonglong),
                        ("ullAvailVirtual", ctypes.c_ulonglong),
                        ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]

        st = MS()
        st.dwLength = ctypes.sizeof(MS)
        ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(st))
        if st.ullTotalPhys:
            # Three quarters. The rest is the operating system, this process,
            # and whatever else the machine is doing.
            return int(st.ullTotalPhys * 0.75)
    except Exception:
        pass
    return 32 * (1 << 30)



def run(las_paths, out_dir, preset="clear", workers=None, cell_m=None,
        make_g1=False, keep_ground=False, rvt_script=None, redo=False,
        manifest=None, verbose=True, memory_gb=None):
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

    budget = _budget_bytes(memory_gb)
    charge = {q: os.path.getsize(q) * MEMORY_FACTOR for q in todo}
    # Threads per worker come from how many tiles will actually be in flight,
    # which memory decides, not from the worker count. The median charge is used
    # rather than the largest: a run sized on its biggest tile would give every
    # worker far too many threads for the hundreds of ordinary ones.
    med = float(np.median(list(charge.values()))) if charge else budget
    concurrent = int(max(1, min(workers, budget // max(med, 1))))
    threads = int(max(1, min(16, (os.cpu_count() or 2) // concurrent)))
    if verbose:
        print("memory budget %.0f GB; largest tile charged %.1f GB; "
              "expect %d tiles at once, %d threads each of %d cores"
              % (budget / (1 << 30), max(charge.values()) / (1 << 30),
                 concurrent, threads, os.cpu_count() or 0), flush=True)

    t0, n_ok, n_bad, done_n = time.time(), 0, 0, 0
    pending, running, committed = list(todo), {}, 0
    total = len(pending)
    with cf.ProcessPoolExecutor(max_workers=workers,
                                max_tasks_per_child=1) as pool:
        while pending or running:
            # Start whatever fits. A tile may always run alone, however large,
            # so an oversized one cannot wedge the queue behind a budget it
            # could never satisfy.
            while pending and len(running) < workers:
                nxt = pending[0]
                if running and committed + charge[nxt] > budget:
                    break
                pending.pop(0)
                committed += charge[nxt]
                fut = pool.submit(_one, (nxt, out_dir, preset, cell_m, make_g1,
                                         keep_ground, rvt_script, threads))
                running[fut] = nxt
            finished, _ = cf.wait(running, return_when=cf.FIRST_COMPLETED)
            for fut in finished:
                path = running.pop(fut)
                committed -= charge[path]
                try:
                    rec = fut.result()
                except BaseException as exc:
                    rec = {"tile": _tile_key(path), "las": str(path),
                           "status": "failed",
                           "error": "worker died: %s" % exc}
                records[rec["tile"]] = rec
                _save(manifest, records)
                done_n += 1
                if rec["status"] == "done":
                    n_ok += 1
                    msg = ("%d pts, %.2f/m2, cell %.2f m, %d cells, %.0f%% idw,"
                           " %.0f s"
                           % (rec["ground_points"], rec["density"],
                              rec["cell_m"], rec["cells_after_trim"],
                              100 * rec.get("fallback_fraction", 0.0),
                              rec["seconds"]))
                else:
                    n_bad += 1
                    msg = "FAILED " + str(rec.get("error", ""))[:80]
                if verbose:
                    elapsed = time.time() - t0
                    eta = elapsed / done_n * (total - done_n)
                    print("[%4d/%d] %-34s %s  (%d running, %.0f GB, "
                          "elapsed %.1f h, eta %.1f h)"
                          % (done_n, total, rec["tile"], msg, len(running),
                             committed / (1 << 30), elapsed / 3600,
                             eta / 3600), flush=True)

    if verbose:
        fin = [r for r in records.values() if r.get("status") == "done"]
        if fin:
            fb = sorted(r.get("fallback_fraction", 0.0) for r in fin)
            print("\ninverse-distance fallback over %d tiles: median %.1f%%, p90 %.1f%%, max %.1f%%"
                  % (len(fb), 100 * fb[len(fb) // 2],
                     100 * fb[int(0.9 * (len(fb) - 1))], 100 * fb[-1]))
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
    p.add_argument("--memory-gb", type=float, default=None,
                   help="memory the run may commit at once "
                        "(default: three quarters of installed RAM)")
    p.add_argument("--redo", action="store_true",
                   help="reprocess tiles already recorded as done")
    a = p.parse_args(argv)

    paths = []
    for item in a.inputs:
        if os.path.isdir(item):
            paths += sorted(glob.glob(os.path.join(item, "**", "*.la[sz]"),
                                      recursive=True))
        else:
            # recursive, so a ** in the pattern means what it looks like it
            # means. Without it "root/**/*.las" quietly matches nothing and the
            # run exits having done no work.
            hits = sorted(glob.glob(item, recursive=True))
            paths += hits if hits else [item]
    paths = [q for q in dict.fromkeys(paths) if os.path.isfile(q)]
    if not paths:
        p.error("no LAS or LAZ files matched")

    run(paths, a.out, preset=a.preset, workers=a.workers, cell_m=a.cell,
        make_g1=a.g1, keep_ground=a.keep_ground, rvt_script=a.rvt_script,
        redo=a.redo, memory_gb=a.memory_gb)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
