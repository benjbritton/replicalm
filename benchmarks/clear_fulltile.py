r"""Clear over three whole G-LiHT tiles, end to end, as the tool ships it.

Whole transects, not clipped windows, because the edge treatment and the nodata
convention are part of what has to be looked at and a crop shows neither. The
DEM comes from `pipeline.process`, which is the path the application runs --
classify, cleanup, krige, fill holes, trim the coverage edge by the search
radius, write with the nodata convention -- rather than from those stages
reassembled here, which is how an earlier version of this script omitted
finalise entirely.

Each tile is then compared against the archetype's own 0.5 m raster and its six
visualizations are built by the published recipe. Nothing in the chain is
optional or run by hand afterwards.
"""
import argparse
import json
import os
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import SRC, TESTS
sys.path.insert(0, SRC)
from replicalm import calibrate, pipeline
from replicalm.config import PRESETS, BASELINE, verify_baseline

# Taken from tests/clips/index.json, which is the authoritative record of
# which tile each test result refers to. Tile ids repeat across regions --
# there are two l0s395 and three l0s444 -- so matching on the id alone picks
# whichever sorts first, which is how the wrong l0s395 got processed once.
# References are the edge-fixed archive, the same rasters every earlier
# measurement was scored against.
TILES = {
    "l0s395": (r"D:\GLiHT_LAS_orig\Yuc_South\South_Glas\South_Glas\AMIGACarb_Yuc_South_GLAS_Apr2013_l0s395.las",
               r"D:\_Archive_EdgeFixed\Yuc_South\South_Glas\South_Glas\South_GLAS_l0s395_DEM_0p5m_v1.tif"),
    "l8s431": (r"D:\GLiHT_LAS_orig\Yuc_Centro\Centro_NFI\Centro_NFI\AMIGACarb_Yuc_Centro_NFI_Apr2013_l8s431.las",
               r"D:\_Archive_EdgeFixed\Yuc_Centro\Centro_NFI\Centro_NFI\Centro_NFI_l8s431_DEM_0p5m_v1.tif"),
    "l0s444": (r"D:\GLiHT_LAS_orig\Yuc_Centro\Centro_Glas\Centro_Glas\AMIGACarb_Yuc_Centro_GLAS_Apr2013_l0s444.las",
               r"D:\_Archive_EdgeFixed\Yuc_Centro\Centro_Glas\Centro_Glas\Centro_GLAS_l0s444_DEM_0p5m_v1.tif"),
}
OUT = os.path.join(TESTS, "clear_fulltile")
# RVT lives in the g1 environment, not in the one that runs PDAL, so the image
# step is a subprocess against that interpreter.
G1_PY = os.environ.get("REPLICALM_G1_PYTHON",
                       r"C:\Users\benja\anaconda3\envs\g1\python.exe")
RVT = os.environ.get("REPLICALM_RVT_SCRIPT", r"C:\g1\tools\GLiHT_rvt.py")


def visualize(dem_path, tile_dir, tile):
    """The six products, from the finalised DEM, by the published recipe."""
    staged = os.path.join(tile_dir, "_rvt_in", "tile")
    os.makedirs(staged, exist_ok=True)
    link = os.path.join(staged, "%s_clear_DEM_0p5m_v1.tif" % tile)
    if not os.path.exists(link):
        import shutil
        shutil.copy2(dem_path, link)
    rvt_root = os.path.join(tile_dir, "rvt")
    # The g1 interpreter must not inherit this one's environment. PYTHONPATH
    # points at Replicalm's src, and GDAL_DATA and PROJ_LIB point into the
    # ArcGIS environment; the child picks them up and its own GDAL fails. The
    # step then returns exit 1 and, because this captured the output without
    # printing it, said nothing about why.
    env = {k: v for k, v in os.environ.items()
           if k not in ("PYTHONPATH", "PYTHONHOME", "GDAL_DATA", "PROJ_LIB",
                        "GDAL_DRIVER_PATH", "PROJ_DATA")}
    # Clearing those is not enough: PATH still leads with the ArcGIS binaries,
    # so the g1 interpreter loads their GDAL DLLs and fails with
    # "No module named '_gdal'". Its own directories have to come first.
    g1_home = os.path.dirname(G1_PY)
    env["PATH"] = os.pathsep.join(
        [g1_home, os.path.join(g1_home, "Library", "bin"),
         os.path.join(g1_home, "Library", "usr", "bin"),
         os.path.join(g1_home, "Scripts"), env.get("PATH", "")])
    proc = subprocess.run(
        [G1_PY, RVT, "--res", "0p5m", "--no-metadata", "--workers", "1",
         "--dem-root", os.path.join(tile_dir, "_rvt_in"),
         "--rvt-root", rvt_root],
        capture_output=True, text=True, env=env)
    if proc.returncode != 0:
        tail = (proc.stderr or proc.stdout or "").strip().splitlines()[-4:]
        print("  [visualize  ] FAILED: %s" % " / ".join(tail), flush=True)
    made = []
    for root, _dirs, files in os.walk(rvt_root):
        made += [os.path.join(root, f) for f in files if f.endswith(".tif")]
    return rvt_root, sorted(made), proc.returncode


MANIFEST = os.path.join(OUT, "clear_fulltile.json")


def save(rows):
    """Merge this run's rows into the manifest, keyed by tile.

    Rebuilt from empty, this file records only whatever the last invocation
    happened to process. That is how l0s395's row came to describe North2: a
    three-tile run overwrote a corrected single-tile one, and the record then
    named a LAS that was never processed into the DEM sitting beside it. Rows
    for tiles this run did not touch are preserved.
    """
    keep = {}
    if os.path.exists(MANIFEST):
        try:
            with open(MANIFEST, encoding="utf-8") as fh:
                keep = {r["tile"]: r for r in json.load(fh)}
        except (ValueError, KeyError):
            keep = {}          # unreadable is not a reason to lose this run
    for r in rows:
        keep[r["tile"]] = r
    with open(MANIFEST, "w", encoding="utf-8") as fh:
        json.dump([keep[k] for k in sorted(keep)], fh, indent=1, default=float)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--tiles", nargs="*", default=sorted(TILES))
    ap.add_argument("--on-grid", action="store_true",
                    help="also build a second DEM on the archive's own lattice "
                         "and score against that, since a production raster and "
                         "an archive raster never share a grid")
    a = ap.parse_args(argv)

    cfg = PRESETS["clear"]
    verify_baseline(cfg, profile="clear")
    os.makedirs(OUT, exist_ok=True)
    print("Clear, locked %s: threshold %.2f m, %d neighbours, min %d, cleanup "
          "at %.2f m over %.2f m"
          % (BASELINE["locked"], cfg.passes[0].threshold_m, cfg.max_points,
             cfg.min_points, cfg.clean_height_m, cfg.clean_patch_m), flush=True)

    rows = []
    for tile in a.tiles:
        las, ref = TILES[tile]
        tile_dir = os.path.join(OUT, tile)
        os.makedirs(tile_dir, exist_ok=True)
        t0 = time.time()
        print("\n=== %s  (%.0f MB) ===" % (tile, os.path.getsize(las) / 1e6),
              flush=True)

        def show(stage, message, fraction=None):
            print("  [%-11s] %s" % (stage, message), flush=True)

        out = pipeline.process(las, tile_dir, cfg=cfg, make_g1=False,
                               keep_ground=bool(a.on_grid), progress=show)
        rec = {"tile": tile, "las": las, "reference": ref,
               "ground_points": out["ground_points"],
               "density": round(out["density"], 3),
               "cell_m": round(out["cell_m"], 4),
               "search_radius_m": round(out["search_radius_m"], 3),
               "fallback_fraction": round(out.get("fallback_fraction", 0.0), 5),
               "variogram_range_m": out.get("variogram_range_m"),
               "holes": out.get("holes"), "filled_cells": out.get("filled_cells"),
               "erode_cells": out.get("erode_cells"),
               "cells_after_trim": out.get("cells_after_trim"),
               "dem": out["dem"]}
        # The production raster derives its cell from measured density and snaps
        # its origin, so it never shares a lattice with the archive's, and
        # compare_to_reference refuses rather than score a half-cell offset.
        # A second candidate is built on the archive's grid for the comparison;
        # classification is reused, so it costs the interpolation only.
        scored = out["dem"]
        if a.on_grid:
            grid_dir = os.path.join(tile_dir, "on_archive_grid")
            os.makedirs(grid_dir, exist_ok=True)
            # reuse_ground looks beside its own output, so the classified ground
            # has to be there. A hard link costs nothing and no bytes; without
            # it the second pass silently reclassifies and the whole point of
            # reusing -- that classification does not depend on the lattice --
            # is lost.
            stem = os.path.splitext(os.path.basename(las))[0]
            src_g = os.path.join(tile_dir, stem + "_ground.laz")
            dst_g = os.path.join(grid_dir, stem + "_ground.laz")
            if os.path.exists(src_g) and not os.path.exists(dst_g):
                try:
                    os.link(src_g, dst_g)
                except OSError:
                    import shutil
                    shutil.copy2(src_g, dst_g)
            og = pipeline.process(las, grid_dir, cfg=cfg, make_g1=False,
                                  keep_ground=True, reuse_ground=True,
                                  on_grid=ref, progress=show)
            scored = og["dem"]
            rec["dem_on_archive_grid"] = scored
            rec["cell_m_on_archive_grid"] = round(og["cell_m"], 6)
        rec["scored"] = scored
        try:
            sc = calibrate.compare_to_reference(scored, ref)
            rec.update({k: sc[k] for k in sc})
            rec["ratio"] = sc["complexity_candidate"] / sc["complexity_reference"]
            print("  [compare    ] rmse %.4f m, bias %+.4f m, texture %.2fx"
                  % (sc["rmse_m"], sc["bias_m"], rec["ratio"]), flush=True)
        except Exception as exc:
            rec["compare_error"] = str(exc)[:200]
            print("  [compare    ] FAILED %s" % str(exc)[:90], flush=True)

        rvt_root, made, code = visualize(out["dem"], tile_dir, tile)
        rec["rvt_root"] = rvt_root
        rec["visualizations"] = [os.path.basename(m) for m in made]
        print("  [visualize  ] %d products in %s (exit %d)"
              % (len(made), rvt_root, code), flush=True)
        rec["seconds"] = round(time.time() - t0, 1)
        print("  [done       ] %.0f s, fallback %.2f%%, %d cells kept"
              % (rec["seconds"], 100 * rec["fallback_fraction"],
                 rec["cells_after_trim"] or 0), flush=True)
        rows.append(rec)
        save(rows)

    print("\n%-8s %10s %7s %6s %7s %9s %9s %7s %8s"
          % ("tile", "ground", "pts/m2", "cell", "idw %", "rmse m", "bias m",
             "cplx", "sec"))
    for r in rows:
        print("%-8s %10d %7.2f %6.2f %6.2f%% %9.4f %+9.4f %6.2fx %8.0f"
              % (r["tile"], r["ground_points"], r["density"], r["cell_m"],
                 100 * r["fallback_fraction"], r.get("rmse_m", float("nan")),
                 r.get("bias_m", float("nan")), r.get("ratio", float("nan")),
                 r["seconds"]))
    print("\nwrote %s" % os.path.join(OUT, "clear_fulltile.json"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
