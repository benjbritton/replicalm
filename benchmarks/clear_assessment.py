r"""Clear at its 9/23 configuration, with the de-duplicated covariance solve.

Run before committing days of machine time to the full corpus. Three questions:

  does it still fall back?   The solve now de-duplicates coincident locations and
                             builds the system from covariances, so the inverse
                             distance path should be idle. Recorded per tile.

  did the surfaces move?     Against the archetype's own rasters, so the answer is
                             in the same units the replication claim is made in.

  what does it cost?         Wall time per tile, to size the full run from
                             measurement rather than from extrapolation.

Clear is the configuration as locked on 2026-09-23: SMRF one pass, slope 0.1584,
elevation threshold 0.50 m, working cell 0.50 m, 16 neighbours, minimum 3,
cleanup at 0.20 m over a 0.75 m patch. Nothing here varies it.
"""
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import SRC, TESTS
sys.path.insert(0, SRC)
import numpy as np
from replicalm import (calibrate, classify, finalise, grid as G, interpolate,
                       kriging as K)
from replicalm.config import PRESETS, BASELINE, verify_baseline

OUT = os.path.join(TESTS, "clear_assessment")
os.makedirs(OUT, exist_ok=True)
with open(os.path.join(TESTS, "clips", "index.json"), encoding="utf-8") as fh:
    clips = {c["tile"]: c for c in json.load(fh)}
RADIUS = {"l0s395": 5.0, "l8s431": 20.0, "l0s444": 5.0}

cfg = PRESETS["clear"]
verify_baseline(cfg, profile="clear")
print("Clear, baseline locked %s: threshold %.2f m, %d neighbours, minimum %d, "
      "cleanup %s at %.2f m"
      % (BASELINE["locked"], cfg.passes[0].threshold_m, cfg.max_points,
         cfg.min_points, cfg.clean_vegetation, cfg.clean_height_m))
print("\n%-8s %9s %8s %7s %8s %9s %9s %8s %7s"
      % ("tile", "ground", "dropped", "idw %", "clear %", "rmse m", "bias m",
         "cplx", "sec"))
print("-" * 84)

rows = []
for tile in ("l0s395", "l8s431", "l0s444"):
    c = clips[tile]
    t0 = time.time()
    g = G.grid_from_raster(c["clipped_reference"])
    wkt = interpolate.source_srs(c["clipped_las"])
    las_p = os.path.join(OUT, tile + "_clear.las")
    classify.classify_tile(c["clipped_las"], las_p, cfg, verbose=False)
    arr, _ = classify.read_points(las_p)

    from replicalm import cleanup
    before = int((arr["Classification"] == 2).sum())
    arr, report = cleanup.clear(arr, patch=cfg.clean_patch_m,
                                height=cfg.clean_height_m,
                                percentile=cfg.clean_percentile, verbose=False)
    gnd = arr[arr["Classification"] == 2]
    x, y, z = gnd["X"], gnd["Y"], gnd["Z"]
    # how many coincident locations the solve will drop
    dropped = len(x) - len(np.unique(np.c_[x, y], axis=0))

    v = K.fit_variogram(x, y, z, model="spherical", max_lag=RADIUS[tile])
    dem, info = K.krige_grid(x, y, z, g, radius=RADIUS[tile],
                             max_points=cfg.max_points,
                             min_points=cfg.min_points, variogram=v,
                             verbose=False)
    # Finalise, as the shipped pipeline does: fill the holes kriging leaves
    # inside coverage, trim the edge back by the search radius so the surface
    # does not extend past the returns that support it, and write with the
    # nodata convention. Skipping it -- which an earlier version of this script
    # did -- leaves a raw kriged grid whose edges reach beyond the data and
    # whose gaps are untreated, and renders visualizations from that.
    dem_p = os.path.join(OUT, tile + "_clear.tif")
    fin = finalise.finalise(dem, g, wkt, dem_p, radius_m=RADIUS[tile],
                            nodata_in=-9999.0, erode_factor=cfg.erode_factor,
                            verbose=False)
    sc = calibrate.compare_to_reference(dem_p, c["clipped_reference"])
    secs = time.time() - t0
    ratio = sc["complexity_candidate"] / sc["complexity_reference"]
    rows.append({"tile": tile, "ground_points": int(len(x)),
                 "coincident_dropped": int(dropped),
                 "clear_removed_pct": round(100 * report["fraction"], 2),
                 "fallback_fraction": info["fallback_fraction"],
                 "variogram_range_m": float(info["params"][2]),
                 "seconds": round(secs, 1), "ratio": ratio,
                 "filled_cells": fin["filled_cells"], "holes": fin["holes"],
                 "erode_cells": fin["erode_cells"],
                 "cells_after_trim": fin["cells_after_trim"], **sc})
    print("%-8s %9d %8d %6.1f%% %7.2f%% %9.4f %+9.4f %7.2fx %7.0f"
          % (tile, len(x), dropped, 100 * info["fallback_fraction"],
             100 * report["fraction"], sc["rmse_m"], sc["bias_m"], ratio, secs),
          flush=True)
    print("         finalise: %d holes filled (%d cells), edge trimmed %d cells, "
          "%d cells kept" % (fin["holes"], fin["filled_cells"],
                             fin["erode_cells"], fin["cells_after_trim"]),
          flush=True)

with open(os.path.join(OUT, "clear_assessment.json"), "w", encoding="utf-8") as fh:
    json.dump(rows, fh, indent=1, default=float)
print("\nwrote %s" % os.path.join(OUT, "clear_assessment.json"))

fb = [r["fallback_fraction"] for r in rows]
print("\nfallback across the test set: max %.2f%%" % (100 * max(fb)))
print("for comparison, the 9/23 solve on l0s444 fell back on 89.0%% of cells")
