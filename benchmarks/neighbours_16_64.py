r"""Does the neighbour count matter? The reference grids with 64, we use 16.

The workflow being replicated searches 20 m and takes up to 64 points per
estimate (GLiHT_Methods_Materials, section 5.2, after Estrada-Belli et al.
2025). Replicalm's locked baseline takes 16, justified in the repository's own
table as binding before the radius does on well-covered ground -- which is true,
and is not the same as being what the source specifies.

More neighbours means a flatter weighting and a more averaged estimate, so the
difference should appear in terrain complexity before it appears in RMSE. This
holds everything else fixed -- same classified cloud, same variogram, same grid,
same radius -- and varies only max_points.

Minimum points also differs (reference 1, baseline 3) and is tested alongside,
since a minimum of 3 leaves cells empty that the reference would have filled.
"""
import json
import os
import sys
import time
from dataclasses import replace

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import SRC, TESTS
sys.path.insert(0, SRC)
import numpy as np
from replicalm import calibrate, classify, grid as G, kriging as K, interpolate
from replicalm.config import PRESETS

OUT = os.path.join(TESTS, "neighbours")
os.makedirs(OUT, exist_ok=True)
with open(os.path.join(TESTS, "clips", "index.json"), encoding="utf-8") as fh:
    clips = {c["tile"]: c for c in json.load(fh)}
RADIUS = {"l0s395": 5.0, "l8s431": 20.0, "l0s444": 5.0}

ARMS = [("16 pts, min 3", 16, 3),      # the locked baseline
        ("64 pts, min 3", 64, 3),      # the reference's neighbour count
        ("64 pts, min 1", 64, 1)]      # the reference's pair, both parameters

print("%-8s %-14s %9s %9s %9s %8s %9s"
      % ("tile", "arm", "rmse m", "bias m", "mae m", "cplx", "filled"))
print("-" * 76)
rows = []
for tile in ("l0s395", "l8s431", "l0s444"):
    c = clips[tile]
    g = G.grid_from_raster(c["clipped_reference"])
    wkt = interpolate.source_srs(c["clipped_las"])
    cfg = PRESETS["baseline"]
    las_p = os.path.join(OUT, tile + "_baseline.las")
    if not os.path.exists(las_p):
        classify.classify_tile(c["clipped_las"], las_p, cfg, verbose=False)
    arr, _ = classify.read_points(las_p)
    gnd = arr[arr["Classification"] == 2]
    x, y, z = gnd["X"], gnd["Y"], gnd["Z"]
    # One variogram per tile, shared by every arm: the neighbour count is the
    # only thing allowed to vary.
    v = K.fit_variogram(x, y, z, model="spherical", max_lag=RADIUS[tile])
    for name, maxp, minp in ARMS:
        t0 = time.time()
        dem, vv = K.krige_grid(x, y, z, g, radius=RADIUS[tile], max_points=maxp,
                               min_points=minp, variogram=v, verbose=False)
        dem_p = os.path.join(OUT, "%s_%s.tif"
                             % (tile, name.replace(" ", "").replace(",", "")))
        K.write_geotiff(dem, g, wkt, dem_p)
        sc = calibrate.compare_to_reference(dem_p, c["clipped_reference"])
        filled = float((dem != -9999.0).mean())
        ratio = sc["complexity_candidate"] / sc["complexity_reference"]
        rows.append({"tile": tile, "arm": name, "max_points": maxp,
                     "min_points": minp, "filled_fraction": round(filled, 5),
                     "seconds": round(time.time() - t0, 1),
                     "ratio": ratio, **sc})
        print("%-8s %-14s %9.4f %+9.4f %9.4f %7.2fx %8.4f"
              % (tile, name, sc["rmse_m"], sc["bias_m"], sc["mae_m"], ratio,
                 filled), flush=True)
    print()

with open(os.path.join(OUT, "neighbours.json"), "w", encoding="utf-8") as fh:
    json.dump(rows, fh, indent=1, default=float)
print("wrote %s (%d arms)" % (os.path.join(OUT, "neighbours.json"), len(rows)))

for tile in ("l0s395", "l8s431", "l0s444"):
    r = {x["arm"]: x for x in rows if x["tile"] == tile}
    if "16 pts, min 3" in r and "64 pts, min 3" in r:
        a, b = r["16 pts, min 3"], r["64 pts, min 3"]
        print("%-8s 16 -> 64: rmse %+.4f m, bias %+.4f m, complexity %+.2fx"
              % (tile, b["rmse_m"] - a["rmse_m"], b["bias_m"] - a["bias_m"],
                 b["ratio"] - a["ratio"]))
