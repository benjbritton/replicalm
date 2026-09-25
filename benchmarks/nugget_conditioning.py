r"""Does a nugget let 64 neighbours actually krige?

At the source's 64 points per estimate, 90 to 99% of cells in the production run
fell back to inverse distance: 64 returns at 9 per square metre all sit within
about 1.5 m, the 65x65 system is built from near-identical semivariances, and it
loses rank. The surfaces were weighted averages wearing the name of kriging.

This is not a property of kriging with 64 neighbours. It is a property of
kriging with 64 neighbours and no nugget. A nugget adds a constant to the
diagonal, which is what conditions the system -- and the reference implementation
fits one per tile, which is the thing the Ka'Kabish incident was about.

So the question is whether the reference's neighbourhood becomes workable with
the reference's nugget, or whether 16 remains the only count this implementation
can actually solve. Everything is held fixed except the nugget floor and the
neighbour count: same classified cloud, same grid, same radius, same model.

Reported per arm: the share of cells that fell back, and whether the resulting
surface is closer to or further from the commercial product.
"""
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import SRC, TESTS
sys.path.insert(0, SRC)
import numpy as np
from replicalm import calibrate, classify, grid as G, kriging as K, interpolate
from replicalm.config import PRESETS

OUT = os.path.join(TESTS, "nugget_conditioning")
os.makedirs(OUT, exist_ok=True)
with open(os.path.join(TESTS, "clips", "index.json"), encoding="utf-8") as fh:
    clips = {c["tile"]: c for c in json.load(fh)}
RADIUS = {"l0s395": 5.0, "l8s431": 20.0, "l0s444": 5.0}

# (neighbours, nugget as a fraction of the sill)
ARMS = [(16, 0.00), (64, 0.00), (64, 0.01), (64, 0.02), (64, 0.05), (64, 0.10)]

print("%-8s %6s %8s %8s %9s %9s %8s %7s"
      % ("tile", "pts", "nugget", "idw %", "rmse m", "bias m", "cplx", "sec"))
print("-" * 74)
rows = []
for tile in ("l0s395", "l8s431", "l0s444"):
    c = clips[tile]
    g = G.grid_from_raster(c["clipped_reference"])
    wkt = interpolate.source_srs(c["clipped_las"])
    las_p = os.path.join(OUT, tile + "_baseline.las")
    if not os.path.exists(las_p):
        classify.classify_tile(c["clipped_las"], las_p, PRESETS["baseline"],
                               verbose=False)
    arr, _ = classify.read_points(las_p)
    gnd = arr[arr["Classification"] == 2]
    x, y, z = gnd["X"], gnd["Y"], gnd["Z"]
    for maxp, nug in ARMS:
        t0 = time.time()
        v = K.fit_variogram(x, y, z, model="spherical", max_lag=RADIUS[tile],
                            min_nugget_fraction=nug)
        dem, info = K.krige_grid(x, y, z, g, radius=RADIUS[tile],
                                 max_points=maxp, min_points=1, variogram=v,
                                 verbose=False)
        dem_p = os.path.join(OUT, "%s_%d_%03d.tif" % (tile, maxp, nug * 100))
        K.write_geotiff(dem, g, wkt, dem_p)
        sc = calibrate.compare_to_reference(dem_p, c["clipped_reference"])
        rows.append({"tile": tile, "max_points": maxp, "nugget_fraction": nug,
                     "fallback_fraction": info["fallback_fraction"],
                     "variogram_range_m": float(info["params"][2]),
                     "seconds": round(time.time() - t0, 1),
                     "ratio": sc["complexity_candidate"] / sc["complexity_reference"],
                     **sc})
        print("%-8s %6d %8.2f %8.1f %9.4f %+9.4f %7.2fx %7.0f"
              % (tile, maxp, nug, 100 * info["fallback_fraction"], sc["rmse_m"],
                 sc["bias_m"], rows[-1]["ratio"], rows[-1]["seconds"]),
              flush=True)
    print()

with open(os.path.join(OUT, "nugget_conditioning.json"), "w",
          encoding="utf-8") as fh:
    json.dump(rows, fh, indent=1, default=float)
print("wrote %s (%d arms)"
      % (os.path.join(OUT, "nugget_conditioning.json"), len(rows)))

print("\nfallback at 64 neighbours, by nugget:")
for nug in (0.00, 0.01, 0.02, 0.05, 0.10):
    v = [r["fallback_fraction"] for r in rows
         if r["max_points"] == 64 and r["nugget_fraction"] == nug]
    if v:
        print("  %.2f -> %5.1f%% median across %d tiles"
              % (nug, 100 * float(np.median(v)), len(v)))
