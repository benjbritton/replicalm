r"""All three calibration tiles, radius derived from each tile's own variogram.

The question is whether the correlation range varies enough between tiles that
deriving the radius matters, or whether one fixed value would have done.
"""
import json, os, sys, time
from dataclasses import replace
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import ROOT, SRC, DEM_ROOT, LAS_ROOT, TESTS
sys.path.insert(0, SRC)
from replicalm import calibrate, classify, grid as G, kriging as K, interpolate
from replicalm.config import PRESETS

with open(TESTS + r"\\clips\index.json", encoding="utf-8") as fh:
    clips = json.load(fh)

base = PRESETS["ncalm"]
cfg = replace(base, passes=[replace(p, slope=0.10, threshold_m=1.0)
                            for p in base.passes])
OUT = TESTS + r"\\three_tiles"
os.makedirs(OUT, exist_ok=True)

print("%-8s %10s %8s %8s %9s %8s %9s %9s %8s"
      % ("tile", "ground", "pts/m2", "range m", "radius m", "rmse m",
         "cplx ref", "cplx cand", "filled%"))
print("-" * 92)
rows = []
for c in clips:
    t0 = time.time()
    tile = c["tile"]
    g = G.grid_from_raster(c["clipped_reference"])
    wkt = interpolate.source_srs(c["clipped_las"])
    gnd_path = os.path.join(OUT, "%s_ground.las" % tile)
    classify.classify_tile(c["clipped_las"], gnd_path, cfg, verbose=False)
    arr, _ = classify.read_points(gnd_path)
    gnd = arr[arr["Classification"] == 2]
    area = (g.bounds[2] - g.bounds[0]) * (g.bounds[3] - g.bounds[1])

    v = K.fit_variogram(gnd["X"], gnd["Y"], gnd["Z"], model="spherical",
                        max_lag=5.0)
    radius = K.radius_from_variogram(v)
    dem, _ = K.krige_grid(gnd["X"], gnd["Y"], gnd["Z"], g, radius=radius,
                          max_points=16, variogram=v, verbose=False)
    p = os.path.join(OUT, "%s_dem.tif" % tile)
    K.write_geotiff(dem, g, wkt, p)
    sc = calibrate.compare_to_reference(p, c["clipped_reference"])
    rows.append({"tile": tile, "region": c["region"],
                 "ground_points": int(len(gnd)),
                 "density": len(gnd) / area,
                 "variogram_range_m": v["params"][2],
                 "radius_m": radius, "seconds": round(time.time() - t0, 1), **sc})
    print("%-8s %10d %8.2f %8.2f %9.2f %8.3f %9.4f %9.4f %7.1f%%"
          % (tile, len(gnd), len(gnd) / area, v["params"][2], radius,
             sc["rmse_m"], sc["complexity_reference"],
             sc["complexity_candidate"], 100 * sc["coverage_candidate"]))

with open(os.path.join(OUT, "three_tiles.json"), "w", encoding="utf-8") as fh:
    json.dump(rows, fh, indent=1, default=float)

rr = [r["variogram_range_m"] for r in rows]
print("\ncorrelation range across tiles: %.2f to %.2f m (%.1fx spread)"
      % (min(rr), max(rr), max(rr) / min(rr)))
print("rmse: %s" % ", ".join("%s %.3f" % (r["tile"], r["rmse_m"]) for r in rows))
