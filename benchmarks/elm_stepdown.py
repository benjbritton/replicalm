r"""How far can ELM be pushed on l0s444 before it reaches the 3 m structure?

The residual there is a third noise over a real signal with a 3.0 m range. That
range is the yardstick: while it holds, the filter is taking noise. If it moves
or disappears, the filter has started clipping the morphology the signal is made
of, and the nugget falling is no longer good news.
"""
import json, os, sys, time
from dataclasses import replace
sys.path.insert(0, r"C:\Replicalm\src")
import numpy as np
from replicalm import calibrate, classify, grid as G, kriging as K, residual, interpolate
from replicalm.config import PRESETS, GroundPass

with open(r"C:\Replicalm\tests\clips\index.json", encoding="utf-8") as fh:
    clips = {c["tile"]: c for c in json.load(fh)}
TILE = sys.argv[1] if len(sys.argv) > 1 else "l0s444"
OUT = r"C:\Replicalm\tests\elm_step_" + TILE
os.makedirs(OUT, exist_ok=True)
c = clips[TILE]
g = G.grid_from_raster(c["clipped_reference"])
wkt = interpolate.source_srs(c["clipped_las"])
area = (g.bounds[2] - g.bounds[0]) * (g.bounds[3] - g.bounds[1])

BEST_FILTER = {
    "l0s444": dict(algorithm="smrf", slope=0.05, threshold_m=0.25),
    "l0s395": dict(algorithm="csf", csf_rigidness=2,
                   csf_threshold_m=0.5),
    "l8s431": dict(algorithm="smrf", slope=0.05, threshold_m=0.5),
}
RADIUS = {"l0s444": 5.0, "l0s395": 5.0, "l8s431": 20.0}
THRESHOLDS = (None, 1.0, 0.5, 0.3, 0.25, 0.1)
print("baseline: nugget/sill 0.36, range 3.0 m, complexity ratio 0.96\n")
print("%-9s %10s %9s %8s %9s %9s %8s %9s"
      % ("elm thr", "ground", "removed", "rmse m", "cplx", "ratio",
         "nug/sill", "range m"))
print("-" * 82)
rows = []
base_n = None
for thr in THRESHOLDS:
    on = thr is not None
    cfg = replace(PRESETS["ncalm"],
                  passes=[GroundPass(**BEST_FILTER[TILE])],
                  remove_low_noise=on, remove_outliers=on,
                  elm_threshold_m=thr if on else 1.0)
    tag = "thr%s" % ("off" if not on else str(thr).replace(".", ""))
    las_p = os.path.join(OUT, "%s.las" % tag)
    dem_p = os.path.join(OUT, "%s.tif" % tag)
    try:
        classify.classify_tile(c["clipped_las"], las_p, cfg, verbose=False)
        arr, _ = classify.read_points(las_p)
        gnd = arr[arr["Classification"] == 2]
        if base_n is None:
            base_n = len(gnd)
        v = K.fit_variogram(gnd["X"], gnd["Y"], gnd["Z"], model="spherical",
                            max_lag=RADIUS[TILE])
        dem, _ = K.krige_grid(gnd["X"], gnd["Y"], gnd["Z"], g,
                              radius=RADIUS[TILE],
                              max_points=16, variogram=v, verbose=False)
        K.write_geotiff(dem, g, wkt, dem_p)
        sc = calibrate.compare_to_reference(dem_p, c["clipped_reference"])
        rs = residual.compare(dem_p, c["clipped_reference"], max_lag_m=15.0)
        ratio = sc["complexity_candidate"] / sc["complexity_reference"]
        rows.append({"elm_threshold": thr, "ground": int(len(gnd)),
                     "ratio": ratio, **rs})
        print("%-9s %10d %9d %8.3f %9.4f %9.2fx %8.2f %9s"
              % ("off" if not on else "%.2f" % thr, len(gnd),
                 base_n - len(gnd), sc["rmse_m"], sc["complexity_candidate"],
                 ratio, rs["nugget_to_sill"],
                 "%.1f" % rs["range_m"] if rs["range_m"] else "none"))
    except Exception as e:
        print("%-9s FAILED %s" % (str(thr), str(e)[:52]))

with open(os.path.join(OUT, "stepdown.json"), "w", encoding="utf-8") as fh:
    json.dump(rows, fh, indent=1, default=float)
print()
print("the range holding near 3.0 m means the filter is still taking noise;")
print("a range that moves or vanishes means it has reached the morphology.")
