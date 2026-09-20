r"""Which did the work: the tighter search, or the tighter variogram?

The previous sweep refitted the variogram to a lag equal to the search radius,
so the two moved together and the 40% improvement at 5 m could belong to either.
This crosses them: every combination of fitting lag and search radius, on one
classification, so each factor can be read on its own.
"""
import itertools, json, os, sys, time
from dataclasses import replace
sys.path.insert(0, r"C:\Replicalm\src")
import numpy as np
from replicalm import calibrate, classify, grid as G, kriging as K, interpolate
from replicalm.config import PRESETS

with open(r"C:\Replicalm\tests\clips\index.json", encoding="utf-8") as fh:
    clips = {c["tile"]: c for c in json.load(fh)}
c = clips["l0s395"]
OUT = r"C:\Replicalm\tests\sweep_decouple"
os.makedirs(OUT, exist_ok=True)

base = PRESETS["ncalm"]
cfg = replace(base, passes=[replace(p, slope=0.10, threshold_m=1.0)
                            for p in base.passes])
g = G.grid_from_raster(c["clipped_reference"])
wkt = interpolate.source_srs(c["clipped_las"])

gnd_path = os.path.join(OUT, "ground.las")
classify.classify_tile(c["clipped_las"], gnd_path, cfg, verbose=False)
arr, _ = classify.read_points(gnd_path)
gnd = arr[arr["Classification"] == 2]
x, y, z = gnd["X"], gnd["Y"], gnd["Z"]
print("ground points: %d\n" % len(x))

FIT_LAGS = (5.0, 20.0)
RADII = (5.0, 20.0)

vgs = {}
for L in FIT_LAGS:
    v = K.fit_variogram(x, y, z, model="spherical", max_lag=L)
    vgs[L] = v
    nug, sill, rng = v["params"]
    print("variogram fitted to %4.0f m: nugget %.4f  sill %.4f  range %6.2f m"
          % (L, nug, sill, rng))

print("\n%9s %8s %9s %11s %11s" % ("fit lag", "radius", "rmse m",
                                   "cplx cand", "cplx loss"))
print("-" * 54)
rows = []
for L, R in itertools.product(FIT_LAGS, RADII):
    tag = "fit%g_r%g" % (L, R)
    p = os.path.join(OUT, "dem_%s.tif" % tag)
    dem, _ = K.krige_grid(x, y, z, g, radius=R, max_points=16,
                          variogram=vgs[L], verbose=False)
    K.write_geotiff(dem, g, wkt, p)
    sc = calibrate.compare_to_reference(p, c["clipped_reference"])
    rows.append({"fit_lag": L, "radius": R, **sc})
    print("%9.0f %8.0f %9.3f %11.4f %11.3f"
          % (L, R, sc["rmse_m"], sc["complexity_candidate"],
             sc["complexity_loss"]))

print("\nreference complexity %.4f" % rows[0]["complexity_reference"])
print("\neffect of each factor, holding the other:")
for L in FIT_LAGS:
    a = [r for r in rows if r["fit_lag"] == L]
    a.sort(key=lambda r: r["radius"])
    print("  fit lag %4.0f m: radius 5 -> 20 changes rmse %.3f -> %.3f"
          % (L, a[0]["rmse_m"], a[-1]["rmse_m"]))
for R in RADII:
    a = [r for r in rows if r["radius"] == R]
    a.sort(key=lambda r: r["fit_lag"])
    print("  radius  %4.0f m: fit 5 -> 20 changes rmse %.3f -> %.3f"
          % (R, a[0]["rmse_m"], a[-1]["rmse_m"]))

with open(os.path.join(OUT, "decouple.json"), "w", encoding="utf-8") as fh:
    json.dump(rows, fh, indent=1, default=float)
