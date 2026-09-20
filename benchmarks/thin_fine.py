r"""Does declustered thinning to Liao's optimum clear the kriging bottleneck?

Four things decide it, and they can disagree:

  1 does the fallback rate fall from 88% toward zero
  2 does the variogram fit and the matrix solve
  3 does elevation error stay near 0.040 m and complexity near 0.96x
  4 does the variance raster show coherent structure rather than a fallback mask

Liao's figure for a metre-scale DEM is 2.08 points per square metre. This tile
carries ten. Whether that figure transfers from Sichuan landslides to Maya karst
is what is being tested, not assumed.
"""
import json, os, sys, time
from dataclasses import replace
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import ROOT, SRC, DEM_ROOT, LAS_ROOT, TESTS
sys.path.insert(0, SRC)
import numpy as np
from replicalm import calibrate, classify, grid as G, kriging as K, residual, interpolate
from replicalm.config import PRESETS, GroundPass

OUT = TESTS + r"\\thin_fine"
os.makedirs(OUT, exist_ok=True)
with open(TESTS + r"\\clips\index.json", encoding="utf-8") as fh:
    clips = {c["tile"]: c for c in json.load(fh)}
c = clips["l0s444"]
g = G.grid_from_raster(c["clipped_reference"])
wkt = interpolate.source_srs(c["clipped_las"])
area = (g.bounds[2] - g.bounds[0]) * (g.bounds[3] - g.bounds[1])

cfg = replace(PRESETS["ncalm"],
              passes=[GroundPass(algorithm="smrf", slope=0.05,
                                 threshold_m=0.25)],
              remove_low_noise=False, remove_outliers=False)
gnd_path = os.path.join(TESTS + r"\\thin", "ground.las")
if not os.path.exists(gnd_path):
    classify.classify_tile(c["clipped_las"], gnd_path, cfg, verbose=False)
arr, _ = classify.read_points(gnd_path)
gnd = arr[arr["Classification"] == 2]
X, Y, Z = gnd["X"], gnd["Y"], gnd["Z"]
print("full set: %d ground points, %.2f pts/m2\n" % (len(Z), len(Z) / area))

# Targets chosen so the ACHIEVED densities span 2.66 to 8.54, the gap
# where the fallback is already solved and the complexity ratio has not
# yet fallen away. Occupancy drops as cells shrink, because the returns
# are clustered: a target of 5 yielded 2.66, a target of 1 yielded 0.83.
# These are set high to compensate.
TARGETS = [None, 22.0, 17.0, 13.0, 10.0, 8.0, 6.0]
print("%-9s %9s %8s %8s %8s %8s %8s %9s %9s"
      % ("target", "points", "pts/m2", "range m", "radius", "fallback",
         "rmse m", "cplx ratio", "var cover"))
print("-" * 92)
rows = []
for t in TARGETS:
    if t is None:
        sel = np.arange(len(Z)); label = "full"
    else:
        sel = K.thin_declustered(X, Y, t); label = "%.2f" % t
    x, y, z = X[sel], Y[sel], Z[sel]
    dens = len(z) / area
    tag = "thin_%s" % label.replace(".", "")
    dem_p = os.path.join(OUT, "%s.tif" % tag)
    var_p = os.path.join(OUT, "%s_var.tif" % tag)
    try:
        v = K.fit_variogram(x, y, z, model="spherical", max_lag=10.0)
        rng_m = v["params"][2]
        radius = K.radius_from_variogram(v)
        dem, info = K.krige_grid(x, y, z, g, radius=radius, max_points=16,
                                 variogram=v, verbose=False)
        K.write_geotiff(dem, g, wkt, dem_p)
        var = info["variance"]
        K.write_geotiff(var, g, wkt, var_p)
        sc = calibrate.compare_to_reference(dem_p, c["clipped_reference"])
        ratio = sc["complexity_candidate"] / sc["complexity_reference"]
        vmask = np.isfinite(var) & (var > -9000)
        rows.append({"target": t, "points": int(len(z)), "density": dens,
                     "range_m": rng_m, "radius_m": radius,
                     "fallback": info["fallback_fraction"],
                     "rmse_m": sc["rmse_m"], "ratio": ratio,
                     "variance_coverage": float(vmask.mean()),
                     "variance_mean": float(var[vmask].mean()) if vmask.any() else None})
        print("%-9s %9d %8.2f %8.2f %8.1f %7.1f%% %8.3f %9.2fx %8.1f%%"
              % (label, len(z), dens, rng_m, radius,
                 100 * info["fallback_fraction"], sc["rmse_m"], ratio,
                 100 * vmask.mean()))
    except Exception as e:
        print("%-9s FAILED %s" % (label, str(e)[:56]))

with open(os.path.join(OUT, "thin.json"), "w", encoding="utf-8") as fh:
    json.dump(rows, fh, indent=1, default=float)

print()
print("criterion 4: spatial structure of the variance raster")
for r in rows:
    if r.get("variance_coverage", 0) > 0.5:
        p = os.path.join(OUT, "thin_%s_var.tif"
                         % ("full" if r["target"] is None
                            else ("%.2f" % r["target"]).replace(".", "")))
        try:
            d = residual.compare(p, p, max_lag_m=15.0)
        except Exception:
            pass
print("  variance rasters written beside each DEM for inspection")
