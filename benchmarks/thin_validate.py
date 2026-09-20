r"""Is 4.4 points per square metre a ceiling or a target?

The fine ladder on l0s444 found 4.43 pts/m2 to be its best operating point: at
its native 8.54 the kriging system degenerated and 89% of cells fell back to
inverse distance, and below about 3 the complexity ratio started to sag. That
result came from one tile, and a tile that was anomalously dense.

l8s431 and l0s395 sit natively at 3.46 and 4.25 pts/m2 -- already inside the
corridor. So they cannot confirm the figure by being thinned to it; they can only
say whether thinning toward it is ever the right move, or whether 4.4 is simply
the point below which l0s444 stopped being pathological.

The two outcomes that matter:

  native is best          4.4 is a ceiling. Thin down to it, never up to it, and
                          a tile already below it needs nothing done.
  thinning still helps    4.4 is a target, and declustering is doing something
                          beyond relieving the matrix -- worth knowing, because
                          it would mean the even spacing matters and not just
                          the count.

Each tile is classified with its own best filter, found in the classification
sweep, not with l0s444's.
"""
import json, os, sys
from dataclasses import replace
sys.path.insert(0, r"C:\Replicalm\src")
import numpy as np
from replicalm import calibrate, classify, grid as G, kriging as K, interpolate
from replicalm.config import PRESETS, GroundPass

# The best classifier per tile, from class_sweep. These differ by tile and the
# difference is large, so reusing one tile's settings on another would measure
# the classifier rather than the thinning.
BEST = {
    "l8s431": dict(algorithm="smrf", slope=0.05, threshold_m=0.5),
    "l0s395": dict(algorithm="csf", csf_rigidness=2, csf_threshold_m=0.5),
    "l0s444": dict(algorithm="smrf", slope=0.05, threshold_m=0.25),
}

TILE = sys.argv[1] if len(sys.argv) > 1 else "l8s431"
OUT = r"C:\Replicalm\tests\thin_val_" + TILE
os.makedirs(OUT, exist_ok=True)
with open(r"C:\Replicalm\tests\clips\index.json", encoding="utf-8") as fh:
    clips = {c["tile"]: c for c in json.load(fh)}
c = clips[TILE]
g = G.grid_from_raster(c["clipped_reference"])
wkt = interpolate.source_srs(c["clipped_las"])
area = (g.bounds[2] - g.bounds[0]) * (g.bounds[3] - g.bounds[1])

cfg = replace(PRESETS["ncalm"], passes=[GroundPass(**BEST[TILE])],
              remove_low_noise=False, remove_outliers=False)
gnd_path = os.path.join(OUT, "%s_ground.las" % TILE)
if not os.path.exists(gnd_path):
    classify.classify_tile(c["clipped_las"], gnd_path, cfg, verbose=False)
arr, _ = classify.read_points(gnd_path)
gnd = arr[arr["Classification"] == 2]
X, Y, Z = gnd["X"], gnd["Y"], gnd["Z"]
native = len(Z) / area
print("%s, %s: %d ground points, %.2f pts/m2 native\n"
      % (TILE, BEST[TILE]["algorithm"], len(Z), native))

# Targets bracket the native density from both sides where that is possible, so
# a tile at 4.25 is asked both whether it wants thinning and whether it would
# have preferred more points than it has.
TARGETS = [None, 4.43, 3.5, 3.0, 2.5, 2.08]
print("%-9s %9s %8s %8s %8s %8s %8s %10s %9s"
      % ("target", "points", "pts/m2", "range m", "radius", "fallback",
         "rmse m", "cplx ratio", "var cover"))
print("-" * 92)
rows = []
for t in TARGETS:
    if t is None:
        sel = np.arange(len(Z)); label = "full"
    elif t >= native:
        print("%-9s skipped: above this tile's native %.2f pts/m2" % (t, native))
        continue
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
        rows.append({"tile": TILE, "target": t, "points": int(len(z)),
                     "density": dens, "native_density": native,
                     "range_m": rng_m, "radius_m": radius,
                     "fallback": info["fallback_fraction"],
                     "rmse_m": sc["rmse_m"], "ratio": ratio,
                     "variance_coverage": float(vmask.mean())})
        print("%-9s %9d %8.2f %8.2f %8.1f %7.1f%% %8.3f %9.2fx %8.1f%%"
              % (label, len(z), dens, rng_m, radius,
                 100 * info["fallback_fraction"], sc["rmse_m"], ratio,
                 100 * vmask.mean()))
    except Exception as e:
        print("%-9s FAILED %s" % (label, str(e)[:56]))

with open(os.path.join(OUT, "thin.json"), "w", encoding="utf-8") as fh:
    json.dump(rows, fh, indent=1, default=float)
print("\nwrote %s" % os.path.join(OUT, "thin.json"))
