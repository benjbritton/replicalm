r"""What the noise filters cost on flat ground, and which settings cost least.

THE OBSERVATION
---------------
Switching the noise filters on fixed the steep ground and put a fine speckle on
the flat. Measured tile-wide on l0s395: cells more than half a metre from the
reference ran 0.09% below 5 degrees with the filters on against 0.01% with them
off, and 0.17% against 0.04% between 5 and 10 degrees. Small, but visible in the
G1 as scattered dark pits across plaza floors.

The filters are not optional -- they are what carried the steep ground, and the
earlier attribution of that to the second SMRF pass was wrong. So the question
is not whether to run them but at what settings.

WHY A FLAT WINDOW
-----------------
This is the mirror image of the relief harness. Speckle is an artifact of flat
ground: a spurious pit is only visible where the surface around it is smooth,
and on a 30 degree flank it is indistinguishable from terrain. A window chosen
for maximum coverage -- the old selection rule -- is exactly right here, and is
why the coverage rule is kept rather than replaced.

WHAT IS MEASURED
----------------
  pits        cells whose elevation departs from the median of their 5x5
              neighbourhood by more than 0.3 m, on ground below 5 degrees.
              This is the speckle directly, and unlike RMSE it does not average
              away against a mostly-correct surface.
  flat_over   share of cells below 5 degrees more than 0.5 m from the reference
  density     ground returns kept per square metre. A filter that removes the
              speckle by removing the ground has not helped.

ARMS
----
ELM threshold against the statistical outlier multiplier, plus each filter on
its own, against both off as the baseline.
"""
import json, os, sys, time
from dataclasses import replace
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import ROOT, SRC, DEM_ROOT, LAS_ROOT, TESTS
sys.path.insert(0, SRC)
import numpy as np
from osgeo import gdal
from scipy import ndimage as ndi
from replicalm import classify, grid as G, kriging as K, interpolate
from replicalm.config import PRESETS
gdal.UseExceptions()

# the coverage-selected clip: 95.2% covered, 0.02% of cells above 20 degrees
LAS = TESTS + r"\\pilot_clips\l2s443_clip.las"
REF = TESTS + r"\\pilot_clips\l2s443_reference.tif"
OUT = TESTS + r"\\flatsweep"
os.makedirs(OUT, exist_ok=True)

ARMS = [
    ("off",            dict(remove_low_noise=False, remove_outliers=False)),
    ("elm1.0_out2.5",  dict(elm_threshold_m=1.0, outlier_multiplier=2.5)),
    ("elm0.5_out2.5",  dict(elm_threshold_m=0.5, outlier_multiplier=2.5)),
    ("elm2.0_out2.5",  dict(elm_threshold_m=2.0, outlier_multiplier=2.5)),
    ("elm1.0_out2.0",  dict(elm_threshold_m=1.0, outlier_multiplier=2.0)),
    ("elm1.0_out3.5",  dict(elm_threshold_m=1.0, outlier_multiplier=3.5)),
    ("elm_only",       dict(elm_threshold_m=1.0, remove_outliers=False)),
    ("outlier_only",   dict(remove_low_noise=False, outlier_multiplier=2.5)),
]

g = G.grid_from_raster(REF)
wkt = interpolate.source_srs(LAS)
area = (g.bounds[2] - g.bounds[0]) * (g.bounds[3] - g.bounds[1])
rds = gdal.Open(REF); rb = rds.GetRasterBand(1)
ref = rb.ReadAsArray().astype("f8"); nod = rb.GetNoDataValue()
rm = np.isfinite(ref) & (ref != 0.0)
if nod is not None:
    rm &= ref != nod
gyy, gxx = np.gradient(np.where(rm, ref, np.nan), g.dy, g.cell)
slope = np.degrees(np.arctan(np.hypot(gxx, gyy)))
flat = rm & (slope < 5)
print("window %d x %d, %.1f%% valid, %.2f%% of cells below 5 deg"
      % (g.width, g.height, 100 * rm.mean(), 100 * flat.sum() / rm.sum()))


def pits(z, m, thresh=0.3):
    """Cells departing from their local median: the speckle, counted."""
    zz = np.where(m, z, np.nan)
    med = ndi.median_filter(np.nan_to_num(zz, nan=0.0), size=5)
    d = np.abs(zz - med)
    k = m & np.isfinite(d)
    return float((d[k] > thresh).mean())


ref_pits = pits(ref, flat)
print("reference speckle on flat ground: %.3f%%\n" % (100 * ref_pits))
print("%-15s %9s %9s %9s %10s %8s" %
      ("arm", "gnd/m2", "radius", "pits", "flat>0.5", "rmse"))
print("-" * 66)
rows = []
for name, params in ARMS:
    t0 = time.time()
    base = dict(remove_low_noise=True, remove_outliers=True)
    base.update(params)
    cfg = replace(PRESETS["ncalm"], **base)
    gp = os.path.join(OUT, "%s.las" % name)
    try:
        if not os.path.exists(gp):
            classify.classify_tile(LAS, gp, cfg, verbose=False)
        arr, _ = classify.read_points(gp)
        gnd = arr[arr["Classification"] == 2]
        x, y, z = (gnd["X"].astype("f8"), gnd["Y"].astype("f8"),
                   gnd["Z"].astype("f8"))
        dens = len(z) / area
        radius = K.radius_for_density(dens, cfg.max_points,
                                      cfg.search_radius_factor,
                                      cfg.search_radius_floor_m,
                                      cfg.search_radius_ceiling_m)
        v = K.fit_variogram(x, y, z, model="spherical", max_lag=10.0)
        dem, info = K.krige_grid(x, y, z, g, radius=radius,
                                 max_points=cfg.max_points, variogram=v,
                                 verbose=False)
        K.write_geotiff(dem, g, wkt, os.path.join(OUT, "%s.tif" % name))
        m = rm & (dem > -9000)
        fk = flat & (dem > -9000)
        e = np.abs(dem - ref)
        row = {"arm": name, "params": base, "density": dens, "radius": radius,
               "pits": pits(dem, fk), "flat_over": float((e[fk] > 0.5).mean()),
               "rmse": float(np.sqrt((e[m] ** 2).mean())),
               "fallback": float(info["fallback_fraction"]),
               "seconds": round(time.time() - t0, 1)}
        rows.append(row)
        print("%-15s %9.2f %9.2f %8.3f%% %9.3f%% %8.4f"
              % (name, dens, radius, 100 * row["pits"],
                 100 * row["flat_over"], row["rmse"]))
    except Exception as ex:
        print("%-15s FAILED %s" % (name, str(ex)[:46]))

with open(os.path.join(OUT, "flatsweep.json"), "w", encoding="utf-8") as fh:
    json.dump({"reference_pits": ref_pits, "arms": rows}, fh, indent=1,
              default=float)
print("\nreference speckle %.3f%%; an arm at or below it adds none of its own"
      % (100 * ref_pits))
