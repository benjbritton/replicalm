r"""Terracing on trench walls: is it SMRF's working grid?

WHAT IS LEFT
------------
Switching the classifier from CSF to SMRF removed the bead artifact and cut
steep-ground error from 7.88% to 1.20% of cells at 20-30 degrees. What remains
is concentric banding on the steepest faces -- the walls of looter trenches and
the sharpest mound edges -- which reads as terracing rather than as smoothing.

THE SUSPECT
-----------
filters.smrf works on its own grid, `cell`, which the config leaves at 1.0 m
while this tile is rendered at 0.5 m. SMRF opens a minimum-elevation surface on
that grid and admits a return as ground if it sits within `threshold` of it.

On flat ground a 1 m cell spans a few centimetres of real elevation and the
threshold means what it says. On a trench wall it spans a metre or more, so
whether a return is admitted depends on where in the cell it fell, not on
whether it is ground. Admitted and rejected returns then alternate down the
face in bands the width of the working grid -- which is what terracing is.

If that is the mechanism, halving the cell should halve the band width and
quartering it should quarter it, with the steep-band error falling as it does.
If the banding is indifferent to `cell`, it is coming from the interpolation or
from the returns themselves and this is the wrong knob.

`threshold` is swept alongside because the two interact: a tighter threshold on
a coarse grid rejects most of the face, while on a fine grid it can afford to be
tight. Scored only on the steepest ground, since that is where the artifact is
and tile-wide figures are dominated by flat terrain that all arms render
identically.
"""
import json, os, sys, time
from dataclasses import replace
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import ROOT, SRC, DEM_ROOT, LAS_ROOT, TESTS
sys.path.insert(0, SRC)
import numpy as np
from osgeo import gdal
from replicalm import classify, clip, grid as G, kriging as K, interpolate
from replicalm.config import PRESETS, GroundPass
gdal.UseExceptions()

LAS = (LAS_ROOT + r"\\Yuc_South\South_Glas\South_Glas"
       r"\AMIGACarb_Yuc_South_GLAS_Apr2013_l0s395.las")
REF = (DEM_ROOT + r"\\Yuc_South\South_Glas\South_Glas"
       r"\South_GLAS_l0s395_DEM_0p5m_v1.tif")
OUT = ROOT + r"\\render_smrf\trench"
os.makedirs(OUT, exist_ok=True)
X0, Y0, N = 1500, 5250, 500        # the worst window, 3.48% of cells over 0.5 m
RADIUS = 5.0

ARMS = [
    ("c1.00_t0.50", dict(cell_m=1.00, threshold_m=0.50)),
    ("c0.50_t0.50", dict(cell_m=0.50, threshold_m=0.50)),
    ("c0.50_t0.25", dict(cell_m=0.50, threshold_m=0.25)),
    ("c0.25_t0.25", dict(cell_m=0.25, threshold_m=0.25)),
    ("c0.25_t0.15", dict(cell_m=0.25, threshold_m=0.15)),
    ("c0.50_t0.15", dict(cell_m=0.50, threshold_m=0.15)),
]

full = G.grid_from_raster(REF)
sub = G.Grid(full.origin_x + X0*full.cell, full.origin_y - Y0*full.dy,
             full.cell, N, N, cell_y=full.cell_y)
mnx, mny, mxx, mxy = sub.bounds
wkt = interpolate.source_srs(LAS)
win = os.path.join(OUT, "window.las")
if not os.path.exists(win):
    clip.clip_las(LAS, win, {"minx": mnx-25, "miny": mny-25,
                             "maxx": mxx+25, "maxy": mxy+25}, verbose=True)

rds = gdal.Open(REF); rb = rds.GetRasterBand(1)
ref = rb.ReadAsArray(X0, Y0, N, N).astype("f8"); nod = rb.GetNoDataValue()
rm = np.isfinite(ref) & (ref != 0.0)
if nod is not None: rm &= ref != nod
gyy, gxx = np.gradient(np.where(rm, ref, np.nan), full.dy, full.cell)
slope = np.degrees(np.arctan(np.hypot(gxx, gyy)))


def banding(z, m, lo=20.0):
    """A terracing measure: how much of the steep surface is locally flat.

    Terracing replaces a continuous face with steps, so it puts near-zero
    gradients where the reference has a steady one. The share of steep cells
    whose own gradient has collapsed is therefore a direct measure of it, and
    it does not care about the sign or size of the elevation error.
    """
    zz = np.where(m, z, np.nan)
    gy, gx = np.gradient(zz, full.dy, full.cell)
    s = np.degrees(np.arctan(np.hypot(gx, gy)))
    k = m & (slope >= lo) & np.isfinite(s)
    if k.sum() < 200:
        return float("nan")
    return float((s[k] < 0.25 * slope[k]).mean())


print("%-13s %8s %9s %9s %9s %9s %9s" %
      ("arm", "ground", "rmse", ">=20 deg", ">=30 deg", "terrace", "sec"))
print("-" * 74)
rows = []
base = banding(ref, rm)
for name, params in ARMS:
    t0 = time.time()
    cfg = replace(PRESETS["ncalm"],
                  passes=[GroundPass(algorithm="smrf", slope=0.05, **params)],
                  remove_low_noise=False, remove_outliers=False)
    gp = os.path.join(OUT, "%s.las" % name)
    try:
        if not os.path.exists(gp):
            classify.classify_tile(win, gp, cfg, verbose=False)
        arr, _ = classify.read_points(gp)
        g = arr[arr["Classification"] == 2]
        x, y, z = g["X"].astype("f8"), g["Y"].astype("f8"), g["Z"].astype("f8")
        v = K.fit_variogram(x, y, z, model="spherical", max_lag=10.0)
        dem, _ = K.krige_grid(x, y, z, sub, radius=RADIUS, max_points=16,
                              variogram=v, verbose=False)
        K.write_geotiff(dem, sub, wkt, os.path.join(OUT, "%s.tif" % name))
        m = rm & (dem > -9000) & np.isfinite(slope)
        e = np.abs(dem - ref)
        b20 = m & (slope >= 20); b30 = m & (slope >= 30)
        row = {"arm": name, "params": params, "ground": int(len(z)),
               "rmse": float(np.sqrt((e[m]**2).mean())),
               "over20": float((e[b20] > 0.5).mean()),
               "over30": float((e[b30] > 0.5).mean()),
               "terrace": banding(dem, m)}
        rows.append(row)
        print("%-13s %8d %9.4f %8.2f%% %8.2f%% %8.2f%% %8.0f"
              % (name, len(z), row["rmse"], 100*row["over20"],
                 100*row["over30"], 100*row["terrace"], time.time()-t0))
    except Exception as ex:
        print("%-13s FAILED %s" % (name, str(ex)[:50]))

with open(os.path.join(OUT, "trench.json"), "w", encoding="utf-8") as fh:
    json.dump({"reference_terrace": base, "arms": rows}, fh, indent=1)
print("\nterrace = share of cells above 20 deg whose own gradient has collapsed")
print("to under a quarter of the reference's. The reference scores %.2f%%."
      % (100 * base))
