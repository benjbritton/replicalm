r"""Why do steep flanks lose their ground returns, and which filter keeps them?

The window at 2100,0 has 74 returns within 1.5 m of the cells that fail, and
zero of them classified as ground. The surface is measured; the classifier is
discarding it. Ground fraction falls from 45% on flat ground to 28% above 30
degrees, and to nil on the flanks themselves.

CSF rigidness is the first suspect. config.py documents the scale as 1 steep,
2 relief, 3 flat, and this tile was processed at 2 because that setting won the
400 m calibration clip -- a clip that is almost entirely flat. A cloth too stiff
to drape over a 30 degree flank leaves those returns further below it than the
0.5 m threshold allows, so they are never ground.

Each arm is scored on what actually matters here: how much of the STEEP ground
gets returns at all, not the tile-wide RMSE that flat terrain dominates.
"""
import json, os, sys, time
from dataclasses import replace
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import ROOT, SRC, DEM_ROOT, LAS_ROOT, TESTS
sys.path.insert(0, SRC)
import numpy as np
from scipy.spatial import cKDTree
from osgeo import gdal
from replicalm import classify, clip, grid as G, kriging as K, interpolate
from replicalm.config import PRESETS, GroundPass
gdal.UseExceptions()

LAS = (LAS_ROOT + r"\\Yuc_South\South_Glas\South_Glas"
       r"\AMIGACarb_Yuc_South_GLAS_Apr2013_l0s395.las")
REF = (DEM_ROOT + r"\\Yuc_South\South_Glas\South_Glas"
       r"\South_GLAS_l0s395_DEM_0p5m_v1.tif")
OUT = ROOT + r"\\render\flank"
os.makedirs(OUT, exist_ok=True)
X0, Y0, NW, NH = 2100, 0, 700, 700
RADIUS = 5.0        # a middle radius: the fixed-cell ladder found 3 to 14.4 m
                    # indistinguishable, so this is not the variable under test

ARMS = [
    ("csf_r2_t050", dict(algorithm="csf", csf_rigidness=2, csf_threshold_m=0.5)),
    ("csf_r1_t050", dict(algorithm="csf", csf_rigidness=1, csf_threshold_m=0.5)),
    ("csf_r1_t030", dict(algorithm="csf", csf_rigidness=1, csf_threshold_m=0.3)),
    ("smrf_s05_t050", dict(algorithm="smrf", slope=0.05, threshold_m=0.5)),
    ("smrf_s20_t050", dict(algorithm="smrf", slope=0.20, threshold_m=0.5)),
    ("smrf_s40_t050", dict(algorithm="smrf", slope=0.40, threshold_m=0.5)),
]

full = G.grid_from_raster(REF)
sub = G.Grid(full.origin_x + X0 * full.cell, full.origin_y - Y0 * full.dy,
             full.cell, NW, NH, cell_y=full.cell_y)
mnx, mny, mxx, mxy = sub.bounds
wkt = interpolate.source_srs(LAS)
win = os.path.join(OUT, "window.las")
if not os.path.exists(win):
    clip.clip_las(LAS, win, {"minx": mnx - 25, "miny": mny - 25,
                             "maxx": mxx + 25, "maxy": mxy + 25}, verbose=True)

rds = gdal.Open(REF); rb = rds.GetRasterBand(1)
ref = rb.ReadAsArray(X0, Y0, NW, NH).astype("f8"); nod = rb.GetNoDataValue()
rm = np.isfinite(ref) & (ref != 0.0)
if nod is not None: rm &= ref != nod
gyy, gxx = np.gradient(np.where(rm, ref, np.nan), full.dy, full.cell)
slope = np.degrees(np.arctan(np.hypot(gxx, gyy)))
steep = rm & (slope >= 20)
gx_, gy_ = sub.cell_centres(); flat = np.c_[gx_.ravel(), gy_.ravel()]

print("%-14s %8s %9s %10s %10s %9s %9s" %
      ("arm", "ground", "gnd/m2", "steep>=3pt", "filled", "rmse", "steep>0.5"))
print("-" * 76)
rows = []
for name, params in ARMS:
    t0 = time.time()
    cfg = replace(PRESETS["ncalm"], passes=[GroundPass(**params)],
                  remove_low_noise=False, remove_outliers=False)
    gp = os.path.join(OUT, "%s.las" % name)
    try:
        if not os.path.exists(gp):
            classify.classify_tile(win, gp, cfg, verbose=False)
        arr, _ = classify.read_points(gp)
        g = arr[arr["Classification"] == 2]
        x, y, z = g["X"].astype("f8"), g["Y"].astype("f8"), g["Z"].astype("f8")
        area = (mxx - mnx) * (mxy - mny)
        tg = cKDTree(np.c_[x, y])
        ng = np.array(tg.query_ball_point(flat, 1.5,
                                          return_length=True)).reshape(NH, NW)
        supported = float((ng[steep] >= 3).mean())

        v = K.fit_variogram(x, y, z, model="spherical", max_lag=10.0)
        dem, info = K.krige_grid(x, y, z, sub, radius=RADIUS, max_points=16,
                                 variogram=v, verbose=False)
        K.write_geotiff(dem, sub, wkt, os.path.join(OUT, "%s.tif" % name))
        m = rm & (dem > -9000) & np.isfinite(slope)
        e = np.abs(dem - ref)
        st = m & (slope >= 20)
        rows.append({"arm": name, "ground": int(len(z)),
                     "density": len(z) / area, "steep_supported": supported,
                     "filled": float(m.sum() / rm.sum()),
                     "rmse": float(np.sqrt((e[m] ** 2).mean())),
                     "steep_over_half": float((e[st] > 0.5).mean())})
        print("%-14s %8d %9.2f %9.1f%% %9.1f%% %9.3f %8.2f%%"
              % (name, len(z), len(z) / area, 100 * supported,
                 100 * m.sum() / rm.sum(), rows[-1]["rmse"],
                 100 * rows[-1]["steep_over_half"]))
    except Exception as e:
        print("%-14s FAILED %s" % (name, str(e)[:50]))

with open(os.path.join(OUT, "flank.json"), "w", encoding="utf-8") as fh:
    json.dump(rows, fh, indent=1)
print("\nsteep>=3pt is the share of cells above 20 degrees with at least three")
print("ground returns within 1.5 m -- the condition the flanks were failing.")
