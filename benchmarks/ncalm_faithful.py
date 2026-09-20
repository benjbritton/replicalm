r"""Run the NCALM method as documented, and see whether the deviations are the problem.

WHY
---
The archive G1 renders these trench walls and mound edges without the errors our
output shows, from the same point cloud, on the same 0.5 m grid. Whatever the
original workflow did, it works here and ours does not, so the burden is on the
deviations rather than on the terrain.

Four of those deviations are mine, introduced during calibration on 400 m clips
that turned out to be almost entirely flat:

  two ground passes -> one     SMRF has no refinement mode, so the second pass
                               was dropped rather than approximated. The source
                               runs Classify Ground twice, at 9 then 12 degrees.
  noise filters off             remove_low_noise and remove_outliers were
                               disabled because they changed nothing measurable
                               on the flat clips.
  20 m search radius -> 5 m     the source states 20 m. I replaced it with a
                               value derived from a variogram that, as the refit
                               showed, has no range to derive anything from.
  1 m SMRF working grid         kept at the config default while rendering at
                               0.5 m, which the sweep has just shown produces
                               terracing on steep faces.

THE ARMS
--------
  ncalm_asdoc      PRESETS["ncalm"] untouched: two passes at 9 and 12 degrees,
                   3 m threshold, 70 m window, 1 m cell, noise filters on,
                   20 m kriging radius. The method as written.
  ncalm_r20        the same, but with the search radius I have been using, to
                   separate the classification deviations from the radius one.
  ncalm_2pass_fine two passes at the fine working grid the sweep prefers, noise
                   filters on, 20 m radius: the source's structure with the one
                   parameter the source could not have specified, since
                   TerraScan has no SMRF cell.
  current          what the tile was last rendered with, for reference.

Scored on steep ground and for terracing, on the window where the artifact is.
"""
import json, os, sys, time
from dataclasses import replace
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import ROOT, SRC, DEM_ROOT, LAS_ROOT, TESTS
sys.path.insert(0, SRC)
import numpy as np
from osgeo import gdal
from replicalm import classify, grid as G, kriging as K, interpolate
from replicalm.config import PRESETS, GroundPass, NCALM_TERRASCAN
gdal.UseExceptions()

LAS = (LAS_ROOT + r"\\Yuc_South\South_Glas\South_Glas"
       r"\AMIGACarb_Yuc_South_GLAS_Apr2013_l0s395.las")
REF = (DEM_ROOT + r"\\Yuc_South\South_Glas\South_Glas"
       r"\South_GLAS_l0s395_DEM_0p5m_v1.tif")
OUT = ROOT + r"\\render_smrf\ncalm"
WIN = ROOT + r"\\render_smrf\trench\window.las"
os.makedirs(OUT, exist_ok=True)
X0, Y0, N = 1500, 5250, 500

full = G.grid_from_raster(REF)
sub = G.Grid(full.origin_x + X0*full.cell, full.origin_y - Y0*full.dy,
             full.cell, N, N, cell_y=full.cell_y)
wkt = interpolate.source_srs(LAS)
rds = gdal.Open(REF); rb = rds.GetRasterBand(1)
ref = rb.ReadAsArray(X0, Y0, N, N).astype("f8"); nod = rb.GetNoDataValue()
rm = np.isfinite(ref) & (ref != 0.0)
if nod is not None: rm &= ref != nod
gyy, gxx = np.gradient(np.where(rm, ref, np.nan), full.dy, full.cell)
slope = np.degrees(np.arctan(np.hypot(gxx, gyy)))


def terrace(z, m, lo=20.0):
    zz = np.where(m, z, np.nan)
    gy, gx = np.gradient(zz, full.dy, full.cell)
    s = np.degrees(np.arctan(np.hypot(gx, gy)))
    k = m & (slope >= lo) & np.isfinite(s)
    return float((s[k] < 0.25 * slope[k]).mean()) if k.sum() > 200 else float("nan")


def fine_passes():
    """The source's two passes, on a working grid it could not have specified."""
    import math
    return [GroundPass(algorithm="smrf", slope=round(math.tan(math.radians(9.0)), 4),
                       threshold_m=0.5, cell_m=0.5, window_m=70.0,
                       stands_in_for="TerraScan pass 1"),
            GroundPass(algorithm="smrf", slope=round(math.tan(math.radians(12.0)), 4),
                       threshold_m=0.5, cell_m=0.5, window_m=70.0, reuse_ground=True,
                       stands_in_for="TerraScan pass 2")]


ARMS = [
    ("ncalm_asdoc", PRESETS["ncalm"], 20.0),
    ("ncalm_r5", PRESETS["ncalm"], 5.0),
    ("ncalm_2pass_fine", replace(PRESETS["ncalm"], passes=fine_passes()), 20.0),
    ("ncalm_2pass_fine_r5", replace(PRESETS["ncalm"], passes=fine_passes()), 5.0),
    ("current_1pass", replace(PRESETS["ncalm"],
        passes=[GroundPass(algorithm="smrf", slope=0.05, threshold_m=0.5,
                           cell_m=0.5)],
        remove_low_noise=False, remove_outliers=False), 5.0),
]

print("source radius %.0f m, source passes: %s"
      % (NCALM_TERRASCAN["kriging_search_radius_m"],
         "9 deg then 12 deg, 3 m, 70 m window"))
print("\n%-20s %6s %8s %9s %9s %9s %9s" %
      ("arm", "radius", "ground", "rmse", ">=20 deg", ">=30 deg", "terrace"))
print("-" * 80)
rows = []
for name, cfg, radius in ARMS:
    t0 = time.time()
    gp = os.path.join(OUT, "%s.las" % name)
    try:
        if not os.path.exists(gp):
            classify.classify_tile(WIN, gp, cfg, verbose=False)
        arr, _ = classify.read_points(gp)
        g = arr[arr["Classification"] == 2]
        if len(g) < 500:
            print("%-20s %6.1f only %d ground points" % (name, radius, len(g)))
            continue
        x, y, z = g["X"].astype("f8"), g["Y"].astype("f8"), g["Z"].astype("f8")
        v = K.fit_variogram(x, y, z, model="spherical", max_lag=10.0)
        dem, _ = K.krige_grid(x, y, z, sub, radius=radius, max_points=16,
                              variogram=v, verbose=False)
        K.write_geotiff(dem, sub, wkt, os.path.join(OUT, "%s.tif" % name))
        m = rm & (dem > -9000) & np.isfinite(slope)
        e = np.abs(dem - ref)
        b20 = m & (slope >= 20); b30 = m & (slope >= 30)
        row = {"arm": name, "radius": radius, "ground": int(len(z)),
               "rmse": float(np.sqrt((e[m]**2).mean())),
               "over20": float((e[b20] > 0.5).mean()),
               "over30": float((e[b30] > 0.5).mean()),
               "terrace": terrace(dem, m)}
        rows.append(row)
        print("%-20s %6.1f %8d %9.4f %8.2f%% %8.2f%% %8.2f%%"
              % (name, radius, len(z), row["rmse"], 100*row["over20"],
                 100*row["over30"], 100*row["terrace"]))
    except Exception as ex:
        print("%-20s %6.1f FAILED %s" % (name, radius, str(ex)[:44]))

with open(os.path.join(OUT, "ncalm.json"), "w", encoding="utf-8") as fh:
    json.dump({"reference_terrace": terrace(ref, rm), "arms": rows}, fh, indent=1)
print("\nreference terracing on this window: %.2f%%" % (100 * terrace(ref, rm)))
