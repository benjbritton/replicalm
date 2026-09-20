r"""What cell size does this point cloud actually support?

Judged with no reference DEM. The archive is 0.5 m, and scoring against it
would make 0.5 m correct by construction, which is the question rather than the
answer.

THE THREE MEASURES
------------------
  fidelity      buffered leave-one-out on ground returns. Resolution-
                independent: it predicts at return locations, so a 0.25 m grid
                and a 1.0 m grid are asked the same question.
  penetration   returns lying below the surface, which none should. Computed
                over every return in the window, not only the classified ones.
  sharpness     elevation change over a fixed one-metre baseline, sampled
                bilinearly so cell size does not flatter a finer grid.

ARMS
----
Each cell size is run on the raw ground classification and on the Clear one --
demote a ground return standing more than 0.20 m above the tenth percentile
within 0.75 m -- because resolution and clutter are coupled. A finer grid gives
each cell less evidence, so residual clutter prints more sharply: 0.25 m on a
contaminated surface should look worse than 0.5 m does, and 0.25 m on a clean
one may look better.

WHAT WOULD SETTLE IT
--------------------
Mean ground-return spacing on this window is about 0.49 m, so 0.5 m is roughly
matched to the data, 0.33 m is 1.5x oversampled and 0.25 m is 2x. Oversampling
adds no terrain information but renders slope and sky-view factor with less
quantisation, which is a display gain rather than an accuracy one. If fidelity
and penetration are flat across cell sizes while sharpness keeps rising, that
is what is happening, and the choice becomes one of rendering rather than of
measurement.
"""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import ROOT, SRC, DEM_ROOT
sys.path.insert(0, SRC)
import numpy as np
from osgeo import gdal
from replicalm import classify, cleanup, evaluate, grid as G, kriging as K
gdal.UseExceptions()

X0, Y0, N = 1500, 5250, 500                 # the Pixoyal window
CELLS = [1.00, 0.50, 0.33, 0.25]
REF = os.path.join(DEM_ROOT, r"Yuc_South\South_Glas\South_Glas"
                             r"\South_GLAS_l0s395_DEM_0p5m_v1.tif")
OUR_GND = os.path.join(ROOT, r"render_ncalm\work\l0s395_00000_04000_gnd.las")
ALL_LAS = os.path.join(ROOT, r"render_ncalm\cleanup\trench_all.las")
OUT = os.path.join(ROOT, "render_ncalm", "deep")
os.makedirs(OUT, exist_ok=True)

rds = gdal.Open(REF)
gt = rds.GetGeoTransform(); cy, cx = abs(gt[5]), abs(gt[1])
ox, oy = gt[0] + X0 * cx, gt[3] - Y0 * cy
minx, maxy = ox, oy
maxx, miny = ox + N * cx, oy - N * cy
wkt = rds.GetProjection()

arr, _ = classify.read_points(OUR_GND)
gnd = arr[arr["Classification"] == 2]
k = ((gnd["X"] >= minx - 25) & (gnd["X"] <= maxx + 25) &
     (gnd["Y"] >= miny - 25) & (gnd["Y"] <= maxy + 25))
gnd = gnd[k]
x = gnd["X"].astype("f8"); y = gnd["Y"].astype("f8"); z = gnd["Z"].astype("f8")
area = (maxx - minx) * (maxy - miny)
density = ((x >= minx) & (x <= maxx) & (y >= miny) & (y <= maxy)).sum() / area
print("%d ground returns, %.2f per m2 inside the window, mean spacing %.2f m"
      % (len(z), density, 1 / np.sqrt(density)))

above = cleanup.height_above_floor(x, y, z, patch=0.75, percentile=10.0)
clear_keep = above <= 0.20
print("Clear keeps %.2f%% of them\n" % (100 * clear_keep.mean()))

allr, _ = classify.read_points(ALL_LAS)
ax = allr["X"].astype("f8"); ay = allr["Y"].astype("f8"); az = allr["Z"].astype("f8")
print("%d returns of every class, for the penetration test\n" % len(az))

print("%-8s %-7s %9s %9s %9s %10s %9s %9s"
      % ("cell", "arm", "loo med", "loo mae", "loo p90", "below %",
         "sharp p90", "sharp p99"))
print("-" * 82)
rows = []
for cell in CELLS:
    g = G.snap_outward(minx, miny, maxx, maxy, cell)
    for arm, keep in (("raw", np.ones(len(z), bool)), ("clear", clear_keep)):
        xx, yy, zz = x[keep], y[keep], z[keep]
        v = K.fit_variogram(xx, yy, zz, model="spherical", max_lag=10.0)
        dem, info = K.krige_grid(xx, yy, zz, g, radius=20.0, max_points=16,
                                 variogram=v, verbose=False)
        tif = os.path.join(OUT, "l0s395_%s_%dcm_DEM.tif"
                           % (arm, int(round(cell * 100))))
        K.write_geotiff(dem, g, wkt, tif)
        # fidelity is scored only on returns Clear trusts, in both arms, so the
        # two are judged against the same targets
        loo = evaluate.buffered_loo(xx, yy, zz, v, radius=20.0, max_points=16,
                                    buffer=0.75, n_samples=1500, seed=7,
                                    trust=(above[keep] <= 0.20))
        below = evaluate.returns_below(dem, g, ax, ay, az)
        sharp = evaluate.sharpness(dem, g, baseline_m=1.0)
        rows.append({"cell": cell, "arm": arm, "loo": loo, "below": below,
                     "sharp": sharp, "fallback": info["fallback_fraction"],
                     "dem": tif})
        print("%-8.2f %-7s %9.4f %9.4f %9.4f %9.2f%% %9.4f %9.4f"
              % (cell, arm, loo["median"], loo["mae"], loo["p90"],
                 100 * below["below_fraction"], sharp["p90"], sharp["p99"]))

with open(os.path.join(OUT, "resolution.json"), "w", encoding="utf-8") as fh:
    json.dump({"density": float(density), "arms": rows}, fh, indent=1,
              default=float)
print("\nmean ground spacing %.2f m; fidelity is resolution-independent, so a "
      "flat\nloo column across cell sizes means finer grids add rendering, not "
      "information." % (1 / np.sqrt(density)))
