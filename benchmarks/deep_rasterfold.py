r"""Fidelity with the grid inside the measurement, so cell size can be judged.

WHY THE PREVIOUS TEST COULD NOT ANSWER THIS
-------------------------------------------
`evaluate.buffered_loo` predicts point to point and never touches the raster,
so it returned identical figures at 1.00, 0.50, 0.33 and 0.25 m -- correct for
comparing classifications, useless for comparing resolutions. And the
returns-below-surface measure improves as cells shrink for a reason that has
nothing to do with quality: a coarse cell averages over more ground, so returns
at the low end of a cell fall beneath its single value.

Neither could see the grid. This one puts it in the path.

THE TEST
--------
Withhold whole spatial blocks, build the DEM at a given cell size from what is
left, then sample THAT RASTER at the withheld returns. The residual is the
difference between a rasterised surface and real measurements that contributed
nothing to it.

BLOCK SIZE IS THE THING TO GET RIGHT
------------------------------------
An earlier attempt used 25 m blocks against a 20 m search radius, which left
block centres beyond any training point: it measured extrapolation across gaps
and reported half-metre errors that said nothing about the surface. Blocks here
are 5 m, so a withheld return is at most about 2.5 m from training data and the
surface is interpolating, as it does in use, while still owing nothing to the
point being predicted.

SCORED ON THE SAME TARGETS IN BOTH ARMS
---------------------------------------
Only returns that Clear trusts are used as targets, in the raw arm as well.
Otherwise the raw surface would be graded on predicting the clutter it contains,
which it does well by construction.
"""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import ROOT, SRC, DEM_ROOT
sys.path.insert(0, SRC)
import numpy as np
from osgeo import gdal
from replicalm import classify, cleanup, grid as G, kriging as K
gdal.UseExceptions()

X0, Y0, N = 1500, 5250, 500
CELLS = [1.00, 0.50, 0.33, 0.25]
BLOCK_M = 5.0
FOLDS = 3
SEED = 20260920
REF = os.path.join(DEM_ROOT, r"Yuc_South\South_Glas\South_Glas"
                             r"\South_GLAS_l0s395_DEM_0p5m_v1.tif")
OUR_GND = os.path.join(ROOT, r"render_ncalm\work\l0s395_00000_04000_gnd.las")
OUT = os.path.join(ROOT, "render_ncalm", "deep")
os.makedirs(OUT, exist_ok=True)

rds = gdal.Open(REF)
gt = rds.GetGeoTransform(); cy, cx = abs(gt[5]), abs(gt[1])
minx, maxy = gt[0] + X0 * cx, gt[3] - Y0 * cy
maxx, miny = minx + N * cx, maxy - N * cy

arr, _ = classify.read_points(OUR_GND)
gnd = arr[arr["Classification"] == 2]
k = ((gnd["X"] >= minx - 25) & (gnd["X"] <= maxx + 25) &
     (gnd["Y"] >= miny - 25) & (gnd["Y"] <= maxy + 25))
gnd = gnd[k]
x = gnd["X"].astype("f8"); y = gnd["Y"].astype("f8"); z = gnd["Z"].astype("f8")
above = cleanup.height_above_floor(x, y, z, patch=0.75, percentile=10.0)
trust = above <= 0.20
inside = (x >= minx) & (x <= maxx) & (y >= miny) & (y <= maxy)
print("%d ground returns, %d trusted, %d inside the window"
      % (len(z), trust.sum(), inside.sum()))

rng = np.random.default_rng(SEED)
bx = np.floor((x - minx) / BLOCK_M).astype(int)
by = np.floor((maxy - y) / BLOCK_M).astype(int)
bid = by * (bx.max() + 2) + bx
blocks = np.unique(bid)
assign = rng.integers(0, FOLDS, len(blocks))
fold_of = dict(zip(blocks.tolist(), assign.tolist()))
fold = np.array([fold_of[b] for b in bid])
print("%d blocks of %.0f m across %d folds\n" % (len(blocks), BLOCK_M, FOLDS))

print("%-8s %-7s %9s %9s %9s %9s %8s"
      % ("cell", "arm", "tested", "median", "mae", "p90", "sec"))
print("-" * 66)
rows = []
for cell in CELLS:
    g = G.snap_outward(minx, miny, maxx, maxy, cell)
    for arm, keep_rule in (("raw", np.ones(len(z), bool)), ("clear", trust)):
        t0 = time.time()
        residuals = []
        for f in range(FOLDS):
            train = keep_rule & (fold != f)
            test = trust & inside & (fold == f)
            if train.sum() < 5000 or test.sum() < 200:
                continue
            v = K.fit_variogram(x[train], y[train], z[train],
                                model="spherical", max_lag=10.0)
            dem, _i = K.krige_grid(x[train], y[train], z[train], g,
                                   radius=20.0, max_points=16, variogram=v,
                                   verbose=False)
            col = np.clip(((x[test] - g.origin_x) / g.cell).astype(int),
                          0, g.width - 1)
            row = np.clip(((g.origin_y - y[test]) / g.dy).astype(int),
                          0, g.height - 1)
            pred = dem[row, col]
            ok = pred > -9000
            residuals.append(np.abs(pred[ok] - z[test][ok]))
        if not residuals:
            print("%-8.2f %-7s no usable fold" % (cell, arm)); continue
        r = np.concatenate(residuals)
        rows.append({"cell": cell, "arm": arm, "tested": int(r.size),
                     "median": float(np.median(r)), "mae": float(r.mean()),
                     "p90": float(np.percentile(r, 90)),
                     "rmse": float(np.sqrt((r ** 2).mean()))})
        print("%-8.2f %-7s %9d %9.4f %9.4f %9.4f %8.0f"
              % (cell, arm, r.size, rows[-1]["median"], rows[-1]["mae"],
                 rows[-1]["p90"], time.time() - t0))

with open(os.path.join(OUT, "rasterfold.json"), "w", encoding="utf-8") as fh:
    json.dump({"block_m": BLOCK_M, "folds": FOLDS, "seed": SEED,
               "arms": rows}, fh, indent=1, default=float)
print("\nresidual between the rasterised surface and returns that did not "
      "build it.\nIf this falls as cells shrink, finer grids carry real "
      "information; if it is flat,\nthey carry rendering only.")
