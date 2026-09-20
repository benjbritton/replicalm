r"""Does the surface fit real measurements better, independently of the archive?

WHY THIS IS NEEDED
------------------
Every metric used so far -- roughness, RMSE, speck count -- is measured against
NCALM's output. If our surface were genuinely more accurate than theirs, all
three would get worse, because the archive would be the thing in error. Those
metrics measure agreement and have been standing in for quality.

This one does not involve the archive. Ground returns are held out, the surface
is kriged without them, and the residual is measured where those returns
actually are. It asks how well the surface predicts real measurements it never
saw.

HOW THE HOLDOUT IS DONE
-----------------------
In spatial blocks, not at random. Returns cluster along scan lines, so a random
holdout leaves a near neighbour of every held-out point in the training set and
the surface interpolates almost exactly -- an optimistic answer that says
nothing. Whole blocks are withheld instead.

WHAT IS EVALUATED AGAINST
-------------------------
Consensus ground: returns that both our classifier and the delivered cloud's
TerraScan pass call ground. Evaluating against our own ground alone would
reward a surface for fitting the clutter we are trying to remove, and the whole
question is whether removing it helps.

THE ARCHIVE COMPARISON IS NOT FAIR, AND IS REPORTED ANYWAY
----------------------------------------------------------
The archive DEM is sampled at the same held-out points, but NCALM built that
surface using them -- it is in-sample for them and out-of-sample for us. Their
figure is therefore a lower bound on their error and should be read as a
landmark, not as a result. Comparing our two arms with each other is the part
that is sound.
"""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import ROOT, SRC, DEM_ROOT
sys.path.insert(0, SRC)
import numpy as np
from osgeo import gdal
from scipy.spatial import cKDTree
from replicalm import classify, grid as G, kriging as K
gdal.UseExceptions()

PATCH, HEIGHT, PCT = 0.75, 0.20, 10.0
BLOCK_M = 25.0          # holdout block size
HOLDOUT = 0.20          # share of blocks withheld
SEED = 20260920
X0, Y0, N = 1500, 5250, 500
OUT = os.path.join(ROOT, "render_ncalm", "cleanup")
REF = os.path.join(DEM_ROOT, r"Yuc_South\South_Glas\South_Glas"
                             r"\South_GLAS_l0s395_DEM_0p5m_v1.tif")
ALL_LAS = os.path.join(OUT, "trench_all.las")
OUR_GND = os.path.join(ROOT, r"render_ncalm\work\l0s395_00000_04000_gnd.las")

rds = gdal.Open(REF)
band = rds.GetRasterBand(1)
ref = band.ReadAsArray(X0, Y0, N, N).astype("f8")
nod = band.GetNoDataValue()
gt = rds.GetGeoTransform(); cy, cx = abs(gt[5]), abs(gt[1])
ox, oy = gt[0] + X0 * cx, gt[3] - Y0 * cy
fullg = G.grid_from_raster(REF)
g = G.Grid(fullg.origin_x + X0 * fullg.cell, fullg.origin_y - Y0 * fullg.dy,
           fullg.cell, N, N, cell_y=fullg.cell_y)
valid = np.isfinite(ref) & (ref != 0) & (np.abs(ref) < 1e6)
if nod is not None:
    valid &= ref != nod

orig, _ = classify.read_points(ALL_LAS)
ours, _ = classify.read_points(OUR_GND)
og = ours[ours["Classification"] == 2]
mnx, mny, mxx, mxy = g.bounds
k = ((og["X"] >= mnx - 25) & (og["X"] <= mxx + 25) &
     (og["Y"] >= mny - 25) & (og["Y"] <= mxy + 25))
og = og[k]
x = og["X"].astype("f8"); y = og["Y"].astype("f8"); z = og["Z"].astype("f8")

t = cKDTree(np.c_[orig["X"], orig["Y"], orig["Z"]])
dist, idx = t.query(np.c_[x, y, z], k=1)
ncalm_cls = np.where(dist < 0.01, orig["Classification"][idx], 0)
consensus = ncalm_cls == 2
print("%d of our ground returns; %d also ground for NCALM (%.1f%%)"
      % (len(z), consensus.sum(), 100 * consensus.mean()))

tree = cKDTree(np.c_[x, y])
nb = tree.query_ball_point(np.c_[x, y], PATCH)
floor = np.array([np.percentile(z[i], PCT) if len(i) >= 3 else z[j]
                  for j, i in enumerate(nb)])
clutter = (z - floor) > HEIGHT
print("rule would demote %.2f%% of them\n" % (100 * clutter.mean()))

rng = np.random.default_rng(SEED)
bx = np.floor((x - x.min()) / BLOCK_M).astype(int)
by = np.floor((y - y.min()) / BLOCK_M).astype(int)
bid = by * (bx.max() + 1) + bx
blocks = np.unique(bid)
held_blocks = rng.choice(blocks, max(1, int(len(blocks) * HOLDOUT)),
                         replace=False)
held = np.isin(bid, held_blocks)
test = held & consensus & ~clutter        # evaluate on clean consensus ground
print("%d blocks, %d withheld; %d evaluation points\n"
      % (len(blocks), len(held_blocks), test.sum()))

results = {}
for name, use in (("ours, no rule", ~held),
                  ("ours, rule applied", ~held & ~clutter)):
    v = K.fit_variogram(x[use], y[use], z[use], model="spherical", max_lag=10.0)
    dem, info = K.krige_grid(x[use], y[use], z[use], g, radius=20.0,
                             max_points=16, variogram=v, verbose=False)
    col = np.clip(((x[test] - ox) / cx).astype(int), 0, N - 1)
    row = np.clip(((oy - y[test]) / cy).astype(int), 0, N - 1)
    pred = dem[row, col]
    ok = pred > -9000
    res = np.abs(pred[ok] - z[test][ok])
    results[name] = {"n": int(ok.sum()), "mae": float(res.mean()),
                     "rmse": float(np.sqrt((res ** 2).mean())),
                     "p90": float(np.percentile(res, 90)),
                     "median": float(np.median(res))}

col = np.clip(((x[test] - ox) / cx).astype(int), 0, N - 1)
row = np.clip(((oy - y[test]) / cy).astype(int), 0, N - 1)
av = ref[row, col]
ok = valid[row, col]
res = np.abs(av[ok] - z[test][ok])
results["archive (in-sample for NCALM)"] = {
    "n": int(ok.sum()), "mae": float(res.mean()),
    "rmse": float(np.sqrt((res ** 2).mean())),
    "p90": float(np.percentile(res, 90)), "median": float(np.median(res))}

print("residual at held-out ground returns, in metres")
print("%-32s %8s %9s %9s %9s" % ("surface", "points", "median", "mae", "rmse"))
print("-" * 72)
for name, r in results.items():
    print("%-32s %8d %9.4f %9.4f %9.4f"
          % (name, r["n"], r["median"], r["mae"], r["rmse"]))
with open(os.path.join(OUT, "crossval.json"), "w", encoding="utf-8") as fh:
    json.dump({"block_m": BLOCK_M, "holdout": HOLDOUT, "seed": SEED,
               "results": results}, fh, indent=1)
print("\nthe archive row is in-sample for NCALM and out-of-sample for us;")
print("it is a landmark, not a fair comparison.")
