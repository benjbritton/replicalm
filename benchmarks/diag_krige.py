r"""Why does the kriged surface leave the data range, and why are cells empty?"""
import sys
sys.path.insert(0, r"C:\Replicalm\src")
import numpy as np
from scipy.spatial import cKDTree
from replicalm import grid as G, kriging as K, classify

LAS = r"C:\Replicalm\tests\out\l2s505_ground.las"
g = G.grid_for_las(LAS, cell=1.0)
arr, _ = classify.read_points(LAS)
x, y, z = arr["X"], arr["Y"], arr["Z"]
print("data z: %.2f - %.2f, mean %.2f" % (z.min(), z.max(), z.mean()))

tree = cKDTree(np.c_[x, y])
gx, gy = g.cell_centres()
flat = np.c_[gx.ravel(), gy.ravel()]

for k, radius in ((32, 20.0), (32, 5.0), (16, 20.0)):
    d, _ = tree.query(flat, k=k, distance_upper_bound=radius)
    nn = np.isfinite(d).sum(axis=1)
    print("k=%-3d r=%-5.1f  cells with >=3 neighbours: %5.1f%%   "
          "median neighbours %d   nearest-point distance median %.2f m"
          % (k, radius, 100 * (nn >= 3).mean(), int(np.median(nn)),
             float(np.median(d[:, 0][np.isfinite(d[:, 0])]))))

print()
v = K.fit_variogram(x, y, z, model="linear", max_lag=20.0)
nug, slope, _ = v["params"]
print("linear variogram: nugget %.5f, slope %.5f per m" % (nug, slope))
print("  semivariance at 1 m %.4f, at 20 m %.4f"
      % (nug + slope * 1, nug + slope * 20))
print("  data variance %.4f" % z.var())
print("  -> the model explains %.1f%% of the variance at 20 m"
      % (100 * (nug + slope * 20) / z.var()))

dem, _ = K.krige_grid(x, y, z, g, radius=20.0, max_points=32, variogram=v,
                      verbose=False)
m = dem > -9000
out_lo = (dem[m] < z.min()).sum()
out_hi = (dem[m] > z.max()).sum()
print("\nkriged cells outside the data range: %d below, %d above, of %d (%.1f%%)"
      % (out_lo, out_hi, m.sum(), 100 * (out_lo + out_hi) / m.sum()))
worst = dem[m].min()
print("worst undershoot: %.2f m, which is %.2f m below the lowest return"
      % (worst, z.min() - worst))

# is the overshoot associated with few neighbours?
d, _ = tree.query(flat, k=32, distance_upper_bound=20.0)
nn = np.isfinite(d).sum(axis=1).reshape(dem.shape)
bad = m & ((dem < z.min()) | (dem > z.max()))
if bad.any():
    print("cells outside range have median %d neighbours; good cells %d"
          % (int(np.median(nn[bad])), int(np.median(nn[m & ~bad]))))
