r"""Which variogram model gives a well-conditioned solve on this terrain?"""
import sys
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import ROOT, SRC, DEM_ROOT, LAS_ROOT, TESTS
sys.path.insert(0, SRC)
import numpy as np
from replicalm import grid as G, kriging as K, classify

LAS = TESTS + r"\\out\l2s505_ground.las"
g = G.grid_for_las(LAS, cell=1.0)
arr, _ = classify.read_points(LAS)
x, y, z = arr["X"], arr["Y"], arr["Z"]
print("data z %.2f - %.2f, variance %.3f\n" % (z.min(), z.max(), z.var()))

print("%-14s %-6s %-34s %8s %8s %7s %s"
      % ("model", "maxlag", "params", "zmin", "zmax", "out%", "filled%"))
print("-" * 96)
for model in ("linear", "spherical", "exponential"):
    for maxlag in (20.0, 60.0):
        try:
            v = K.fit_variogram(x, y, z, model=model, max_lag=maxlag, n_lags=14)
        except Exception as e:
            print("%-14s %-6.0f fit failed: %s" % (model, maxlag, e)); continue
        dem, _ = K.krige_grid(x, y, z, g, radius=20.0, max_points=32,
                              variogram=v, verbose=False)
        m = dem > -9000
        out = ((dem[m] < z.min()) | (dem[m] > z.max())).sum()
        p = tuple(round(float(q), 4) if q is not None else None for q in v["params"])
        print("%-14s %-6.0f %-34s %8.2f %8.2f %6.2f%% %7.1f%%"
              % (model, maxlag, str(p), dem[m].min(), dem[m].max(),
                 100 * out / m.sum(), 100 * m.mean()))
