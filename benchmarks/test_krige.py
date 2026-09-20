r"""Krige the classified tile onto a declared grid, and test reproducibility.

The point of the declared grid is that a rerun gives the same answer. That is
checked here rather than asserted, and it is the property CloudCompare could
not provide.
"""
import json, sys, time
sys.path.insert(0, r"C:\Replicalm\src")
import numpy as np
from replicalm import grid as G, kriging as K, classify, interpolate

LAS = r"C:\Replicalm\tests\out\l2s505_ground.las"
OUT = r"C:\Replicalm\tests\out\l2s505_dem_krige.tif"

g = G.grid_for_las(LAS, cell=1.0)
print("declared grid:", g.describe())
print("  geotransform:", tuple(round(v, 3) for v in g.geotransform))

arr, _ = classify.read_points(LAS)
x, y, z = arr["X"], arr["Y"], arr["Z"]
print("\nground points: %d over %.0f m2 = %.2f pts/m2"
      % (len(z), (g.bounds[2]-g.bounds[0])*(g.bounds[3]-g.bounds[1]),
         len(z) / ((g.bounds[2]-g.bounds[0])*(g.bounds[3]-g.bounds[1]))))

print("\nfitting variogram:")
v = K.fit_variogram(x, y, z, model="spherical", max_lag=20.0)
print("  model %s  params %s  from %d pairs"
      % (v["model"], tuple(round(p, 5) if p else p for p in v["params"]),
         v["pairs_used"]))

t0 = time.time()
dem, _ = K.krige_grid(x, y, z, g, radius=20.0, max_points=32, variogram=v)
print("  %.1f s" % (time.time() - t0))

wkt = interpolate.source_srs(LAS)
K.write_geotiff(dem, g, wkt, OUT)
print("\nwrote", OUT)
for k, val in interpolate.raster_stats(OUT).items():
    print("  %-16s %s" % (k, val))

print("\nreproducibility: kriging the same input twice")
dem2, _ = K.krige_grid(x, y, z, g, radius=20.0, max_points=32, variogram=v,
                       verbose=False)
same = np.array_equal(np.nan_to_num(dem), np.nan_to_num(dem2))
print("  bit-identical:", same)

print("\ngrid alignment: a neighbouring tile would share cell edges")
g2 = G.Grid(g.origin_x + 96.0, g.origin_y, 1.0, 96, 316)
print("  aligns_with:", g.aligns_with(g2))
