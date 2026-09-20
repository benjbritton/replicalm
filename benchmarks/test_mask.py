r"""The nodata convention: black outside, never black inside."""
import sys
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import ROOT, SRC, DEM_ROOT, LAS_ROOT, TESTS
sys.path.insert(0, SRC)
import numpy as np
from replicalm import grid as G, kriging as K, classify, interpolate, mask

LAS = TESTS + r"\\out\l2s505_ground.las"
OUT = TESTS + r"\\out\l2s505_dem_masked.tif"

g = G.grid_for_las(LAS, cell=1.0)
arr, _ = classify.read_points(LAS)
gnd = arr[arr["Classification"] == 2]
dem, v = K.krige_grid(gnd["X"], gnd["Y"], gnd["Z"], g, radius=20.0,
                      verbose=False)
valid = mask.coverage_mask(dem)
print("kriged: %d of %d cells carry a value (%.1f%%)"
      % (valid.sum(), valid.size, 100 * valid.mean()))

wkt = interpolate.source_srs(LAS)
p, rep = mask.write_masked_geotiff(dem, g, wkt, OUT, erode_cells=1)
print("\nwrote", p)
for k in ("cells_valid", "valid_fraction", "zeroed_outside", "lifted_inside",
          "lift_step", "bands"):
    print("  %-16s %s" % (k, rep[k]))

print("\nreading it back and checking the convention holds on disk:")
from osgeo import gdal
gdal.UseExceptions()
d = gdal.Open(p)
z = d.GetRasterBand(1).ReadAsArray()
m = d.GetRasterBand(2).ReadAsArray().astype(bool)
print("  band 1: %s" % d.GetRasterBand(1).GetDescription())
print("  band 2: %s" % d.GetRasterBand(2).GetDescription())
print("  zeros inside the mask : %d" % int((m & (z == 0)).sum()))
print("  non-zeros outside     : %d" % int((~m & (z != 0)).sum()))
print("  elevation range where valid: %.2f - %.2f m" % (z[m].min(), z[m].max()))
mask.verify(z, m)
print("  verify() passed")

print("\nan artificial zero inside the mask must be caught:")
z2 = z.copy(); z2[m][0] = 0
zz = z.copy(); idx = np.argwhere(m)[0]; zz[tuple(idx)] = 0
try:
    mask.verify(zz, m)
    print("   FAILED: it was not caught")
except mask.MaskError as e:
    print("   caught:", str(e)[:78])
