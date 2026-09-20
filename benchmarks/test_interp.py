r"""Rasterise the classified tile and confirm the result is a real GeoTIFF."""
import os, sys
sys.path.insert(0, r"C:\Replicalm\src")
from replicalm import interpolate

LAS = r"C:\Replicalm\tests\out\l2s505_ground.las"
OUT = r"C:\Replicalm\tests\out\l2s505_dem.tif"

print("source CRS from LAS header:")
wkt = interpolate.source_srs(LAS)
print("  ", (wkt or "(none)")[:78])

print("\nrasterising, kriging fill, 1 m cells:")
p = interpolate.rasterize_cloudcompare(LAS, OUT, cell_m=1.0, fill="KRIGING")

print("\nresult:")
for k, v in interpolate.raster_stats(p).items():
    print("  %-16s %s" % (k, v))

print("\nrejection check: an unknown fill must fail loudly, not silently:")
try:
    interpolate.rasterize_cloudcompare(
        LAS, os.path.join(os.path.dirname(OUT), "bad.tif"), fill="DELAUNAY")
    print("   FAILED: it accepted a bad fill strategy")
except interpolate.InterpolateError as e:
    print("   raised as it should:", str(e)[:70])
