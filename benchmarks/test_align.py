r"""Confirm the alignment check passes on good output and would catch a shift."""
import sys
sys.path.insert(0, r"C:\Replicalm\src")
from replicalm import interpolate

LAS = r"C:\Replicalm\tests\out\l2s505_ground.las"
OUT = r"C:\Replicalm\tests\out\l2s505_dem.tif"

print("rasterise with the check wired in:")
interpolate.rasterize_cloudcompare(LAS, OUT, cell_m=1.0, fill="KRIGING")

print("\nalignment:")
for k, v in interpolate.verify_alignment(OUT, LAS).items():
    print("  %-14s %s" % (k, round(v, 3) if isinstance(v, float) else v))

print("\nwould it catch a shift? tightening tolerance to 0.001 cells:")
try:
    interpolate.verify_alignment(OUT, LAS, tolerance_cells=0.001)
    print("   no -- check is too loose to detect anything")
except interpolate.InterpolateError as e:
    print("   yes:", str(e)[:88])
