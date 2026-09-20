r"""Does an explicit -GLOBAL_SHIFT change where the raster lands?

Three runs on identical input: CloudCompare's automatic shift, an explicit
shift chosen by us, and no shift at all. If all three place the grid at the
same world coordinates, the shift is being undone on export and the output
position is safe whatever mode is used. If they differ, the mode matters and we
must pin it.
"""
import glob, os, subprocess, sys
sys.path.insert(0, r"C:\Replicalm\src")
from osgeo import gdal
gdal.UseExceptions()

CC = r"C:\Program Files\CloudCompare\CloudCompare.exe"
LAS = r"C:\Replicalm\tests\out\l2s505_ground.las"
D = os.path.dirname(LAS)

MODES = [("auto (default)", []),
         ("AUTO explicit", ["-GLOBAL_SHIFT", "AUTO"]),
         ("explicit value", ["-GLOBAL_SHIFT", "-715000", "-2020000", "0"]),
         ("none", ["-GLOBAL_SHIFT", "0", "0", "0"])]

print("%-16s %-9s %14s %15s %s" % ("mode", "accepted", "origin x", "origin y", "cells"))
print("-" * 74)
for name, extra in MODES:
    before = set(glob.glob(os.path.join(D, "*.tif")))
    cmd = [CC, "-SILENT", "-O"] + extra + [LAS, "-RASTERIZE",
           "-GRID_STEP", "1.0", "-PROJ", "AVG",
           "-EMPTY_FILL", "KRIGING", "-OUTPUT_RASTER_Z"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    made = sorted(set(glob.glob(os.path.join(D, "*.tif"))) - before)
    msg = (r.stdout or "") + (r.stderr or "")
    bad = "nknown" in msg or "nvalid" in msg or "rror" in msg
    if not made:
        print("%-16s %-9s %s" % (name, "no output", msg.strip()[:40]))
        continue
    d = gdal.Open(made[-1])
    gt = d.GetGeoTransform()
    print("%-16s %-9s %14.2f %15.2f %d x %d"
          % (name, "warned" if bad else "ok", gt[0], gt[3],
             d.RasterXSize, d.RasterYSize))
    d = None
    os.remove(made[-1])
