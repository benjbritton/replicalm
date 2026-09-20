r"""Is the 1-2 cm a displaced surface, or only a differently-labelled origin?

If the elevation values are identical cell for cell and only the geotransform
differs, nothing moved: the grid was described slightly differently. If the
values differ, the surface itself was resampled and the error is real.
"""
import glob, os, subprocess
import numpy as np
from osgeo import gdal
gdal.UseExceptions()

CC = r"C:\Program Files\CloudCompare\CloudCompare.exe"
LAS = r"C:\Replicalm\tests\out\l2s505_ground.las"
D = os.path.dirname(LAS)


def run(extra, tag):
    before = set(glob.glob(os.path.join(D, "*.tif")))
    cmd = [CC, "-SILENT", "-O"] + extra + [LAS, "-RASTERIZE", "-GRID_STEP",
           "1.0", "-PROJ", "AVG", "-EMPTY_FILL", "KRIGING", "-OUTPUT_RASTER_Z"]
    subprocess.run(cmd, capture_output=True, text=True)
    made = sorted(set(glob.glob(os.path.join(D, "*.tif"))) - before)
    out = os.path.join(D, "prec_%s.tif" % tag)
    os.replace(made[-1], out)
    return out


a = run([], "default")
b = run(["-GLOBAL_SHIFT", "AUTO"], "autoexp")

da, db = gdal.Open(a), gdal.Open(b)
ga, gb = da.GetGeoTransform(), db.GetGeoTransform()
Aa = da.GetRasterBand(1).ReadAsArray().astype("f8")
Ab = db.GetRasterBand(1).ReadAsArray().astype("f8")

print("origins   default %.4f %.4f" % (ga[0], ga[3]))
print("          explicit %.4f %.4f" % (gb[0], gb[3]))
print("          delta    %.4f %.4f m" % (ga[0] - gb[0], ga[3] - gb[3]))
print("shape     %s vs %s" % (Aa.shape, Ab.shape))

if Aa.shape == Ab.shape:
    m = np.isfinite(Aa) & np.isfinite(Ab)
    d = np.abs(Aa[m] - Ab[m])
    print("\ncell-for-cell elevation difference:")
    print("  identical cells : %.2f%%" % (100 * (d < 1e-9).mean()))
    print("  max difference  : %.6f m" % d.max())
    print("  mean difference : %.6f m" % d.mean())
    print("  rms difference  : %.6f m" % np.sqrt((d ** 2).mean()))
    verdict = ("only the origin label differs; the surface is identical"
               if d.max() < 1e-6 else
               "the surface itself differs, so the grid was resampled")
    print("\n  ->", verdict)

# what a canonically-aligned grid would look like
print("\nfor reference, a grid snapped to whole metres would start at")
print("  x %.2f   y %.2f" % (np.floor(ga[0]), np.ceil(ga[3])))
print("  offset from default origin: %.2f  %.2f m"
      % (ga[0] - np.floor(ga[0]), np.ceil(ga[3]) - ga[3]))
