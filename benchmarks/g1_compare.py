r"""Put the Replicalm G1 beside the archive G1, at the same cells.

Two products are written:

  overview   the whole 1.4 x 8.0 km strip, downsampled, both versions side by
             side. Useful for spotting anything structural -- block seams, a
             coverage difference, a tonal drift along the strip.
  details    full-resolution crops at a few places, unscaled, so the question
             "does this resolve the same structures" can be answered by eye
             rather than by a statistic.

Crops are placed where the elevation difference is SMALL, not where it is
large. A crop chosen at the worst disagreement would show the pipeline at its
worst and say nothing about the typical case; the detail views are here to show
what the method normally produces, and the overview already carries the tail.
One crop is placed at a block boundary on purpose, since a blocked driver has to
be shown not to seam.
"""
import os, sys
import numpy as np
from osgeo import gdal
gdal.UseExceptions()

CAND = (r"C:\Replicalm\render\rvt\Yuc_South\South_Glas\South_Glas"
        r"\South_GLAS_l0s395_G1_0p5m_v1.tif")
ARCH = (r"D:\_Archive_EdgeFixed\Yuc_South\South_Glas\South_Glas"
        r"\South_GLAS_l0s395_G1_0p5m_v1.tif")
CDEM = (r"C:\Replicalm\render\dem\Yuc_South\South_Glas\South_Glas"
        r"\South_GLAS_l0s395_DEM_0p5m_v1.tif")
RDEM = (r"D:\_Archive_EdgeFixed\Yuc_South\South_Glas\South_Glas"
        r"\South_GLAS_l0s395_DEM_0p5m_v1.tif")
OUT = r"C:\Replicalm\render\compare"
os.makedirs(OUT, exist_ok=True)
CROP = 700              # detail crop size in cells, 350 m at 0.5 m
BLOCK_ROWS = 2000       # where the driver's block boundaries fall


def rgb(path, x, y, w, h, bw=None, bh=None):
    d = gdal.Open(path)
    w = min(w, d.RasterXSize - x); h = min(h, d.RasterYSize - y)
    kw = {}
    if bw: kw = {"buf_xsize": bw, "buf_ysize": bh}
    n = min(3, d.RasterCount)
    a = np.dstack([d.GetRasterBand(i + 1).ReadAsArray(x, y, w, h, **kw)
                   for i in range(n)])
    if n == 1:
        a = np.repeat(a, 3, axis=2)
    return a.astype("u1")


def save(arr, path):
    from PIL import Image
    Image.fromarray(arr).save(path)
    print("  %s  %d x %d" % (path, arr.shape[1], arr.shape[0]))


def pad_to(a, h):
    if a.shape[0] >= h:
        return a[:h]
    return np.vstack([a, np.zeros((h - a.shape[0],) + a.shape[1:], "u1")])


if not os.path.exists(CAND):
    sys.exit("the Replicalm G1 is not written yet: %s" % CAND)

d = gdal.Open(CAND)
W, H = d.RasterXSize, d.RasterYSize
print("G1 %d x %d" % (W, H))

# --- overview, both strips side by side -------------------------------------
sc = 8
ow, oh = W // sc, H // sc
a = rgb(CAND, 0, 0, W, H, ow, oh)
b = rgb(ARCH, 0, 0, W, H, ow, oh)
gap = np.zeros((oh, 12, 3), "u1"); gap[:, :] = 60
save(np.hstack([b, gap, a]), os.path.join(OUT, "overview_archive_left_replicalm_right.png"))

# --- find quiet places to crop ----------------------------------------------
# the dataset must outlive the band: GDAL frees the band with its dataset, and
# a temporary gdal.Open(...) expression leaves the band pointing at freed memory
cds, rds = gdal.Open(CDEM), gdal.Open(RDEM)
ca = cds.GetRasterBand(1)
cb = rds.GetRasterBand(1)
nod = cb.GetNoDataValue()
best = []
for y in range(0, H - CROP, CROP):
    for x in range(0, W - CROP, CROP):
        u = ca.ReadAsArray(x, y, CROP, CROP).astype("f8")
        v = cb.ReadAsArray(x, y, CROP, CROP).astype("f8")
        m = (u > -9000) & np.isfinite(v) & (v != 0.0)
        if nod is not None:
            m &= v != nod
        if m.mean() < 0.9:
            continue
        diff = np.abs(u[m] - v[m])
        best.append((float(np.median(diff)), float((diff > 0.5).mean()),
                     x, y, float(v[m].std())))
best.sort(key=lambda t: t[1])
# the most textured quiet windows first: flat ground shows nothing either way
quiet = sorted(best[:max(6, len(best) // 4)], key=lambda t: -t[4])[:3]

seam_y = ((H // 2) // BLOCK_ROWS) * BLOCK_ROWS      # a real block boundary
seam = None
for med, bad, x, y, sd in best:
    if abs((y + CROP // 2) - seam_y) < CROP:
        seam = (med, bad, x, seam_y - CROP // 2, sd)
        break

picks = [("detail%d" % i, p) for i, p in enumerate(quiet, 1)]
if seam:
    picks.append(("seam_at_row_%d" % seam_y, seam))

for name, (med, bad, x, y, sd) in picks:
    y = max(0, y)
    left = rgb(ARCH, x, y, CROP, CROP)
    right = rgb(CAND, x, y, CROP, CROP)
    h = max(left.shape[0], right.shape[0])
    left, right = pad_to(left, h), pad_to(right, h)
    gap = np.zeros((h, 12, 3), "u1"); gap[:, :] = 60
    print("  %-18s at %5d,%-6d  median diff %.3f m, %.2f%% over 0.5 m, "
          "relief sd %.2f m" % (name, x, y, med, 100 * bad, sd))
    save(np.hstack([left, gap, right]),
         os.path.join(OUT, "%s_archive_left_replicalm_right.png" % name))

print("\nleft is the archive G1 (TerraScan), right is Replicalm")
