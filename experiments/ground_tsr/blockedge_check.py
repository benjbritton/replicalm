"""Is the coverage boundary the data's edge, or my solver blocks' edge?

`surfer_grid_blocked` fills a whole 150 m core wherever it managed to solve,
so every cell of that core carries a value whether or not a ground return is
anywhere near it. If that is what is happening, the "active data area" is a
union of 150 m squares and the polygon traced around it is a staircase of
block edges, not the footprint of the survey.

Two tests. How far each active cell sits from the nearest ground return,
ours against the archive's; and whether the boundary cells pile up at
multiples of the 150 m block pitch.
"""
import json
import os
import sys

import numpy as np
import rasterio
from scipy import ndimage

sys.path.insert(0, r"C:\Replicalm\src")
from replicalm import classify
from replicalm.preprocess import dedupe_xy

DEM = (r"C:\Replicalm\test100526\centro454_nofix50"
       r"\AMIGACarb_Yuc_Centro_GLAS_Apr2013_l0s454_DEM.tif")
GND = (r"C:\Replicalm\test100526\centro454_nofix50"
       r"\AMIGACarb_Yuc_Centro_GLAS_Apr2013_l0s454_ground_cleaned.laz")
ARCH = (r"D:\_Archive_EdgeFixed\Yuc_Centro\Centro_Glas\Centro_Glas"
        r"\Centro_GLAS_l0s454_DEM_0p5m_v1.tif")
CORE_M = 150.0

with rasterio.open(DEM) as r:
    Z = r.read(1).astype("f8")
    at, cell, nd = r.transform, abs(r.transform.a), r.nodata
ours = np.isfinite(Z) & (Z != 0.0)
if nd is not None:
    ours &= Z != nd
with rasterio.open(ARCH) as r:
    A = r.read(1).astype("f8")
    andv = r.nodata
if andv is not None:
    A[A == andv] = np.nan
A[A > 1e30] = np.nan
A[A == 0] = np.nan
arch = np.isfinite(A)
H, W = Z.shape
print("ours %d active cells, archive %d, lattice %d"
      % (int(ours.sum()), int(arch.sum()), ours.size))

p = classify._pdal().Pipeline(json.dumps({"pipeline": [
    {"type": "readers.las", "filename": GND, "nosrs": True}]}))
p.execute()
a = p.arrays[0]
gx, gy, gz, _ = dedupe_xy(a["X"].astype("f8"), a["Y"].astype("f8"),
                          a["Z"].astype("f8"), keep="mean")
del a, p
print("ground returns in the delivered cloud: %d" % len(gz))

col = np.floor((gx - at.c) / cell).astype(np.int64)
row = np.floor((at.f - gy) / cell).astype(np.int64)
ok = (col >= 0) & (col < W) & (row >= 0) & (row < H)
occ = np.zeros((H, W), bool)
occ[row[ok], col[ok]] = True
print("cells containing at least one ground return: %d (%.2f%% of lattice)"
      % (int(occ.sum()), 100.0 * occ.mean()))

dist = ndimage.distance_transform_edt(~occ, sampling=cell)
print("")
print("DISTANCE FROM AN ACTIVE CELL TO THE NEAREST GROUND RETURN")
print("%-10s %9s %9s %9s %9s %9s %9s"
      % ("mask", "p50", "p90", "p99", "p99.9", "max", ">20 m"))
for nm, m in (("ours", ours), ("archive", arch)):
    d = dist[m]
    print("%-10s %9.2f %9.2f %9.2f %9.2f %9.2f %8.2f%%"
          % (nm, *[float(np.percentile(d, q)) for q in (50, 90, 99, 99.9)],
             float(d.max()), 100.0 * float((d > 20.0).mean())))

print("")
print("ACTIVE CELLS FARTHER THAN A GIVEN DISTANCE FROM ANY RETURN")
print("%-10s %12s %12s" % ("threshold", "ours", "archive"))
for t in (2.0, 5.0, 10.0, 20.0, 50.0, 100.0):
    print("%-10s %12d %12d"
          % ("%.0f m" % t, int((ours & (dist > t)).sum()),
             int((arch & (dist > t)).sum())))

# block alignment of the boundary
pitch = int(round(CORE_M / cell))
edge = ours ^ ndimage.binary_erosion(ours)
ys, xs = np.nonzero(edge)
print("")
print("BOUNDARY ALIGNMENT, block pitch %d cells (%.1f m)" % (pitch, CORE_M))
for nm, idx in (("column", xs), ("row", ys)):
    r_ = idx % pitch
    near = float(np.mean((r_ < 3) | (r_ > pitch - 4)))
    print("  %-7s boundary cells within 3 cells of a block line: %.1f%% "
          "(chance %.1f%%)" % (nm, 100.0 * near, 100.0 * 6.0 / pitch))

eh = np.bincount(xs % pitch, minlength=pitch)
top = np.argsort(eh)[::-1][:6]
print("  most-occupied column offsets within a block: %s"
      % ", ".join("%d (%d cells)" % (int(t), int(eh[t])) for t in top))
print("  mean cells per offset: %.0f" % eh.mean())
print("")
if float(np.mean(((xs % pitch) < 3) | ((xs % pitch) > pitch - 4))) > 0.25:
    print("VERDICT: the boundary follows the solver's block lines")
else:
    print("VERDICT: the boundary does not align to block lines")
