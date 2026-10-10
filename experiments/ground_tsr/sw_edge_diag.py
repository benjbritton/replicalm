"""What is happening at the south-west edge of centro454?

Two reported symptoms: seam that the strip did not remove, and a straight
vertical cut that took active data. Both are testable against the mask.

The across-track strip removes 21 m inward from the first and last valid cell
of each row. That treats a boundary running north-south. Where the transect
ends, the boundary runs east-west, and a per-row strip cannot touch it -- so
the along-track ends should show untreated seam by construction.

A straight vertical cut is the other side of the same mechanism: if one row
reaches further west than its neighbours, the strip starts from that row's own
outlying cell, and the next row starts 20 m east of it. Stacked, that is a
staircase, and a single long outlier makes one long vertical face.
"""
import numpy as np
import rasterio
from scipy import ndimage

E2E = (r"C:\Replicalm\test100526\centro454_e2e"
       r"\AMIGACarb_Yuc_Centro_GLAS_Apr2013_l0s454_DEM.tif")
PRE = (r"C:\Replicalm\test100526\centro454_nofix50"
       r"\AMIGACarb_Yuc_Centro_GLAS_Apr2013_l0s454_DEM.tif")
ARCH = (r"D:\_Archive_EdgeFixed\Yuc_Centro\Centro_Glas\Centro_Glas"
        r"\Centro_GLAS_l0s454_DEM_0p5m_v1.tif")


def rd(p):
    with rasterio.open(p) as r:
        z = r.read(1).astype("f8")
        nd, cc, tr = r.nodata, abs(r.transform.a), r.transform
    if nd is not None:
        z[z == nd] = np.nan
    z[z > 1e30] = np.nan
    z[z == 0] = np.nan
    return np.isfinite(z), cc, tr


M, cell, tr = rd(E2E)
P, _, _ = rd(PRE)
A, _, _ = rd(ARCH)
H, W = M.shape
print("grid %d x %d at %.4f m" % (H, W, cell))

rows = np.flatnonzero(M.any(axis=1))
r0, r1 = int(rows[0]), int(rows[-1])
print("data spans rows %d to %d" % (r0, r1))

# orientation of the boundary, row by row: how fast does the west edge move?
west = np.full(H, -1)
east = np.full(H, -1)
for r in rows:
    c = np.flatnonzero(M[r])
    west[r] = c[0]
    east[r] = c[-1]
dw = np.diff(west[rows])
print("")
print("WEST BOUNDARY STEP BETWEEN ADJACENT ROWS (cells)")
print("  median |step| %.1f, p95 %.1f, p99 %.1f, max %d"
      % (float(np.median(np.abs(dw))), float(np.percentile(np.abs(dw), 95)),
         float(np.percentile(np.abs(dw), 99)), int(np.abs(dw).max())))
big = np.flatnonzero(np.abs(dw) > 40)
print("  steps over 40 cells (20 m): %d" % len(big))
for k in big[:8]:
    rr = int(rows[k])
    print("    row %5d: west edge %d -> %d  (%+d cells, %+.1f m)"
          % (rr, west[rr], west[rows[k + 1]], dw[k], dw[k] * cell))

# the southern end: is the boundary there along-track?
south = rows[-1]
tail = rows[rows > south - 400]
print("")
print("SOUTHERN END, last 200 m of rows")
print("  rows %d to %d; west edge varies %d to %d, east edge %d to %d"
      % (int(tail[0]), int(tail[-1]), int(west[tail].min()),
         int(west[tail].max()), int(east[tail].min()), int(east[tail].max())))
widths = (east[tail] - west[tail] + 1)
print("  row width: first %d, last %d cells" % (int(widths[0]),
                                                int(widths[-1])))
untreated = int((P[tail[-1]].sum() - M[tail[-1]].sum()))
print("  cells removed on the very last row: %d (a strip would remove ~84)"
      % untreated)

# how much was removed at the southern tip versus the western flank?
band = np.zeros_like(M)
band[south - 200:south + 1, :] = True
flank = np.zeros_like(M)
flank[(r0 + r1) // 2 - 100:(r0 + r1) // 2 + 100, :] = True
for nm, b in (("south tip, 100 m", band), ("mid flank, 100 m", flank)):
    pre_n = int((P & b).sum())
    post_n = int((M & b).sum())
    print("  %-18s pre %8d, post %8d, removed %7d (%.1f%%)"
          % (nm, pre_n, post_n, pre_n - post_n,
             100.0 * (pre_n - post_n) / max(pre_n, 1)))

# straight vertical faces in the final mask
print("")
print("STRAIGHT VERTICAL FACES on the west boundary")
vals, counts = np.unique(west[rows], return_counts=True)
top = np.argsort(counts)[::-1][:6]
for t in top:
    col = int(vals[t])
    n = int(counts[t])
    if n < 20:
        continue
    rr = rows[west[rows] == col]
    x = tr.c + (col + 0.5) * cell
    print("  column %5d (x=%.1f) is the west edge for %5d rows, "
          "rows %d-%d, %.1f m tall"
          % (col, x, n, int(rr.min()), int(rr.max()),
             (rr.max() - rr.min() + 1) * cell))

# does the archive have the same faces?
westA = np.full(H, -1)
rowsA = np.flatnonzero(A.any(axis=1))
for r in rowsA:
    c = np.flatnonzero(A[r])
    westA[r] = c[0]
valsA, countsA = np.unique(westA[rowsA], return_counts=True)
print("  archive's longest vertical face: %d rows" % int(countsA.max()))
print("  ours: %d rows" % int(counts.max()))
