"""How wide is the defective band at the along-track ends?

The sides carry a 21 m stitching seam, which is where adjacent flight passes
were merged with differing classification. The ends have no adjacent pass, so
there is no seam there. What they have is the end-of-pass taper: the scan
narrows, ground returns thin, and the solver has support on one side only.

That is a different defect and need not be 21 m wide, so it is measured:
error against the archive and ground-return density, both as a function of
distance from the last row of data.
"""
import numpy as np
import rasterio
from scipy import ndimage

PRE = (r"C:\Replicalm\test100526\centro454_nofix50"
       r"\AMIGACarb_Yuc_Centro_GLAS_Apr2013_l0s454_DEM.tif")
ARCH = (r"D:\_Archive_EdgeFixed\Yuc_Centro\Centro_Glas\Centro_Glas"
        r"\Centro_GLAS_l0s454_DEM_0p5m_v1.tif")


def rd(p):
    with rasterio.open(p) as r:
        z = r.read(1).astype("f8")
        nd, cc = r.nodata, abs(r.transform.a)
    if nd is not None:
        z[z == nd] = np.nan
    z[z > 1e30] = np.nan
    z[z == 0] = np.nan
    return z, cc


Z, cell = rd(PRE)
A, _ = rd(ARCH)
mz, ma = np.isfinite(Z), np.isfinite(A)
H, W = Z.shape
rz = np.flatnonzero(mz.any(axis=1))
ra = np.flatnonzero(ma.any(axis=1))
print("our data spans rows %d-%d; the archive spans %d-%d"
      % (int(rz[0]), int(rz[-1]), int(ra[0]), int(ra[-1])))
print("so the archive is %d rows (%.1f m) short at the north end and "
      "%d rows (%.1f m) at the south"
      % (int(ra[0] - rz[0]), (ra[0] - rz[0]) * cell,
         int(rz[-1] - ra[-1]), (rz[-1] - ra[-1]) * cell))

both = mz & ma
print("")
print("AGAINST THE ARCHIVE, BY DISTANCE FROM THE NEAREST ALONG-TRACK END")
print("%-14s %10s %9s %9s %10s" % ("band", "cells", "RMSE", "MAE", "row width"))
for lo, hi in ((0, 10), (10, 20), (20, 30), (30, 50), (50, 75),
               (75, 100), (100, 200), (200, 400)):
    r_lo, r_hi = int(lo / cell), int(hi / cell)
    sel = np.zeros((H, W), bool)
    n0, n1 = int(rz[0]) + r_lo, int(rz[0]) + r_hi
    s1, s0 = int(rz[-1]) - r_hi, int(rz[-1]) - r_lo
    sel[n0:n1, :] = True
    sel[s1:s0, :] = True
    m = sel & both
    if m.sum() < 500:
        continue
    e = (Z - A)[m]
    wid = float(np.median(mz[n0:n1].sum(axis=1).tolist()
                          + mz[s1:s0].sum(axis=1).tolist()))
    print("%-14s %10d %9.4f %9.4f %10.0f"
          % ("%d-%d m" % (lo, hi), int(m.sum()),
             float(np.sqrt(np.mean(e ** 2))), float(np.mean(np.abs(e))), wid))

print("")
print("ROW WIDTH APPROACHING EACH END, our coverage")
print("%-10s %12s %12s" % ("from end", "north width", "south width"))
for d in (0, 5, 10, 20, 30, 50, 75, 100, 200, 400):
    k = int(d / cell)
    n = int(mz[min(int(rz[0]) + k, H - 1)].sum())
    s = int(mz[max(int(rz[-1]) - k, 0)].sum())
    print("%-10s %12d %12d" % ("%d m" % d, n, s))

mid = int((rz[0] + rz[-1]) / 2)
full = int(np.median([mz[r].sum() for r in range(mid - 200, mid + 200)]))
print("")
print("mid-transect row width is %d cells (%.0f m)" % (full, full * cell))
for nm, idx, step in (("north", int(rz[0]), 1), ("south", int(rz[-1]), -1)):
    d = None
    for k in range(0, 1200):
        r = idx + step * k
        if 0 <= r < H and mz[r].sum() >= 0.9 * full:
            d = k * cell
            break
    print("  %s end reaches 90%% of full width at %.1f m" % (nm, d))
