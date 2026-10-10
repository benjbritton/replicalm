"""Where does l0s454's error actually sit?

RMSE 0.70 m against MAE 0.087 m is not a surface that is uniformly worse; it
is a surface that is fine almost everywhere and catastrophically wrong in a few
places. On the other two tiles the same measurement gave 0.079 and 0.055 m, so
either this tile is genuinely different or something is wrong with the inputs.

This locates the bad cells rather than guessing: the error distribution, the
coverage of each surface, whether the worst cells cluster, and whether the
extreme values look like terrain or like a sentinel that was not caught.
"""
import numpy as np
import rasterio
from scipy import ndimage

TS = (r"D:\_Archive_EdgeFixed\Yuc_Centro\Centro_Glas\Centro_Glas"
      r"\Centro_GLAS_l0s454_DEM_0p5m_v1.tif")
SM = r"C:\Replicalm\test100526\centro454_050\dem_smrf_050.tif"
TR = r"C:\Replicalm\test100526\centro454_050\dem_tsr_050.tif"


def read(p):
    with rasterio.open(p) as r:
        z = r.read(1).astype("f8")
        nd = r.nodata
        tr, cc = r.transform, abs(r.transform.a)
    print("  %-22s nodata=%s  raw min %.4f  raw max %.4f"
          % (p.rsplit("\\", 1)[-1], nd, float(np.nanmin(z)),
             float(np.nanmax(z))))
    if nd is not None:
        z[z == nd] = np.nan
    z[z > 1e30] = np.nan
    z[z == 0.0] = np.nan
    return z, tr, cc


print("RAW FILE RANGES, before any masking")
A, at, ac = read(TS)
S, _, _ = read(SM)
T, _, _ = read(TR)

print("")
print("AFTER MASKING")
for nm, z in (("T/S", A), ("smrf", S), ("tsr", T)):
    v = z[np.isfinite(z)]
    print("  %-5s %9d cells, min %9.3f, p1 %8.3f, p50 %8.3f, p99 %8.3f, "
          "max %9.3f" % (nm, len(v), v.min(), np.percentile(v, 1),
                         np.percentile(v, 50), np.percentile(v, 99),
                         v.max()))

common = np.isfinite(A) & np.isfinite(S) & np.isfinite(T)
e = (S - A)[common]
print("")
print("SMRF MINUS T/S on %d common cells" % int(common.sum()))
for q in (0.1, 1, 5, 25, 50, 75, 95, 99, 99.9):
    print("  p%-5s %+9.4f m" % (q, float(np.percentile(e, q))))
print("  mean %+.4f, RMSE %.4f, MAE %.4f" % (e.mean(),
                                             np.sqrt((e ** 2).mean()),
                                             np.abs(e).mean()))
for t in (0.25, 0.5, 1.0, 2.0, 5.0, 10.0):
    n = int((np.abs(e) > t).sum())
    print("  |e| > %5.2f m: %8d cells  %6.3f%%  carrying %5.1f%% of the "
          "squared error" % (t, n, 100.0 * n / len(e),
                             100.0 * (e[np.abs(e) > t] ** 2).sum()
                             / (e ** 2).sum()))

bad = np.zeros(A.shape, bool)
bad[common] = np.abs(S - A)[common] > 1.0
lab, n = ndimage.label(bad)
print("")
print("CLUSTERING of cells over 1 m: %d connected patches" % n)
if n:
    sz = ndimage.sum(bad, lab, np.arange(1, n + 1))
    order = np.argsort(sz)[::-1][:6]
    print("  largest patches, cells and location")
    for k in order:
        ys, xs = np.where(lab == (k + 1))
        cx = at.c + (xs.mean() + 0.5) * ac
        cy = at.f - (ys.mean() + 0.5) * ac
        rows = "rows %d-%d" % (ys.min(), ys.max())
        print("    %7d cells at %.1f, %.1f  (%s)  mean error %+.3f m"
              % (int(sz[k]), cx, cy, rows,
                 float(np.mean((S - A)[lab == (k + 1)]))))
    print("  top 6 patches hold %.1f%% of all over-1 m cells"
          % (100.0 * sz[order].sum() / sz.sum()))

print("")
print("EDGE PROXIMITY of the bad cells")
dist = ndimage.distance_transform_edt(common) * ac
print("  all common cells      median distance to edge %.1f m"
      % float(np.median(dist[common])))
if bad.any():
    print("  cells over 1 m        median distance to edge %.1f m"
          % float(np.median(dist[bad])))
    print("  fraction of bad cells within 25 m of an edge %.1f%%"
          % (100.0 * float((dist[bad] < 25).mean())))
