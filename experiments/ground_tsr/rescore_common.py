"""Rescore every surface on the cells all of them filled.

The TIN build left 31.3 percent of block cells outside a convex hull and
finished with 11,412,389 scored cells against the biharmonic's 12,278,982.
Those are different populations, and the cells a triangulation cannot reach
are the sparse ones, which are also the hard ones -- so a TIN scored only
where it succeeded is scored on an easier subset. Comparing the engines
requires the intersection.
"""
import json
import os

import numpy as np
import rasterio
from scipy import ndimage

ARCH = (r"D:\_Archive_EdgeFixed\Yuc_Centro\Centro_Glas\Centro_Glas"
        r"\Centro_GLAS_l0s441_DEM_0p5m_v1.tif")
RUNS = [
    ("tin_c2c8", r"C:\Replicalm\test100526\terrascan_tin_surface"
                 r"\dem_tin_c2c8.tif"),
    ("biharm_c2c8", r"C:\Replicalm\test100526\terrascan_ceiling"
                    r"\dem_ts_c2c8.tif"),
    ("smrf", r"C:\Replicalm\test100526\cellsize_validation"
             r"\dem_B_cell1.00.tif"),
    ("mcc", r"C:\Replicalm\test100526\mcc_fulltile\dem_mcc.tif"),
]
OUT = r"C:\Replicalm\test100526\terrascan_tin_surface"
PROM_T = 0.25

with rasterio.open(ARCH) as r:
    A = r.read(1).astype("f4")
    at, ac = r.transform, abs(r.transform.a)
    nd = r.nodata
if nd is not None:
    A[A == nd] = np.nan
A[A > 1e30] = np.nan
A[A == 0] = np.nan
Am = np.isfinite(A)


def nc_smooth(z, mask, scale):
    s = scale / ac / 2.0
    num = ndimage.gaussian_filter(np.where(mask, z, 0.0).astype("f4"), s)
    den = ndimage.gaussian_filter(mask.astype("f4"), s)
    return np.where(den > 1e-6, num / np.maximum(den, 1e-6), np.nan)


sm2 = nc_smooth(A, Am, 2.0)
sm5 = nc_smooth(A, Am, 5.0)
sm10 = nc_smooth(A, Am, 10.0)
prom = (sm2 - sm10).astype("f4")
curv = (ndimage.laplace(np.nan_to_num(sm2)) / (ac * ac)).astype("f4")
gy_, gx_ = np.gradient(np.nan_to_num(sm5), ac)
slope_r = np.degrees(np.arctan(np.hypot(gx_, gy_))).astype("f4")
CLS = np.zeros(A.shape, np.int8)
CLS[Am & (prom > PROM_T)] = 1
CLS[Am & (prom < -PROM_T)] = 2
sharp = np.abs(curv) > np.nanpercentile(np.abs(curv[Am]), 90)
CLS[(CLS == 0) & Am & (slope_r > 5) & sharp & (curv < 0)] = 3
NM = {0: "planar", 1: "mound", 2: "depression", 3: "shoulder"}
del sm2, sm5, sm10, curv, slope_r, sharp, prom

Z = {}
for tag, path in RUNS:
    with rasterio.open(path) as r:
        z = r.read(1).astype("f4")
    z[z == 0.0] = np.nan
    Z[tag] = z
    print("%-13s %11d valid cells" % (tag, int((Am & np.isfinite(z)).sum())))

common = Am.copy()
for tag in Z:
    common &= np.isfinite(Z[tag])
print("")
print("cells valid in the archive and all four surfaces: %d" % int(common.sum()))
print("that is %.2f%% of the archive and %.2f%% of the biharmonic's own set"
      % (100 * common.sum() / Am.sum(),
         100 * common.sum() / (Am & np.isfinite(Z["biharm_c2c8"])).sum()))

SHARE = {}
for kk in (0, 1, 2, 3):
    SHARE[NM[kk]] = float((CLS[common] == kk).mean())
print("")
print("landform shares on the common cells: " + ", ".join(
    "%s %.2f%%" % (n, 100 * SHARE[n]) for n in SHARE))


def detail(z, mask, scale):
    s = scale / ac / 2.0
    num = ndimage.gaussian_filter(np.where(mask, z, 0.0).astype("f4"), s)
    den = ndimage.gaussian_filter(mask.astype("f4"), s)
    sm = np.where(den > 1e-6, num / np.maximum(den, 1e-6), 0.0)
    return np.where(mask, z - sm, np.nan)


inner = ndimage.binary_erosion(common, structure=np.ones((3, 3)),
                               iterations=int(5.0 / ac))
k = max(3, int(round(1.0 / ac)) | 1)
dA = detail(A, common, 1.5)
dAm = ndimage.median_filter(np.nan_to_num(dA, nan=0.0), size=k)
sel0 = inner & np.isfinite(dA)
arch_coh = float(np.sqrt(np.mean(dAm[sel0] ** 2)))
arch_spk = float(np.sqrt(np.mean((dA - dAm)[sel0] ** 2)))

res = []
for tag, _ in RUNS:
    z = Z[tag]
    e = (z - A)[common]
    dB = detail(z, common, 1.5)
    dBm = ndimage.median_filter(np.nan_to_num(dB, nan=0.0), size=k)
    sel = sel0 & np.isfinite(dB)
    by, bias = {}, {}
    for kk in (0, 1, 2, 3):
        s2 = common & (CLS == kk)
        by[NM[kk]] = float(np.sqrt(np.mean((z - A)[s2] ** 2)))
        bias[NM[kk]] = float(np.mean((z - A)[s2]))
    res.append(dict(tag=tag, rmse=float(np.sqrt(np.mean(e ** 2))),
                    mae=float(np.mean(np.abs(e))), bias=float(np.mean(e)),
                    p99abs=float(np.percentile(np.abs(e), 99)),
                    coherent=float(np.sqrt(np.mean(dBm[sel] ** 2))),
                    speckle=float(np.sqrt(np.mean((dB - dBm)[sel] ** 2))),
                    rough=float(np.sqrt(np.mean(dB[sel] ** 2))
                                / max(np.sqrt(np.mean(dA[sel] ** 2)), 1e-9)),
                    rmse_by_form=by, bias_by_form=bias))

print("")
print("ON THE COMMON CELLS")
print("%-13s %9s %8s %9s %9s %9s %9s %7s"
      % ("run", "RMSE", "MAE", "bias", "p99|e|", "coherent", "speckle",
         "rough"))
print("-" * 82)
for r in res:
    print("%-13s %9.4f %8.4f %9.4f %9.4f %9.4f %9.4f %7.3f"
          % (r["tag"], r["rmse"], r["mae"], r["bias"], r["p99abs"],
             r["coherent"], r["speckle"], r["rough"]))
print("archive own texture on these cells: coherent %.4f, speckle %.4f"
      % (arch_coh, arch_spk))

print("")
print("RMSE BY LANDFORM, common cells")
print("%-13s %10s %10s %12s %10s" % ("run", "planar", "mound", "depression",
                                     "shoulder"))
for r in res:
    b = r["rmse_by_form"]
    print("%-13s " % r["tag"] + " ".join(
        "%10.4f" % b[n] for n in ("planar", "mound", "depression",
                                  "shoulder")))
print("")
print("MEAN SIGNED ERROR BY LANDFORM, common cells")
print("%-13s %10s %10s %12s %10s" % ("run", "planar", "mound", "depression",
                                     "shoulder"))
for r in res:
    b = r["bias_by_form"]
    print("%-13s " % r["tag"] + " ".join(
        "%+10.4f" % b[n] for n in ("planar", "mound", "depression",
                                   "shoulder")))
print("")
print("SHARE OF ERROR VARIANCE, area-weighted on the common cells")
print("%-13s %10s %10s %12s %10s %12s"
      % ("run", "planar", "mound", "depression", "shoulder", "non-planar"))
for r in res:
    b = r["rmse_by_form"]
    v = {n: SHARE[n] * b[n] ** 2 for n in SHARE}
    t = sum(v.values())
    print("%-13s " % r["tag"] + " ".join(
        "%9.1f%%" % (100 * v[n] / t) for n in ("planar", "mound",
                                               "depression", "shoulder"))
        + "%12s" % ("%.1f%%" % (100 * sum(
            v[n] for n in ("mound", "depression", "shoulder")) / t)))

with open(os.path.join(OUT, "metrics_common.json"), "w") as fh:
    json.dump({"cells": int(common.sum()), "shares": SHARE,
               "arch_coherent": arch_coh, "arch_speckle": arch_spk,
               "runs": res}, fh, indent=2)
print("")
print("wrote %s" % os.path.join(OUT, "metrics_common.json"))
