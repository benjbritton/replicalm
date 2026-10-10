"""l0s454 scored away from the coverage boundary.

The diagnostic established that 0.94 percent of cells exceed 1 m and carry 99.1
percent of the squared error, that those cells sit a median 4.6 m from the edge
of the common area against 82.9 m for cells in general, and that 91.4 percent
of them are within 25 m of an edge. So the tile-wide RMSE is measuring the trim,
not the surface.

This rescores with the boundary excluded, at several setbacks, so the figure
can be compared with l0s441 and South 395 -- and so the setback at which the
artifact stops dominating is a measured number rather than a choice.
"""
import json

import numpy as np
import rasterio
from scipy import ndimage

TS = (r"D:\_Archive_EdgeFixed\Yuc_Centro\Centro_Glas\Centro_Glas"
      r"\Centro_GLAS_l0s454_DEM_0p5m_v1.tif")
SM = r"C:\Replicalm\test100526\centro454_050\dem_smrf_050.tif"
TR = r"C:\Replicalm\test100526\centro454_050\dem_tsr_050.tif"
OUT = r"C:\Replicalm\test100526\centro454_050\interior_metrics.json"
PROM_T = 0.25


def read(p):
    with rasterio.open(p) as r:
        z = r.read(1).astype("f8")
        nd, tr, cc = r.nodata, r.transform, abs(r.transform.a)
    if nd is not None:
        z[z == nd] = np.nan
    z[z > 1e30] = np.nan
    z[z == 0.0] = np.nan
    return z, tr, cc


A, at, ac = read(TS)
S, _, _ = read(SM)
T, _, _ = read(TR)
common = np.isfinite(A) & np.isfinite(S) & np.isfinite(T)
dist = ndimage.distance_transform_edt(common) * ac

num = ndimage.gaussian_filter(np.where(common, A, 0.0), 2.0 / ac / 2.0)
den = ndimage.gaussian_filter(common.astype("f8"), 2.0 / ac / 2.0)
s2 = np.where(den > 1e-6, num / np.maximum(den, 1e-6), np.nan)
num = ndimage.gaussian_filter(np.where(common, A, 0.0), 10.0 / ac / 2.0)
den = ndimage.gaussian_filter(common.astype("f8"), 10.0 / ac / 2.0)
s10 = np.where(den > 1e-6, num / np.maximum(den, 1e-6), np.nan)
prom = (s2 - s10)
CLS = np.zeros(A.shape, np.int8)
CLS[common & (prom > PROM_T)] = 1
CLS[common & (prom < -PROM_T)] = 2
NM = {0: "planar", 1: "mound", 2: "depression"}

res = {"cell_m": ac, "setbacks": {}}
print("%-9s %10s %8s %9s %9s %9s %9s %9s"
      % ("setback", "cells", "kept", "smrf RMSE", "tsr RMSE", "smrf MAE",
         "tsr MAE", "tsr-smrf"))
print("-" * 82)
for sb in (0.0, 5.0, 10.0, 25.0, 50.0):
    m = common & (dist > sb)
    if m.sum() < 1000:
        continue
    es, et = (S - A)[m], (T - A)[m]
    d = (T - S)[m]
    row = {"cells": int(m.sum()),
           "kept_fraction": float(m.sum() / common.sum()),
           "smrf_rmse": float(np.sqrt(np.mean(es ** 2))),
           "tsr_rmse": float(np.sqrt(np.mean(et ** 2))),
           "smrf_mae": float(np.mean(np.abs(es))),
           "tsr_mae": float(np.mean(np.abs(et))),
           "smrf_bias": float(np.mean(es)),
           "tsr_bias": float(np.mean(et)),
           "tsr_minus_smrf_rmse": float(np.sqrt(np.mean(d ** 2))),
           "by_form": {}}
    for k in (0, 1, 2):
        mk = m & (CLS == k)
        if mk.sum() > 100:
            row["by_form"][NM[k]] = {
                "smrf": float(np.sqrt(np.mean((S - A)[mk] ** 2))),
                "tsr": float(np.sqrt(np.mean((T - A)[mk] ** 2)))}
    res["setbacks"]["%.0f" % sb] = row
    print("%-9s %10d %7.1f%% %9.4f %9.4f %9.4f %9.4f %9.4f"
          % ("%.0f m" % sb, row["cells"], 100 * row["kept_fraction"],
             row["smrf_rmse"], row["tsr_rmse"], row["smrf_mae"],
             row["tsr_mae"], row["tsr_minus_smrf_rmse"]))

print("")
print("BY LANDFORM at a 25 m setback")
r25 = res["setbacks"]["25"]
print("%-6s %10s %10s %12s" % ("run", "planar", "mound", "depression"))
for tag in ("smrf", "tsr"):
    print("%-6s " % tag + " ".join(
        "%10.4f" % r25["by_form"][n][tag] if n in r25["by_form"] else
        "%10s" % "-" for n in ("planar", "mound", "depression")))
print("")
print("bias at 25 m: smrf %+.4f, tsr %+.4f"
      % (r25["smrf_bias"], r25["tsr_bias"]))
with open(OUT, "w") as fh:
    json.dump(res, fh, indent=2)
print("")
print("wrote %s" % OUT)
