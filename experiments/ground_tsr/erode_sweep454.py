"""How deep should the erosion be, measured rather than inherited?

21 m is the figure from mask_dem_edges.py, where it was applied to a DEM that
Surfer had already bounded. This interpolator reaches much further out -- 67.4
percent of the lattice against the archive's 43 -- so eroding 21 m from our
boundary lands somewhere quite different from eroding 21 m from Surfer's.

Two criteria, reported together because they may not agree: the depth at which
our coverage matches the archive's, and the depth at which the elevation at
the new boundary stops being extrapolation.
"""
import json
import os

import numpy as np
import rasterio
from scipy import ndimage

SRC = (r"C:\Replicalm\test100526\centro454_nofix50"
       r"\AMIGACarb_Yuc_Centro_GLAS_Apr2013_l0s454_DEM.tif")
ARCH = (r"D:\_Archive_EdgeFixed\Yuc_Centro\Centro_Glas\Centro_Glas"
        r"\Centro_GLAS_l0s454_DEM_0p5m_v1.tif")
OUT = r"C:\Replicalm\test100526\centro454_polyedge"

with rasterio.open(SRC) as r:
    Z = r.read(1).astype("f8")
    cell, nd = abs(r.transform.a), r.nodata
valid = np.isfinite(Z) & (Z != 0.0)
if nd is not None:
    valid &= Z != nd
closed = ndimage.binary_fill_holes(valid)

with rasterio.open(ARCH) as r:
    A = r.read(1).astype("f8")
    andv = r.nodata
if andv is not None:
    A[A == andv] = np.nan
A[A > 1e30] = np.nan
A[A == 0] = np.nan
Am = np.isfinite(A)
n_arch = int(Am.sum())

# distance inward from the untrimmed boundary: one transform serves every depth
din = ndimage.distance_transform_edt(closed) * cell
print("untrimmed active %d cells, archive %d, lattice %d"
      % (int(closed.sum()), n_arch, closed.size))
print("")
print("%-8s %11s %10s %12s %11s %11s %9s"
      % ("erode", "cells", "vs arch", "RMSE at the", "RMSE whole",
         "MAE whole", "bias"))
print("%-8s %11s %10s %12s %11s %11s %9s"
      % ("", "", "", "new edge", "area", "area", ""))
print("-" * 80)
rows = []
for d in (0, 10, 21, 30, 40, 50, 60, 80, 100):
    m = closed & (din > d)
    if m.sum() < 10000:
        continue
    both = m & Am & valid
    if both.sum() < 10000:
        continue
    e = (Z - A)[both]
    # the band just inside the new boundary, 5 m wide: is it still extrapolating?
    band = both & (din <= d + 5.0)
    eb = (Z - A)[band] if band.sum() > 1000 else np.array([np.nan])
    row = {"erode_m": d, "cells": int(m.sum()),
           "vs_archive_pct": 100.0 * (m.sum() / n_arch - 1),
           "edge_band_rmse": float(np.sqrt(np.nanmean(eb ** 2))),
           "rmse": float(np.sqrt(np.mean(e ** 2))),
           "mae": float(np.mean(np.abs(e))), "bias": float(np.mean(e))}
    rows.append(row)
    print("%-8s %11d %9.1f%% %12.4f %11.4f %11.4f %+9.4f"
          % ("%d m" % d, row["cells"], row["vs_archive_pct"],
             row["edge_band_rmse"], row["rmse"], row["mae"], row["bias"]))

match = min(rows, key=lambda r: abs(r["vs_archive_pct"]))
clean = next((r for r in rows if r["edge_band_rmse"] < 0.10), None)
print("")
print("coverage matches the archive nearest at %d m (%+.1f%%)"
      % (match["erode_m"], match["vs_archive_pct"]))
if clean:
    print("the new boundary stops extrapolating at %d m "
          "(edge band RMSE %.4f m)" % (clean["erode_m"],
                                       clean["edge_band_rmse"]))
else:
    print("no tested depth brings the edge band under 0.10 m")
with open(os.path.join(OUT, "erode_sweep.json"), "w") as fh:
    json.dump(rows, fh, indent=2)
print("wrote %s" % os.path.join(OUT, "erode_sweep.json"))
