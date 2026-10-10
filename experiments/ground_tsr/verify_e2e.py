"""Does the pipeline's own output match what the scratch harness produced?

The harness result is only evidence if the production path reproduces it.
This compares the end-to-end DEM against the archive and against the scratch
version cell for cell.
"""
import json
import os

import numpy as np
import rasterio
from scipy import ndimage

E2E = (r"C:\Replicalm\test100526\centro454_e2e"
       r"\AMIGACarb_Yuc_Centro_GLAS_Apr2013_l0s454_DEM.tif")
HARNESS = (r"C:\Replicalm\test100526\centro454_autoedge"
           r"\dem_autoedge_20m_strip21m.tif")
ARCH = (r"D:\_Archive_EdgeFixed\Yuc_Centro\Centro_Glas\Centro_Glas"
        r"\Centro_GLAS_l0s454_DEM_0p5m_v1.tif")
CFG = (r"C:\Replicalm\test100526\centro454_e2e"
       r"\AMIGACarb_Yuc_Centro_GLAS_Apr2013_l0s454_config.json")


def rd(p):
    with rasterio.open(p) as r:
        z = r.read(1).astype("f8")
        nd, cc = r.nodata, abs(r.transform.a)
    if nd is not None:
        z[z == nd] = np.nan
    z[z > 1e30] = np.nan
    z[z == 0] = np.nan
    return z, cc


Z, cell = rd(E2E)
Hh, _ = rd(HARNESS)
A, _ = rd(ARCH)
mz, mh, ma = np.isfinite(Z), np.isfinite(Hh), np.isfinite(A)
print("pipeline %d cells, harness %d, archive %d"
      % (int(mz.sum()), int(mh.sum()), int(ma.sum())))
print("pipeline vs archive %+.2f%%, vs harness %+d cells"
      % (100.0 * (mz.sum() / ma.sum() - 1), int(mz.sum() - mh.sum())))
both = mz & mh
d = (Z - Hh)[both]
print("where both have data: %d cells, identical %d (%.4f%%), max |diff| %.6f"
      % (len(d), int((d == 0).sum()), 100.0 * (d == 0).mean(),
         float(np.abs(d).max())))


def shp(m):
    f = ndimage.binary_fill_holes(m)
    per = int((f ^ ndimage.binary_erosion(f)).sum())
    ar = int(f.sum())
    _l, k = ndimage.label(m)
    runs, rows = 0, int(m.any(axis=1).sum())
    for r in np.flatnonzero(m.any(axis=1))[::15]:
        dd = np.diff(np.concatenate([[0], m[r].astype(np.int8), [0]]))
        runs += int((dd == 1).sum())
    return (ar, per, per / max(2.0 * np.sqrt(np.pi * ar), 1e-9), int(k),
            runs / max(len(np.flatnonzero(m.any(axis=1))[::15]), 1))


print("")
print("%-12s %10s %10s %11s %7s %10s" % ("mask", "area", "perimeter",
                                         "raggedness", "pieces", "runs/row"))
for nm, m in (("archive", ma), ("pipeline", mz)):
    ar, per, rg, k, rr = shp(m)
    print("%-12s %10d %10d %11.2f %7d %10.2f" % (nm, ar, per, rg, k, rr))

common = mz & ma
dist = ndimage.distance_transform_edt(common) * cell
print("")
print("%-10s %10s %9s %9s %9s" % ("setback", "cells", "RMSE", "MAE", "bias"))
for sb in (0.0, 25.0, 50.0):
    m = common & (dist > sb)
    e = (Z - A)[m]
    print("%-10s %10d %9.4f %9.4f %+9.4f"
          % ("%.0f m" % sb, int(m.sum()), float(np.sqrt(np.mean(e ** 2))),
             float(np.mean(np.abs(e))), float(np.mean(e))))

if os.path.exists(CFG):
    with open(CFG) as fh:
        c = json.load(fh)
    rep = c.get("report", c)
    print("")
    print("RUN REPORT")
    for k in ("footprint_radius_m", "edge_strip_m", "edge_axis",
              "cells_after_trim", "interior_cells_restored"):
        if k in rep:
            print("  %-26s %s" % (k, rep[k]))
    fp = rep.get("footprint")
    if isinstance(fp, dict):
        for k in sorted(fp):
            print("  footprint.%-16s %s" % (k, fp[k]))
