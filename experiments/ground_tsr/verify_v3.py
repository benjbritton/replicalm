"""Verify the v3 run: dimensions, ends, frame, vectors and the report."""
import json
import os

import numpy as np
import rasterio
from scipy import ndimage

D = r"C:\Replicalm\test100526\centro454_e2e_v3"
BASE = "AMIGACarb_Yuc_Centro_GLAS_Apr2013_l0s454"
V2 = (r"C:\Replicalm\test100526\centro454_e2e"
      r"\AMIGACarb_Yuc_Centro_GLAS_Apr2013_l0s454_DEM.tif")
ARCH = (r"D:\_Archive_EdgeFixed\Yuc_Centro\Centro_Glas\Centro_Glas"
        r"\Centro_GLAS_l0s454_DEM_0p5m_v1.tif")


def rd(p):
    with rasterio.open(p) as r:
        z = r.read(1).astype("f8")
        nd, cc, tr = r.nodata, abs(r.transform.a), r.transform
        sz = (r.height, r.width)
    if nd is not None:
        z[z == nd] = np.nan
    z[z > 1e30] = np.nan
    z[z == 0] = np.nan
    return z, cc, tr, sz


Z, cell, tr, sz = rd(os.path.join(D, BASE + "_DEM.tif"))
P, _, _, szp = rd(V2)
A, _, _, sza = rd(ARCH)
M, mp, ma = np.isfinite(Z), np.isfinite(P), np.isfinite(A)
print("DIMENSIONS  v3 %s, previous %s, archive %s  identical: %s"
      % (sz, szp, sza, sz == szp == sza))
print("cells: v3 %d, previous %d, archive %d"
      % (int(M.sum()), int(mp.sum()), int(ma.sum())))
print("v3 vs archive %+.2f%%" % (100.0 * (M.sum() / ma.sum() - 1)))

H, W = sz
rows = np.flatnonzero(M.any(axis=1))
strip = int(round(21.0 / cell))
print("")
print("ENDS: data now spans rows %d-%d of %d" % (rows[0], rows[-1], H))
ra = np.flatnonzero(ma.any(axis=1))
print("  archive spans %d-%d; ours is inset %d at the north, %d at the south"
      % (ra[0], ra[-1], int(rows[0]), int(H - 1 - rows[-1])))
print("  expected inset %d rows (21 m): %s"
      % (strip, rows[0] >= strip - 2 and (H - 1 - rows[-1]) >= strip - 2))

west = np.array([np.flatnonzero(M[r])[0] for r in rows])
print("")
print("FRAME: rows whose west edge is at column 0 (kept, not stripped): %d"
      % int((west == 0).sum()))
print("  rows whose west edge sits exactly at the strip width: %d"
      % int((west == strip).sum()))
recovered = int((mp & ~ma).sum() - (M & ~ma).sum())
print("  cells beyond the archive that v3 keeps vs the previous run: %+d"
      % -recovered)

common = M & ma
dist = ndimage.distance_transform_edt(common) * cell
print("")
print("%-10s %10s %9s %9s %9s" % ("setback", "cells", "RMSE", "MAE", "bias"))
for sb in (0.0, 25.0, 50.0):
    m = common & (dist > sb)
    e = (Z - A)[m]
    print("%-10s %10d %9.4f %9.4f %+9.4f"
          % ("%.0f m" % sb, int(m.sum()), float(np.sqrt(np.mean(e ** 2))),
             float(np.mean(np.abs(e))), float(np.mean(e))))


def shp(m):
    f = ndimage.binary_fill_holes(m)
    per = int((f ^ ndimage.binary_erosion(f)).sum())
    ar = int(f.sum())
    _l, k = ndimage.label(m)
    return ar, per, per / max(2.0 * np.sqrt(np.pi * ar), 1e-9), int(k)


print("")
print("%-12s %10s %10s %11s %7s" % ("mask", "area", "perimeter",
                                    "raggedness", "pieces"))
for nm, m in (("archive", ma), ("v3", M)):
    a_, p_, r_, k_ = shp(m)
    print("%-12s %10d %10d %11.2f %7d" % (nm, a_, p_, r_, k_))

gj = os.path.join(D, BASE + "_footprint.geojson")
if os.path.exists(gj):
    with open(gj) as fh:
        fc = json.load(fh)
    print("")
    print("VECTORS: %s" % os.path.basename(gj))
    for f in fc["features"]:
        pr = f["properties"]
        print("  %-16s %10.3f ha, %d ring(s)"
              % (pr["kind"], pr["area_ha"],
                 len(f["geometry"]["coordinates"])))

cfg = os.path.join(D, BASE + "_config.json")
if os.path.exists(cfg):
    with open(cfg) as fh:
        c = json.load(fh)
    rep = c.get("report", c)
    print("")
    print("RUN REPORT")
    for k in ("footprint_radius_m", "edge_strip_m", "edge_axis"):
        print("  %-22s %s" % (k, rep.get(k)))
    fp = rep.get("footprint")
    if isinstance(fp, dict):
        for k in sorted(fp):
            if k != "vectors":
                print("  footprint.%-12s %s" % (k, fp[k]))
        if isinstance(fp.get("vectors"), dict):
            v = fp["vectors"]
            print("  footprint.vectors     %d feature(s), %d active, %d void"
                  % (v["features"], v["active"], v["voids"]))
