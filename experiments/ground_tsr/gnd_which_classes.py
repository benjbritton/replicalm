"""Which delivered classes lie on the archive surface, and where they came from.

Two measurements, both direct.

1. Residual of each delivered point against the archive DEM. The points that
   generated a TIN surface must lie on it: their residual distribution is
   centred on zero and narrow, limited only by the 0.5 m grid resampling of a
   triangulated surface. Points that were present but not used have no such
   constraint. This distinguishes the sets by measurement instead of by
   convention.

2. Whether the delivered ground classes came from the flightline's own class 2
   or from its unclassified returns. The original LAS is already carried at
   class 2 / class 1 by the provider, so a coordinate match says whether the
   delivery re-labelled an existing ground set or classified the full cloud.
   Coordinates are matched on the coarser of the two scales, 0.01 m.
"""
import json
import sys

import numpy as np
import rasterio

sys.path.insert(0, r"C:\Replicalm\src")
from replicalm import classify

GND = (r"D:\_LAZ\Yuc_Centro\Centro_Glas\Centro_Glas"
       r"\Centro_GLAS_l0s441_GND_v1.laz")
RAW = (r"D:\GLiHT_LAS_orig\Yuc_Centro\Centro_Glas\Centro_Glas"
       r"\AMIGACarb_Yuc_Centro_GLAS_Apr2013_l0s441.las")
ARCH = (r"D:\_Archive_EdgeFixed\Yuc_Centro\Centro_Glas\Centro_Glas"
        r"\Centro_GLAS_l0s441_DEM_0p5m_v1.tif")


def load(path, want_cls=False):
    p = classify._pdal().Pipeline(json.dumps({"pipeline": [
        {"type": "readers.las", "filename": path, "nosrs": True}]}))
    p.execute()
    a = p.arrays[0]
    out = [a["X"].astype("f8"), a["Y"].astype("f8"), a["Z"].astype("f8")]
    if want_cls:
        out.append(a["Classification"].astype(np.int32))
    return out


with rasterio.open(ARCH) as r:
    A = r.read(1).astype("f8")
    at, ac = r.transform, abs(r.transform.a)
    HT, WT = r.height, r.width
    nd = r.nodata
if nd is not None:
    A[A == nd] = np.nan
A[A > 1e30] = np.nan
A[A == 0] = np.nan

gx, gy, gz, gc = load(GND, want_cls=True)
print("delivery: %d points, classes %s"
      % (len(gz), ", ".join("%d=%d" % (v, int((gc == v).sum()))
                            for v in np.unique(gc))), flush=True)

col = np.floor((gx - at.c) / ac).astype(np.int64)
row = np.floor((at.f - gy) / ac).astype(np.int64)
ins = (col >= 0) & (col < WT) & (row >= 0) & (row < HT)
res = np.full(len(gz), np.nan)
res[ins] = gz[ins] - A[row[ins], col[ins]]

print("")
print("RESIDUAL AGAINST THE ARCHIVE DEM, point elevation minus surface")
print("%-8s %10s %9s %9s %9s %9s %9s %9s"
      % ("class", "n in grid", "median", "mean", "IQR", "|r|<0.02",
         "|r|<0.05", "|r|<0.25"))
for v in list(np.unique(gc)) + ["2+8"]:
    m = ins & (np.isin(gc, [2, 8]) if v == "2+8" else (gc == v))
    m &= np.isfinite(res)
    if m.sum() < 1000:
        continue
    r_ = res[m]
    q1, q3 = np.percentile(r_, [25, 75])
    print("%-8s %10d %9.4f %9.4f %9.4f %8.2f%% %8.2f%% %8.2f%%"
          % (v, int(m.sum()), float(np.median(r_)), float(r_.mean()),
             float(q3 - q1), 100 * float((np.abs(r_) < 0.02).mean()),
             100 * float((np.abs(r_) < 0.05).mean()),
             100 * float((np.abs(r_) < 0.25).mean())))

print("")
print("ARCHIVE CELL COVERAGE, valid cells containing at least one point")
Am = np.isfinite(A)
for v, lbl in ((2, "class 2"), (8, "class 8"), (None, "class 2 + 8")):
    m = ins & (np.isin(gc, [2, 8]) if v is None else (gc == v))
    hit = np.zeros((HT, WT), bool)
    hit[row[m], col[m]] = True
    print("  %-12s %10d cells of %d valid  %6.2f%%"
          % (lbl, int((hit & Am).sum()), int(Am.sum()),
             100 * float((hit & Am).sum()) / int(Am.sum())))

print("")
print("PROVENANCE OF THE DELIVERED POINTS IN THE FLIGHTLINE LAS")
rx, ry, rz, rc = load(RAW, want_cls=True)
print("  flightline: %d points, classes %s"
      % (len(rz), ", ".join("%d=%d" % (v, int((rc == v).sum()))
                            for v in np.unique(rc))), flush=True)


def key(x, y, z):
    xi = np.round(x * 100).astype(np.int64)
    yi = np.round(y * 100).astype(np.int64)
    zi = np.round(z * 100).astype(np.int64)
    return (xi * 2000003 + yi) * 2000003 + zi


kr = key(rx, ry, rz)
o = np.argsort(kr, kind="stable")
kr_s, rc_s = kr[o], rc[o]
kg = key(gx, gy, gz)
pos = np.searchsorted(kr_s, kg)
pos_c = np.clip(pos, 0, len(kr_s) - 1)
found = kr_s[pos_c] == kg
src = np.where(found, rc_s[pos_c], -1)
print("  matched on 0.01 m coordinates: %d of %d  %.2f%%"
      % (int(found.sum()), len(kg), 100 * float(found.mean())))
for v in np.unique(gc):
    m = gc == v
    print("    delivered class %-3d n=%-10d" % (v, int(m.sum())), end="")
    for s in (-1, 1, 2):
        lbl = "unmatched" if s == -1 else "flightline class %d" % s
        print("  %s %.2f%%" % (lbl, 100 * float((src[m] == s).mean())),
              end="")
    print("")
nmatch = int(found.sum())
print("  flightline class 2 not present in the delivery: %d of %d  %.2f%%"
      % (int((rc == 2).sum()) - nmatch, int((rc == 2).sum()),
         100 * (1 - nmatch / max(int((rc == 2).sum()), 1))))
