"""Rejection of the full delivered ground population, class 2 and class 8.

Every rejection figure reported so far was measured against class 2 alone,
8,448,955 of the 17,705,532 delivered ground points, and the half used is the
half sitting systematically below the surface. This recomputes them against
the union, which is the only population with a reproducible definition, so
the SMRF and MCC numbers can be set beside the TIN prototype's.

Both cached classifications are reused, so the ground sets are the same ones
that produced the full-tile DEMs. The two 600 m cores match the prototype's
geometry exactly, which is what makes those rows comparable.
"""
import json
import os
import sys

import numpy as np
import rasterio
from rasterio.warp import transform as warp_transform
from scipy import ndimage
from scipy.spatial import cKDTree

sys.path.insert(0, r"C:\Replicalm\src")
from replicalm import classify
from replicalm.preprocess import dedupe_xy

GND = (r"D:\_LAZ\Yuc_Centro\Centro_Glas\Centro_Glas"
       r"\Centro_GLAS_l0s441_GND_v1.laz")
ARCH = (r"D:\_Archive_EdgeFixed\Yuc_Centro\Centro_Glas\Centro_Glas"
        r"\Centro_GLAS_l0s441_DEM_0p5m_v1.tif")
SMRF_G = (r"C:\Replicalm\test100526\cellsize_validation"
          r"\ground_B_cell1.00.npy")
MCC_G = r"C:\Replicalm\test100526\mcc_fulltile\ground_mcc.npy"
OUT = r"C:\Replicalm\test100526\axelsson"
PROM_T, NEAR, CORE_M = 0.25, 0.25, 600.0
FORENSIC = (("mound", -88.9673900, 19.2242934),
            ("depression", -88.9664570, 19.2223460))

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


def cls_at(x, y):
    c = np.clip(((x - at.c) / ac).astype(np.int64), 0, CLS.shape[1] - 1)
    r = np.clip(((at.f - y) / ac).astype(np.int64), 0, CLS.shape[0] - 1)
    return CLS[r, c]


def ref(expr):
    p = classify._pdal().Pipeline(json.dumps({"pipeline": [
        {"type": "readers.las", "filename": GND, "nosrs": True},
        {"type": "filters.expression", "expression": expr}]}))
    p.execute()
    a = p.arrays[0]
    x, y, z, _ = dedupe_xy(a["X"].astype("f8"), a["Y"].astype("f8"),
                           a["Z"].astype("f8"), keep="lowest")
    return x, y, z


REFS = {
    "class 2 only": "Classification == 2",
    "class 2 + 8": "Classification == 2 || Classification == 8",
}
MODELS = {}
for tag, path in (("smrf", SMRF_G), ("mcc", MCC_G)):
    if os.path.exists(path):
        g = np.load(path)
        MODELS[tag] = g
        print("%s ground set: %d points" % (tag, len(g)), flush=True)

rows = []
for rname, expr in REFS.items():
    tx, ty, tz = ref(expr)
    tcls = cls_at(tx, ty)
    print("")
    print("reference %s: %d points, " % (rname, len(tx))
          + ", ".join("%s %.2f%%" % (NM[k], 100 * float((tcls == k).mean()))
                      for k in (0, 1, 2, 3)), flush=True)
    for tag, g in MODELS.items():
        tr = cKDTree(np.column_stack([g[:, 0], g[:, 1]]))
        d, _ = tr.query(np.column_stack([tx, ty]), workers=-1)
        rj = d > NEAR
        rec = {"reference": rname, "model": tag, "scope": "full tile",
               "n_ref": int(len(tx)), "overall": float(rj.mean())}
        for kk in (0, 1, 2, 3):
            m = tcls == kk
            rec[NM[kk]] = float(rj[m].mean()) if m.sum() > 100 else None
        rows.append(rec)
        print("  %-5s overall %.2f%%, " % (tag, 100 * rec["overall"])
              + ", ".join("%s %.2f%%" % (NM[kk], 100 * rec[NM[kk]])
                          for kk in (0, 1, 2, 3)
                          if rec[NM[kk]] is not None), flush=True)

tx, ty, tz = ref(REFS["class 2 + 8"])
tcls = cls_at(tx, ty)
print("")
print("THE TWO 600 m CORES, matching the TIN prototype geometry")
for nm, lo, la in FORENSIC:
    xs, ys = warp_transform("EPSG:4326", "EPSG:32616", [lo], [la])
    cx, cy = xs[0], ys[0]
    mt = ((tx >= cx - CORE_M / 2) & (tx <= cx + CORE_M / 2)
          & (ty >= cy - CORE_M / 2) & (ty <= cy + CORE_M / 2))
    sub = tcls[mt]
    print("  core %s: %d delivered ground points" % (nm, int(mt.sum())))
    for tag, g in MODELS.items():
        k = ((g[:, 0] >= cx - CORE_M) & (g[:, 0] <= cx + CORE_M)
             & (g[:, 1] >= cy - CORE_M) & (g[:, 1] <= cy + CORE_M))
        tr = cKDTree(np.column_stack([g[k, 0], g[k, 1]]))
        d, _ = tr.query(np.column_stack([tx[mt], ty[mt]]), workers=-1)
        rj = d > NEAR
        rec = {"reference": "class 2 + 8", "model": tag,
               "scope": "core %s" % nm, "n_ref": int(mt.sum()),
               "overall": float(rj.mean())}
        for kk in (0, 1, 2, 3):
            m = sub == kk
            rec[NM[kk]] = float(rj[m].mean()) if m.sum() > 100 else None
        rows.append(rec)
        print("    %-5s overall %.2f%%, " % (tag, 100 * rec["overall"])
              + ", ".join("%s %.2f%%" % (NM[kk], 100 * rec[NM[kk]])
                          for kk in (0, 1, 2, 3)
                          if rec[NM[kk]] is not None), flush=True)

with open(os.path.join(OUT, "rejection_union.json"), "w") as fh:
    json.dump(rows, fh, indent=2)
print("")
print("wrote %s" % os.path.join(OUT, "rejection_union.json"))
