"""The delivered ground population through a TIN surface, full tile.

Replication test for the surface stage alone. The input is the reference's own
ground set -- class 2 and class 8 together, 17,705,532 points -- so the
classifier is held at the reference and the only thing varying is the surface
engine. The archive was produced by a TIN interpolator through these points,
so this is the paradigm-matched build.

`terrascan_tin.tin_grid` is Delaunay plus Clough-Tocher with no tuning
parameters. It is run in the same 300-cell blocks with a 40-cell buffer as the
biharmonic ceiling run, through the same support mask and edge strip, scored
by the same metrics, so the two rows differ in the interpolator and nothing
else.

Clough-Tocher returns NaN outside the convex hull of each block's points
rather than a fallback estimate, so hull misses are counted and reported
instead of being hidden by an extrapolation.
"""
import json
import os
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from scipy import ndimage

sys.path.insert(0, r"C:\Replicalm\src")
from replicalm import classify, finalise
from replicalm import terrascan_tin as TIN
from replicalm.grid import Grid
from replicalm.preprocess import dedupe_xy

GND = (r"D:\_LAZ\Yuc_Centro\Centro_Glas\Centro_Glas"
       r"\Centro_GLAS_l0s441_GND_v1.laz")
ARCH = (r"D:\_Archive_EdgeFixed\Yuc_Centro\Centro_Glas\Centro_Glas"
        r"\Centro_GLAS_l0s441_DEM_0p5m_v1.tif")
PREV = r"C:\Replicalm\test100526\cellsize_validation"
MCCDIR = r"C:\Replicalm\test100526\mcc_fulltile"
CEILDIR = r"C:\Replicalm\test100526\terrascan_ceiling"
OUTDIR = r"C:\Replicalm\test100526\terrascan_tin_surface"
CORE, BBUF = 300, 40
PROM_T = 0.25
os.makedirs(OUTDIR, exist_ok=True)

with rasterio.open(ARCH) as r:
    A = r.read(1).astype("f4")
    at, ac = r.transform, abs(r.transform.a)
    HT, WT = r.height, r.width
    wkt = r.crs.to_wkt()
    nd = r.nodata
if nd is not None:
    A[A == nd] = np.nan
A[A > 1e30] = np.nan
A[A == 0] = np.nan
Am = np.isfinite(A)
print("tile %d x %d, archive extent %.2f%%" % (HT, WT, 100 * Am.mean()),
      flush=True)


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

p = classify._pdal().Pipeline(json.dumps({"pipeline": [
    {"type": "readers.las", "filename": GND, "nosrs": True},
    {"type": "filters.expression",
     "expression": "Classification == 2 || Classification == 8"}]}))
p.execute()
a = p.arrays[0]
cl = a["Classification"].astype(np.int32)
G = np.column_stack([a["X"].astype("f8"), a["Y"].astype("f8"),
                     a["Z"].astype("f8")])
print("delivered ground: %d points, class 2 = %d, class 8 = %d"
      % (len(G), int((cl == 2).sum()), int((cl == 8).sum())), flush=True)
del a, p, cl


def build(g, tag):
    f = os.path.join(OUTDIR, "dem_%s.tif" % tag)
    if os.path.exists(f):
        with rasterio.open(f) as r:
            z = r.read(1).astype("f4")
        z[z == 0.0] = np.nan
        return z, {}
    gx, gy, gz, _r0 = dedupe_xy(g[:, 0], g[:, 1], g[:, 2], keep="lowest")
    print("  dedupe: %d -> %d points" % (len(g), len(gz)), flush=True)
    out = np.full((HT, WT), -9999.0, dtype="f4")
    NBI = (HT + CORE - 1) // CORE
    NBJ = (WT + CORE - 1) // CORE
    bi = np.clip(np.floor((at.f - gy) / ac).astype(np.int64) // CORE, 0,
                 NBI - 1)
    bj = np.clip(np.floor((gx - at.c) / ac).astype(np.int64) // CORE, 0,
                 NBJ - 1)
    bid = bi * NBJ + bj
    o = np.argsort(bid, kind="stable")
    gx, gy, gz, bid = gx[o], gy[o], gz[o], bid[o]
    starts = np.searchsorted(bid, np.arange(NBI * NBJ + 1))
    t0 = time.time()
    tally = {"blocks": 0, "failed": 0, "hull_miss": 0, "cells": 0,
             "triangles": 0}
    for i in range(NBI):
        for j in range(NBJ):
            r0, c0 = i * CORE, j * CORE
            h = min(CORE, HT - r0)
            w = min(CORE, WT - c0)
            sel = []
            for ii in range(max(0, i - 1), min(NBI, i + 2)):
                for jj in range(max(0, j - 1), min(NBJ, j + 2)):
                    k = ii * NBJ + jj
                    if starts[k + 1] > starts[k]:
                        sel.append(np.arange(starts[k], starts[k + 1]))
            if not sel:
                continue
            idx = np.concatenate(sel)
            ox, oy = at.c + c0 * ac, at.f - r0 * ac
            m = ((gx[idx] >= ox - BBUF * ac)
                 & (gx[idx] <= ox + (w + BBUF) * ac)
                 & (gy[idx] <= oy + BBUF * ac)
                 & (gy[idx] >= oy - (h + BBUF) * ac))
            if m.sum() < 20:
                continue
            k = idx[m]
            gr = Grid(origin_x=ox, origin_y=oy, cell=ac, width=w, height=h)
            tally["blocks"] += 1
            try:
                z, info = TIN.tin_grid(gx[k], gy[k], gz[k], gr, verbose=False)
                zz = np.where(np.isfinite(z), z, -9999.0)
                out[r0:r0 + h, c0:c0 + w] = zz.astype("f4")
                tally["hull_miss"] += int(info["outside_hull"])
                tally["cells"] += int(z.size)
                tally["triangles"] += int(info["triangles"])
            except Exception:
                tally["failed"] += 1
        if i % 10 == 0:
            print("    %s dem row %d/%d, %.1f min"
                  % (tag, i, NBI, (time.time() - t0) / 60.0), flush=True)
    print("  blocks %d, failed %d, cells outside a hull %d of %d (%.3f%%), "
          "%d triangles" % (tally["blocks"], tally["failed"],
                            tally["hull_miss"], tally["cells"],
                            100.0 * tally["hull_miss"]
                            / max(tally["cells"], 1), tally["triangles"]),
          flush=True)
    tile = Grid(origin_x=at.c, origin_y=at.f, cell=ac, width=WT, height=HT)
    sup = finalise.support_mask(gx, gy, tile, 4.0)
    rep = finalise.finalise(out, tile, wkt, f, radius_m=4.0, nodata_in=-9999.0,
                            support=sup, strip_m=21.0, axis="x", verbose=False)
    print("  %s finalised: %d cells, %d interior restored"
          % (tag, rep["cells_after_trim"], rep["interior_cells_restored"]),
          flush=True)
    with rasterio.open(f) as r:
        z = r.read(1).astype("f4")
    z[z == 0.0] = np.nan
    return z, tally


def detail(z, mask, scale):
    s = scale / ac / 2.0
    num = ndimage.gaussian_filter(np.where(mask, z, 0.0).astype("f4"), s)
    den = ndimage.gaussian_filter(mask.astype("f4"), s)
    sm = np.where(den > 1e-6, num / np.maximum(den, 1e-6), 0.0)
    return np.where(mask, z - sm, np.nan)


def metrics(z, tag):
    ok = Am & np.isfinite(z)
    e = (z - A)[ok]
    inner = ndimage.binary_erosion(ok, structure=np.ones((3, 3)),
                                   iterations=int(5.0 / ac))
    k = max(3, int(round(1.0 / ac)) | 1)
    dB = detail(z, ok, 1.5)
    dA = detail(A, ok, 1.5)
    dBm = ndimage.median_filter(np.nan_to_num(dB, nan=0.0), size=k)
    dAm = ndimage.median_filter(np.nan_to_num(dA, nan=0.0), size=k)
    sel = inner & np.isfinite(dB) & np.isfinite(dA)
    zz = np.nan_to_num(z, nan=0.0)
    msk = ok[:, 1:] & ok[:, :-1]
    nug = float(0.5 * np.mean((zz[:, 1:][msk] - zz[:, :-1][msk]) ** 2))
    by, bias = {}, {}
    for kk in (0, 1, 2, 3):
        s2 = ok & (CLS == kk)
        if s2.sum() > 100:
            by[NM[kk]] = float(np.sqrt(np.mean((z - A)[s2] ** 2)))
            bias[NM[kk]] = float(np.mean((z - A)[s2]))
        else:
            by[NM[kk]] = None
            bias[NM[kk]] = None
    return dict(tag=tag, cells=int(ok.sum()),
                rmse=float(np.sqrt(np.mean(e ** 2))),
                mae=float(np.mean(np.abs(e))),
                bias=float(np.mean(e)),
                p99abs=float(np.percentile(np.abs(e), 99)),
                coherent=float(np.sqrt(np.mean(dBm[sel] ** 2))),
                speckle=float(np.sqrt(np.mean((dB - dBm)[sel] ** 2))),
                rough=float(np.sqrt(np.mean(dB[sel] ** 2))
                            / max(np.sqrt(np.mean(dA[sel] ** 2)), 1e-9)),
                arch_coherent=float(np.sqrt(np.mean(dAm[sel] ** 2))),
                arch_speckle=float(np.sqrt(np.mean((dA - dAm)[sel] ** 2))),
                nugget=nug, rmse_by_form=by, bias_by_form=bias)


Z, tally = build(G, "tin_c2c8")
res = [metrics(Z, "tin_c2c8")]
for tag, path in (("biharm_c2c8", os.path.join(CEILDIR, "dem_ts_c2c8.tif")),
                  ("smrf", os.path.join(PREV, "dem_B_cell1.00.tif")),
                  ("mcc", os.path.join(MCCDIR, "dem_mcc.tif"))):
    if not os.path.exists(path):
        continue
    with rasterio.open(path) as r:
        z = r.read(1).astype("f4")
    z[z == 0.0] = np.nan
    res.append(metrics(z, tag))

print("")
print("%-13s %11s %9s %8s %9s %9s %9s %9s %7s %10s"
      % ("run", "cells", "RMSE", "MAE", "bias", "p99|e|", "coherent",
         "speckle", "rough", "nugget"))
print("-" * 116)
for r in res:
    print("%-13s %11d %9.4f %8.4f %9.4f %9.4f %9.4f %9.4f %7.3f %10.6f"
          % (r["tag"], r["cells"], r["rmse"], r["mae"], r["bias"],
             r["p99abs"], r["coherent"], r["speckle"], r["rough"],
             r["nugget"]))
print("archive own texture: coherent %.4f, speckle %.4f"
      % (res[0]["arch_coherent"], res[0]["arch_speckle"]))

SHARE = {"planar": 0.8879, "mound": 0.0507, "depression": 0.0478,
         "shoulder": 0.0136}
print("")
print("RMSE BY LANDFORM")
print("%-13s %10s %10s %12s %10s" % ("run", "planar", "mound", "depression",
                                     "shoulder"))
for r in res:
    b = r["rmse_by_form"]
    print("%-13s " % r["tag"] + " ".join(
        "%10s" % ("%.4f" % b[n] if b[n] is not None else "-")
        for n in ("planar", "mound", "depression", "shoulder")))
print("")
print("MEAN SIGNED ERROR BY LANDFORM")
print("%-13s %10s %10s %12s %10s" % ("run", "planar", "mound", "depression",
                                     "shoulder"))
for r in res:
    b = r["bias_by_form"]
    print("%-13s " % r["tag"] + " ".join(
        "%10s" % ("%+.4f" % b[n] if b[n] is not None else "-")
        for n in ("planar", "mound", "depression", "shoulder")))
print("")
print("SHARE OF ERROR VARIANCE BY LANDFORM, area-weighted")
print("%-13s %10s %10s %12s %10s %12s"
      % ("run", "planar", "mound", "depression", "shoulder", "non-planar"))
for r in res:
    b = r["rmse_by_form"]
    v = {n: SHARE[n] * b[n] ** 2 for n in SHARE if b[n] is not None}
    t = sum(v.values())
    print("%-13s " % r["tag"] + " ".join(
        "%10s" % ("%.1f%%" % (100 * v[n] / t) if n in v else "-")
        for n in ("planar", "mound", "depression", "shoulder"))
        + "%12s" % ("%.1f%%" % (100 * sum(
            v[n] for n in ("mound", "depression", "shoulder") if n in v) / t)))

SH = 8


def shrink(x, k=SH):
    h, w = x.shape
    x = x[:h // k * k, :w // k * k]
    with np.errstate(invalid="ignore"):
        return np.nanmean(x.reshape(h // k, k, w // k, k), axis=(1, 3))


def hs(z):
    zf = np.nan_to_num(z, nan=float(np.nanmedian(z)))
    gy3, gx3 = np.gradient(zf, ac)
    s = np.arctan(np.hypot(gx3, gy3))
    al = np.arctan2(-gx3, gy3)
    AA, E = np.radians(360 - 315 + 90), np.radians(45)
    o = np.clip(np.sin(E) * np.cos(s) + np.cos(E) * np.sin(s)
                * np.cos(AA - al), 0, 1)
    return np.where(np.isfinite(z), o, np.nan)


with rasterio.open(os.path.join(CEILDIR, "dem_ts_c2c8.tif")) as r:
    ZB = r.read(1).astype("f4")
ZB[ZB == 0.0] = np.nan
d1, d2 = A - Z, A - ZB
v = float(np.nanpercentile(np.abs(d2[Am]), 98))
panels = [("archive", shrink(hs(A)), "gray", 0, 1),
          ("TIN, class 2 + 8", shrink(hs(Z)), "gray", 0, 1),
          ("biharmonic, class 2 + 8", shrink(hs(ZB)), "gray", 0, 1),
          ("archive - TIN, \u00b1%.2f m" % v, shrink(d1), "RdBu_r", -v, v),
          ("archive - biharmonic, \u00b1%.2f m" % v, shrink(d2), "RdBu_r",
           -v, v)]
for cm in ("gray", "RdBu_r"):
    matplotlib.colormaps[cm].set_bad("white")
fig, axes = plt.subplots(1, 5, figsize=(21.0, 17.0))
for ax, (t, img, cm, lo, hi) in zip(axes, panels):
    ax.imshow(img, cmap=cm, vmin=lo, vmax=hi, interpolation="nearest")
    ax.set_title(t, fontsize=11)
    ax.set_xticks([])
    ax.set_yticks([])
fig.suptitle("Centro GLAS l0s441 \u2014 the delivered ground population "
             "through a TIN surface against the biharmonic surface",
             fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.97])
pp = os.path.join(OUTDIR, "ts_tin_surface.png")
fig.savefig(pp, dpi=200)
print("")
print("wrote %s" % pp, flush=True)
with open(os.path.join(OUTDIR, "metrics_tin.json"), "w") as fh:
    json.dump({"metrics": res, "tally": tally}, fh, indent=2)
print("wrote %s" % os.path.join(OUTDIR, "metrics_tin.json"), flush=True)
