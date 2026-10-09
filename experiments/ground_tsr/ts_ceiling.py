"""Full-tile DEM from the delivered ground points, class 2 and class 8 together.

This sets the ceiling. The input is the reference's own ground set, so whatever
residual remains against the archive is attributable to the surface engine and
the grid handling, not to classification. Comparing it against the SMRF row
(0.0789 m) and the MCC row (0.0899 m) from the same harness decomposes the
error into a classification part and an interpolation part.

The harness is held identical to those two runs: the same 300-cell blocks with
a 40-cell buffer, the same global biharmonic engine at pool 0.25, the same 4 m
support mask and 21 m across-track edge strip, the same landform classes from
the archive's own prominence, the same metrics.

Both delivered classes are used. Class 2 alone sits a median 16 mm below the
archive surface and class 8 a median 12 mm above it, while the union is
centred to 3 mm; class 2 alone also occupies only 35.75 percent of the valid
archive cells against the union's 48.32 percent. Restricting the input to
class 2 would therefore measure a different surface than the one delivered.

`dedupe_xy` runs with round_to None, so it removes only bit-identical repeats
rather than thinning to one point per cell. The pooling inside the engine is a
Gaussian-weighted mean of the constraints in a cell, so neither population is
preferentially taken.
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
from replicalm import surfer_algoA as ALGOA
from replicalm.grid import Grid
from replicalm.preprocess import dedupe_xy

GND = (r"D:\_LAZ\Yuc_Centro\Centro_Glas\Centro_Glas"
       r"\Centro_GLAS_l0s441_GND_v1.laz")
ARCH = (r"D:\_Archive_EdgeFixed\Yuc_Centro\Centro_Glas\Centro_Glas"
        r"\Centro_GLAS_l0s441_DEM_0p5m_v1.tif")
PREV = r"C:\Replicalm\test100526\cellsize_validation"
MCCDIR = r"C:\Replicalm\test100526\mcc_fulltile"
OUTDIR = r"C:\Replicalm\test100526\terrascan_ceiling"
CORE, BBUF, POOL = 300, 40, 0.25
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
        return z
    gx, gy, gz, rep0 = dedupe_xy(g[:, 0], g[:, 1], g[:, 2], keep="lowest")
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
            try:
                z, _i = ALGOA.surfer_grid(gx[k], gy[k], gz[k], gr,
                                          pool_sigma_m=POOL, verbose=False)
                out[r0:r0 + h, c0:c0 + w] = z.astype("f4")
            except Exception:
                pass
        if i % 10 == 0:
            print("    %s dem row %d/%d, %.1f min"
                  % (tag, i, NBI, (time.time() - t0) / 60.0), flush=True)
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
    return z


def detail(z, mask, scale):
    s = scale / ac / 2.0
    num = ndimage.gaussian_filter(np.where(mask, z, 0.0).astype("f4"), s)
    den = ndimage.gaussian_filter(mask.astype("f4"), s)
    sm = np.where(den > 1e-6, num / np.maximum(den, 1e-6), 0.0)
    return np.where(mask, z - sm, np.nan)


def metrics(z, n_ground, tag):
    ok = Am & np.isfinite(z)
    e = (z - A)[ok]
    zf = np.nan_to_num(z, nan=float(np.nanmean(z)))
    gy2, gx2 = np.gradient(zf, ac)
    sl = np.degrees(np.arctan(np.hypot(gx2, gy2)))[ok]
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
    return dict(tag=tag, n_ground=int(n_ground), cells=int(ok.sum()),
                rmse=float(np.sqrt(np.mean(e ** 2))),
                mae=float(np.mean(np.abs(e))),
                bias=float(np.mean(e)),
                p95=float(np.percentile(sl, 95)),
                p99=float(np.percentile(sl, 99)),
                coherent=float(np.sqrt(np.mean(dBm[sel] ** 2))),
                speckle=float(np.sqrt(np.mean((dB - dBm)[sel] ** 2))),
                rough=float(np.sqrt(np.mean(dB[sel] ** 2))
                            / max(np.sqrt(np.mean(dA[sel] ** 2)), 1e-9)),
                arch_coherent=float(np.sqrt(np.mean(dAm[sel] ** 2))),
                arch_speckle=float(np.sqrt(np.mean((dA - dAm)[sel] ** 2))),
                nugget=nug, rmse_by_form=by, bias_by_form=bias)


Z = build(G, "ts_c2c8")
res = [metrics(Z, len(G), "ts_c2c8")]

for tag, path in (("smrf", os.path.join(PREV, "dem_B_cell1.00.tif")),
                  ("mcc", os.path.join(MCCDIR, "dem_mcc.tif"))):
    if not os.path.exists(path):
        continue
    with rasterio.open(path) as r:
        z = r.read(1).astype("f4")
    z[z == 0.0] = np.nan
    res.append(metrics(z, 0, tag))

print("")
print("%-9s %12s %11s %9s %8s %9s %9s %9s %7s %10s"
      % ("run", "ground pts", "cells", "RMSE", "MAE", "bias", "coherent",
         "speckle", "rough", "nugget"))
print("-" * 112)
for r in res:
    print("%-9s %12d %11d %9.4f %8.4f %9.4f %9.4f %9.4f %7.3f %10.6f"
          % (r["tag"], r["n_ground"], r["cells"], r["rmse"], r["mae"],
             r["bias"], r["coherent"], r["speckle"], r["rough"], r["nugget"]))
print("archive own texture: coherent %.4f, speckle %.4f"
      % (res[0]["arch_coherent"], res[0]["arch_speckle"]))

print("")
print("RMSE BY LANDFORM")
print("%-9s %10s %10s %12s %10s" % ("run", "planar", "mound", "depression",
                                    "shoulder"))
for r in res:
    b = r["rmse_by_form"]
    print("%-9s " % r["tag"] + " ".join(
        "%10s" % ("%.4f" % b[n] if b[n] is not None else "-")
        for n in ("planar", "mound", "depression", "shoulder")))
print("")
print("MEAN SIGNED ERROR BY LANDFORM, surface minus archive")
print("%-9s %10s %10s %12s %10s" % ("run", "planar", "mound", "depression",
                                    "shoulder"))
for r in res:
    b = r["bias_by_form"]
    print("%-9s " % r["tag"] + " ".join(
        "%10s" % ("%+.4f" % b[n] if b[n] is not None else "-")
        for n in ("planar", "mound", "depression", "shoulder")))

if len(res) >= 2:
    c = res[0]
    print("")
    for other in res[1:]:
        print("%s RMSE %.4f, ceiling %.4f -> interpolation accounts for "
              "%.1f%% of it" % (other["tag"], other["rmse"], c["rmse"],
                                100 * c["rmse"] / other["rmse"]))

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


d = A - Z
v = float(np.nanpercentile(np.abs(d[Am]), 98))
panels = [("archive", shrink(hs(A)), "gray", 0, 1),
          ("delivered class 2 + 8, biharmonic", shrink(hs(Z)), "gray", 0, 1),
          ("archive - this surface, \u00b1%.2f m" % v, shrink(d), "RdBu_r",
           -v, v)]
for cm in ("gray", "RdBu_r"):
    matplotlib.colormaps[cm].set_bad("white")
fig, axes = plt.subplots(1, 3, figsize=(13.0, 17.0))
for ax, (t, img, cm, lo, hi) in zip(axes, panels):
    ax.imshow(img, cmap=cm, vmin=lo, vmax=hi, interpolation="nearest")
    ax.set_title(t, fontsize=11)
    ax.set_xticks([])
    ax.set_yticks([])
fig.suptitle("Centro GLAS l0s441 \u2014 the delivered ground set through the "
             "Replicalm surface engine", fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.97])
pp = os.path.join(OUTDIR, "ts_ceiling.png")
fig.savefig(pp, dpi=200)
print("")
print("wrote %s" % pp, flush=True)
with open(os.path.join(OUTDIR, "metrics_ceiling.json"), "w") as fh:
    json.dump(res, fh, indent=2)
print("wrote %s" % os.path.join(OUTDIR, "metrics_ceiling.json"), flush=True)
