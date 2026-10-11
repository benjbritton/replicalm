"""Full tile through ground_tsr at its calibrated defaults.

Block means and a full-tile RMSE are not the same measurement, so the held-out
block figure is not the deliverable. This runs the whole of l0s441 from the raw
flightline, with the delivered labels nowhere in the execution path, through
the same 400 m tiling, the same biharmonic engine, the same support mask and
edge strip as every other full-tile run in this series.

Scored against two surfaces. The reference is the delivered class 2 + 8
population through the same biharmonic engine, which isolates classification:
an exact reproduction would read zero. The archive is the scale the 0.052 m
benchmark was quoted on. Both are reported on the cells every surface filled,
because the comparison is otherwise between different cell populations.
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
from replicalm import ground_tsr as TSR
from replicalm import surfer_algoA as ALGOA
from replicalm.grid import Grid
from replicalm.preprocess import dedupe_xy

RAW = (r"D:\GLiHT_LAS_orig\Yuc_Centro\Centro_Glas\Centro_Glas"
       r"\AMIGACarb_Yuc_Centro_GLAS_Apr2013_l0s441.las")
REF = r"C:\Replicalm\test100526\terrascan_ceiling\dem_ts_c2c8.tif"
ARCH = (r"D:\_Archive_EdgeFixed\Yuc_Centro\Centro_Glas\Centro_Glas"
        r"\Centro_GLAS_l0s441_DEM_0p5m_v1.tif")
SMRF_DEM = (r"C:\Replicalm\test100526\cellsize_validation"
            r"\dem_B_cell1.00.tif")
OUT = r"C:\Replicalm\test100526\ground_replication"
TILE_M, TBUF_M = 400.0, 40.0
CORE, BBUF, POOL = 300, 40, 0.25
PROM_T = 0.25
os.makedirs(OUT, exist_ok=True)


def read(path):
    with rasterio.open(path) as r:
        z = r.read(1).astype("f4")
        nd = r.nodata
        tr, cc, wkt = r.transform, abs(r.transform.a), r.crs.to_wkt()
        h, w = r.height, r.width
    if nd is not None:
        z[z == nd] = np.nan
    z[z > 1e30] = np.nan
    z[z == 0.0] = np.nan
    return z, tr, cc, wkt, h, w


R, at, ac, wkt, HT, WT = read(REF)
A = read(ARCH)[0]
Am = np.isfinite(A)
Rm = np.isfinite(R)
print("tile %d x %d, archive %.2f%%, reference %.2f%%"
      % (HT, WT, 100 * Am.mean(), 100 * Rm.mean()), flush=True)


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

x0, x1 = at.c, at.c + WT * ac
y1, y0 = at.f, at.f - HT * ac
GPATH = os.path.join(OUT, "ground_tsr.npy")

if os.path.exists(GPATH):
    G = np.load(GPATH)
    print("cached ground: %d returns" % len(G), flush=True)
else:
    t0 = time.time()
    pl = classify._pdal().Pipeline(json.dumps({"pipeline": [
        {"type": "readers.las", "filename": RAW}]}))
    pl.execute()
    a = pl.arrays[0]
    RX = a["X"].astype("f8")
    RY = a["Y"].astype("f8")
    RZ = a["Z"].astype("f8")
    del a, pl
    print("read %d raw returns in %.0f s" % (len(RZ), time.time() - t0),
          flush=True)
    params = TSR.GroundTSRParams()
    out, reps = [], []
    nx = int(np.ceil((x1 - x0) / TILE_M))
    ny = int(np.ceil((y1 - y0) / TILE_M))
    t0 = time.time()
    for i in range(ny):
        for j in range(nx):
            cy1 = y1 - i * TILE_M
            cy0 = max(y0, cy1 - TILE_M)
            cx0 = x0 + j * TILE_M
            cx1 = min(x1, cx0 + TILE_M)
            m = ((RX >= cx0 - TBUF_M) & (RX <= cx1 + TBUF_M)
                 & (RY >= cy0 - TBUF_M) & (RY <= cy1 + TBUF_M))
            if m.sum() < 2000:
                continue
            bx, by, bz = RX[m], RY[m], RZ[m]
            try:
                mask, rep = TSR.classify_ground(bx, by, bz, params=params)
            except TSR.GroundTSRError as exc:
                print("    tile %d,%d skipped: %s" % (i, j, exc), flush=True)
                continue
            reps.append(rep)
            gx, gy, gz = bx[mask], by[mask], bz[mask]
            k = (gx >= cx0) & (gx < cx1) & (gy >= cy0) & (gy < cy1)
            if k.any():
                out.append(np.column_stack([gx[k], gy[k], gz[k]]))
        print("  row %d/%d, %.1f min" % (i + 1, ny,
                                         (time.time() - t0) / 60.0),
              flush=True)
    G = np.vstack(out)
    np.save(GPATH, G)
    dd = [r["density_per_m2"] for r in reps]
    cc2 = [r["smrf_cell_m"] for r in reps]
    print("ground_tsr: %d returns in %.1f min over %d tiles"
          % (len(G), (time.time() - t0) / 60.0, len(reps)), flush=True)
    print("  density across tiles: min %.1f, median %.1f, max %.1f /m2"
          % (np.min(dd), np.median(dd), np.max(dd)))
    print("  working cell chosen: min %.2f, median %.2f, max %.2f m"
          % (np.min(cc2), np.median(cc2), np.max(cc2)))
    with open(os.path.join(OUT, "tsr_tile_reports.json"), "w") as fh:
        json.dump(reps, fh, indent=2)
    del RX, RY, RZ

DEM = os.path.join(OUT, "dem_tsr.tif")
if os.path.exists(DEM):
    Z = read(DEM)[0]
else:
    gx, gy, gz, _ = dedupe_xy(G[:, 0], G[:, 1], G[:, 2], keep="lowest")
    print("dedupe: %d -> %d" % (len(G), len(gz)), flush=True)
    outg = np.full((HT, WT), -9999.0, dtype="f4")
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
                outg[r0:r0 + h, c0:c0 + w] = z.astype("f4")
            except Exception:
                pass
        if i % 10 == 0:
            print("    dem row %d/%d, %.1f min"
                  % (i, NBI, (time.time() - t0) / 60.0), flush=True)
    tile = Grid(origin_x=at.c, origin_y=at.f, cell=ac, width=WT, height=HT)
    sup = finalise.support_mask(gx, gy, tile, 4.0)
    rp = finalise.finalise(outg, tile, wkt, DEM, radius_m=4.0,
                           nodata_in=-9999.0, support=sup, strip_m=21.0,
                           axis="x", verbose=False)
    print("  finalised: %d cells, %d interior restored"
          % (rp["cells_after_trim"], rp["interior_cells_restored"]),
          flush=True)
    Z = read(DEM)[0]

S = read(SMRF_DEM)[0]
common = Am & Rm & np.isfinite(Z) & np.isfinite(S)
print("")
print("common cells: %d (%.2f%% of the reference)"
      % (int(common.sum()), 100 * common.sum() / Rm.sum()))
SHARE = {NM[k]: float((CLS[common] == k).mean()) for k in (0, 1, 2, 3)}


def detail(z, mask, scale):
    s = scale / ac / 2.0
    num = ndimage.gaussian_filter(np.where(mask, z, 0.0).astype("f4"), s)
    den = ndimage.gaussian_filter(mask.astype("f4"), s)
    sm = np.where(den > 1e-6, num / np.maximum(den, 1e-6), 0.0)
    return np.where(mask, z - sm, np.nan)


inner = ndimage.binary_erosion(common, structure=np.ones((3, 3)),
                               iterations=int(5.0 / ac))
kk = max(3, int(round(1.0 / ac)) | 1)
dA = detail(A, common, 1.5)
dAm = ndimage.median_filter(np.nan_to_num(dA, nan=0.0), size=kk)
sel0 = inner & np.isfinite(dA)

rows = []
for tag, z in (("ground_tsr", Z), ("smrf baseline", S), ("reference", R)):
    e = (z - R)[common]
    ea = (z - A)[common]
    dB = detail(z, common, 1.5)
    dBm = ndimage.median_filter(np.nan_to_num(dB, nan=0.0), size=kk)
    sel = sel0 & np.isfinite(dB)
    by = {NM[k]: float(np.sqrt(np.mean((z - A)[common & (CLS == k)] ** 2)))
          for k in (0, 1, 2, 3)}
    rows.append(dict(tag=tag,
                     rmse_ref=float(np.sqrt(np.mean(e ** 2))),
                     rmse_arch=float(np.sqrt(np.mean(ea ** 2))),
                     mae_arch=float(np.mean(np.abs(ea))),
                     bias_arch=float(np.mean(ea)),
                     p99_arch=float(np.percentile(np.abs(ea), 99)),
                     coherent=float(np.sqrt(np.mean(dBm[sel] ** 2))),
                     speckle=float(np.sqrt(np.mean((dB - dBm)[sel] ** 2))),
                     rough=float(np.sqrt(np.mean(dB[sel] ** 2))
                                 / max(np.sqrt(np.mean(dA[sel] ** 2)), 1e-9)),
                     rmse_arch_by_form=by))

print("")
print("%-14s %10s %11s %9s %9s %9s %9s %9s %7s"
      % ("run", "RMSE ref", "RMSE arch", "MAE", "bias", "p99|e|", "coherent",
         "speckle", "rough"))
print("-" * 100)
for r in rows:
    print("%-14s %10.4f %11.4f %9.4f %9.4f %9.4f %9.4f %9.4f %7.3f"
          % (r["tag"], r["rmse_ref"], r["rmse_arch"], r["mae_arch"],
             r["bias_arch"], r["p99_arch"], r["coherent"], r["speckle"],
             r["rough"]))
print("archive own texture: coherent %.4f, speckle %.4f"
      % (float(np.sqrt(np.mean(dAm[sel0] ** 2))),
         float(np.sqrt(np.mean((dA - dAm)[sel0] ** 2)))))

print("")
print("RMSE AGAINST THE ARCHIVE BY LANDFORM")
print("%-14s %10s %10s %12s %10s" % ("run", "planar", "mound", "depression",
                                     "shoulder"))
for r in rows:
    b = r["rmse_arch_by_form"]
    print("%-14s " % r["tag"] + " ".join(
        "%10.4f" % b[n] for n in ("planar", "mound", "depression",
                                  "shoulder")))

bench = [r for r in rows if r["tag"] == "reference"][0]["rmse_arch"]
mine = [r for r in rows if r["tag"] == "ground_tsr"][0]["rmse_arch"]
base = [r for r in rows if r["tag"] == "smrf baseline"][0]["rmse_arch"]
print("")
print("benchmark (delivered union through this engine) %.4f m" % bench)
print("ground_tsr from raw returns only               %.4f m  (%+.1f%% "
      "of benchmark)" % (mine, 100 * (mine / bench - 1)))
print("smrf baseline                                  %.4f m  (%+.1f%%)"
      % (base, 100 * (base / bench - 1)))
print("ground_tsr against the baseline: %+.1f%%" % (100 * (mine / base - 1)))

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


d1, d2 = A - Z, A - S
v = float(np.nanpercentile(np.abs(d2[common]), 98))
panels = [("archive", shrink(hs(A)), "gray", 0, 1),
          ("ground_tsr, raw input only", shrink(hs(Z)), "gray", 0, 1),
          ("smrf baseline", shrink(hs(S)), "gray", 0, 1),
          ("archive - ground_tsr, \u00b1%.2f m" % v, shrink(d1), "RdBu_r",
           -v, v),
          ("archive - smrf, \u00b1%.2f m" % v, shrink(d2), "RdBu_r", -v, v)]
for cm in ("gray", "RdBu_r"):
    matplotlib.colormaps[cm].set_bad("white")
fig, axes = plt.subplots(1, 5, figsize=(21.0, 17.0))
for ax, (t, img, cm, lo, hi) in zip(axes, panels):
    ax.imshow(img, cmap=cm, vmin=lo, vmax=hi, interpolation="nearest")
    ax.set_title(t, fontsize=11)
    ax.set_xticks([])
    ax.set_yticks([])
fig.suptitle("Centro GLAS l0s441 \u2014 ground_tsr from raw returns against "
             "the locked baseline", fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.97])
pp = os.path.join(OUT, "tsr_fulltile.png")
fig.savefig(pp, dpi=200)
print("")
print("wrote %s" % pp, flush=True)
with open(os.path.join(OUT, "tsr_fulltile_metrics.json"), "w") as fh:
    json.dump({"cells": int(common.sum()), "shares": SHARE, "runs": rows},
              fh, indent=2)
print("wrote %s" % os.path.join(OUT, "tsr_fulltile_metrics.json"))
