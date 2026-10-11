"""Full-tile MCC against the SMRF baseline on l0s441, cell 1.0 m.

Everything downstream of classification is held fixed: the same original LAS,
the same statistical outlier filter, the same 400 m tiling with a 40 m buffer,
the same global biharmonic interpolation at pool 0.25, the same 4 m support
mask and 21 m across-track edge strip. The SMRF side is the cached ground set
from the cell-size validation, so the two runs differ in the classifier and
nothing else.

One asymmetry from the block comparison is removed here. In that run SMRF saw
the outlier filter and MCC saw raw returns, because MCC was called on the
point array directly. Here both get the denoised input, which is what the
pipeline actually does.

MCC is implemented locally. Evans & Hudak use a thin-plate spline as the
curvature reference and this uses a Gaussian-smoothed mean surface at the same
three scales. It is a close relative of their method, not their method, and
the tolerance and scale schedule are the published defaults rather than
anything tuned here.
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
from scipy.spatial import cKDTree

sys.path.insert(0, r"C:\Replicalm\src")
from replicalm import classify, finalise
from replicalm import surfer_algoA as ALGOA
from replicalm.config import PRESETS
from replicalm.grid import Grid
from replicalm.preprocess import dedupe_xy

LAS = (r"D:\GLiHT_LAS_orig\Yuc_Centro\Centro_Glas\Centro_Glas"
       r"\AMIGACarb_Yuc_Centro_GLAS_Apr2013_l0s441.las")
TS = (r"D:\_LAZ\Yuc_Centro\Centro_Glas\Centro_Glas"
      r"\Centro_GLAS_l0s441_GND_v1.laz")
ARCH = (r"D:\_Archive_EdgeFixed\Yuc_Centro\Centro_Glas\Centro_Glas"
        r"\Centro_GLAS_l0s441_DEM_0p5m_v1.tif")
PREV = r"C:\Replicalm\test100526\cellsize_validation"
OUTDIR = r"C:\Replicalm\test100526\mcc_fulltile"
TILE_M, TBUF_M = 400.0, 40.0
CORE, BBUF, POOL = 300, 40, 0.25
PROM_T, NEAR = 0.25, 0.25
os.makedirs(OUTDIR, exist_ok=True)

with rasterio.open(ARCH) as r:
    A = r.read(1).astype("f4")
    at, ac, crs = r.transform, abs(r.transform.a), r.crs
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

cfg = PRESETS["clear"]
p0 = cfg.passes[0]
x0, x1 = at.c, at.c + WT * ac
y1, y0 = at.f, at.f - HT * ac


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
print("landform cells: " + ", ".join(
    "%s %.2f%%" % (NM[k], 100 * float((CLS[Am] == k).mean()))
    for k in (0, 1, 2, 3)), flush=True)
del sm2, sm5, sm10, curv, slope_r, sharp, prom


def cls_at(x, y):
    c = np.clip(((x - at.c) / ac).astype(np.int64), 0, CLS.shape[1] - 1)
    r = np.clip(((at.f - y) / ac).astype(np.int64), 0, CLS.shape[0] - 1)
    return CLS[r, c]


t0 = time.time()
pl = classify._pdal().Pipeline(json.dumps(
    {"pipeline": [{"type": "readers.las", "filename": LAS}]}))
pl.execute()
arr = pl.arrays[0]
RX = arr["X"].astype("f8")
RY = arr["Y"].astype("f8")
RZ = arr["Z"].astype("f8")
del arr, pl
print("read %d raw returns in %.0f s" % (len(RZ), time.time() - t0),
      flush=True)


def denoise(px, py, pz):
    a = np.zeros(len(pz), dtype=[("X", "f8"), ("Y", "f8"), ("Z", "f8"),
                                 ("Classification", "u1")])
    a["X"], a["Y"], a["Z"] = px, py, pz
    stg = [{"type": "filters.assign", "value": "Classification = 0"},
           {"type": "filters.outlier", "method": "statistical",
            "mean_k": cfg.outlier_neighbours,
            "multiplier": cfg.outlier_multiplier, "class": cfg.noise_class},
           {"type": "filters.expression",
            "expression": "Classification != %d" % cfg.noise_class}]
    p = classify._pdal().Pipeline(json.dumps({"pipeline": stg}), arrays=[a])
    p.execute()
    g = p.arrays[0]
    return g["X"].astype("f8"), g["Y"].astype("f8"), g["Z"].astype("f8")


def mcc(x, y, z, scales=(0.5, 1.0, 1.5), tol=0.30, iters=4):
    """Multiscale curvature tolerance, after Evans & Hudak 2007."""
    keep = np.ones(len(z), bool)
    ox, oy = x.min(), y.max()
    for lam in scales:
        for _ in range(iters):
            w = int((x.max() - ox) / lam) + 2
            h = int((oy - y.min()) / lam) + 2
            c = np.clip(((x - ox) / lam).astype(int), 0, w - 1)
            r = np.clip(((oy - y) / lam).astype(int), 0, h - 1)
            S = np.zeros((h, w))
            N = np.zeros((h, w))
            np.add.at(S, (r[keep], c[keep]), z[keep])
            np.add.at(N, (r[keep], c[keep]), 1.0)
            empty = N < 1
            with np.errstate(invalid="ignore"):
                M = np.where(empty, np.nan, S / np.maximum(N, 1))
            idx = ndimage.distance_transform_edt(
                empty, return_distances=False, return_indices=True)
            M = np.where(empty, M[tuple(idx)], M)
            M = ndimage.gaussian_filter(np.nan_to_num(
                M, nan=float(np.nanmean(M))), 1.0)
            drop = keep & (z - M[r, c] > tol)
            if drop.sum() == 0:
                break
            keep &= ~drop
    return keep


def classify_mcc():
    f = os.path.join(OUTDIR, "ground_mcc.npy")
    if os.path.exists(f):
        g = np.load(f)
        print("mcc cached: %d ground returns" % len(g), flush=True)
        return g
    out = []
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
            if m.sum() < 500:
                continue
            px, py, pz = denoise(RX[m], RY[m], RZ[m])
            if len(pz) < 500:
                continue
            k = mcc(px, py, pz)
            gx, gy, gz = px[k], py[k], pz[k]
            k = (gx >= cx0) & (gx < cx1) & (gy >= cy0) & (gy < cy1)
            if k.any():
                out.append(np.column_stack([gx[k], gy[k], gz[k]]))
        print("  mcc row %d/%d, %.1f min" % (i + 1, ny,
                                             (time.time() - t0) / 60.0),
              flush=True)
    g = np.vstack(out)
    np.save(f, g)
    print("mcc: %d ground returns in %.1f min"
          % (len(g), (time.time() - t0) / 60.0), flush=True)
    return g


def build(g, tag):
    f = os.path.join(OUTDIR, "dem_%s.tif" % tag)
    if os.path.exists(f):
        with rasterio.open(f) as r:
            z = r.read(1).astype("f4")
        z[z == 0.0] = np.nan
        return z
    gx, gy, gz, _d = dedupe_xy(g[:, 0], g[:, 1], g[:, 2], keep="lowest")
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
    by = {}
    for kk in (0, 1, 2, 3):
        s2 = ok & (CLS == kk)
        by[NM[kk]] = (float(np.sqrt(np.mean((z - A)[s2] ** 2)))
                      if s2.sum() > 100 else None)
    return dict(tag=tag, n_ground=int(n_ground), cells=int(ok.sum()),
                rmse=float(np.sqrt(np.mean(e ** 2))),
                mae=float(np.mean(np.abs(e))),
                p95=float(np.percentile(sl, 95)),
                p99=float(np.percentile(sl, 99)),
                coherent=float(np.sqrt(np.mean(dBm[sel] ** 2))),
                speckle=float(np.sqrt(np.mean((dB - dBm)[sel] ** 2))),
                rough=float(np.sqrt(np.mean(dB[sel] ** 2))
                            / max(np.sqrt(np.mean(dA[sel] ** 2)), 1e-9)),
                arch_coherent=float(np.sqrt(np.mean(dAm[sel] ** 2))),
                arch_speckle=float(np.sqrt(np.mean((dA - dAm)[sel] ** 2))),
                nugget=nug, rmse_by_form=by)


def rejection(g, tcls, tx, ty, tag):
    tr = cKDTree(np.column_stack([g[:, 0], g[:, 1]]))
    d, _ = tr.query(np.column_stack([tx, ty]), workers=-1)
    rej = d > NEAR
    out = {"overall": float(rej.mean()), "n_terrascan": int(len(tx))}
    for kk in (0, 1, 2, 3):
        m = tcls == kk
        out[NM[kk]] = float(rej[m].mean()) if m.sum() > 100 else None
    print("  %s rejection: overall %.2f%%, " % (tag, 100 * out["overall"])
          + ", ".join("%s %.2f%%" % (NM[kk], 100 * out[NM[kk]])
                      for kk in (0, 1, 2, 3) if out[NM[kk]] is not None),
          flush=True)
    return out


pl = classify._pdal().Pipeline(json.dumps({"pipeline": [
    {"type": "readers.las", "filename": TS, "nosrs": True},
    {"type": "filters.expression", "expression": "Classification == 2"}]}))
pl.execute()
a = pl.arrays[0]
tx, ty, tz, _ = dedupe_xy(a["X"].astype("f8"), a["Y"].astype("f8"),
                          a["Z"].astype("f8"), keep="lowest")
del a, pl
tcls = cls_at(tx, ty)
print("TerraScan class 2: %d points" % len(tx), flush=True)

g_mcc = classify_mcc()
del RX, RY, RZ
g_smrf = np.load(os.path.join(PREV, "ground_B_cell1.00.npy"))
print("smrf cached: %d ground returns" % len(g_smrf), flush=True)

res, dems, rej = [], {}, {}
for tag, g in (("smrf", g_smrf), ("mcc", g_mcc)):
    rej[tag] = rejection(g, tcls, tx, ty, tag)
    if tag == "smrf":
        with rasterio.open(os.path.join(PREV, "dem_B_cell1.00.tif")) as r:
            z = r.read(1).astype("f4")
        z[z == 0.0] = np.nan
    else:
        z = build(g, tag)
    dems[tag] = z
    res.append(metrics(z, len(g), tag))

print("")
print("%-6s %12s %11s %9s %8s %7s %7s %9s %9s %7s %10s"
      % ("run", "class-2", "cells", "RMSE", "MAE", "p95", "p99", "coherent",
         "speckle", "rough", "nugget"))
print("-" * 110)
for r in res:
    print("%-6s %12d %11d %9.4f %8.4f %7.1f %7.1f %9.4f %9.4f %7.3f %10.6f"
          % (r["tag"], r["n_ground"], r["cells"], r["rmse"], r["mae"],
             r["p95"], r["p99"], r["coherent"], r["speckle"], r["rough"],
             r["nugget"]))
print("archive own texture: coherent %.4f, speckle %.4f"
      % (res[0]["arch_coherent"], res[0]["arch_speckle"]))
print("")
print("RMSE BY LANDFORM")
print("%-6s %10s %10s %12s %10s" % ("run", "planar", "mound", "depression",
                                    "shoulder"))
for r in res:
    b = r["rmse_by_form"]
    print("%-6s " % r["tag"] + " ".join(
        "%10s" % ("%.4f" % b[n] if b[n] is not None else "-")
        for n in ("planar", "mound", "depression", "shoulder")))
print("")
print("TERRASCAN GROUND REJECTED, TILE-WIDE")
print("%-6s %10s %10s %10s %12s %10s" % ("run", "overall", "planar", "mound",
                                         "depression", "shoulder"))
for tag in ("smrf", "mcc"):
    o = rej[tag]
    print("%-6s " % tag + " ".join(
        "%10s" % ("%.2f%%" % (100 * o[n]) if o[n] is not None else "-")
        for n in ("overall", "planar", "mound", "depression", "shoulder")))
a_, b_ = res
print("")
print("mcc against smrf: class-2 %+.1f%%, RMSE %+.1f%%, MAE %+.1f%%, "
      "mound rejection %.2f%% -> %.2f%%"
      % (100 * (b_["n_ground"] / a_["n_ground"] - 1),
         100 * (b_["rmse"] / a_["rmse"] - 1),
         100 * (b_["mae"] / a_["mae"] - 1),
         100 * rej["smrf"]["mound"], 100 * rej["mcc"]["mound"]))

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


dS, dM = A - dems["smrf"], A - dems["mcc"]
v = float(np.nanpercentile(np.abs(dS[Am]), 98))
panels = [("archive", shrink(hs(A)), "gray", 0, 1),
          ("SMRF, cell 1.00 m", shrink(hs(dems["smrf"])), "gray", 0, 1),
          ("MCC, cell 1.00 m", shrink(hs(dems["mcc"])), "gray", 0, 1),
          ("archive - SMRF, \u00b1%.2f m" % v, shrink(dS), "RdBu_r", -v, v),
          ("archive - MCC, \u00b1%.2f m" % v, shrink(dM), "RdBu_r", -v, v)]
for cm in ("gray", "RdBu_r"):
    matplotlib.colormaps[cm].set_bad("white")
fig, axes = plt.subplots(1, 5, figsize=(21.0, 17.0))
for ax, (t, img, cm, lo, hi) in zip(axes, panels):
    ax.imshow(img, cmap=cm, vmin=lo, vmax=hi, interpolation="nearest")
    ax.set_title(t, fontsize=11)
    ax.set_xticks([])
    ax.set_yticks([])
fig.suptitle("Centro GLAS l0s441, full tile \u2014 SMRF against MCC, "
             "cell 1.00 m, everything downstream unchanged", fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.97])
p = os.path.join(OUTDIR, "mcc_fulltile.png")
fig.savefig(p, dpi=200)
print("wrote %s" % p, flush=True)
with open(os.path.join(OUTDIR, "metrics_mcc.json"), "w") as fh:
    json.dump({"metrics": res, "rejection": rej}, fh, indent=2)
print("wrote %s" % os.path.join(OUTDIR, "metrics_mcc.json"), flush=True)
