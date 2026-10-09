"""Why the block-level gain reversed at full-tile scale.

On 240 m blocks the adjudicated routine beat the opening alone by 8 percent.
Over the whole tile it lost to the locked baseline by 4 percent. Two things
changed between those two measurements and only one run is needed to tell them
apart, because the classifier can now emit its stage masks.

  the working cell   calibration blocks ran at 19 to 24 returns per square
                     metre, so density adaptation left the cell at 1.0 m. The
                     full tile spans 5.0 to 24.4, and the adaptation moved the
                     cell between 0.90 and 2.00 m -- off the value validated
                     full-tile. Pinning it with adapt_density False tests this.
  the adjudication   if the confident set alone, at the pinned cell, already
                     matches the baseline, then what costs the 4 percent is
                     the gate rather than the cell.

Three surfaces come out of one classification pass over identical input: the
confident set, the candidate set, and the adjudicated result. Scoring them
against each other on the same points removes the possibility that something
else varied between runs.
"""
import json
import os
import sys
import time

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
VARIANTS = ("conf", "cand", "adj")
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
Am, Rm = np.isfinite(A), np.isfinite(R)


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
paths = {v: os.path.join(OUT, "ground_pin_%s.npy" % v) for v in VARIANTS}

if all(os.path.exists(p) for p in paths.values()):
    GS = {v: np.load(paths[v]) for v in VARIANTS}
    for v in VARIANTS:
        print("cached %s: %d returns" % (v, len(GS[v])), flush=True)
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
    acc = {v: [] for v in VARIANTS}
    reps = []
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
                # pinned: the cell stays at the value validated full-tile
                mk, rep = TSR.classify_ground(bx, by, bz, params=params,
                                              adapt_density=False,
                                              keep_masks=True)
            except TSR.GroundTSRError as exc:
                print("    tile %d,%d skipped: %s" % (i, j, exc), flush=True)
                continue
            core = ((bx >= cx0) & (bx < cx1) & (by >= cy0) & (by < cy1))
            for v, mm in (("conf", rep.pop("confident_mask")),
                          ("cand", rep.pop("candidate_mask")),
                          ("adj", mk)):
                k = mm & core
                if k.any():
                    acc[v].append(np.column_stack([bx[k], by[k], bz[k]]))
            reps.append(rep)
        print("  row %d/%d, %.1f min"
              % (i + 1, ny, (time.time() - t0) / 60.0), flush=True)
    GS = {}
    for v in VARIANTS:
        GS[v] = np.vstack(acc[v])
        np.save(paths[v], GS[v])
        print("  %s: %d returns" % (v, len(GS[v])), flush=True)
    dd = [r["density_per_m2"] for r in reps]
    print("density across %d tiles: min %.1f, median %.1f, max %.1f /m2"
          % (len(reps), np.min(dd), np.median(dd), np.max(dd)))
    print("working cell pinned at %.2f m everywhere" % params.smrf_cell_m)
    del RX, RY, RZ


def build(G, tag):
    f = os.path.join(OUT, "dem_pin_%s.tif" % tag)
    if os.path.exists(f):
        return read(f)[0]
    gx, gy, gz, _ = dedupe_xy(G[:, 0], G[:, 1], G[:, 2], keep="lowest")
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
        if i % 15 == 0:
            print("    %s dem row %d/%d, %.1f min"
                  % (tag, i, NBI, (time.time() - t0) / 60.0), flush=True)
    tile = Grid(origin_x=at.c, origin_y=at.f, cell=ac, width=WT, height=HT)
    sup = finalise.support_mask(gx, gy, tile, 4.0)
    finalise.finalise(out, tile, wkt, f, radius_m=4.0, nodata_in=-9999.0,
                      support=sup, strip_m=21.0, axis="x", verbose=False)
    return read(f)[0]


Z = {v: build(GS[v], v) for v in VARIANTS}
S = read(SMRF_DEM)[0]
common = Am & Rm & np.isfinite(S)
for v in VARIANTS:
    common &= np.isfinite(Z[v])
print("")
print("common cells: %d" % int(common.sum()))


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
items = [("confident (pinned)", Z["conf"], len(GS["conf"])),
         ("candidates (pinned)", Z["cand"], len(GS["cand"])),
         ("adjudicated (pinned)", Z["adj"], len(GS["adj"])),
         ("smrf baseline", S, 17702632),
         ("reference", R, 15313221)]
for tag, z, n in items:
    e = (z - R)[common]
    ea = (z - A)[common]
    dB = detail(z, common, 1.5)
    dBm = ndimage.median_filter(np.nan_to_num(dB, nan=0.0), size=kk)
    sel = sel0 & np.isfinite(dB)
    by = {NM[k]: float(np.sqrt(np.mean((z - A)[common & (CLS == k)] ** 2)))
          for k in (0, 1, 2, 3)}
    rows.append(dict(tag=tag, points=int(n),
                     rmse_ref=float(np.sqrt(np.mean(e ** 2))),
                     rmse_arch=float(np.sqrt(np.mean(ea ** 2))),
                     mae_arch=float(np.mean(np.abs(ea))),
                     bias_arch=float(np.mean(ea)),
                     speckle=float(np.sqrt(np.mean((dB - dBm)[sel] ** 2))),
                     rough=float(np.sqrt(np.mean(dB[sel] ** 2))
                                 / max(np.sqrt(np.mean(dA[sel] ** 2)), 1e-9)),
                     rmse_arch_by_form=by))

print("")
print("%-22s %11s %10s %11s %9s %9s %9s %7s"
      % ("run", "points", "RMSE ref", "RMSE arch", "MAE", "bias", "speckle",
         "rough"))
print("-" * 98)
for r in rows:
    print("%-22s %11d %10.4f %11.4f %9.4f %9.4f %9.4f %7.3f"
          % (r["tag"], r["points"], r["rmse_ref"], r["rmse_arch"],
             r["mae_arch"], r["bias_arch"], r["speckle"], r["rough"]))

print("")
print("RMSE AGAINST THE ARCHIVE BY LANDFORM")
print("%-22s %10s %10s %12s %10s"
      % ("run", "planar", "mound", "depression", "shoulder"))
for r in rows:
    b = r["rmse_arch_by_form"]
    print("%-22s " % r["tag"] + " ".join(
        "%10.4f" % b[n] for n in ("planar", "mound", "depression",
                                  "shoulder")))

d = {r["tag"]: r["rmse_arch"] for r in rows}
print("")
print("the two candidate causes")
print("  cell adaptation: adjudicated was 0.0817 adapted, %.4f pinned"
      % d["adjudicated (pinned)"])
print("  the gate:        confident %.4f, adjudicated %.4f, baseline %.4f"
      % (d["confident (pinned)"], d["adjudicated (pinned)"],
         d["smrf baseline"]))
with open(os.path.join(OUT, "tsr_pinned_metrics.json"), "w") as fh:
    json.dump(rows, fh, indent=2)
print("")
print("wrote %s" % os.path.join(OUT, "tsr_pinned_metrics.json"))
