"""South GLAS l0s395: SMRF, ground_tsr and the delivered T/S surface.

A tile that has never been processed here. Every parameter in both classifiers
was fixed on Centro GLAS l0s441, so nothing on this tile was fitted and the
comparison is a generalisation test rather than a second look at the
calibration data.

"T/S" is read as the delivered TerraScan/Surfer product -- the 0.5 m DEM in
the archive tree. It is rendered from that DEM through the same G1 code as the
other two rather than from the delivered G1 image, because a shared stretch is
what makes the three panels comparable and the delivered image carries its own.

Held identical to the l0s441 work: 400 m classification tiles with a 40 m
buffer, the global biharmonic engine at pool 0.25, a 4 m support mask, a 21 m
across-track edge strip, vertical exaggeration 1, sun elevation 35 degrees,
8 hillshade directions, 16 sky-view directions at radius 10 cells, and one
1-99 percentile stretch per layer shared across all three surfaces.

ground_tsr runs with adapt_density False. Its density rule was measured to
cost 4 percent at full-tile scale on l0s441 and has not been replaced, so the
working cell is pinned at the value validated there. That is a known open
fault, not a tuning choice made for this tile.
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
from replicalm import visualize as VIS
from replicalm.config import PRESETS
from replicalm.grid import Grid
from replicalm.preprocess import dedupe_xy

STEM = "South_GLAS_l0s395"
RAW = (r"D:\GLiHT_LAS_orig\Yuc_South\South_Glas\South_Glas"
       r"\AMIGACarb_Yuc_South_GLAS_Apr2013_l0s395.las")
GND = (r"D:\_LAZ\Yuc_South\South_Glas\South_Glas"
       r"\South_GLAS_l0s395_GND_v1.laz")
TSDEM = (r"D:\_Archive_EdgeFixed\Yuc_South\South_Glas\South_Glas"
         r"\South_GLAS_l0s395_DEM_0p5m_v1.tif")
OUT = r"C:\Replicalm\test100526\south395"
SCRATCH = (r"C:\Users\benja\AppData\Local\Temp\claude\C--m2"
           r"\c9449ab5-87f2-48b7-862d-4c77b170476d\scratchpad\g1_cache_s395")
TILE_M, TBUF_M = 400.0, 40.0
CORE, BBUF, POOL = 300, 40, 0.25
BLOCK, HALO, N_STRIPS = 2048, 48, 4
VE_FACTOR, PROM_T = 1.0, 0.25
os.makedirs(OUT, exist_ok=True)
os.makedirs(SCRATCH, exist_ok=True)


def read(path):
    with rasterio.open(path) as r:
        z = r.read(1).astype("f4")
        nd, prof = r.nodata, r.profile
        tr, cc, wkt = r.transform, abs(r.transform.a), r.crs.to_wkt()
    if nd is not None:
        z[z == nd] = np.nan
    z[z > 1e30] = np.nan
    z[z == 0.0] = np.nan
    return z, tr, cc, wkt, prof


TS, at, ac, wkt, prof = read(TSDEM)
HT, WT = TS.shape
TSm = np.isfinite(TS)
print("%s: %d x %d at %.3f m, T/S covers %.2f%%"
      % (STEM, HT, WT, ac, 100 * TSm.mean()), flush=True)
x0, x1 = at.c, at.c + WT * ac
y1, y0 = at.f, at.f - HT * ac
axis = finalise.infer_transect_axis(
    Grid(origin_x=at.c, origin_y=at.f, cell=ac, width=WT, height=HT), STEM)
print("edge strip axis inferred from the tile name: %s" % axis, flush=True)

cfg = PRESETS["clear"]
p0 = cfg.passes[0]
print("SMRF: window %.0f, slope %.3f, threshold %.2f, scalar %.1f, cell %.2f"
      % (p0.window_m, p0.slope, p0.threshold_m, p0.scalar, p0.cell_m),
      flush=True)

RXYZ = os.path.join(SCRATCH, "raw.npy")
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


def smrf_block(px, py, pz):
    arr = np.zeros(len(pz), dtype=[("X", "f8"), ("Y", "f8"), ("Z", "f8"),
                                   ("Classification", "u1")])
    arr["X"], arr["Y"], arr["Z"] = px, py, pz
    stg = [{"type": "filters.assign", "value": "Classification = 0"},
           {"type": "filters.outlier", "method": "statistical",
            "mean_k": cfg.outlier_neighbours,
            "multiplier": cfg.outlier_multiplier, "class": cfg.noise_class},
           {"type": "filters.smrf", "window": p0.window_m, "slope": p0.slope,
            "threshold": p0.threshold_m, "scalar": p0.scalar,
            "cell": p0.cell_m,
            "ignore": "Classification[%d:%d]" % (cfg.noise_class,
                                                 cfg.noise_class)},
           {"type": "filters.expression", "expression": "Classification == 2"}]
    p = classify._pdal().Pipeline(json.dumps({"pipeline": stg}), arrays=[arr])
    p.execute()
    g = p.arrays[0]
    return g["X"].astype("f8"), g["Y"].astype("f8"), g["Z"].astype("f8")


def classify_all(which):
    f = os.path.join(OUT, "ground_%s.npy" % which)
    if os.path.exists(f):
        g = np.load(f)
        print("%s cached: %d returns" % (which, len(g)), flush=True)
        return g
    out, reps = [], []
    params = TSR.GroundTSRParams()
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
            if which == "smrf":
                gx, gy, gz = smrf_block(bx, by, bz)
            else:
                try:
                    mask, rep = TSR.classify_ground(bx, by, bz, params=params,
                                                    adapt_density=False)
                except TSR.GroundTSRError as exc:
                    print("    tile %d,%d skipped: %s" % (i, j, exc),
                          flush=True)
                    continue
                reps.append(rep)
                gx, gy, gz = bx[mask], by[mask], bz[mask]
            k = (gx >= cx0) & (gx < cx1) & (gy >= cy0) & (gy < cy1)
            if k.any():
                out.append(np.column_stack([gx[k], gy[k], gz[k]]))
        print("  %s row %d/%d, %.1f min"
              % (which, i + 1, ny, (time.time() - t0) / 60.0), flush=True)
    g = np.vstack(out)
    np.save(f, g)
    print("%s: %d returns in %.1f min" % (which, len(g),
                                          (time.time() - t0) / 60.0),
          flush=True)
    if reps:
        dd = [r["density_per_m2"] for r in reps]
        print("  density across %d tiles: min %.1f, median %.1f, max %.1f /m2"
              % (len(reps), np.min(dd), np.median(dd), np.max(dd)), flush=True)
        with open(os.path.join(OUT, "tsr_reports.json"), "w") as fh:
            json.dump(reps, fh, indent=2)
    return g


def build(G, tag):
    f = os.path.join(OUT, "dem_%s.tif" % tag)
    if os.path.exists(f):
        return read(f)[0]
    gx, gy, gz, _ = dedupe_xy(G[:, 0], G[:, 1], G[:, 2], keep="lowest")
    print("  %s dedupe: %d -> %d" % (tag, len(G), len(gz)), flush=True)
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
    rp = finalise.finalise(out, tile, wkt, f, radius_m=4.0, nodata_in=-9999.0,
                           support=sup, strip_m=21.0, axis=axis,
                           verbose=False)
    print("  %s finalised: %d cells, %d interior restored"
          % (tag, rp["cells_after_trim"], rp["interior_cells_restored"]),
          flush=True)
    return read(f)[0]


DEM = {"ts": TS}
for which in ("smrf", "tsr"):
    DEM[which] = build(classify_all(which), which)
del RX, RY, RZ

common = TSm.copy()
for k in DEM:
    common &= np.isfinite(DEM[k])
print("")
print("cells present in all three: %d (%.2f%% of the tile)"
      % (int(common.sum()), 100 * common.mean()), flush=True)


# ------------------------------------------------------------------- G1
def layer_blocked(dem, name, path):
    out = np.memmap(path, dtype="float32", mode="w+", shape=dem.shape)
    for i in range(int(np.ceil(HT / BLOCK))):
        for j in range(int(np.ceil(WT / BLOCK))):
            r0, c0 = i * BLOCK, j * BLOCK
            r1, c1 = min(HT, r0 + BLOCK), min(WT, c0 + BLOCK)
            hr0, hc0 = max(0, r0 - HALO), max(0, c0 - HALO)
            hr1, hc1 = min(HT, r1 + HALO), min(WT, c1 + HALO)
            sub = dem[hr0:hr1, hc0:hc1]
            if not np.isfinite(sub).any():
                out[r0:r1, c0:c1] = np.nan
                continue
            fill = float(np.nanmedian(sub))
            d = np.where(np.isfinite(sub), sub, fill).astype("f8")
            raw = VIS._centre_into(VIS.compute_layer(name, d, ac), sub.shape)
            raw = np.where(np.isfinite(sub), raw, np.nan)
            out[r0:r1, c0:c1] = raw[r0 - hr0:r0 - hr0 + (r1 - r0),
                                    c0 - hc0:c0 - hc0 + (c1 - c0)]
    out.flush()
    return out


print("")
print("computing G1 layers")
cache, samples = {}, {nm: [] for nm in VIS.LAYER_NAMES}
for tag in DEM:
    cache[tag] = {}
    for nm in VIS.LAYER_NAMES:
        t0 = time.time()
        mp = os.path.join(SCRATCH, "%s_%s.f32" % (tag, nm))
        arr = layer_blocked(DEM[tag], nm, mp)
        cache[tag][nm] = mp
        s = arr[::7, ::7].ravel()
        samples[nm].append(s[np.isfinite(s)])
        print("  %-5s %-8s %5.1f s" % (tag, nm, time.time() - t0), flush=True)
        del arr
shared = {}
for nm in VIS.LAYER_NAMES:
    allv = np.concatenate(samples[nm])
    shared[nm] = (float(np.percentile(allv, 1)),
                  float(np.percentile(allv, 99)))
    print("  shared stretch %-8s %.4f to %.4f" % (nm, shared[nm][0],
                                                  shared[nm][1]))
del samples

NAMES = {"ts": "TS", "smrf": "SMRF", "tsr": "tsr"}
G1 = {}
for tag in DEM:
    canvas = None
    for nm, transp, mode, br, co, ga in VIS.G1_LAYERS:
        lo_, hi_ = shared[nm]
        arr = np.memmap(cache[tag][nm], dtype="float32", mode="r",
                        shape=(HT, WT))
        v = np.clip((np.asarray(arr) - lo_) / max(hi_ - lo_, 1e-9), 0.0, 1.0)
        del arr
        v = VIS.adjust(v.astype("f4"), br, co, ga)
        canvas = v if canvas is None else VIS.blend(
            canvas, v, mode, 1.0 - transp / 100.0)
    ok = np.isfinite(DEM[tag])
    G1[tag] = np.where(ok, canvas, np.nan)
    b = np.clip(np.nan_to_num(G1[tag], nan=0.0) * 255.0, 0,
                255).astype("uint8")
    b = np.where(ok & (b == 0), 1, b).astype("uint8")
    rgb = np.dstack([b, b, b])
    rgb[~ok] = 0
    p = dict(prof)
    p.update(count=3, dtype="uint8", nodata=None, compress="lzw", tiled=True,
             blockxsize=512, blockysize=512, predictor=2,
             BIGTIFF="IF_SAFER")
    p.pop("photometric", None)
    tifp = os.path.join(OUT, "G1_%s_%s.tif" % (STEM, NAMES[tag]))
    with rasterio.open(tifp, "w", **p) as dst:
        for k in range(3):
            dst.write(rgb[:, :, k], k + 1)
        dst.update_tags(VERTICAL_EXAGGERATION="%.1f" % VE_FACTOR,
                        HS_SUN_ELEVATION="%d" % VIS.HS_SUN_ELEVATION,
                        HS_DIRECTIONS="%d" % VIS.HS_DIRECTIONS,
                        SVF_DIRECTIONS="%d" % VIS.SVF_DIRECTIONS,
                        SVF_RADIUS_CELLS="%d" % VIS.SVF_RADIUS_CELLS,
                        STRETCH="shared 1-99 percentile over TS, SMRF, tsr")
    print("  wrote %s" % tifp, flush=True)
    from PIL import Image
    Image.MAX_IMAGE_PIXELS = None
    rgb[~ok] = 255
    step = int(np.ceil(HT / N_STRIPS))
    for k in range(N_STRIPS):
        r0, r1 = k * step, min(HT, (k + 1) * step)
        if r0 >= r1:
            break
        pp = os.path.join(OUT, "G1_%s_%s_strip%d_rows%d-%d.png"
                          % (STEM, NAMES[tag], k + 1, r0, r1))
        Image.fromarray(rgb[r0:r1]).save(pp)
        print("    %s" % pp, flush=True)
    del canvas, b, rgb


# --------------------------------------------------------------- metrics
def detail(z, mask, scale):
    s = scale / ac / 2.0
    num = ndimage.gaussian_filter(np.where(mask, z, 0.0).astype("f4"), s)
    den = ndimage.gaussian_filter(mask.astype("f4"), s)
    sm = np.where(den > 1e-6, num / np.maximum(den, 1e-6), 0.0)
    return np.where(mask, z - sm, np.nan)


sm2 = detail(TS, common, 0.0) if False else None
num = ndimage.gaussian_filter(np.where(common, TS, 0.0), 2.0 / ac / 2.0)
den = ndimage.gaussian_filter(common.astype("f4"), 2.0 / ac / 2.0)
s2 = np.where(den > 1e-6, num / np.maximum(den, 1e-6), np.nan)
num = ndimage.gaussian_filter(np.where(common, TS, 0.0), 10.0 / ac / 2.0)
den = ndimage.gaussian_filter(common.astype("f4"), 10.0 / ac / 2.0)
s10 = np.where(den > 1e-6, num / np.maximum(den, 1e-6), np.nan)
prom = (s2 - s10).astype("f4")
CLS = np.zeros(TS.shape, np.int8)
CLS[common & (prom > PROM_T)] = 1
CLS[common & (prom < -PROM_T)] = 2
NM = {0: "planar", 1: "mound", 2: "depression"}
del s2, s10, prom, num, den

kk = max(3, int(round(1.0 / ac)) | 1)
inner = ndimage.binary_erosion(common, structure=np.ones((3, 3)),
                               iterations=int(5.0 / ac))
dT = detail(TS, common, 1.5)
dTm = ndimage.median_filter(np.nan_to_num(dT, nan=0.0), size=kk)
sel0 = inner & np.isfinite(dT)
met = {"tile": STEM, "cells": int(common.sum()), "cell_m": ac,
       "edge_axis": axis, "vertical_exaggeration": VE_FACTOR,
       "shared_stretch": shared, "surfaces": {},
       "ts_texture": {"coherent": float(np.sqrt(np.mean(dTm[sel0] ** 2))),
                      "speckle": float(np.sqrt(np.mean(
                          (dT - dTm)[sel0] ** 2)))}}
for tag in ("smrf", "tsr"):
    z = DEM[tag]
    e = (z - TS)[common]
    dB = detail(z, common, 1.5)
    dBm = ndimage.median_filter(np.nan_to_num(dB, nan=0.0), size=kk)
    sel = sel0 & np.isfinite(dB)
    met["surfaces"][tag] = {
        "rmse_vs_ts": float(np.sqrt(np.mean(e ** 2))),
        "mae_vs_ts": float(np.mean(np.abs(e))),
        "bias_vs_ts": float(np.mean(e)),
        "p99_abs": float(np.percentile(np.abs(e), 99)),
        "coherent": float(np.sqrt(np.mean(dBm[sel] ** 2))),
        "speckle": float(np.sqrt(np.mean((dB - dBm)[sel] ** 2))),
        "roughness_ratio": float(np.sqrt(np.mean(dB[sel] ** 2))
                                 / max(np.sqrt(np.mean(dT[sel] ** 2)), 1e-9)),
        "rmse_by_form": {NM[k]: float(np.sqrt(np.mean(
            (z - TS)[common & (CLS == k)] ** 2))) for k in (0, 1, 2)}}
d = (DEM["tsr"] - DEM["smrf"])[common]
met["tsr_minus_smrf"] = {"rmse": float(np.sqrt(np.mean(d ** 2))),
                         "bias": float(np.mean(d)),
                         "p99_abs": float(np.percentile(np.abs(d), 99))}

print("")
print("%-6s %10s %9s %9s %9s %9s %9s %8s"
      % ("run", "RMSE v TS", "MAE", "bias", "p99|e|", "coherent", "speckle",
         "rough"))
for tag in ("smrf", "tsr"):
    s = met["surfaces"][tag]
    print("%-6s %10.4f %9.4f %+9.4f %9.4f %9.4f %9.4f %8.3f"
          % (tag, s["rmse_vs_ts"], s["mae_vs_ts"], s["bias_vs_ts"],
             s["p99_abs"], s["coherent"], s["speckle"],
             s["roughness_ratio"]))
print("T/S own texture: coherent %.4f, speckle %.4f"
      % (met["ts_texture"]["coherent"], met["ts_texture"]["speckle"]))
print("")
print("%-6s %10s %10s %12s" % ("run", "planar", "mound", "depression"))
for tag in ("smrf", "tsr"):
    b = met["surfaces"][tag]["rmse_by_form"]
    print("%-6s %10.4f %10.4f %12.4f"
          % (tag, b["planar"], b["mound"], b["depression"]))
print("")
print("tsr minus smrf: RMSE %.4f, bias %+.4f, p99 %.4f"
      % (met["tsr_minus_smrf"]["rmse"], met["tsr_minus_smrf"]["bias"],
         met["tsr_minus_smrf"]["p99_abs"]))

with open(os.path.join(OUT, "south395_metrics.json"), "w") as fh:
    json.dump(met, fh, indent=2)
print("")
print("wrote %s" % os.path.join(OUT, "south395_metrics.json"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

K = 10


def shrink(x, k=K):
    h, w = x.shape
    x = x[:h // k * k, :w // k * k]
    with np.errstate(invalid="ignore"):
        return np.nanmean(x.reshape(h // k, k, w // k, k), axis=(1, 3))


matplotlib.colormaps["gray"].set_bad("white")
fig, axes = plt.subplots(1, 3, figsize=(16.0, 18.0))
for ax, tag in zip(axes, ("smrf", "tsr", "ts")):
    ax.imshow(shrink(G1[tag]), cmap="gray", vmin=0, vmax=1,
              interpolation="nearest")
    ax.set_title(NAMES[tag], fontsize=12)
    ax.set_xticks([])
    ax.set_yticks([])
fig.suptitle("%s \u2014 G1 at a shared stretch. Overview only; inspect the "
             "GeoTIFFs and strips at native 0.5 m." % STEM, fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.97])
pp = os.path.join(OUT, "south395_overview.png")
fig.savefig(pp, dpi=200)
print("wrote %s" % pp)
