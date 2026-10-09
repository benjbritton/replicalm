"""Full-tile G1 package at native resolution, for close visual inspection.

Three surfaces -- the delivered class 2 + 8 reference, SMRF, ground_tsr -- and
three difference rasters, all on the same 15158 by 3157 grid at 0.5 m, one set
of hillshade parameters, one vertical exaggeration of 1, and one colour ramp
per product type.

Two things are done differently from `visualize.g1`, both required for the
panels to be comparable rather than merely adjacent.

Each G1 layer is normalised against a 1st-99th percentile range shared by all
three surfaces, not by its own min and max. Per-surface stretching would give
each panel its own contrast and any difference seen between them could be the
stretch rather than the terrain.

The layers are computed in blocks with a 48-cell halo and the halo is then
discarded. The horizon scans behind the openness and sky-view layers reach 10
cells, and the multiple hillshade 1, so a 48-cell halo makes the blocked result
identical to an unblocked one while keeping peak memory bounded -- a whole-tile
SVF at this size would otherwise allocate 16 directions over 47.8 million cells
at once.

Output is georeferenced GeoTIFF at native resolution, which is what a GIS needs
for close inspection, plus PNG strips at native resolution for viewing without
one. Nothing is downsampled except the single overview contact sheet.
"""
import json
import os
import sys
import time

import numpy as np
import rasterio
from rasterio.enums import Compression
from scipy import ndimage

sys.path.insert(0, r"C:\Replicalm\src")
from replicalm import visualize as VIS

REF = r"C:\Replicalm\test100526\terrascan_ceiling\dem_ts_c2c8.tif"
SMRF = r"C:\Replicalm\test100526\cellsize_validation\dem_B_cell1.00.tif"
TSRD = r"C:\Replicalm\test100526\ground_replication\dem_pin_adj.tif"
ARCH = (r"D:\_Archive_EdgeFixed\Yuc_Centro\Centro_Glas\Centro_Glas"
        r"\Centro_GLAS_l0s441_DEM_0p5m_v1.tif")
OUT = r"C:\Replicalm\test100526\g1_fulltile"
SCRATCH = (r"C:\Users\benja\AppData\Local\Temp\claude\C--m2"
           r"\c9449ab5-87f2-48b7-862d-4c77b170476d\scratchpad\g1_cache")
BLOCK, HALO = 2048, 48
VE_FACTOR = 1.0
N_STRIPS = 4
PROM_T = 0.25
SURFACES = [("reference", REF), ("smrf", SMRF), ("tsr", TSRD)]
os.makedirs(OUT, exist_ok=True)
os.makedirs(SCRATCH, exist_ok=True)


def read(path):
    with rasterio.open(path) as r:
        z = r.read(1).astype("f4")
        nd, prof = r.nodata, r.profile
        tr, cc, crs = r.transform, abs(r.transform.a), r.crs
    if nd is not None:
        z[z == nd] = np.nan
    z[z > 1e30] = np.nan
    z[z == 0.0] = np.nan
    return z, tr, cc, crs, prof


DEM = {}
for tag, path in SURFACES:
    DEM[tag], at, ac, crs, prof = read(path)
A = read(ARCH)[0]
HT, WT = DEM["reference"].shape
print("tile %d x %d at %.3f m, vertical exaggeration %.1f"
      % (HT, WT, ac, VE_FACTOR), flush=True)
common = np.isfinite(A)
for tag in DEM:
    common &= np.isfinite(DEM[tag])
print("cells present in the archive and all three surfaces: %d (%.2f%%)"
      % (int(common.sum()), 100 * common.mean()), flush=True)


# --------------------------------------------------------------- G1 layers
def layer_blocked(dem, name, memmap_path):
    """One rvt layer over the whole grid, computed in haloed blocks."""
    out = np.memmap(memmap_path, dtype="float32", mode="w+", shape=dem.shape)
    nbi = int(np.ceil(HT / BLOCK))
    nbj = int(np.ceil(WT / BLOCK))
    for i in range(nbi):
        for j in range(nbj):
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
            raw = VIS.compute_layer(name, d, ac)
            raw = VIS._centre_into(raw, sub.shape)
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
        print("  %-10s %-8s %5.1f s" % (tag, nm, time.time() - t0),
              flush=True)
        del arr

shared = {}
for nm in VIS.LAYER_NAMES:
    allv = np.concatenate(samples[nm])
    shared[nm] = (float(np.percentile(allv, 1)),
                  float(np.percentile(allv, 99)))
    print("  shared stretch %-8s %.4f to %.4f" % (nm, shared[nm][0],
                                                  shared[nm][1]))
del samples


def write_rgb(rgb, path):
    p = dict(prof)
    p.update(count=3, dtype="uint8", nodata=None, compress="lzw",
             tiled=True, blockxsize=512, blockysize=512, predictor=2,
             BIGTIFF="IF_SAFER")
    p.pop("photometric", None)
    with rasterio.open(path, "w", **p) as dst:
        for b in range(3):
            dst.write(rgb[:, :, b], b + 1)
        dst.update_tags(VERTICAL_EXAGGERATION="%.1f" % VE_FACTOR,
                        HS_SUN_ELEVATION="%d" % VIS.HS_SUN_ELEVATION,
                        HS_DIRECTIONS="%d" % VIS.HS_DIRECTIONS,
                        SVF_DIRECTIONS="%d" % VIS.SVF_DIRECTIONS,
                        SVF_RADIUS_CELLS="%d" % VIS.SVF_RADIUS_CELLS,
                        STRETCH="shared 1-99 percentile across all surfaces")
    return path


print("")
print("blending and writing G1")
G1 = {}
for tag in DEM:
    canvas = None
    for nm, transp, mode, br, co, ga in VIS.G1_LAYERS:
        lo_, hi_ = shared[nm]
        arr = np.memmap(cache[tag][nm], dtype="float32", mode="r",
                        shape=(HT, WT))
        a = np.clip((np.asarray(arr) - lo_) / max(hi_ - lo_, 1e-9), 0.0, 1.0)
        del arr
        a = VIS.adjust(a.astype("f4"), br, co, ga)
        canvas = a if canvas is None else VIS.blend(
            canvas, a, mode, 1.0 - transp / 100.0)
    ok = np.isfinite(DEM[tag])
    g = np.where(ok, canvas, np.nan)
    G1[tag] = g
    v = np.clip(np.nan_to_num(g, nan=0.0) * 255.0, 0, 255).astype("uint8")
    # a legitimate G1 value can land on zero; lift valid pixels off pure
    # black so nodata stays distinguishable from deep shadow
    v = np.where(ok & (v == 0), 1, v).astype("uint8")
    rgb = np.dstack([v, v, v])
    rgb[~ok] = 0
    p = write_rgb(rgb, os.path.join(OUT, "G1_%s.tif" % tag))
    print("  wrote %s" % p, flush=True)
    del canvas, v, rgb

# ---------------------------------------------------------------- differences
DIFFS = [("smrf_minus_reference", "smrf", "reference"),
         ("tsr_minus_reference", "tsr", "reference"),
         ("tsr_minus_smrf", "tsr", "smrf")]
D = {}
for name, a_, b_ in DIFFS:
    D[name] = (DEM[a_] - DEM[b_]).astype("f4")
    p = dict(prof)
    p.update(count=1, dtype="float32", nodata=-9999.0, compress="lzw",
             tiled=True, blockxsize=512, blockysize=512,
             BIGTIFF="IF_SAFER")
    with rasterio.open(os.path.join(OUT, "diff_%s.tif" % name), "w",
                       **p) as dst:
        dst.write(np.where(np.isfinite(D[name]), D[name], -9999.0), 1)
    print("  wrote %s" % os.path.join(OUT, "diff_%s.tif" % name), flush=True)

sc = np.concatenate([np.abs(D[n][common]) for n, _, _ in DIFFS])
SHARED_V = float(np.percentile(sc[np.isfinite(sc)], 98))
OWN_V = {n: float(np.percentile(np.abs(D[n][common])[
    np.isfinite(np.abs(D[n][common]))], 98)) for n, _, _ in DIFFS}
print("")
print("shared difference scale +/- %.4f m" % SHARED_V)
for n in OWN_V:
    print("  own scale %-22s +/- %.4f m" % (n, OWN_V[n]))


def diverging_rgb(d, v):
    """Blue-white-red by hand, so no plotting backend sees 48 million cells."""
    t = np.clip(d / max(v, 1e-9), -1.0, 1.0)
    out = np.zeros(d.shape + (3,), dtype="uint8")
    neg, pos = t < 0, t >= 0
    s = np.abs(t)
    r = np.where(neg, 1.0 - s, 1.0)
    g = 1.0 - s
    b = np.where(neg, 1.0, 1.0 - s)
    out[..., 0] = np.clip(r * 255, 0, 255)
    out[..., 1] = np.clip(g * 255, 0, 255)
    out[..., 2] = np.clip(b * 255, 0, 255)
    bad = ~np.isfinite(d)
    out[bad] = 255
    return out


def write_png_strips(img, stem, n_strips=N_STRIPS):
    from PIL import Image
    Image.MAX_IMAGE_PIXELS = None
    h = img.shape[0]
    step = int(np.ceil(h / n_strips))
    paths = []
    for k in range(n_strips):
        r0, r1 = k * step, min(h, (k + 1) * step)
        if r0 >= r1:
            break
        p = os.path.join(OUT, "%s_strip%d_rows%d-%d.png"
                         % (stem, k + 1, r0, r1))
        Image.fromarray(img[r0:r1]).save(p, optimize=False)
        paths.append(p)
    return paths


print("")
print("writing PNG strips at native resolution")
for tag in DEM:
    ok = np.isfinite(DEM[tag])
    v = np.clip(np.nan_to_num(G1[tag], nan=0.0) * 255.0, 0, 255).astype("uint8")
    v = np.where(ok & (v == 0), 1, v).astype("uint8")
    rgb = np.dstack([v, v, v])
    rgb[~ok] = 255
    for p in write_png_strips(rgb, "G1_%s" % tag):
        print("  %s" % p, flush=True)
    del v, rgb
for name, _, _ in DIFFS:
    for lbl, vv in (("shared", SHARED_V), ("own", OWN_V[name])):
        img = diverging_rgb(D[name], vv)
        for p in write_png_strips(img, "diff_%s_%s%.3fm"
                                  % (name, lbl, vv)):
            print("  %s" % p, flush=True)
        del img

# ------------------------------------------------------------------- metrics
sm2 = ndimage.gaussian_filter(np.where(common, A, 0.0), 2.0 / ac / 2.0)
sm2d = ndimage.gaussian_filter(common.astype("f4"), 2.0 / ac / 2.0)
sm2 = np.where(sm2d > 1e-6, sm2 / np.maximum(sm2d, 1e-6), np.nan)
sm10 = ndimage.gaussian_filter(np.where(common, A, 0.0), 10.0 / ac / 2.0)
sm10d = ndimage.gaussian_filter(common.astype("f4"), 10.0 / ac / 2.0)
sm10 = np.where(sm10d > 1e-6, sm10 / np.maximum(sm10d, 1e-6), np.nan)
prom = (sm2 - sm10).astype("f4")
CLS = np.zeros(A.shape, np.int8)
CLS[common & (prom > PROM_T)] = 1
CLS[common & (prom < -PROM_T)] = 2
NM = {0: "planar", 1: "mound", 2: "depression"}
del sm2, sm10, sm2d, sm10d, prom


def detail(z, mask, scale):
    s = scale / ac / 2.0
    num = ndimage.gaussian_filter(np.where(mask, z, 0.0).astype("f4"), s)
    den = ndimage.gaussian_filter(mask.astype("f4"), s)
    sm = np.where(den > 1e-6, num / np.maximum(den, 1e-6), 0.0)
    return np.where(mask, z - sm, np.nan)


kk = max(3, int(round(1.0 / ac)) | 1)
inner = ndimage.binary_erosion(common, structure=np.ones((3, 3)),
                               iterations=int(5.0 / ac))
dA = detail(A, common, 1.5)
dAm = ndimage.median_filter(np.nan_to_num(dA, nan=0.0), size=kk)
sel0 = inner & np.isfinite(dA)
met = {"cells": int(common.sum()), "cell_m": ac,
       "vertical_exaggeration": VE_FACTOR,
       "hillshade": {"sun_elevation": VIS.HS_SUN_ELEVATION,
                     "directions": VIS.HS_DIRECTIONS,
                     "svf_directions": VIS.SVF_DIRECTIONS,
                     "svf_radius_cells": VIS.SVF_RADIUS_CELLS},
       "shared_difference_scale_m": SHARED_V, "own_scales_m": OWN_V,
       "surfaces": {}, "differences": {},
       "archive_texture": {}}
for tag in DEM:
    z = DEM[tag]
    e = (z - A)[common]
    dB = detail(z, common, 1.5)
    dBm = ndimage.median_filter(np.nan_to_num(dB, nan=0.0), size=kk)
    sel = sel0 & np.isfinite(dB)
    met["surfaces"][tag] = {
        "rmse_vs_archive": float(np.sqrt(np.mean(e ** 2))),
        "mae_vs_archive": float(np.mean(np.abs(e))),
        "bias_vs_archive": float(np.mean(e)),
        "coherent": float(np.sqrt(np.mean(dBm[sel] ** 2))),
        "speckle": float(np.sqrt(np.mean((dB - dBm)[sel] ** 2))),
        "roughness_ratio": float(np.sqrt(np.mean(dB[sel] ** 2))
                                 / max(np.sqrt(np.mean(dA[sel] ** 2)), 1e-9)),
        "rmse_vs_archive_by_form": {
            NM[k]: float(np.sqrt(np.mean((z - A)[common & (CLS == k)] ** 2)))
            for k in (0, 1, 2)}}
met["archive_texture"] = {
    "coherent": float(np.sqrt(np.mean(dAm[sel0] ** 2))),
    "speckle": float(np.sqrt(np.mean((dA - dAm)[sel0] ** 2)))}
for name, _, _ in DIFFS:
    d = D[name][common]
    d = d[np.isfinite(d)]
    met["differences"][name] = {
        "rmse": float(np.sqrt(np.mean(d ** 2))),
        "mae": float(np.mean(np.abs(d))), "bias": float(np.mean(d)),
        "p99_abs": float(np.percentile(np.abs(d), 99)),
        "max_abs": float(np.max(np.abs(d))),
        "fraction_over_0.10m": float((np.abs(d) > 0.10).mean()),
        "fraction_over_0.25m": float((np.abs(d) > 0.25).mean()),
        "by_form": {NM[k]: float(np.sqrt(np.mean(
            D[name][common & (CLS == k)] ** 2))) for k in (0, 1, 2)}}

print("")
print("%-11s %11s %10s %10s %10s %10s %8s" % ("surface", "RMSE arch", "MAE",
                                              "bias", "coherent", "speckle",
                                              "rough"))
for tag in ("reference", "smrf", "tsr"):
    s = met["surfaces"][tag]
    print("%-11s %11.4f %10.4f %10.4f %10.4f %10.4f %8.3f"
          % (tag, s["rmse_vs_archive"], s["mae_vs_archive"],
             s["bias_vs_archive"], s["coherent"], s["speckle"],
             s["roughness_ratio"]))
print("archive own texture: coherent %.4f, speckle %.4f"
      % (met["archive_texture"]["coherent"],
         met["archive_texture"]["speckle"]))
print("")
print("%-24s %9s %9s %10s %9s %12s %12s"
      % ("difference", "RMSE", "MAE", "bias", "p99|d|", ">0.10 m", ">0.25 m"))
for name, _, _ in DIFFS:
    d = met["differences"][name]
    print("%-24s %9.4f %9.4f %+10.4f %9.4f %11.2f%% %11.2f%%"
          % (name, d["rmse"], d["mae"], d["bias"], d["p99_abs"],
             100 * d["fraction_over_0.10m"],
             100 * d["fraction_over_0.25m"]))
print("")
print("%-24s %10s %10s %12s" % ("difference RMSE by form", "planar", "mound",
                                "depression"))
for name, _, _ in DIFFS:
    b = met["differences"][name]["by_form"]
    print("%-24s %10.4f %10.4f %12.4f"
          % (name, b["planar"], b["mound"], b["depression"]))

with open(os.path.join(OUT, "g1_fulltile_metrics.json"), "w") as fh:
    json.dump(met, fh, indent=2)
print("")
print("wrote %s" % os.path.join(OUT, "g1_fulltile_metrics.json"))

# overview contact sheet, the only downsampled product
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

K = 10


def shrink(x, k=K):
    h, w = x.shape
    x = x[:h // k * k, :w // k * k]
    with np.errstate(invalid="ignore"):
        return np.nanmean(x.reshape(h // k, k, w // k, k), axis=(1, 3))


for cm in ("gray", "RdBu_r"):
    matplotlib.colormaps[cm].set_bad("white")
panels = [("TerraScan 2+8 reference", shrink(G1["reference"]), "gray", 0, 1),
          ("SMRF", shrink(G1["smrf"]), "gray", 0, 1),
          ("ground_tsr", shrink(G1["tsr"]), "gray", 0, 1)]
panels += [(n.replace("_", " ") + "\n\u00b1%.3f m" % SHARED_V,
            shrink(D[n]), "RdBu_r", -SHARED_V, SHARED_V)
           for n, _, _ in DIFFS]
fig, axes = plt.subplots(1, 6, figsize=(26.0, 18.0))
for ax, (t, img, cm, lo_, hi_) in zip(axes, panels):
    ax.imshow(img, cmap=cm, vmin=lo_, vmax=hi_, interpolation="nearest")
    ax.set_title(t, fontsize=10)
    ax.set_xticks([])
    ax.set_yticks([])
fig.suptitle("Centro GLAS l0s441 full tile \u2014 G1 at a shared stretch, and "
             "the three difference rasters at a shared scale. Overview only; "
             "inspect the GeoTIFFs and strips at native 0.5 m.", fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.97])
pp = os.path.join(OUT, "overview_contact_sheet.png")
fig.savefig(pp, dpi=200)
print("wrote %s" % pp)
