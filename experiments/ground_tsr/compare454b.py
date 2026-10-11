"""Local spline against global biharmonic, both through the pipeline.

Three changes separate the two runs and they are not independent, so the
comparison is of configurations rather than of any single edit:

  spline run   local tensioned spline, dedupe keeping the lowest of each
               coincident pair, coverage from a 4 m disc union
  biharmonic   global minimum curvature in 150 m cores, dedupe keeping the
               mean, coverage from one span per row

The texture ratios are the measurement that matters here. A surface that is
merely cleaner loses speckle and keeps coherent detail; one that is smeared
loses both. The archive is the reference for each.
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from scipy import ndimage

ARCH_DEM = (r"D:\_Archive_EdgeFixed\Yuc_Centro\Centro_Glas\Centro_Glas"
            r"\Centro_GLAS_l0s454_DEM_0p5m_v1.tif")
ARCH_G1 = (r"D:\_Archive_EdgeFixed\Yuc_Centro\Centro_Glas\Centro_Glas"
           r"\Centro_GLAS_l0s454_G1_0p5m_v1.tif")
RUNS = [
    ("spline", r"C:\Replicalm\test100526\centro454_clear50"),
    ("biharmonic", r"C:\Replicalm\test100526\centro454_bih50"),
]
BASE = "AMIGACarb_Yuc_Centro_GLAS_Apr2013_l0s454"
OUT = r"C:\Replicalm\test100526\centro454_bih50"


def rd(p, band=1):
    with rasterio.open(p) as r:
        z = r.read(band).astype("f8")
        nd, cc = r.nodata, abs(r.transform.a)
    if nd is not None:
        z[z == nd] = np.nan
    z[z > 1e30] = np.nan
    if band == 1 and z.dtype.kind == "f":
        z[z == 0.0] = np.nan
    return z, cc


A, ac = rd(ARCH_DEM)
Am = np.isfinite(A)
gA, _ = rd(ARCH_G1)
D, G = {}, {}
for tag, d in RUNS:
    D[tag] = rd(os.path.join(d, BASE + "_DEM.tif"))[0]
    with rasterio.open(os.path.join(d, "G1", BASE + "_G1.tif")) as r:
        G[tag] = r.read(1).astype("f4")
    print("%-11s %9d valid cells" % (tag, int(np.isfinite(D[tag]).sum())))
print("%-11s %9d valid cells" % ("archive", int(Am.sum())))


def detail(z, mask, scale):
    s = scale / ac / 2.0
    num = ndimage.gaussian_filter(np.where(mask, z, 0.0), s)
    den = ndimage.gaussian_filter(mask.astype("f8"), s)
    sm = np.where(den > 1e-6, num / np.maximum(den, 1e-6), 0.0)
    return np.where(mask, z - sm, np.nan)


def shape(mask):
    fill = ndimage.binary_fill_holes(mask)
    per = int((fill ^ ndimage.binary_erosion(fill)).sum())
    area = int(fill.sum())
    _l, n = ndimage.label(mask)
    runs = 0
    rows = int(mask.any(axis=1).sum())
    for r in np.flatnonzero(mask.any(axis=1)):
        d = np.diff(np.concatenate([[0], mask[r].astype(np.int8), [0]]))
        runs += int((d == 1).sum())
    return {"perimeter": per, "area": area,
            "raggedness": float(per / max(2.0 * np.sqrt(np.pi * area), 1e-9)),
            "pieces": int(n), "runs_per_row": float(runs / max(rows, 1))}


res = {"archive_edge": shape(Am), "runs": {}}
print("")
print("EDGE SHAPE")
print("%-11s %9s %10s %11s %8s %13s" % ("run", "area", "perimeter",
                                        "raggedness", "pieces", "runs/row"))
e = res["archive_edge"]
print("%-11s %9d %10d %11.2f %8d %13.2f"
      % ("archive", e["area"], e["perimeter"], e["raggedness"], e["pieces"],
         e["runs_per_row"]))
for tag, _ in RUNS:
    s = shape(np.isfinite(D[tag]))
    res["runs"].setdefault(tag, {})["edge"] = s
    print("%-11s %9d %10d %11.2f %8d %13.2f"
          % (tag, s["area"], s["perimeter"], s["raggedness"], s["pieces"],
             s["runs_per_row"]))

print("")
print("ELEVATION AGAINST THE ARCHIVE, by setback from the common boundary")
print("%-11s %9s %9s %9s %9s %9s"
      % ("run", "0 m RMSE", "25 m", "50 m", "25 m MAE", "25 m bias"))
for tag, _ in RUNS:
    common = Am & np.isfinite(D[tag])
    dist = ndimage.distance_transform_edt(common) * ac
    row = {}
    for sb in (0.0, 25.0, 50.0):
        m = common & (dist > sb)
        er = (D[tag] - A)[m]
        row["%.0f" % sb] = {"rmse": float(np.sqrt(np.mean(er ** 2))),
                            "mae": float(np.mean(np.abs(er))),
                            "bias": float(np.mean(er)),
                            "cells": int(m.sum())}
    res["runs"][tag]["elevation"] = row
    print("%-11s %9.4f %9.4f %9.4f %9.4f %+9.4f"
          % (tag, row["0"]["rmse"], row["25"]["rmse"], row["50"]["rmse"],
             row["25"]["mae"], row["25"]["bias"]))

print("")
print("TEXTURE AGAINST THE ARCHIVE, 25 m inside the boundary")
print("%-11s %8s %14s %14s" % ("run", "scale", "coherent ours/arch",
                               "speckle ours/arch"))
kk = max(3, int(round(1.0 / ac)) | 1)
for tag, _ in RUNS:
    common = Am & np.isfinite(D[tag])
    dist = ndimage.distance_transform_edt(common) * ac
    inner = common & (dist > 25.0)
    tex = {}
    for sc in (0.75, 1.5, 3.0, 6.0):
        dA_ = detail(A, common, sc)
        dO_ = detail(D[tag], common, sc)
        mA = ndimage.median_filter(np.nan_to_num(dA_), size=kk)
        mO = ndimage.median_filter(np.nan_to_num(dO_), size=kk)
        s = inner & np.isfinite(dA_) & np.isfinite(dO_)
        ca = float(np.sqrt(np.mean(mA[s] ** 2)))
        co = float(np.sqrt(np.mean(mO[s] ** 2)))
        sa = float(np.sqrt(np.mean((dA_ - mA)[s] ** 2)))
        so = float(np.sqrt(np.mean((dO_ - mO)[s] ** 2)))
        tex["%.2f" % sc] = {"coherent_ratio": co / max(ca, 1e-9),
                            "speckle_ratio": so / max(sa, 1e-9),
                            "arch_coherent": ca, "our_coherent": co}
        print("%-11s %8s %14.3f %14.3f"
              % (tag if sc == 0.75 else "", "%.2f m" % sc,
                 co / max(ca, 1e-9), so / max(sa, 1e-9)))
    res["runs"][tag]["texture"] = tex

with open(os.path.join(OUT, "spline_vs_biharmonic.json"), "w") as fh:
    json.dump(res, fh, indent=2)
print("")
print("wrote %s" % os.path.join(OUT, "spline_vs_biharmonic.json"))

PANELS = [("local spline", G["spline"]),
          ("global biharmonic", G["biharmonic"]),
          ("archive EdgeFixed", gA)]
H, W = gA.shape
win = 500
var = ndimage.uniform_filter((gA.astype("f4") -
                              ndimage.uniform_filter(gA.astype("f4"), 9)) ** 2,
                             65)
var = np.where(Am, var, -1)
picks, taken = [], np.zeros(gA.shape, bool)
for idx in np.argsort(var.ravel())[::-1][:400000]:
    r, c = divmod(int(idx), W)
    if r < win or r > H - win or c < win or c > W - win or taken[r, c]:
        continue
    taken[max(0, r - 900):r + 900, max(0, c - 900):c + 900] = True
    picks.append((r, c))
    if len(picks) == 3:
        break
for n, (r, c) in enumerate(picks, 1):
    sl = (slice(r - win // 2, r + win // 2), slice(c - win // 2, c + win // 2))
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 6.2))
    for ax, (t, img) in zip(axes, PANELS):
        ax.imshow(img[sl], cmap="gray", vmin=0, vmax=255,
                  interpolation="nearest")
        ax.set_title(t, fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
    fig.suptitle("l0s454 crop %d \u2014 250 m across at 0.5 m, 1:1 pixels"
                 % n, fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    p = os.path.join(OUT, "bih_crop%d.png" % n)
    fig.savefig(p, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("  wrote %s" % p)

ys, xs = np.where(Am ^ ndimage.binary_erosion(Am))
j = len(ys) // 2
r = int(np.clip(ys[j], win, H - win))
c = int(np.clip(xs[j], win, W - win))
sl = (slice(r - win // 2, r + win // 2), slice(c - win // 2, c + win // 2))
fig, axes = plt.subplots(1, 3, figsize=(16.5, 6.2))
for ax, (t, img) in zip(axes, PANELS):
    ax.imshow(img[sl], cmap="gray", vmin=0, vmax=255, interpolation="nearest")
    ax.set_title(t, fontsize=9)
    ax.set_xticks([])
    ax.set_yticks([])
fig.suptitle("l0s454 coverage edge \u2014 250 m across at 0.5 m", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.94])
p = os.path.join(OUT, "bih_edge.png")
fig.savefig(p, dpi=200, bbox_inches="tight")
plt.close(fig)
print("  wrote %s" % p)
