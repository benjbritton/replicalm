"""Three panels for l0s454 and the measurements behind them.

Panels, all on the archive's own 0.5 m lattice and all shown exactly as the
files are -- no restretching here at all:

  left    my earlier scratch render: SMRF ground with no cleanup, and a
          stretch shared across three surfaces. Kept in the figure because it
          is what was wrong, not because it is a candidate.
  middle  the Clear profile through replicalm.pipeline.process, G1 by
          visualize.g1 with its own per-layer per-tile normalisation.
  right   the delivered archive G1.

The predictions being tested are specific, so the measurements are too.
Smearing and cleaning look alike by eye and separate numerically: a cleaner
surface loses speckle and keeps coherent detail, a smeared one loses both.
Lost ground shows up as coverage, not as texture. Edge quality is the shape of
the no-data boundary, which is measurable as the ratio of its perimeter to the
perimeter of an equally large smooth region.
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
CLEAR_DEM = (r"C:\Replicalm\test100526\centro454_clear50"
             r"\AMIGACarb_Yuc_Centro_GLAS_Apr2013_l0s454_DEM.tif")
CLEAR_G1 = (r"C:\Replicalm\test100526\centro454_clear50\G1"
            r"\AMIGACarb_Yuc_Centro_GLAS_Apr2013_l0s454_G1.tif")
OLD_G1 = (r"C:\Replicalm\test100526\centro454_050"
          r"\G1_Centro_GLAS_l0s454_SMRF_0p5m.tif")
OUT = r"C:\Replicalm\test100526\centro454_clear50"


def dem(p):
    with rasterio.open(p) as r:
        z = r.read(1).astype("f8")
        nd, tr, cc = r.nodata, r.transform, abs(r.transform.a)
    if nd is not None:
        z[z == nd] = np.nan
    z[z > 1e30] = np.nan
    z[z == 0.0] = np.nan
    return z, tr, cc


def g1(p):
    with rasterio.open(p) as r:
        b = r.read(1).astype("f4")
        h, w = r.height, r.width
    return b, (h, w)


A, at, ac = dem(ARCH_DEM)
C, _, _ = dem(CLEAR_DEM)
gA, shA = g1(ARCH_G1)
gC, shC = g1(CLEAR_G1)
gO, shO = g1(OLD_G1)
print("grids: archive G1 %s, clear G1 %s, old G1 %s" % (shA, shC, shO))
Am, Cm = np.isfinite(A), np.isfinite(C)
common = Am & Cm
print("archive cells %d, clear cells %d, common %d"
      % (int(Am.sum()), int(Cm.sum()), int(common.sum())))

res = {"cells": {"archive": int(Am.sum()), "clear": int(Cm.sum()),
                 "common": int(common.sum())}}

# ---- coverage: did the cleanup or the trim cost ground?
only_a = Am & ~Cm
only_c = Cm & ~Am
print("")
print("COVERAGE")
print("  archive has, we do not   %8d cells  %5.2f%% of the archive"
      % (int(only_a.sum()), 100.0 * only_a.sum() / Am.sum()))
print("  we have, archive does not %8d cells  %5.2f%% of ours"
      % (int(only_c.sum()), 100.0 * only_c.sum() / Cm.sum()))
res["coverage"] = {"archive_only": int(only_a.sum()),
                   "ours_only": int(only_c.sum())}

# ---- edge shape
def edge_shape(mask, name):
    fill = ndimage.binary_fill_holes(mask)
    per = int((fill ^ ndimage.binary_erosion(fill)).sum())
    area = int(fill.sum())
    circ = per / max(2.0 * np.sqrt(np.pi * area), 1e-9)
    holes = int(fill.sum() - mask.sum())
    lab, n = ndimage.label(~mask & fill)
    print("  %-8s perimeter %8d cells, area %9d, raggedness %6.2f, "
          "interior holes %6d in %d patches"
          % (name, per, area, circ, holes, n))
    return {"perimeter": per, "area": area, "raggedness": float(circ),
            "interior_hole_cells": holes, "interior_hole_patches": int(n)}


print("")
print("EDGE SHAPE  (raggedness 1.0 = a perfect circle; higher = more "
      "convoluted boundary for the same area)")
res["edge"] = {"archive": edge_shape(Am, "archive"),
               "clear": edge_shape(Cm, "clear")}

# ---- elevation agreement, with setbacks because the edge dominates
dist = ndimage.distance_transform_edt(common) * ac
print("")
print("ELEVATION AGAINST THE ARCHIVE")
print("%-9s %10s %8s %9s %9s %9s" % ("setback", "cells", "kept", "RMSE",
                                     "MAE", "bias"))
res["elevation"] = {}
for sb in (0.0, 10.0, 25.0, 50.0):
    m = common & (dist > sb)
    if m.sum() < 1000:
        continue
    e = (C - A)[m]
    row = {"cells": int(m.sum()),
           "rmse": float(np.sqrt(np.mean(e ** 2))),
           "mae": float(np.mean(np.abs(e))), "bias": float(np.mean(e))}
    res["elevation"]["%.0f" % sb] = row
    print("%-9s %10d %7.1f%% %9.4f %9.4f %+9.4f"
          % ("%.0f m" % sb, row["cells"], 100.0 * m.sum() / common.sum(),
             row["rmse"], row["mae"], row["bias"]))

# ---- texture: smeared or merely cleaner?
def detail(z, mask, scale):
    s = scale / ac / 2.0
    num = ndimage.gaussian_filter(np.where(mask, z, 0.0), s)
    den = ndimage.gaussian_filter(mask.astype("f8"), s)
    sm = np.where(den > 1e-6, num / np.maximum(den, 1e-6), 0.0)
    return np.where(mask, z - sm, np.nan)


inner = common & (dist > 25.0)
kk = max(3, int(round(1.0 / ac)) | 1)
print("")
print("TEXTURE, 25 m inside the boundary, by detail scale")
print("%-10s %11s %11s %11s %11s" % ("scale", "arch coherent", "our coherent",
                                     "arch speckle", "our speckle"))
res["texture"] = {}
for sc in (0.75, 1.5, 3.0, 6.0):
    dA_ = detail(A, common, sc)
    dC_ = detail(C, common, sc)
    mA = ndimage.median_filter(np.nan_to_num(dA_), size=kk)
    mC = ndimage.median_filter(np.nan_to_num(dC_), size=kk)
    s = inner & np.isfinite(dA_) & np.isfinite(dC_)
    row = {"arch_coherent": float(np.sqrt(np.mean(mA[s] ** 2))),
           "our_coherent": float(np.sqrt(np.mean(mC[s] ** 2))),
           "arch_speckle": float(np.sqrt(np.mean((dA_ - mA)[s] ** 2))),
           "our_speckle": float(np.sqrt(np.mean((dC_ - mC)[s] ** 2)))}
    row["coherent_ratio"] = row["our_coherent"] / max(row["arch_coherent"],
                                                      1e-9)
    row["speckle_ratio"] = row["our_speckle"] / max(row["arch_speckle"], 1e-9)
    res["texture"]["%.2f" % sc] = row
    print("%-10s %11.5f %11.5f %11.5f %11.5f"
          % ("%.2f m" % sc, row["arch_coherent"], row["our_coherent"],
             row["arch_speckle"], row["our_speckle"]))
print("")
print("%-10s %16s %16s" % ("scale", "coherent ours/arch", "speckle ours/arch"))
for k, v in res["texture"].items():
    print("%-10s %16.3f %16.3f" % ("%s m" % k, v["coherent_ratio"],
                                   v["speckle_ratio"]))

with open(os.path.join(OUT, "comparison_metrics.json"), "w") as fh:
    json.dump(res, fh, indent=2)
print("")
print("wrote %s" % os.path.join(OUT, "comparison_metrics.json"))

# ---- panels
H, W = gA.shape
PANELS = [("my earlier render\nno cleanup, shared stretch", gO),
          ("Clear through the pipeline\nstock G1 normalisation", gC),
          ("archive EdgeFixed", gA)]


def save(fig, name):
    p = os.path.join(OUT, name)
    fig.savefig(p, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("  wrote %s" % p)


K = 8
def shrink(x, k=K):
    h, w = x.shape
    x = x[:h // k * k, :w // k * k].astype("f4")
    return x.reshape(h // k, k, w // k, k).mean(axis=(1, 3))


fig, axes = plt.subplots(1, 3, figsize=(15.0, 20.0))
for ax, (t, img) in zip(axes, PANELS):
    ax.imshow(shrink(img), cmap="gray", vmin=0, vmax=255,
              interpolation="nearest")
    ax.set_title(t, fontsize=10)
    ax.set_xticks([])
    ax.set_yticks([])
fig.suptitle("Centro GLAS l0s454, full tile at 0.5 m \u2014 panels exactly as "
             "the files are", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.97])
save(fig, "compare_fulltile.png")

# crops where the archive carries the most local detail
win = 500
var = ndimage.uniform_filter((gA.astype("f4") -
                              ndimage.uniform_filter(gA.astype("f4"), 9)) ** 2,
                             65)
var = np.where(Am, var, -1)
picks, taken = [], np.zeros(gA.shape, bool)
flat = np.argsort(var.ravel())[::-1]
for idx in flat[:400000]:
    r, c = divmod(int(idx), W)
    if r < win or r > H - win or c < win or c > W - win:
        continue
    if taken[r, c]:
        continue
    taken[max(0, r - 900):r + 900, max(0, c - 900):c + 900] = True
    picks.append((r, c))
    if len(picks) == 3:
        break
print("")
print("detail crops at rows/cols: %s" % picks)
for n, (r, c) in enumerate(picks, 1):
    sl = (slice(r - win // 2, r + win // 2), slice(c - win // 2, c + win // 2))
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 6.2))
    for ax, (t, img) in zip(axes, PANELS):
        ax.imshow(img[sl], cmap="gray", vmin=0, vmax=255,
                  interpolation="nearest")
        ax.set_title(t, fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
    fig.suptitle("l0s454 detail crop %d \u2014 %.0f m across at 0.5 m, "
                 "1:1 pixels" % (n, win * ac), fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    save(fig, "compare_crop%d.png" % n)

# an edge crop, since the trim is the stated concern
ys, xs = np.where(Am ^ ndimage.binary_erosion(Am))
j = len(ys) // 2
r, c = int(ys[j]), int(xs[j])
r = int(np.clip(r, win, H - win))
c = int(np.clip(c, win, W - win))
sl = (slice(r - win // 2, r + win // 2), slice(c - win // 2, c + win // 2))
fig, axes = plt.subplots(1, 3, figsize=(16.5, 6.2))
for ax, (t, img) in zip(axes, PANELS):
    ax.imshow(img[sl], cmap="gray", vmin=0, vmax=255,
              interpolation="nearest")
    ax.set_title(t, fontsize=9)
    ax.set_xticks([])
    ax.set_yticks([])
fig.suptitle("l0s454 coverage edge \u2014 %.0f m across at 0.5 m, 1:1 pixels"
             % (win * ac), fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.94])
save(fig, "compare_edge.png")
