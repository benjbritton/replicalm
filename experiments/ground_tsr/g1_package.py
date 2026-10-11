"""G1 inspection package for the forensic mound and depression.

The question this is built to answer is whether the terrain supports the
current account of the failure: that the stage-4 reference surface is too
smooth to represent a mound, so on curved landforms the gate adjudicates
against a surface already wrong by more than the gate width.

Everything is rendered on one spatial extent per site, one elevation scale,
one difference scale, and one set of hillshade parameters.

One deliberate departure from `visualize.g1`. That function normalises each
G1 layer by its own min and max, which is right for a standalone tile and
wrong here: it would give every panel its own stretch and the panels would not
be comparable. Each layer is instead normalised against a range shared by all
four surfaces at a site, so a difference between two panels is a difference in
the terrain rather than in the stretch.

The stage-4 surface is sampled onto the output grid by nearest neighbour,
which is what the gate itself does, so the panel shows the surface the
algorithm actually measures against rather than a tidier version of it.
"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from rasterio.warp import transform as warp_transform
from scipy.spatial import cKDTree

sys.path.insert(0, r"C:\Replicalm\src")
from replicalm import classify
from replicalm import ground_tsr as TSR
from replicalm import visualize as VIS
from replicalm.preprocess import dedupe_xy

RAW = (r"D:\GLiHT_LAS_orig\Yuc_Centro\Centro_Glas\Centro_Glas"
       r"\AMIGACarb_Yuc_Centro_GLAS_Apr2013_l0s441.las")
GND = (r"D:\_LAZ\Yuc_Centro\Centro_Glas\Centro_Glas"
       r"\Centro_GLAS_l0s441_GND_v1.laz")
REF = r"C:\Replicalm\test100526\terrascan_ceiling\dem_ts_c2c8.tif"
SMRF = r"C:\Replicalm\test100526\cellsize_validation\dem_B_cell1.00.tif"
TSRD = r"C:\Replicalm\test100526\ground_replication\dem_pin_adj.tif"
OUT = r"C:\Replicalm\test100526\g1_inspection"
WIN_M, BUF_M, NEAR = 120.0, 60.0, 0.25
SITES = [("mound", -88.9673900, 19.2242934),
         ("depression", -88.9664570, 19.2223460)]
os.makedirs(OUT, exist_ok=True)


def read(path):
    with rasterio.open(path) as r:
        z = r.read(1).astype("f8")
        nd, tr, cc = r.nodata, r.transform, abs(r.transform.a)
    if nd is not None:
        z[z == nd] = np.nan
    z[z > 1e30] = np.nan
    z[z == 0.0] = np.nan
    return z, tr, cc


R, at, ac = read(REF)
S, _, _ = read(SMRF)
T, _, _ = read(TSRD)
HT, WT = R.shape
print("grids loaded, cell %.3f m" % ac, flush=True)

print("reading clouds", flush=True)
pl = classify._pdal().Pipeline(json.dumps({"pipeline": [
    {"type": "readers.las", "filename": RAW}]}))
pl.execute()
a = pl.arrays[0]
RX, RY, RZ = (a["X"].astype("f8"), a["Y"].astype("f8"), a["Z"].astype("f8"))
del a, pl
pl = classify._pdal().Pipeline(json.dumps({"pipeline": [
    {"type": "readers.las", "filename": GND, "nosrs": True},
    {"type": "filters.expression",
     "expression": "Classification == 2 || Classification == 8"}]}))
pl.execute()
a = pl.arrays[0]
TX, TY, TZ, _ = dedupe_xy(a["X"].astype("f8"), a["Y"].astype("f8"),
                          a["Z"].astype("f8"), keep="lowest")
del a, pl
print("raw %d, delivered union %d" % (len(RZ), len(TZ)), flush=True)

PTS = [("source returns", "src"), ("confident set", "conf"),
       ("candidate set", "cand"), ("final accepted set", "adj"),
       ("TerraScan class 2 + 8", "ref")]
SURF = [("stage-4 reference surface", "stage4"), ("SMRF surface", "smrf"),
        ("ground_tsr surface", "tsr"),
        ("TerraScan 2 + 8 surface", "refsurf")]
DIFF = [("stage-4 minus reference", "stage4"), ("SMRF minus reference", "smrf"),
        ("ground_tsr minus reference", "tsr")]
metrics = {}

for name, lo, la in SITES:
    xs, ys = warp_transform("EPSG:4326", "EPSG:32616", [lo], [la])
    cx, cy = xs[0], ys[0]
    n = int(round(WIN_M / ac))
    c0 = int(round((cx - WIN_M / 2 - at.c) / ac))
    r0 = int(round((at.f - (cy + WIN_M / 2)) / ac))
    c0 = int(np.clip(c0, 0, WT - n))
    r0 = int(np.clip(r0, 0, HT - n))
    ext = [at.c + c0 * ac, at.c + (c0 + n) * ac,
           at.f - (r0 + n) * ac, at.f - r0 * ac]
    print("")
    print("%s: window %.0f m at %.1f, %.1f" % (name, WIN_M, cx, cy),
          flush=True)

    mb = ((RX >= ext[0] - BUF_M) & (RX <= ext[1] + BUF_M)
          & (RY >= ext[2] - BUF_M) & (RY <= ext[3] + BUF_M))
    bx, by, bz = RX[mb], RY[mb], RZ[mb]
    p = TSR.GroundTSRParams()
    mask, rep = TSR.classify_ground(bx, by, bz, params=p,
                                    adapt_density=False, keep_masks=True)
    conf = rep.pop("confident_mask")
    cand = rep.pop("candidate_mask")
    print("  returns %d, confident %d, candidates %d, accepted %d"
          % (len(bz), conf.sum(), cand.sum(), mask.sum()), flush=True)

    # stage 4 exactly as the gate sees it
    S4, geom = TSR.reference_surface(bx, by, bz, conf, p)
    cols = np.arange(n)
    rows = np.arange(n)
    gxm, gym = np.meshgrid(at.c + (c0 + cols + 0.5) * ac,
                           at.f - (r0 + rows + 0.5) * ac)
    stage4 = TSR._sample(S4, geom, gxm.ravel(), gym.ravel()).reshape(n, n)

    sl = (slice(r0, r0 + n), slice(c0, c0 + n))
    surf = {"stage4": stage4, "smrf": S[sl], "tsr": T[sl], "refsurf": R[sl]}
    valid = np.isfinite(surf["refsurf"])
    for k in surf:
        surf[k] = np.where(valid, surf[k], np.nan)

    inw = lambda X, Y: ((X >= ext[0]) & (X <= ext[1])
                        & (Y >= ext[2]) & (Y <= ext[3]))
    mt = inw(TX, TY)
    pops = {"src": (bx, by, bz, inw(bx, by)),
            "conf": (bx, by, bz, conf & inw(bx, by)),
            "cand": (bx, by, bz, cand & inw(bx, by)),
            "adj": (bx, by, bz, mask & inw(bx, by)),
            "ref": (TX, TY, TZ, mt)}

    # ---- metrics
    met = {"window_m": WIN_M, "centre_lonlat": [lo, la],
           "centre_utm": [float(cx), float(cy)],
           "cells_scored": int(valid.sum()),
           "density_per_m2": rep["density_per_m2"], "populations": {},
           "surfaces": {}}
    tref = np.column_stack([TX[mt], TY[mt]])
    tree_ref = cKDTree(tref) if mt.sum() > 10 else None
    for key, (X, Y, Z, m) in pops.items():
        d = {"points": int(m.sum()),
             "per_m2": float(m.sum() / (WIN_M * WIN_M))}
        if tree_ref is not None and m.sum() > 10 and key != "ref":
            own = np.column_stack([X[m], Y[m]])
            d1, _ = cKDTree(own).query(tref, workers=-1)
            d2, _ = tree_ref.query(own, workers=-1)
            rec = float((d1 <= NEAR).mean())
            pre = float((d2 <= NEAR).mean())
            d.update(recall=rec, precision=pre,
                     f1=float(2 * pre * rec / max(pre + rec, 1e-9)))
        met["populations"][key] = d
    for key in surf:
        if key == "refsurf":
            continue
        e = (surf[key] - surf["refsurf"])[valid]
        e = e[np.isfinite(e)]
        if len(e) < 10:
            continue
        met["surfaces"][key] = {
            "rmse_vs_reference": float(np.sqrt(np.mean(e ** 2))),
            "mae": float(np.mean(np.abs(e))),
            "bias": float(np.mean(e)),
            "p99_abs": float(np.percentile(np.abs(e), 99)),
            "max_abs": float(np.max(np.abs(e)))}
    metrics[name] = met
    for k, v in met["surfaces"].items():
        print("    %-8s vs reference: RMSE %.4f, bias %+.4f, p99 %.4f, "
              "max %.4f" % (k, v["rmse_vs_reference"], v["bias"],
                            v["p99_abs"], v["max_abs"]), flush=True)

    # ---- G1 with a shared stretch per layer
    layers = {}
    for key in surf:
        dem = surf[key]
        fill = float(np.nanmedian(dem)) if np.isfinite(dem).any() else 0.0
        d = np.where(np.isfinite(dem), dem, fill)
        # rvt trims the border on some layers, so every layer is centred
        # back into the DEM's shape -- the same thing visualize.g1 does
        layers[key] = {nm: VIS._centre_into(VIS.compute_layer(nm, d, ac),
                                            dem.shape)
                       for nm in VIS.LAYER_NAMES}
    shared = {}
    for nm in VIS.LAYER_NAMES:
        allv = np.concatenate([layers[k][nm][np.isfinite(layers[k][nm])]
                               .ravel() for k in layers])
        shared[nm] = (float(np.nanpercentile(allv, 1)),
                      float(np.nanpercentile(allv, 99)))
    g1 = {}
    for key in surf:
        canvas = None
        for nm, transp, mode, br, co, ga in VIS.G1_LAYERS:
            lo_, hi_ = shared[nm]
            arr = np.clip((layers[key][nm] - lo_) / max(hi_ - lo_, 1e-9),
                          0.0, 1.0).astype(np.float32)
            arr = VIS.adjust(arr, br, co, ga)
            canvas = arr if canvas is None else VIS.blend(
                canvas, arr, mode, 1.0 - transp / 100.0)
        g1[key] = np.where(valid, canvas, np.nan)

    zall = np.concatenate([Z[m] for _k, (_X, _Y, Z, m) in pops.items()
                           if m.sum() > 0])
    zlo, zhi = (float(np.percentile(zall, 1)), float(np.percentile(zall, 99)))
    dall = np.concatenate([(surf[k] - surf["refsurf"])[valid].ravel()
                           for k in ("stage4", "smrf", "tsr")])
    dall = dall[np.isfinite(dall)]
    dv = float(np.percentile(np.abs(dall), 98))

    for cm in ("terrain", "RdBu_r", "gray", "viridis"):
        matplotlib.colormaps[cm].set_bad("white")

    fig = plt.figure(figsize=(21.0, 15.5))
    gs = fig.add_gridspec(3, 5, hspace=0.16, wspace=0.06)
    for i, (lbl, key) in enumerate(PTS):
        ax = fig.add_subplot(gs[0, i])
        X, Y, Z, m = pops[key]
        ax.imshow(np.zeros((n, n)), cmap="gray", vmin=0, vmax=1,
                  extent=ext, origin="upper")
        ax.scatter(X[m], Y[m], c=Z[m], s=0.25, cmap="terrain", vmin=zlo,
                   vmax=zhi, linewidths=0, marker=",")
        d = met["populations"][key]
        sub = "%d pts, %.1f/m2" % (d["points"], d["per_m2"])
        if "f1" in d:
            sub += "\nrecall %.3f, precision %.3f" % (d["recall"],
                                                      d["precision"])
        ax.set_title("%s\n%s" % (lbl, sub), fontsize=9)
        ax.set_xlim(ext[0], ext[1])
        ax.set_ylim(ext[2], ext[3])
        ax.set_xticks([])
        ax.set_yticks([])
    for i, (lbl, key) in enumerate(SURF):
        ax = fig.add_subplot(gs[1, i])
        ax.imshow(g1[key], extent=ext, origin="upper",
                  interpolation="nearest")
        ax.set_title("%s\nG1, shared stretch" % lbl, fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
    ax = fig.add_subplot(gs[1, 4])
    ax.imshow(surf["refsurf"], cmap="terrain", vmin=zlo, vmax=zhi,
              extent=ext, origin="upper", interpolation="nearest")
    ax.set_title("reference elevation\n%.1f to %.1f m" % (zlo, zhi),
                 fontsize=9)
    ax.set_xticks([])
    ax.set_yticks([])
    for i, (lbl, key) in enumerate(DIFF):
        ax = fig.add_subplot(gs[2, i])
        im = ax.imshow(surf[key] - surf["refsurf"], cmap="RdBu_r",
                       vmin=-dv, vmax=dv, extent=ext, origin="upper",
                       interpolation="nearest")
        v = met["surfaces"][key]
        ax.set_title("%s\nRMSE %.4f, bias %+.4f, max %.3f"
                     % (lbl, v["rmse_vs_reference"], v["bias"],
                        v["max_abs"]), fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
        if i == 2:
            fig.colorbar(im, ax=ax, fraction=0.046, label="metres")
    ax = fig.add_subplot(gs[2, 3:])
    ax.axis("off")
    lines = ["%s, %.0f m window at %.7f, %.7f" % (name, WIN_M, lo, la),
             "cell %.3f m, %d cells scored, %.1f returns/m2"
             % (ac, met["cells_scored"], met["density_per_m2"]), "",
             "%-22s %9s %9s %9s %9s" % ("surface vs reference", "RMSE",
                                        "bias", "p99|e|", "max|e|")]
    for lbl, key in DIFF:
        v = met["surfaces"][key]
        lines.append("%-22s %9.4f %+9.4f %9.4f %9.4f"
                     % (lbl.replace(" minus reference", ""),
                        v["rmse_vs_reference"], v["bias"], v["p99_abs"],
                        v["max_abs"]))
    lines += ["", "%-22s %10s %8s %10s %9s"
              % ("population", "points", "per m2", "recall", "precision")]
    for lbl, key in PTS:
        d = met["populations"][key]
        lines.append("%-22s %10d %8.1f %10s %9s"
                     % (lbl, d["points"], d["per_m2"],
                        "%.3f" % d["recall"] if "recall" in d else "-",
                        "%.3f" % d["precision"] if "precision" in d else "-"))
    ax.text(0.0, 1.0, "\n".join(lines), va="top", ha="left", fontsize=9,
            family="monospace", transform=ax.transAxes)
    fig.suptitle("Centro GLAS l0s441 \u2014 %s, identical extent, elevation "
                 "scale, difference scale and hillshade parameters" % name,
                 fontsize=13)
    pp = os.path.join(OUT, "g1_%s.png" % name)
    fig.savefig(pp, dpi=170, bbox_inches="tight")
    plt.close(fig)
    print("  wrote %s" % pp, flush=True)

with open(os.path.join(OUT, "g1_metrics.json"), "w") as fh:
    json.dump(metrics, fh, indent=2)
print("")
print("wrote %s" % os.path.join(OUT, "g1_metrics.json"))
