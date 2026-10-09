"""Second calibration pass. Three things the first pass showed were wrong.

The gate floor never bound: every value of gate_floor_above_m gave an
identical score, which means the computed gate was always wider than the
floor. So the floor was not the parameter doing the work, and the direction
the sweep wanted -- lower k_above was always better -- says the gate should be
narrower than anything the first sweep tried. This pass sweeps below 1.0 and
sweeps the slope allowance, which is the term that was holding the gate open.

Precision, not recall, was binding. Opening the gate took recall from 0.985 to
0.995 and made the surface worse every time. Admitted vegetation costs more
than recovered ground returns.

And the gate could only ever add. `accepted |= conf` forced every error the
opening made into the output, so the adjudication could recover ground but
never reject a wrong confident point -- which is the half that matters when
precision binds. Pruning is now switchable and swept.

Blocks are drawn at random from the archive's valid area with a fixed seed,
eight for development and six held out, rather than hand-picked, so the
reported validation number is not a selection of favourable terrain. The two
forensic sites are forced into the development set because they are the terrain
the routine has to handle.
"""
import itertools
import json
import os
import sys
import time

import numpy as np
import rasterio
from rasterio.warp import transform as warp_transform

sys.path.insert(0, r"C:\Replicalm\src")
from replicalm import classify
from replicalm import ground_tsr as TSR
from replicalm import surfer_algoA as ALGOA
from replicalm.grid import Grid
from replicalm.preprocess import dedupe_xy

RAW = (r"D:\GLiHT_LAS_orig\Yuc_Centro\Centro_Glas\Centro_Glas"
       r"\AMIGACarb_Yuc_Centro_GLAS_Apr2013_l0s441.las")
GND = (r"D:\_LAZ\Yuc_Centro\Centro_Glas\Centro_Glas"
       r"\Centro_GLAS_l0s441_GND_v1.laz")
REF = r"C:\Replicalm\test100526\terrascan_ceiling\dem_ts_c2c8.tif"
ARCH = (r"D:\_Archive_EdgeFixed\Yuc_Centro\Centro_Glas\Centro_Glas"
        r"\Centro_GLAS_l0s441_DEM_0p5m_v1.tif")
OUT = r"C:\Replicalm\test100526\ground_replication"
BLOCK_M, BUF_M, POOL = 240.0, 40.0, 0.25
N_DEV, N_VAL = 8, 6
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
A, _, _ = read(ARCH)
HT, WT = R.shape
OKR = np.isfinite(R) & np.isfinite(A)

print("reading clouds", flush=True)
t0 = time.time()
pl = classify._pdal().Pipeline(json.dumps({"pipeline": [
    {"type": "readers.las", "filename": RAW}]}))
pl.execute()
a = pl.arrays[0]
RX = a["X"].astype("f8")
RY = a["Y"].astype("f8")
RZ = a["Z"].astype("f8")
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
print("raw %d, delivered union %d, %.0f s"
      % (len(RZ), len(TZ), time.time() - t0), flush=True)

FORENSIC = [("mound", -88.9673900, 19.2242934),
            ("depression", -88.9664570, 19.2223460)]
blocks = []
for nm, lo, la in FORENSIC:
    xs, ys = warp_transform("EPSG:4326", "EPSG:32616", [lo], [la])
    blocks.append((nm, xs[0], ys[0]))
rng = np.random.default_rng(11)
rr, cc_ = np.nonzero(OKR)
taken = set()
for i in rng.permutation(len(rr)):
    if len(blocks) >= N_DEV + N_VAL:
        break
    cx = at.c + cc_[i] * ac
    cy = at.f - rr[i] * ac
    key = (int(cx // BLOCK_M), int(cy // BLOCK_M))
    if key in taken:
        continue
    n = int(round(BLOCK_M / ac))
    c0 = int(round((cx - BLOCK_M / 2 - at.c) / ac))
    r0 = int(round((at.f - (cy + BLOCK_M / 2)) / ac))
    if c0 < 0 or r0 < 0 or c0 + n > WT or r0 + n > HT:
        continue
    if OKR[r0:r0 + n, c0:c0 + n].mean() < 0.80:
        continue
    taken.add(key)
    blocks.append(("blk%02d" % len(blocks), cx, cy))
DEV, VAL = blocks[:N_DEV], blocks[N_DEV:]
print("development %s" % ", ".join(b[0] for b in DEV))
print("held out    %s" % ", ".join(b[0] for b in VAL))

CACHE = {}


def prep(tag, cx, cy):
    if tag in CACHE:
        return CACHE[tag]
    m = ((RX >= cx - BLOCK_M / 2 - BUF_M) & (RX <= cx + BLOCK_M / 2 + BUF_M)
         & (RY >= cy - BLOCK_M / 2 - BUF_M) & (RY <= cy + BLOCK_M / 2 + BUF_M))
    CACHE[tag] = (RX[m], RY[m], RZ[m])
    return CACHE[tag]


def score(gx, gy, gz, cx, cy):
    n = int(round(BLOCK_M / ac))
    c0 = int(round((cx - BLOCK_M / 2 - at.c) / ac))
    r0 = int(round((at.f - (cy + BLOCK_M / 2)) / ac))
    gr = Grid(origin_x=at.c + c0 * ac, origin_y=at.f - r0 * ac, cell=ac,
              width=n, height=n)
    dx, dy, dz, _ = dedupe_xy(gx, gy, gz, keep="lowest")
    if len(dz) < 40:
        return None
    z, _i = ALGOA.surfer_grid(dx, dy, dz, gr, pool_sigma_m=POOL,
                              verbose=False)
    rr_ = R[r0:r0 + n, c0:c0 + n]
    aa = A[r0:r0 + n, c0:c0 + n]
    ok = np.isfinite(z) & np.isfinite(rr_)
    oa = np.isfinite(z) & np.isfinite(aa)
    if ok.sum() < 100:
        return None
    e = (z - rr_)[ok]
    ea = (z - aa)[oa] if oa.sum() > 100 else np.array([np.nan])
    return dict(rmse_ref=float(np.sqrt(np.mean(e ** 2))),
                mae_ref=float(np.mean(np.abs(e))),
                bias_ref=float(np.mean(e)),
                p99_ref=float(np.percentile(np.abs(e), 99)),
                rmse_arch=float(np.sqrt(np.nanmean(ea ** 2))))


def agreement(gx, gy, cx, cy):
    from scipy.spatial import cKDTree
    mt = ((TX >= cx - BLOCK_M / 2) & (TX <= cx + BLOCK_M / 2)
          & (TY >= cy - BLOCK_M / 2) & (TY <= cy + BLOCK_M / 2))
    mg = ((gx >= cx - BLOCK_M / 2) & (gx <= cx + BLOCK_M / 2)
          & (gy >= cy - BLOCK_M / 2) & (gy <= cy + BLOCK_M / 2))
    if mt.sum() < 100 or mg.sum() < 100:
        return {}
    d1, _ = cKDTree(np.column_stack([gx[mg], gy[mg]])).query(
        np.column_stack([TX[mt], TY[mt]]), workers=-1)
    d2, _ = cKDTree(np.column_stack([TX[mt], TY[mt]])).query(
        np.column_stack([gx[mg], gy[mg]]), workers=-1)
    rec = float((d1 <= 0.25).mean())
    pre = float((d2 <= 0.25).mean())
    return {"recall": rec, "precision": pre,
            "f1": float(2 * pre * rec / max(pre + rec, 1e-9))}


def run_config(blks, over, label):
    base = TSR.GroundTSRParams()
    p = TSR.GroundTSRParams(**{k: over.get(k, getattr(base, k))
                               for k in base.__dataclass_fields__})
    rows = []
    for tag, cx, cy in blks:
        bx, by, bz = prep(tag, cx, cy)
        try:
            mask, rep = TSR.classify_ground(bx, by, bz, params=p)
        except TSR.GroundTSRError:
            continue
        s = score(bx[mask], by[mask], bz[mask], cx, cy)
        if s is None:
            continue
        s.update(agreement(bx[mask], by[mask], cx, cy))
        s["block"] = tag
        s["ground"] = int(mask.sum())
        rows.append(s)
    if not rows:
        return None
    g = lambda k: float(np.mean([r[k] for r in rows if k in r]))
    return {"label": label, "params": over, "n": len(rows),
            "rmse_ref": g("rmse_ref"), "rmse_arch": g("rmse_arch"),
            "mae_ref": g("mae_ref"), "bias_ref": g("bias_ref"),
            "p99_ref": g("p99_ref"), "f1": g("f1"), "recall": g("recall"),
            "precision": g("precision"),
            "worst_ref": float(np.max([r["rmse_ref"] for r in rows])),
            "blocks": rows}


HDR = ("%-32s %9s %11s %9s %8s %9s %9s"
       % ("configuration", "RMSE ref", "RMSE arch", "worst", "F1", "recall",
          "precis"))
results = []
print("")
print("STAGE CONTRIBUTION, %d development blocks" % len(DEV))
print(HDR)
print("-" * 92)
for label, over in [
        ("opening alone", {"gate_k_above": 0.0, "gate_k_below": 0.0,
                           "gate_floor_above_m": 0.0,
                           "gate_floor_below_m": 0.0,
                           "slope_allowance": 0.0, "refine_iterations": 0,
                           "prune_confident": False}),
        ("curvature alone, ungated", {"gate_ceiling_m": 999.0,
                                      "gate_floor_above_m": 999.0,
                                      "gate_floor_below_m": 999.0,
                                      "refine_iterations": 0,
                                      "prune_confident": False}),
        ("adjudicated, additive only", {"prune_confident": False}),
        ("adjudicated, pruning", {"prune_confident": True})]:
    agg = run_config(DEV, over, label)
    if agg is None:
        continue
    results.append(agg)
    print("%-32s %9.4f %11.4f %9.4f %8.3f %9.3f %9.3f"
          % (label[:32], agg["rmse_ref"], agg["rmse_arch"], agg["worst_ref"],
             agg["f1"], agg["recall"], agg["precision"]), flush=True)

print("")
print("GATE SWEEP, k_above / slope_allowance / pruning")
print(HDR)
print("-" * 92)
best, best_rmse = None, 1e9
for ka, sa, pr in itertools.product((0.25, 0.5, 0.75, 1.0),
                                    (0.0, 0.3, 0.6), (True, False)):
    over = {"gate_k_above": ka, "slope_allowance": sa,
            "prune_confident": pr, "gate_k_below": 4.0,
            "gate_floor_above_m": 0.05}
    lbl = "%.2f / %.1f / %s" % (ka, sa, "prune" if pr else "add")
    agg = run_config(DEV, over, lbl)
    if agg is None:
        continue
    results.append(agg)
    mark = ""
    if agg["rmse_ref"] < best_rmse:
        best_rmse, best, mark = agg["rmse_ref"], over, "  <-"
    print("%-32s %9.4f %11.4f %9.4f %8.3f %9.3f %9.3f%s"
          % (lbl, agg["rmse_ref"], agg["rmse_arch"], agg["worst_ref"],
             agg["f1"], agg["recall"], agg["precision"], mark), flush=True)

print("")
print("best on development: %s" % json.dumps(best))
print("  RMSE ref %.4f" % best_rmse)
val = run_config(VAL, best, "held out")
if val:
    print("")
    print("HELD OUT, %d blocks never used for selection" % val["n"])
    print("  RMSE ref %.4f, RMSE arch %.4f, worst %.4f, MAE %.4f, "
          "bias %+.4f" % (val["rmse_ref"], val["rmse_arch"],
                          val["worst_ref"], val["mae_ref"], val["bias_ref"]))
    print("  F1 %.3f, recall %.3f, precision %.3f"
          % (val["f1"], val["recall"], val["precision"]))
    for b in val["blocks"]:
        print("    %-8s RMSE ref %.4f, arch %.4f, ground %d"
              % (b["block"], b["rmse_ref"], b["rmse_arch"], b["ground"]))

with open(os.path.join(OUT, "tsr_calibration2.json"), "w") as fh:
    json.dump({"development": results, "best": best, "held_out": val,
               "dev_blocks": [b[0] for b in DEV],
               "val_blocks": [b[0] for b in VAL]}, fh, indent=2)
print("")
print("wrote %s" % os.path.join(OUT, "tsr_calibration2.json"))
