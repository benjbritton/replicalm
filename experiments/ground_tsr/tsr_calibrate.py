"""Calibrate ground_tsr against the reference surface. Construction, not survey.

The objective is the brief's: DEM RMSE after biharmonic interpolation. Both
sides of the comparison go through the same surface engine, so the score moves
only with classification.

Two references are carried because they answer different questions. The
reference surface -- the delivered class 2 + 8 population through the same
biharmonic engine -- isolates classification, and a routine that reproduced the
delivered population exactly would score zero against it. The archive is the
scale the 0.052 m benchmark was quoted on, so it is reported alongside to say
how close the whole chain is.

Development blocks are disjoint from the validation blocks. Parameters are
chosen on the first set and the chosen configuration is scored once on the
second, so the reported number is not the number that was optimised.
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
os.makedirs(OUT, exist_ok=True)

DEV = [("mound", -88.9673900, 19.2242934),
       ("depression", -88.9664570, 19.2223460),
       ("sparse", -88.9693115, 19.2046148),
       ("sparse2", -88.9691149, 19.1978981)]
VAL = [("val1", -88.9695692, 19.2062540),
       ("val2", -88.9664570, 19.2142000),
       ("val3", -88.9680000, 19.2180000)]


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

print("reading raw cloud", flush=True)
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


def utm(lo, la):
    xs, ys = warp_transform("EPSG:4326", "EPSG:32616", [lo], [la])
    return xs[0], ys[0]


def block(cx, cy, pad):
    m = ((RX >= cx - BLOCK_M / 2 - pad) & (RX <= cx + BLOCK_M / 2 + pad)
         & (RY >= cy - BLOCK_M / 2 - pad) & (RY <= cy + BLOCK_M / 2 + pad))
    return RX[m], RY[m], RZ[m]


def dem_from(gx, gy, gz, cx, cy):
    n = int(round(BLOCK_M / ac))
    c0 = int(round((cx - BLOCK_M / 2 - at.c) / ac))
    r0 = int(round((at.f - (cy + BLOCK_M / 2)) / ac))
    gr = Grid(origin_x=at.c + c0 * ac, origin_y=at.f - r0 * ac, cell=ac,
              width=n, height=n)
    dx, dy, dz, _ = dedupe_xy(gx, gy, gz, keep="lowest")
    z, _i = ALGOA.surfer_grid(dx, dy, dz, gr, pool_sigma_m=POOL,
                             verbose=False)
    return z, r0, c0, n


def score(gx, gy, gz, cx, cy):
    z, r0, c0, n = dem_from(gx, gy, gz, cx, cy)
    rr = R[r0:r0 + n, c0:c0 + n]
    aa = A[r0:r0 + n, c0:c0 + n]
    ok = np.isfinite(z) & np.isfinite(rr)
    oa = np.isfinite(z) & np.isfinite(aa)
    if ok.sum() < 100:
        return None
    e = (z - rr)[ok]
    ea = (z - aa)[oa] if oa.sum() > 100 else np.array([np.nan])
    return dict(cells=int(ok.sum()),
                rmse_ref=float(np.sqrt(np.mean(e ** 2))),
                mae_ref=float(np.mean(np.abs(e))),
                bias_ref=float(np.mean(e)),
                p99_ref=float(np.percentile(np.abs(e), 99)),
                rmse_arch=float(np.sqrt(np.nanmean(ea ** 2))))


def agreement(gx, gy, cx, cy):
    """Diagnostic only, per the brief: never the optimisation target."""
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


CACHE = {}


def prep(tag, lo, la):
    if tag in CACHE:
        return CACHE[tag]
    cx, cy = utm(lo, la)
    bx, by, bz = block(cx, cy, BUF_M)
    CACHE[tag] = (cx, cy, bx, by, bz)
    print("  %-11s %d returns, %.2f/m2"
          % (tag, len(bz), TSR.return_density(bx, by)), flush=True)
    return CACHE[tag]


print("")
print("development blocks")
for t, lo, la in DEV:
    prep(t, lo, la)
print("validation blocks")
for t, lo, la in VAL:
    prep(t, lo, la)


def run_config(blocks, over, label, verbose=False):
    base = TSR.GroundTSRParams()
    p = TSR.GroundTSRParams(**{k: over.get(k, getattr(base, k))
                               for k in base.__dataclass_fields__})
    rows = []
    for t, lo, la in blocks:
        cx, cy, bx, by, bz = prep(t, lo, la)
        try:
            mask, rep = TSR.classify_ground(bx, by, bz, params=p,
                                            verbose=verbose)
        except TSR.GroundTSRError as exc:
            print("    %s failed: %s" % (t, exc), flush=True)
            continue
        s = score(bx[mask], by[mask], bz[mask], cx, cy)
        if s is None:
            continue
        s.update(agreement(bx[mask], by[mask], cx, cy))
        s["block"] = t
        s["ground"] = int(mask.sum())
        s["density"] = rep["density_per_m2"]
        rows.append(s)
    if not rows:
        return None
    agg = {"label": label,
           "rmse_ref": float(np.mean([r["rmse_ref"] for r in rows])),
           "rmse_arch": float(np.mean([r["rmse_arch"] for r in rows])),
           "mae_ref": float(np.mean([r["mae_ref"] for r in rows])),
           "bias_ref": float(np.mean([r["bias_ref"] for r in rows])),
           "f1": float(np.mean([r.get("f1", np.nan) for r in rows])),
           "recall": float(np.mean([r.get("recall", np.nan) for r in rows])),
           "precision": float(np.mean([r.get("precision", np.nan)
                                       for r in rows])),
           "blocks": rows, "params": over}
    return agg


print("")
print("STAGE CONTRIBUTION on the development blocks")
print("%-34s %9s %11s %8s %9s %9s"
      % ("configuration", "RMSE ref", "RMSE arch", "F1", "recall", "precis"))
print("-" * 86)
stages = [
    ("confident set only (opening alone)",
     {"refine_iterations": -1, "gate_k_above": 0.0, "gate_k_below": 0.0,
      "gate_floor_above_m": 0.0, "gate_floor_below_m": 0.0,
      "slope_allowance": 0.0}),
    ("candidates ungated (curvature alone)",
     {"gate_ceiling_m": 999.0, "gate_floor_above_m": 999.0,
      "gate_floor_below_m": 999.0, "refine_iterations": 0}),
    ("adjudicated, no refinement", {"refine_iterations": 0}),
    ("adjudicated, one refinement", {}),
]
results = []
for label, over in stages:
    t0 = time.time()
    agg = run_config(DEV, over, label)
    if agg is None:
        continue
    results.append(agg)
    print("%-34s %9.4f %11.4f %8.3f %9.3f %9.3f"
          % (label[:34], agg["rmse_ref"], agg["rmse_arch"], agg["f1"],
             agg["recall"], agg["precision"]), flush=True)

print("")
print("GATE SWEEP on the development blocks")
print("%-34s %9s %11s %8s %9s %9s"
      % ("k_above / k_below / floor_above", "RMSE ref", "RMSE arch", "F1",
         "recall", "precis"))
print("-" * 86)
best, best_rmse = None, 1e9
for ka, kb, fa in itertools.product((1.0, 1.5, 2.0, 3.0), (2.0, 3.0, 4.0),
                                    (0.05, 0.08, 0.15)):
    over = {"gate_k_above": ka, "gate_k_below": kb,
            "gate_floor_above_m": fa}
    agg = run_config(DEV, over, "ka%.1f_kb%.1f_fa%.2f" % (ka, kb, fa))
    if agg is None:
        continue
    results.append(agg)
    mark = ""
    if agg["rmse_ref"] < best_rmse:
        best_rmse, best, mark = agg["rmse_ref"], over, "  <-"
    print("%-34s %9.4f %11.4f %8.3f %9.3f %9.3f%s"
          % ("%.1f / %.1f / %.2f" % (ka, kb, fa), agg["rmse_ref"],
             agg["rmse_arch"], agg["f1"], agg["recall"], agg["precision"],
             mark), flush=True)

print("")
print("best on development: %s at RMSE ref %.4f" % (best, best_rmse))
val = run_config(VAL, best, "validation")
print("")
print("HELD-OUT VALIDATION with that configuration")
if val:
    print("  RMSE ref %.4f, RMSE arch %.4f, MAE %.4f, bias %+.4f, F1 %.3f"
          % (val["rmse_ref"], val["rmse_arch"], val["mae_ref"],
             val["bias_ref"], val["f1"]))
    for b in val["blocks"]:
        print("    %-11s RMSE ref %.4f, arch %.4f, ground %d"
              % (b["block"], b["rmse_ref"], b["rmse_arch"], b["ground"]))

with open(os.path.join(OUT, "tsr_calibration.json"), "w") as fh:
    json.dump({"development": results, "best": best, "validation": val},
              fh, indent=2)
print("")
print("wrote %s" % os.path.join(OUT, "tsr_calibration.json"))
