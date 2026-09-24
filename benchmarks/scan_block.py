"""Score every 20 m cell in the New River block for the three-rung site selection.

Tiles are 500 m and cells 20 m, so cells align to tile origins exactly and no
cell is split across a tile boundary.

Per cell, from ground returns (classes 2 and 8): a least-squares plane, its
slope and residual scatter. From all returns: height above that plane, which
gives canopy height and the fraction of returns sitting in the 0.15-0.35 m
window the middle rung needs. Strip coverage comes from PointSourceId, counting
only strips holding at least 5% of the cell, so a handful of stray returns from
a neighbouring pass does not read as independent coverage.
"""
import glob, gzip, json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import ROOT, TESTS
sys.path.insert(0, os.path.join(ROOT, "src"))
import pdal
from replicalm.classify import las_reader

# Where the New River block lives is a property of the machine, not of
# the method, so it comes from the environment with this machine as the
# fallback -- the same arrangement as CLOUDCOMPARE in interpolate.py.
D = os.environ.get(
    "REPLICALM_NR_ROOT",
    r"C:\NR_Lidar\Original_LAS_Files_15March2023\Original_LAS_Files_15March2023")
# The full scan is 23 MB and goes to the work directory. What the
# repository keeps is the subset covered by three or more strips -- the
# requirement for any interswath comparison -- plus the denominators the
# full scan supports.
# Results live in the repository now; McLellan gave permission to publish
# measurements derived from the block on 2026-09-24. The environment
# variable still overrides, for working outside a checkout.
RESULTS = os.environ.get(
    "REPLICALM_NR_RESULTS",
    os.path.join(os.path.dirname(os.path.abspath(__file__)),
                 "results", "nr_block"))
os.makedirs(RESULTS, exist_ok=True)
OUT = os.path.join(TESTS, "nr_block", "cells.json")
TRIM = os.path.join(RESULTS, "cells_3strip.json.gz")
SUMMARY = os.path.join(RESULTS, "cells_summary.json")
for _d in (os.path.dirname(OUT), RESULTS):
    os.makedirs(_d, exist_ok=True)
W = 20.0
rows, t0 = [], time.time()
files = sorted(glob.glob(os.path.join(D, "*.las")))
for n, f in enumerate(files, 1):
    try:
        pl = pdal.Pipeline(json.dumps({"pipeline": [las_reader(f)]}))
        pl.execute()
        a = pl.arrays[0]
    except Exception as e:
        print("  skip %s: %s" % (os.path.basename(f), e), flush=True)
        continue
    x, y, z = a["X"], a["Y"], a["Z"]
    cls, sid, sa = a["Classification"], a["PointSourceId"], a["ScanAngleRank"]
    ox, oy = np.floor(x.min() / W) * W, np.floor(y.min() / W) * W
    ix = ((x - ox) // W).astype(np.int32)
    iy = ((y - oy) // W).astype(np.int32)
    key = ix.astype(np.int64) * 100000 + iy
    o = np.argsort(key, kind="stable")
    key, x, y, z, cls, sid, sa = (v[o] for v in (key, x, y, z, cls, sid, sa))
    uk, start = np.unique(key, return_index=True)
    end = np.append(start[1:], len(key))
    for s, e in zip(start, end):
        g = np.isin(cls[s:e], (2, 8))
        ng = int(g.sum())
        if ng < 400:
            continue
        gx, gy, gz = x[s:e][g], y[s:e][g], z[s:e][g]
        mx, my = float(gx.mean()), float(gy.mean())
        A = np.c_[gx - mx, gy - my, np.ones(ng)]
        coef, *_ = np.linalg.lstsq(A, gz, rcond=None)
        res = gz - A @ coef
        h = z[s:e] - (np.c_[x[s:e] - mx, y[s:e] - my, np.ones(e - s)] @ coef)
        ids, counts = np.unique(sid[s:e], return_counts=True)
        keep = ids[counts >= 0.05 * (e - s)]
        angles = [[int(sa[s:e][sid[s:e] == k].min()),
                   int(sa[s:e][sid[s:e] == k].max())] for k in keep]
        rows.append({
            "x": round(mx, 1), "y": round(my, 1),
            "slope_pct": round(float(100 * np.hypot(coef[0], coef[1])), 3),
            "resid_sd_m": round(float(res.std()), 4),
            "ground_per_m2": round(ng / (W * W), 2),
            "returns": int(e - s),
            "canopy_p90_m": round(float(np.percentile(h, 90)), 2),
            "frac_over_2m": round(float(np.mean(h > 2.0)), 4),
            "frac_band_015_035": round(float(np.mean((h >= 0.15) & (h <= 0.35))), 4),
            "strips": [int(k) for k in keep], "scan_angles": angles,
            "tile": os.path.basename(f)[:-4]})
    if n % 12 == 0:
        print("  %3d/%d tiles, %d cells, %.0f s" % (n, len(files), len(rows), time.time() - t0), flush=True)
json.dump(rows, open(OUT, "w"), indent=1)
three = [c for c in rows if len(c["strips"]) >= 3]
with gzip.open(TRIM, "wt", encoding="utf-8", compresslevel=9) as fh:
    json.dump(three, fh, separators=(",", ":"))
forest = [c for c in rows if c["canopy_p90_m"] > 8 and c["frac_over_2m"] > 0.35]
band = np.array([c["frac_band_015_035"] for c in forest]) if forest else np.zeros(1)
json.dump({
    "scan": {"cell_m": int(W), "tiles": len({c["tile"] for c in rows}),
             "cells_scored": len(rows), "min_ground_returns_per_cell": 400},
    "coverage": {"cells_with_3_strips": len(three),
                 "cells_with_2_strips": sum(1 for c in rows if len(c["strips"]) == 2),
                 "cells_with_1_strip": sum(1 for c in rows if len(c["strips"]) <= 1)},
    "populations": {"flat_under_1p5pct_and_sd_under_6cm":
                        sum(1 for c in rows if c["slope_pct"] < 1.5
                            and c["resid_sd_m"] < 0.06),
                    "forest_canopy_over_8m_and_35pct_above_2m": len(forest)},
    "understory_band_over_forest_cells": {
        "median": round(float(np.median(band)), 5),
        "p90": round(float(np.percentile(band, 90)), 5),
        "max": round(float(band.max()), 5)},
}, open(SUMMARY, "w"), indent=1)
print("tracked: %s (%d cells), %s" % (TRIM, len(three), SUMMARY))
three = [c for c in rows if len(c["strips"]) >= 3]
with gzip.open(TRIM, "wt", encoding="utf-8", compresslevel=9) as fh:
    json.dump(three, fh, separators=(",", ":"))
forest = [c for c in rows if c["canopy_p90_m"] > 8 and c["frac_over_2m"] > 0.35]
band = np.array([c["frac_band_015_035"] for c in forest]) if forest else np.zeros(1)
json.dump({
    "scan": {"cell_m": int(W), "tiles": len({c["tile"] for c in rows}),
             "cells_scored": len(rows), "min_ground_returns_per_cell": 400},
    "coverage": {"cells_with_3_strips": len(three),
                 "cells_with_2_strips": sum(1 for c in rows if len(c["strips"]) == 2),
                 "cells_with_1_strip": sum(1 for c in rows if len(c["strips"]) <= 1)},
    "populations": {
        "flat_under_1p5pct_and_sd_under_6cm":
            sum(1 for c in rows if c["slope_pct"] < 1.5 and c["resid_sd_m"] < 0.06),
        "forest_canopy_over_8m_and_35pct_above_2m": len(forest)},
    "understory_band_over_forest_cells": {
        "median": round(float(np.median(band)), 5),
        "p90": round(float(np.percentile(band, 90)), 5),
        "max": round(float(band.max()), 5)},
}, open(SUMMARY, "w"), indent=1)
print("wrote %s: %d cells from %d tiles in %.0f s" % (OUT, len(rows), len(files), time.time() - t0))
print("tracked: %s (%d cells), %s" % (TRIM, len(three), SUMMARY))
