"""Filter behaviour against understory fraction, across the New River block.

Design. Each candidate is a 20 m cell that is flat, under closed canopy, and
covered by three strips. Morphological filters need context well beyond a 20 m
cell -- SMRF's window is 70 m -- so classification runs on a 140 m buffered
neighbourhood and only the central cell is measured.

Two responses, neither referencing a proprietary product:

  divergence  the mean elevation of each filter's ground surface in the central
              cell, differenced against the other filters. How far apart the
              three algorithms place the ground.

  repeatability  each strip classified alone, then the spread of the three
              resulting mean elevations. Independent observations of the same
              ground, so this is precision against physics.

Both are regressed on the cell's understory band fraction.
"""
import gzip, json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import ensure_blas, SRC
sys.path.insert(0, SRC)
ensure_blas()
import pdal
from replicalm.classify import _ground_stage, las_reader
from replicalm.config import GroundPass

# Where the New River block lives is a property of the machine, not of
# the method, so it comes from the environment with this machine as the
# fallback -- the same arrangement as CLOUDCOMPARE in interpolate.py.
D = os.environ.get(
    "REPLICALM_NR_ROOT",
    r"C:\NR_Lidar\Original_LAS_Files_15March2023\Original_LAS_Files_15March2023")
# Results live in the repository now; McLellan gave permission to publish
# measurements derived from the block on 2026-09-24. The environment
# variable still overrides, for working outside a checkout.
RESULTS = os.environ.get(
    "REPLICALM_NR_RESULTS",
    os.path.join(os.path.dirname(os.path.abspath(__file__)),
                 "results", "nr_block"))
os.makedirs(RESULTS, exist_ok=True)
# The tracked scan is gzipped and already restricted to three-strip cells.
CELLS = os.path.join(RESULTS, "cells_3strip.json.gz")
OUT = os.path.join(RESULTS, sys.argv[1] if len(sys.argv) > 1
                   else "band_regression_tuned.json")
HALF, BUF = 10.0, 70.0

# `pmf` is a stage dict rather than a GroundPass because GroundPass exposes only
# three of filters.pmf's knobs and leaves `initial_distance` at PDAL's default of
# 0.15 m. At that value PMF is starved: on the four-rung ladder it recovered no
# ground at all under closed canopy, and on this population its density response
# to understory was -0.06 returns per m2 per 10% band -- which is what a filter
# that classifies almost nothing looks like, not evidence that it resists
# understory. benchmarks/pmf_tune.py swept 36 combinations on four surfaces;
# these bring it within 0.1 points per m2 of SMRF.
#
# `pmf_tight` is the old setting, kept in the same run so the two are compared on
# identical cells by identical code rather than across two runs.
FILTERS = [("smrf", GroundPass(algorithm="smrf", slope=0.05, threshold_m=0.5)),
           ("pmf", {"type": "filters.pmf", "max_window_size": 16.0,
                    "slope": 0.15, "max_distance": 1.5,
                    "initial_distance": 0.5, "cell_size": 1.0}),
           ("pmf_tight", GroundPass(algorithm="pmf", slope=0.05,
                                    threshold_m=0.5)),
           ("csf",  GroundPass(algorithm="csf",  csf_rigidness=3, csf_threshold_m=0.3))]

def classify(arr, gp):
    """Ground points from one filter, run on the array as given."""
    stage = dict(gp) if isinstance(gp, dict) else _ground_stage(gp)
    pl = pdal.Pipeline(json.dumps({"pipeline": [stage]}), arrays=[arr])
    pl.execute()
    out = pl.arrays[0]
    return out[np.isin(out["Classification"], (2,))]

# The tracked scan is gzipped and already restricted to three-strip cells; the
# full uncompressed scan works too if it is to hand.
if CELLS.endswith(".gz"):
    R = json.load(gzip.open(CELLS, "rt", encoding="utf-8"))
else:
    R = json.load(open(CELLS))
cand = [c for c in R if c["canopy_p90_m"] > 8 and c["frac_over_2m"] > 0.35
        and c["frac_band_015_035"] > 0.03 and len(c["strips"]) >= 3
        and c["slope_pct"] < 2.0]
by_tile = {}
for c in cand:
    by_tile.setdefault(c["tile"], []).append(c)
print("%d cells across %d tiles" % (len(cand), len(by_tile)), flush=True)

rows, t0 = [], time.time()
for ti, (tile, cells) in enumerate(sorted(by_tile.items()), 1):
    pl = pdal.Pipeline(json.dumps({"pipeline": [las_reader(os.path.join(D, tile + ".las"))]}))
    pl.execute()
    tile_arr = pl.arrays[0]
    tx, ty = tile_arr["X"], tile_arr["Y"]
    for c in cells:
        nb = ((np.abs(tx - c["x"]) <= BUF) & (np.abs(ty - c["y"]) <= BUF))
        if nb.sum() < 5000:
            continue                      # neighbourhood runs off the tile edge
        sub = tile_arr[nb]
        rec = {"x": c["x"], "y": c["y"], "tile": tile,
               "band": c["frac_band_015_035"], "canopy_p90_m": c["canopy_p90_m"],
               "slope_pct": c["slope_pct"], "ground_per_m2": c["ground_per_m2"],
               "filters": {}}
        for name, gp in FILTERS:
            try:
                g = classify(sub, gp)
            except Exception as e:
                rec["filters"][name] = {"error": str(e)[:120]}
                continue
            inner = (np.abs(g["X"] - c["x"]) <= HALF) & (np.abs(g["Y"] - c["y"]) <= HALF)
            if inner.sum() < 100:
                rec["filters"][name] = {"n": int(inner.sum())}
                continue
            z = g["Z"][inner]
            # per-strip repeatability: same filter, each pass classified alone
            per = []
            for s in c["strips"]:
                m = sub["PointSourceId"] == s
                if m.sum() < 2000:
                    continue
                try:
                    gs = classify(sub[m], gp)
                except Exception:
                    continue
                q = (np.abs(gs["X"] - c["x"]) <= HALF) & (np.abs(gs["Y"] - c["y"]) <= HALF)
                if q.sum() >= 50:
                    per.append(float(np.median(gs["Z"][q])))
            rec["filters"][name] = {
                "n": int(inner.sum()), "density": round(inner.sum() / (4 * HALF * HALF), 2),
                "mean_z": round(float(z.mean()), 4), "median_z": round(float(np.median(z)), 4),
                "sd_z": round(float(z.std()), 4),
                "strip_medians": [round(v, 4) for v in per],
                "strip_spread_m": round(float(np.ptp(per)), 4) if len(per) >= 2 else None}
        rows.append(rec)
    print("  tile %2d/%d  %s  %d cells done  %.0f s"
          % (ti, len(by_tile), tile, len(rows), time.time() - t0), flush=True)
    json.dump(rows, open(OUT, "w"), indent=1)
print("wrote %s: %d cells in %.0f s" % (OUT, len(rows), time.time() - t0), flush=True)
