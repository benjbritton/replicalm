r"""The three-rung ladder: filter divergence measured against a noise floor.

Four surfaces, each covered by three overlapping flight strips, differing in
what stands on them -- bare ground, low understory, closed forest, tall scrub.

NOT THE SAME THREE STRIPS ACROSS RUNGS. Each rung has three passes over it, which
is what repeatability needs, but they are different passes at different places:
bare is strips 1, 2, 3; low understory 3, 4, 7; forest 2, 3, 4; tall scrub
3, 4, 5. So a spread measured on one rung is that rung's own geometry, and the
floors are not strictly interchangeable between rungs. Comparing a filter
divergence against the floor measured on the SAME rung is sound; carrying one
rung's floor to another is not.

Each filter is run on each rung twice over:

  pooled        all three strips together, which is the surface a user gets
  per strip     each pass classified alone, then the spread of the three
                resulting medians -- three independent observations of the same
                ground, so this is precision against physics rather than against
                a commercial product

The point of the bare rung is that it carries no vegetation, so whatever spread
appears there is the measurement and processing noise of the survey itself.
A divergence between two filters that is smaller than that floor is not a
difference. That is the comparison the Optimization claim needs, and it needs no
surveyed control.

Clear's cleanup is measured on the same cells in the same pass: what fraction of
each filter's ground it removes, and how far the surface moves when it does.
Rung 2 carries a third of its returns in the band the cleanup cuts, and the bare
rung carries almost none, so the two together say whether the threshold acts on
vegetation or on everything.

Morphological filters need context well past a 20 m cell -- SMRF's window is
70 m -- so classification runs on a buffered neighbourhood and only the central
cell is measured.

Runs under the PDAL interpreter named by REPLICALM_PDAL_PYTHON (the ArcGIS Pro
one on this machine). The conda `replicalm` env has a broken BLAS: any numpy
matmul or lstsq aborts it outright, exit 127, no traceback.
"""
import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import ensure_blas, ROOT
sys.path.insert(0, os.path.join(ROOT, "src"))

ensure_blas()

import pdal
from replicalm.classify import _ground_stage, las_reader
from replicalm.config import GroundPass, PRESETS

D = os.environ.get(
    "REPLICALM_NR_ROOT",
    r"C:\NR_Lidar\Original_LAS_Files_15March2023\Original_LAS_Files_15March2023")
RESULTS = os.environ.get(
    "REPLICALM_NR_RESULTS",
    os.path.join(os.path.dirname(os.path.abspath(__file__)),
                 "results", "nr_block"))
os.makedirs(RESULTS, exist_ok=True)
OUT = os.path.join(RESULTS, "ladder.json")
HALF, BUF = 10.0, 70.0

# Rungs chosen 2026-09-26 from benchmarks/results/nr_block/band_profile.json.
# The previous low-understory pair was not understory at all: 82% and 75% of its
# returns stood above 5 m, so it was a second forest cell.
RUNGS = [
    ("bare",           320750.0, 1966969.9, "320500_1966500"),
    ("low understory", 320030.5, 1968030.4, "320000_1968000"),
    ("forest",         316450.8, 1970610.0, "316000_1970500"),
    ("tall scrub",     321450.7, 1967388.7, "321000_1967000"),
]

# `pmf_raw` is a stage dict rather than a GroundPass because GroundPass exposes
# only three of filters.pmf's knobs and leaves `initial_distance` at PDAL's
# default of 0.15. At that value PMF returned 0 ground points in the forest cell
# and 32 in tall scrub -- it was starved, not discriminating. benchmarks/
# pmf_tune.py swept 36 combinations on every rung; these recover ground density
# within 0.1 pts/m2 of SMRF on three rungs and a plausible surface on the fourth.
# `pmf_tight` is kept alongside so the starved behaviour stays visible rather
# than being quietly replaced.
FILTERS = [("smrf", GroundPass(algorithm="smrf", slope=0.05, threshold_m=0.5)),
           ("pmf", {"type": "filters.pmf", "max_window_size": 16.0,
                    "slope": 0.15, "max_distance": 1.5,
                    "initial_distance": 0.5, "cell_size": 1.0}),
           ("pmf_tight", GroundPass(algorithm="pmf", slope=0.05,
                                    threshold_m=0.5)),
           ("csf",  GroundPass(algorithm="csf", csf_rigidness=3,
                               csf_threshold_m=0.3))]


def classify(arr, gp):
    stage = dict(gp) if isinstance(gp, dict) else _ground_stage(gp)
    pl = pdal.Pipeline(json.dumps({"pipeline": [stage]}), arrays=[arr])
    pl.execute()
    out = pl.arrays[0]
    return out[np.isin(out["Classification"], (2,))]


def cleanup_effect(g, cx, cy, cfg):
    """What Clear's threshold removes from this filter's ground, in this cell.

    Heights are computed on the whole buffered set so points near the cell edge
    have their neighbours, then restricted to the central cell.
    """
    from replicalm import cleanup
    if len(g) < 200:
        return None
    above = cleanup.height_above_floor(
        g["X"].astype("f8"), g["Y"].astype("f8"), g["Z"].astype("f8"),
        patch=cfg.clean_patch_m, percentile=getattr(cfg, "clean_percentile", 10.0))
    inner = (np.abs(g["X"] - cx) <= HALF) & (np.abs(g["Y"] - cy) <= HALF)
    if inner.sum() < 100:
        return None
    keep = above <= cfg.clean_height_m
    z_before = float(np.median(g["Z"][inner]))
    kept = inner & keep
    if kept.sum() < 50:
        return {"removed_frac": round(float((~keep)[inner].mean()), 4),
                "kept": int(kept.sum())}
    return {"removed_frac": round(float((~keep)[inner].mean()), 4),
            "kept": int(kept.sum()),
            "median_before_m": round(z_before, 4),
            "median_after_m": round(float(np.median(g["Z"][kept])), 4),
            "shift_m": round(float(np.median(g["Z"][kept]) - z_before), 4)}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--rungs", nargs="*", default=None)
    a = ap.parse_args(argv)

    cfg = PRESETS["clear"]
    print("Clear cleanup under test: %.2f m above the %gth percentile of %.2f m"
          % (cfg.clean_height_m, getattr(cfg, "clean_percentile", 10.0),
             cfg.clean_patch_m), flush=True)

    rows, t0 = [], time.time()
    for label, cx, cy, tile in RUNGS:
        if a.rungs and label not in a.rungs:
            continue
        path = os.path.join(D, tile + ".las")
        if not os.path.exists(path):
            print("%-16s no LAS at %s" % (label, path)); continue
        pl = pdal.Pipeline(json.dumps({"pipeline": [las_reader(path)]}))
        pl.execute()
        arr = pl.arrays[0]
        nb = (np.abs(arr["X"] - cx) <= BUF) & (np.abs(arr["Y"] - cy) <= BUF)
        sub = arr[nb]
        ids, counts = np.unique(sub["PointSourceId"], return_counts=True)
        strips = [int(k) for k in ids[counts >= 0.05 * len(sub)]]
        print("\n=== %s === %.1f %.1f  %d returns in the %.0f m neighbourhood, "
              "strips %s  (%.0f s)"
              % (label, cx, cy, len(sub), 2 * BUF, strips, time.time() - t0),
              flush=True)

        rec = {"rung": label, "x": cx, "y": cy, "tile": tile,
               "returns": int(len(sub)), "strips": strips, "filters": {}}
        for name, gp in FILTERS:
            g = classify(sub, gp)
            inner = (np.abs(g["X"] - cx) <= HALF) & (np.abs(g["Y"] - cy) <= HALF)
            if inner.sum() < 100:
                rec["filters"][name] = {"n": int(inner.sum())}
                print("  %-5s only %d ground in cell" % (name, inner.sum()))
                continue
            z = g["Z"][inner]
            per = []
            for s in strips:
                m = sub["PointSourceId"] == s
                if m.sum() < 2000:
                    continue
                gs = classify(sub[m], gp)
                q = ((np.abs(gs["X"] - cx) <= HALF) & (np.abs(gs["Y"] - cy) <= HALF))
                if q.sum() >= 50:
                    per.append(float(np.median(gs["Z"][q])))
            d = {"n": int(inner.sum()),
                 "density": round(inner.sum() / (4 * HALF * HALF), 2),
                 "median_z": round(float(np.median(z)), 4),
                 "mean_z": round(float(z.mean()), 4),
                 "sd_z": round(float(z.std()), 4),
                 "strip_medians": [round(v, 4) for v in per],
                 "strip_spread_m": round(float(np.ptp(per)), 4)
                 if len(per) >= 2 else None,
                 "cleanup": cleanup_effect(g, cx, cy, cfg)}
            rec["filters"][name] = d
            print("  %-5s n %6d  density %6.2f  median %9.3f  spread %s  "
                  "cleanup removes %s"
                  % (name, d["n"], d["density"], d["median_z"],
                     ("%.4f m" % d["strip_spread_m"]) if d["strip_spread_m"]
                     is not None else "  n/a  ",
                     ("%.1f%%" % (100 * d["cleanup"]["removed_frac"]))
                     if d["cleanup"] else "n/a"), flush=True)

        ok = {k: v for k, v in rec["filters"].items() if v.get("median_z") is not None}
        rec["divergence_m"] = {
            "%s-%s" % (p, q): round(ok[p]["median_z"] - ok[q]["median_z"], 4)
            for i, p in enumerate(ok) for q in list(ok)[i + 1:]}
        rec["floor_m"] = round(float(np.mean(
            [v["strip_spread_m"] for v in ok.values()
             if v.get("strip_spread_m") is not None])), 4) if ok else None
        rows.append(rec)
        with open(OUT, "w", encoding="utf-8") as fh:
            json.dump(rows, fh, indent=1)

    print("\n%-16s %6s %10s %10s %10s %10s"
          % ("rung", "filter", "median z", "3-strip", "clean %", "shift m"))
    print("-" * 68)
    for r in rows:
        for name in ("smrf", "pmf", "pmf_tight", "csf"):
            d = r["filters"].get(name, {})
            if "median_z" not in d:
                continue
            c = d.get("cleanup") or {}
            print("%-16s %6s %10.3f %10s %9s%% %10s"
                  % (r["rung"], name, d["median_z"],
                     ("%.4f" % d["strip_spread_m"]) if d["strip_spread_m"]
                     is not None else "n/a",
                     ("%.1f" % (100 * c["removed_frac"])) if c else "n/a",
                     ("%+.4f" % c["shift_m"]) if c.get("shift_m") is not None
                     else "n/a"))

    print("\n%-16s %10s   %s" % ("rung", "floor", "filter divergence (m)"))
    print("-" * 68)
    for r in rows:
        div = "  ".join("%s %+.3f" % (k, v) for k, v in r["divergence_m"].items())
        print("%-16s %10s   %s"
              % (r["rung"], ("%.4f" % r["floor_m"]) if r["floor_m"] else "n/a", div))
    print("\nwrote %s" % OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
