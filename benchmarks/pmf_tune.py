r"""Is PMF's failure on the vegetated rungs a property of PMF or of our settings?

In `ladder.py` PMF returned 0 ground points in the forest cell and 32 in tall
scrub, where SMRF and CSF both returned thousands. That is a failure, not a
measurement, and it cannot be reported as filter behaviour until the parameters
are ruled out.

Replicalm's `GroundPass` exposes only three of PDAL's `filters.pmf` knobs --
`max_window_size`, `slope`, `max_distance` -- and leaves `initial_distance`,
`cell_size` and `exponential` at their defaults. This sweeps the stage directly,
below that interface, so the question is asked of the algorithm rather than of
our wrapper.

What counts as a working setting: ground density in the central cell within
reach of what SMRF and CSF recover there, and a median elevation that is not
obviously riding the canopy. Those two references are read from ladder.json so
the comparison is against measured values, not against judgement.
"""
import argparse
import itertools
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
from replicalm.classify import las_reader

D = os.environ.get(
    "REPLICALM_NR_ROOT",
    r"C:\NR_Lidar\Original_LAS_Files_15March2023\Original_LAS_Files_15March2023")
RESULTS = os.environ.get(
    "REPLICALM_NR_RESULTS",
    os.path.join(os.path.dirname(os.path.abspath(__file__)),
                 "results", "nr_block"))
OUT = os.path.join(RESULTS, "pmf_tune.json")
HALF, BUF = 10.0, 70.0

RUNGS = {
    "forest":     (316450.8, 1970610.0, "316000_1970500"),
    "tall scrub": (321450.7, 1967388.7, "321000_1967000"),
    "low understory": (320030.5, 1968030.4, "320000_1968000"),
    "bare":       (320750.0, 1966969.9, "320500_1966500"),
}

# The current setting is the first row: window 70, slope 0.05, max_distance 0.5,
# and PDAL's own defaults for the rest.
GRID = {
    "max_window_size": [16.0, 33.0, 70.0],
    "slope": [0.05, 0.15],
    "max_distance": [0.5, 1.5, 3.0],
    "initial_distance": [0.15, 0.5],
    "cell_size": [1.0],
}


def run(sub, cx, cy, params):
    stage = dict(params)
    stage["type"] = "filters.pmf"
    pl = pdal.Pipeline(json.dumps({"pipeline": [stage]}), arrays=[sub])
    pl.execute()
    out = pl.arrays[0]
    g = out[out["Classification"] == 2]
    if not len(g):
        return {"n": 0, "density": 0.0}
    inner = (np.abs(g["X"] - cx) <= HALF) & (np.abs(g["Y"] - cy) <= HALF)
    if inner.sum() < 20:
        return {"n": int(inner.sum()), "density": round(inner.sum() / (4 * HALF * HALF), 3)}
    z = g["Z"][inner]
    return {"n": int(inner.sum()),
            "density": round(inner.sum() / (4 * HALF * HALF), 2),
            "median_z": round(float(np.median(z)), 4),
            "sd_z": round(float(z.std()), 4)}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--rungs", nargs="*",
                    default=["forest", "tall scrub"])
    a = ap.parse_args(argv)

    ref = {}
    lp = os.path.join(RESULTS, "ladder.json")
    if os.path.exists(lp):
        for r in json.load(open(lp, encoding="utf-8")):
            ref[r["rung"]] = {k: v.get("median_z") and
                              {"n": v["n"], "density": v["density"],
                               "median_z": v["median_z"]}
                              for k, v in r["filters"].items()}

    keys = list(GRID)
    combos = list(itertools.product(*(GRID[k] for k in keys)))
    print("%d parameter combinations per rung" % len(combos), flush=True)

    rows, t0 = [], time.time()
    for rung in a.rungs:
        cx, cy, tile = RUNGS[rung]
        pl = pdal.Pipeline(json.dumps({"pipeline": [
            las_reader(os.path.join(D, tile + ".las"))]}))
        pl.execute()
        arr = pl.arrays[0]
        nb = (np.abs(arr["X"] - cx) <= BUF) & (np.abs(arr["Y"] - cy) <= BUF)
        sub = arr[nb]
        r = ref.get(rung, {})
        tgt = [v["density"] for k, v in r.items()
               if k in ("smrf", "csf") and v]
        target = float(np.mean(tgt)) if tgt else None
        print("\n=== %s === %d returns; SMRF/CSF recover %s pts/m2 at %s m"
              % (rung, len(sub),
                 ("%.2f" % target) if target else "?",
                 ("%.3f" % r["smrf"]["median_z"]) if r.get("smrf") else "?"),
              flush=True)
        print("%9s %6s %8s %9s %8s %9s %10s"
              % ("window", "slope", "maxdist", "initdist", "n", "dens", "median z"))
        print("-" * 66)
        for vals in combos:
            p = dict(zip(keys, vals))
            try:
                res = run(sub, cx, cy, p)
            except Exception as exc:
                print("%9.0f %6.2f %8.1f %9.2f   FAILED %s"
                      % (p["max_window_size"], p["slope"], p["max_distance"],
                         p["initial_distance"], str(exc)[:40]), flush=True)
                continue
            rows.append({"rung": rung, "params": p, "result": res,
                         "reference": r})
            print("%9.0f %6.2f %8.1f %9.2f %8d %9.2f %10s"
                  % (p["max_window_size"], p["slope"], p["max_distance"],
                     p["initial_distance"], res["n"], res["density"],
                     ("%.3f" % res["median_z"]) if "median_z" in res else "-"),
                  flush=True)
            with open(OUT, "w", encoding="utf-8") as fh:
                json.dump(rows, fh, indent=1)
        print("  %.0f s" % (time.time() - t0), flush=True)

    print("\nwrote %s" % OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
