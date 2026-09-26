r"""Return-height profiles for shortlisted cells, across many bands not one.

The block scan records a single figure per cell, `frac_band_015_035`. That
window straddles Clear's 0.20 m threshold: a return at 0.17 m should be kept and
one at 0.30 m removed, and one fraction cannot tell them apart. So a cell picked
on that number is ambiguous with respect to the threshold it is meant to test.

This reads the clouds for a shortlist and reports the whole height distribution,
per cell and per strip. The datum is unchanged -- a least-squares plane through
the delivered ground returns (classes 2 and 8) in the same 20 m cell, with all
returns measured above it -- so these numbers sit beside the inventory rather
than beneath a different definition.

Per strip as well as pooled, because a rung only carries the interswath
comparison if all three passes see the same vegetation.
"""
import argparse
import glob
import gzip
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import ensure_blas, ROOT
sys.path.insert(0, os.path.join(ROOT, "src"))

ensure_blas()

D = os.environ.get(
    "REPLICALM_NR_ROOT",
    r"C:\NR_Lidar\Original_LAS_Files_15March2023\Original_LAS_Files_15March2023")
RESULTS = os.environ.get(
    "REPLICALM_NR_RESULTS",
    os.path.join(os.path.dirname(os.path.abspath(__file__)),
                 "results", "nr_block"))
W = 20.0

# Finer either side of 0.20 m than anywhere else: that is the only place the
# answer changes. The top bands exist to say what else is standing in the cell.
EDGES = [-0.05, 0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.35,
         0.50, 1.00, 2.00, 5.00, 1e9]

# The ladder as it stands, plus the four candidates the inventory offers.
SHORTLIST = [
    ("current low-understory A", 316509.0, 1970589.0),
    ("current low-understory B", 316509.0, 1970631.0),
    ("Ka'Kabish test surface",   316450.8, 1970610.0),
    ("Cocochan plaza",           320630.1, 1966971.9),
    ("Cocochan bare control",    320750.0, 1966969.9),
    ("candidate 1",              320030.5, 1968030.4),
    ("candidate 2",              319489.0, 1968429.1),
    ("candidate 3",              319950.5, 1968090.5),
    ("candidate 4",              321450.7, 1967388.7),
]


def locate(targets):
    """Attach each target to its scanned cell, so the tile is known."""
    path = os.path.join(RESULTS, "cells_3strip.json.gz")
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        cells = json.load(fh)
    cx = np.array([c["x"] for c in cells])
    cy = np.array([c["y"] for c in cells])
    out = []
    for label, x, y in targets:
        d = np.hypot(cx - x, cy - y)
        j = int(np.argmin(d))
        if d[j] > W:
            print("%-26s no three-strip cell within %.0f m (nearest %.1f m)"
                  % (label, W, d[j]))
            continue
        c = cells[j]
        out.append({"label": label, "asked": (x, y), "cell": c,
                    "offset_m": round(float(d[j]), 1)})
    return out


def profile(tile_path, wanted):
    """Height distribution for every wanted cell inside one tile."""
    import pdal
    from replicalm.classify import las_reader
    pl = pdal.Pipeline(json.dumps({"pipeline": [las_reader(tile_path)]}))
    pl.execute()
    a = pl.arrays[0]
    x, y, z = a["X"], a["Y"], a["Z"]
    cls, sid = a["Classification"], a["PointSourceId"]
    rows = []
    for item in wanted:
        c = item["cell"]
        # the scan reports the mean of the ground returns, not the cell corner,
        # so recover the cell from the same origin convention it used
        ox = np.floor(x.min() / W) * W
        oy = np.floor(y.min() / W) * W
        ix = int((c["x"] - ox) // W)
        iy = int((c["y"] - oy) // W)
        m = (((x - ox) // W).astype(int) == ix) & (((y - oy) // W).astype(int) == iy)
        if m.sum() < 400:
            print("%-26s only %d returns in the cell" % (item["label"], m.sum()))
            continue
        gx, gy, gz = x[m], y[m], z[m]
        g = np.isin(cls[m], (2, 8))
        if g.sum() < 100:
            print("%-26s only %d ground returns" % (item["label"], g.sum()))
            continue
        mx, my = float(gx[g].mean()), float(gy[g].mean())
        A = np.c_[gx[g] - mx, gy[g] - my, np.ones(int(g.sum()))]
        coef, *_ = np.linalg.lstsq(A, gz[g], rcond=None)
        h = gz - (np.c_[gx - mx, gy - my, np.ones(m.sum())] @ coef)
        s = sid[m]
        ids, counts = np.unique(s, return_counts=True)
        strips = [int(k) for k in ids[counts >= 0.05 * m.sum()]]

        hist = np.histogram(h, bins=EDGES)[0] / float(len(h))
        per_strip = {}
        for k in strips:
            hk = h[s == k]
            per_strip[str(k)] = {
                "returns": int(hk.size),
                "band_015_020": round(float(np.mean((hk >= 0.15) & (hk < 0.20))), 4),
                "band_020_035": round(float(np.mean((hk >= 0.20) & (hk < 0.35))), 4),
                "band_035_100": round(float(np.mean((hk >= 0.35) & (hk < 1.00))), 4),
                "p90_m": round(float(np.percentile(hk, 90)), 3)}
        rows.append({
            "label": item["label"], "x": c["x"], "y": c["y"],
            "tile": c["tile"], "offset_from_asked_m": item["offset_m"],
            "slope_pct": c["slope_pct"], "returns": int(m.sum()),
            "ground_returns": int(g.sum()), "strips": strips,
            "edges": EDGES[:-1] + [None],
            "fractions": [round(float(v), 5) for v in hist],
            "per_strip": per_strip})
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(RESULTS, "band_profile.json"))
    a = ap.parse_args(argv)

    located = locate(SHORTLIST)
    by_tile = {}
    for item in located:
        by_tile.setdefault(item["cell"]["tile"], []).append(item)
    print("%d cells across %d tiles" % (len(located), len(by_tile)), flush=True)

    rows, t0 = [], time.time()
    for tile, wanted in sorted(by_tile.items()):
        hits = glob.glob(os.path.join(D, tile + ".las"))
        if not hits:
            print("  %s: no LAS in %s" % (tile, D))
            continue
        print("  %s: %d cells (%.0f s)" % (tile, len(wanted), time.time() - t0),
              flush=True)
        rows += profile(hits[0], wanted)

    order = {label: i for i, (label, _x, _y) in enumerate(SHORTLIST)}
    rows.sort(key=lambda r: order.get(r["label"], 99))

    names = ["<0", "0-.05", ".05-.10", ".10-.15", ".15-.20", ".20-.25",
             ".25-.35", ".35-.50", ".50-1", "1-2", "2-5", ">5"]
    print("")
    print("%-26s %6s %5s " % ("cell", "slope", "strp") +
          " ".join("%7s" % n for n in names))
    print("-" * (39 + 8 * len(names)))
    for r in rows:
        print("%-26s %5.2f%% %5d " % (r["label"], r["slope_pct"], len(r["strips"]))
              + " ".join("%6.2f%%" % (100 * f) for f in r["fractions"]))

    print("")
    print("%-26s %11s %11s %11s   %s"
          % ("cell", "keep .15-.20", "cut .20-.35", "cut .35-1.0", "per-strip .20-.35"))
    print("-" * 96)
    for r in rows:
        f = r["fractions"]
        keep = f[4]
        cut1 = f[5] + f[6]
        cut2 = f[7] + f[8]
        ps = " ".join("%.1f%%" % (100 * v["band_020_035"])
                      for v in r["per_strip"].values())
        print("%-26s %10.2f%% %10.2f%% %10.2f%%   %s"
              % (r["label"], 100 * keep, 100 * cut1, 100 * cut2, ps))

    with open(a.out, "w", encoding="utf-8") as fh:
        json.dump({"cell_m": W, "edges": EDGES, "datum":
                   "least-squares plane through delivered ground classes 2,8",
                   "cells": rows}, fh, indent=1)
    print("")
    print("wrote %s" % a.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
