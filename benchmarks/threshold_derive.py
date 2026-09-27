r"""Derive Clear's cleanup threshold from repetition, instead of inheriting it.

Clear removes classified-ground returns standing more than 0.20 m above the
ground around them. That number came from the source method; nothing in this
project has ever measured it, and the ladder showed it does not do what the
draft claims -- on understory it removes 2.7% and moves the surface 1 cm.

THE IDEA
--------
Three flights cover each test surface. Vegetation is present in every pass, so
presence is not the test. What differs is the height that comes back: a porous
canopy is penetrated to different depths at different scan angles, so its
returns scatter between passes. A rock, a road, a platform riser are opaque and
rigid -- every pass reports the same surface. So inter-pass agreement separates
the two, and the height at which agreement collapses is the threshold.

WHAT IS MEASURED, PER CELL
--------------------------
  h   the cell's height above a plane fitted to ground within PLANE_R of it.
      Local, because scale decides the answer: on a road that is demonstrably
      bare, a plane fitted over 20 m reports 45% of returns standing above
      0.15 m, and the same returns against a 3 m plane report 0%. The ground
      itself ranges 0.66 m over 15 m and 0.14 m over 3 m.

  d   the spread of the three single-strip median elevations. Three
      independent observations of the same ground, so this is precision
      against physics rather than against a commercial product.

THE RUNGS DO THE LABELLING, NOT OUR OWN OUTPUT
----------------------------------------------
Labelling cells clutter or not from the DEM Clear produced would argue in a
circle. Instead each surface was chosen for what stands on it:

  bare            a dirt road; every return is ground, nothing above 0.15 m
  clutter_dense   21% of returns in 0.15-0.35 m, no canopy, 15 m from the road
  clutter_sparse  3.7% in that band -- the sparse case observation 6 describes
  riser_open      16% slope, no canopy: the clean false-positive test, since
                  nothing stands here and anything removed is a mistake
  mound_flank     a 40 degree mound flank under canopy: real architecture, but
                  removals here are ambiguous between vegetation and terrain
  forest          flat ground under 19 m canopy: the hardest case

WHAT COMES OUT
--------------
Not a filter stage -- most G-LiHT transects are single-pass and could not run
one. A number, with the cost of choosing it: at each candidate threshold, what
fraction of each surface is removed. The bare and riser rungs give the false
positives, the clutter rungs the true ones.
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
from replicalm.config import PRESETS

D = os.environ.get(
    "REPLICALM_NR_ROOT",
    r"C:\NR_Lidar\Original_LAS_Files_15March2023\Original_LAS_Files_15March2023")
RESULTS = os.environ.get(
    "REPLICALM_NR_RESULTS",
    os.path.join(os.path.dirname(os.path.abspath(__file__)),
                 "results", "nr_block"))
OUT = os.path.join(RESULTS, "threshold_derive.json")

BUF = 70.0         # context for SMRF, whose window runs to 70 m
PLANE_R = 3.0      # radius the detrending plane is fitted over
MIN_PER_STRIP = 3  # returns a strip needs in a cell to vote on its elevation

# Cell size is set per rung, because it is bounded from both sides and the
# bounds differ. A cell needs enough returns for each strip to report a median
# -- at 1.0 m under canopy that failed, giving three strips in 24 cells of 321
# on the riser and 1 of 217 on the mound -- and it must stay small enough to
# resolve a speck, which observation 6 measures at about one 0.5 m cell.
# Where ground density is high, 1.0 m satisfies both; under canopy and on the
# steep flanks it has to be 2.0 m.
#
# Cell size differing between rungs is fine for the question that matters,
# which is asked WITHIN a rung: does inter-pass agreement fall away as height
# rises? The cost table is computed per return instead of per cell, so it stays
# comparable across rungs regardless.

# Each rung is only homogeneous out to a certain distance, and measuring beyond
# it labels the wrong ground: a 40 m square on the road is mostly verge, which
# is why the first run reported 7% of a demonstrably clean surface as standing
# proud. The half-widths below are what each surface actually holds -- the road
# is clean to about 5 m, the wider surfaces further.
#
# name, easting, northing, class for the cost table, half-width
# name, easting, northing, class, half-width, cell
RUNGS = [
    ("bare",           317124.75, 1970163.18, "should keep",   5.0, 1.0),
    ("clutter_dense",  317130.60, 1970149.50, "should remove", 9.0, 1.0),
    ("clutter_sparse", 316698.41, 1970492.90, "should remove", 9.0, 1.0),
    ("riser_open",     316810.00, 1970329.60, "should keep",  12.0, 2.0),
    # widened to 50 m: at 24 m only 10 cells of 116 had all three strips
    # voting, because per-strip ground recovery on a 40 degree flank is thin.
    # The wider square takes in the crest and the flat below, which are still
    # real ground and still carry the should-keep label, so the slope bands in
    # the report are what to read rather than the rung average.
    ("mound_flank",    316622.34, 1970507.42, "should keep",  25.0, 2.0),
    ("forest",         316450.80, 1970610.00, "unlabelled",   12.0, 2.0),
]
THRESHOLDS = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50]


def tile_for(x, y):
    return "%d_%d" % (int(np.floor(x / 500) * 500), int(np.floor(y / 500) * 500))


def classify(arr, gp):
    pl = pdal.Pipeline(json.dumps({"pipeline": [_ground_stage(gp)]}), arrays=[arr])
    pl.execute()
    out = pl.arrays[0]
    return out[out["Classification"] == 2]


def cell_table(gx, gy, gz, cx, cy, half, cell):
    """Per-cell median elevation and height above a locally fitted plane."""
    ix = np.floor((gx - (cx - half)) / cell).astype(int)
    iy = np.floor((gy - (cy - half)) / cell).astype(int)
    n = int(round(2 * half / cell))
    keep = (ix >= 0) & (ix < n) & (iy >= 0) & (iy < n)
    ix, iy = ix[keep], iy[keep]
    gx, gy, gz = gx[keep], gy[keep], gz[keep]
    key = iy.astype(np.int64) * n + ix
    order = np.argsort(key, kind="stable")
    key, gx, gy, gz = key[order], gx[order], gy[order], gz[order]
    uk, start = np.unique(key, return_index=True)
    end = np.append(start[1:], len(key))

    rows = []
    for k, s, e in zip(uk, start, end):
        if e - s < 2:
            continue
        ex = float(gx[s:e].mean()); ey = float(gy[s:e].mean())
        # plane from every ground return within PLANE_R of this cell
        near = (np.abs(gx - ex) <= PLANE_R) & (np.abs(gy - ey) <= PLANE_R)
        if near.sum() < 30:
            continue
        px, py, pz = gx[near], gy[near], gz[near]
        mx, my = px.mean(), py.mean()
        A = np.c_[px - mx, py - my, np.ones(int(near.sum()))]
        co, *_ = np.linalg.lstsq(A, pz, rcond=None)
        plane_here = co[0] * (ex - mx) + co[1] * (ey - my) + co[2]
        zc = gz[s:e]
        rows.append((int(k), ex, ey,
                     float(np.median(zc)) - plane_here,
                     float(np.percentile(zc, 90)) - plane_here,
                     float(zc.max()) - plane_here,
                     int(e - s), float(100 * np.hypot(co[0], co[1]))))
    return rows


def strip_spreads(sub, strips, gp, cx, cy, half, cell):
    """Per-cell spread of the single-strip surfaces, keyed the same way."""
    n = int(round(2 * half / cell))
    per = {}
    for s in strips:
        m = sub["PointSourceId"] == s
        if m.sum() < 2000:
            continue
        g = classify(sub[m], gp)
        if not len(g):
            continue
        ix = np.floor((g["X"] - (cx - half)) / cell).astype(int)
        iy = np.floor((g["Y"] - (cy - half)) / cell).astype(int)
        ok = (ix >= 0) & (ix < n) & (iy >= 0) & (iy < n)
        key = iy[ok].astype(np.int64) * n + ix[ok]
        z = g["Z"][ok]
        order = np.argsort(key, kind="stable")
        key, z = key[order], z[order]
        uk, st = np.unique(key, return_index=True)
        en = np.append(st[1:], len(key))
        for k, a, b in zip(uk, st, en):
            if b - a < MIN_PER_STRIP:
                continue
            per.setdefault(int(k), []).append(
                (float(np.median(z[a:b])), float(np.percentile(z[a:b], 90))))
    # A median is the most stable statistic a strip can report, so it is the
    # least sensitive to a few raised returns -- which is exactly what a speck
    # is. The spread of the strips' 90th percentiles is kept beside it because
    # porous cover is penetrated to different depths by different passes, and
    # that shows in the top of the distribution before it shows in the middle.
    out = {}
    for k, v in per.items():
        if len(v) < 2:
            continue
        med = [a for a, _b in v]
        p90 = [b for _a, b in v]
        out[k] = (max(med) - min(med), max(p90) - min(p90), len(v))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--rungs", nargs="*", default=None)
    a = ap.parse_args(argv)

    cfg = PRESETS["clear"]
    gp = cfg.passes[0]
    print("ground pass: %s slope %.4f threshold %.2f m window %.0f m"
          % (gp.algorithm, gp.slope, gp.threshold_m, gp.window_m), flush=True)
    print("plane radius %.1f m, %d returns per strip to vote"
          % (PLANE_R, MIN_PER_STRIP), flush=True)

    # merged, not replaced: a single-rung run must not discard the others, the
    # same fault that made clear_fulltile.json describe the wrong tile
    prev = {}
    if os.path.exists(OUT):
        try:
            with open(OUT, encoding="utf-8") as fh:
                prev = {r["rung"]: r for r in json.load(fh)["rungs"]}
        except (ValueError, KeyError):
            prev = {}
    out, t0 = [], time.time()
    for name, cx, cy, label, half, cell in RUNGS:
        if a.rungs and name not in a.rungs:
            continue
        tile = tile_for(cx, cy)
        path = os.path.join(D, tile + ".las")
        if not os.path.exists(path):
            print("%-15s no LAS at %s" % (name, path)); continue
        pl = pdal.Pipeline(json.dumps({"pipeline": [las_reader(path)]}))
        pl.execute()
        arr = pl.arrays[0]
        nb = (np.abs(arr["X"] - cx) <= BUF) & (np.abs(arr["Y"] - cy) <= BUF)
        sub = arr[nb]
        ids, counts = np.unique(sub["PointSourceId"], return_counts=True)
        strips = [int(k) for k in ids[counts >= 0.05 * len(sub)]]

        g = classify(sub, gp)
        gxa = g["X"].astype("f8"); gya = g["Y"].astype("f8")
        gza = g["Z"].astype("f8")
        cells = cell_table(gxa, gya, gza, cx, cy, half, cell)
        spreads = strip_spreads(sub, strips, gp, cx, cy, half, cell)
        # per-return heights above the same local planes, for the cost table
        ret_h = []
        inside = (np.abs(gxa - cx) <= half) & (np.abs(gya - cy) <= half)
        for _k, ex, ey, _h, _h90, _hm, _cnt, _sl in cells:
            near = inside & (np.abs(gxa - ex) <= cell / 2) & \
                   (np.abs(gya - ey) <= cell / 2)
            if not near.any():
                continue
            nb = (np.abs(gxa - ex) <= PLANE_R) & (np.abs(gya - ey) <= PLANE_R)
            px, py, pz = gxa[nb], gya[nb], gza[nb]
            mx, my = px.mean(), py.mean()
            A = np.c_[px - mx, py - my, np.ones(int(nb.sum()))]
            co, *_ = np.linalg.lstsq(A, pz, rcond=None)
            ret_h.append(gza[near] - (co[0] * (gxa[near] - mx)
                                      + co[1] * (gya[near] - my) + co[2]))
        ret_h = np.concatenate(ret_h) if ret_h else np.zeros(0)

        rows = [{"h": h, "h90": h90, "hmax": hmax, "n": cnt, "slope_pct": sl,
                 "d": spreads[k][0] if k in spreads else None,
                 "d90": spreads[k][1] if k in spreads else None,
                 "strips": spreads[k][2] if k in spreads else 0}
                for k, ex, ey, h, h90, hmax, cnt, sl in cells]
        withd = [r for r in rows if r["d"] is not None and r["strips"] >= 3]
        print("\n=== %-15s %s ===  %.0f m square, %d cells, %d with three "
              "strips  (%.0f s)"
              % (name, label, 2 * half, len(rows), len(withd),
                 time.time() - t0), flush=True)
        if withd:
            h90 = np.array([r["h90"] for r in withd])
            dmed = np.array([r["d"] for r in withd])
            dp90 = np.array([r["d90"] for r in withd])
            print("  h of the cell's 90th pct return above the local plane: "
                  "median %+.3f  p90 %+.3f  p99 %+.3f"
                  % (np.median(h90), np.percentile(h90, 90),
                     np.percentile(h90, 99)))
            print("  %9s %7s %11s %11s"
                  % ("h90 band", "cells", "spread med", "spread p90"))
            for lo, hi in ((-9, 0.05), (0.05, 0.10), (0.10, 0.15), (0.15, 0.20),
                           (0.20, 0.30), (0.30, 0.50), (0.50, 9)):
                sel = (h90 >= lo) & (h90 < hi)
                if sel.sum() < 4:
                    continue
                print("  %4.2f-%4.2f %7d %11.3f %11.3f"
                      % (max(lo, 0), min(hi, 1.0), sel.sum(),
                         np.median(dmed[sel]), np.median(dp90[sel])))
        out.append({"rung": name, "label": label, "x": cx, "y": cy,
                    "tile": tile, "strips": strips, "half_m": half,
                    "cell_m": cell, "returns": int(ret_h.size),
                    "removed_frac": {("%.2f" % t): round(
                        float(np.mean(ret_h > t)), 5) for t in THRESHOLDS},
                    "cells": rows})
        with open(OUT, "w", encoding="utf-8") as fh:
            merged = dict(prev)
            for rr in out:
                merged[rr["rung"]] = rr
            order = [n for n, _x, _y, _l, _h, _c in RUNGS]
            json.dump({"plane_radius_m": PLANE_R,
                       "rungs": [merged[k] for k in order if k in merged]}, fh)

    # cost of each candidate threshold
    print("\n%-15s %-13s %s" % ("rung", "label",
                                " ".join("%6.2f" % t for t in THRESHOLDS)))
    print("-" * (30 + 7 * len(THRESHOLDS)))
    for r in out:
        f = r["removed_frac"]
        print("%-15s %-13s %s"
              % (r["rung"], r["label"],
                 " ".join("%5.1f%%" % (100 * f["%.2f" % t]) for t in THRESHOLDS)))
    print("\n(percent of classified-ground RETURNS standing above the "
          "threshold, measured against a %.0f m local plane)" % PLANE_R)
    print("wrote %s" % OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
