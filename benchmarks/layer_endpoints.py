r"""How much do the G1 source layers' endpoints move from tile to tile?

The published recipe normalises each of the four source visualizations to its own
minimum and maximum within each tile, then applies fixed brightness, contrast and
gamma. If those endpoints are stable across terrain, that per-tile normalisation
is nearly harmless and a single post-composite transfer function is enough to fix
the darkness. If they swing -- particularly if single cells are setting them --
then tiles are differently compressed before they are blended, and no curve
applied afterwards can make them comparable.

So this measures, for each layer, across a sample of archive DEMs:

  absolute min and max     what the recipe currently uses as endpoints
  2nd and 98th percentile  what it would use if outliers were excluded
  the gap between them     how much of the range one extreme cell is claiming

Read on the archive's own 0.5 m DEMs rather than ours, because if a fixed stretch
is adopted it should be one both products can share.

Nothing is modified. This reads DEMs and computes layers in memory.
"""
import argparse
import glob
import json
import os
import random
import sys
import time

import numpy as np

DEM_ROOT = os.environ.get("REPLICALM_DEM_ROOT", r"D:\_DEMs")
OUT = os.environ.get("REPLICALM_LAYER_STATS",
                     r"C:\Replicalm\tests\layer_endpoints\layer_endpoints.json")
LAYERS = ("MultiHS", "Slope", "OpnsPos", "SVF")


def compute(vis, dem, cell):
    """The four source layers, by the same calls the published recipe makes."""
    import rvt.vis
    if vis == "MultiHS":
        out = rvt.vis.multi_hillshade(dem=dem, resolution_x=cell,
                                      resolution_y=cell, nr_directions=8,
                                      sun_elevation=35, no_data=np.nan)
        return np.nanmean(out, axis=0) if out.ndim == 3 else out
    if vis == "Slope":
        r = rvt.vis.slope_aspect(dem=dem, resolution_x=cell, resolution_y=cell,
                                 output_units="degree", no_data=np.nan)
        return r["slope"] if isinstance(r, dict) else r
    if vis == "SVF":
        r = rvt.vis.sky_view_factor(dem=dem, resolution=cell, compute_svf=True,
                                    compute_opns=False, compute_asvf=False,
                                    svf_n_dir=16, svf_r_max=10, svf_noise=0,
                                    no_data=np.nan)
        return r["svf"] if isinstance(r, dict) else r
    if vis == "OpnsPos":
        r = rvt.vis.sky_view_factor(dem=dem, resolution=cell, compute_svf=False,
                                    compute_opns=True, compute_asvf=False,
                                    svf_n_dir=16, svf_r_max=10, svf_noise=0,
                                    no_data=np.nan)
        return r["opns"] if isinstance(r, dict) else r
    raise ValueError(vis)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=12, help="tiles to sample")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--max-cells", type=int, default=4_000_000,
                    help="cap per tile; layers are O(cells) and SVF is slow")
    a = ap.parse_args(argv)

    import rasterio
    tiles = sorted(glob.glob(os.path.join(DEM_ROOT, "**", "*0p5m*.tif"),
                             recursive=True))
    if not tiles:
        raise SystemExit("no 0.5 m DEMs found under %s" % DEM_ROOT)
    random.Random(a.seed).shuffle(tiles)
    tiles = tiles[:a.n]
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    print("sampling %d of the archive's 0.5 m DEMs\n" % len(tiles), flush=True)

    rows = []
    for i, path in enumerate(tiles, 1):
        t0 = time.time()
        with rasterio.open(path) as r:
            dem = r.read(1).astype(np.float32)
            nd = r.nodata
            cell = abs(r.transform.a)
        if nd is not None:
            dem[dem == nd] = np.nan
        dem[dem == 0] = np.nan           # the nodata convention in these products
        # crop centrally if very large: the endpoints are a distribution question
        h, w = dem.shape
        if h * w > a.max_cells:
            k = int(np.sqrt(a.max_cells))
            r0, c0 = max(0, (h - k) // 2), max(0, (w - k) // 2)
            dem = dem[r0:r0 + k, c0:c0 + k]
        if np.isfinite(dem).sum() < 10000:
            print("  skip %s (too little data)" % os.path.basename(path)[:44])
            continue

        rec = {"tile": os.path.basename(path), "cells": int(np.isfinite(dem).sum())}
        for lyr in LAYERS:
            try:
                arr = compute(lyr, dem, cell).astype(np.float32)
            except Exception as exc:
                rec[lyr] = {"error": str(exc)[:80]}
                continue
            v = arr[np.isfinite(arr)]
            if not v.size:
                rec[lyr] = {"error": "all nodata"}
                continue
            lo, hi = float(v.min()), float(v.max())
            p2, p98 = (float(np.percentile(v, 2)), float(np.percentile(v, 98)))
            rec[lyr] = {"min": lo, "max": hi, "p2": p2, "p98": p98,
                        # what share of the recipe's range sits outside the
                        # percentile band -- the part one extreme cell claims
                        "tail_share": round(((p2 - lo) + (hi - p98)) / (hi - lo), 4)
                        if hi > lo else 0.0}
        rows.append(rec)
        print("  [%2d/%d] %-46s %5.0f s" % (i, len(tiles),
                                            os.path.basename(path)[:46],
                                            time.time() - t0), flush=True)
        json.dump(rows, open(OUT, "w"), indent=1)

    print("\n%-9s %10s %10s %10s %10s %9s"
          % ("layer", "min", "max", "p2", "p98", "tail %"))
    print("-" * 64)
    for lyr in LAYERS:
        vals = [r[lyr] for r in rows if lyr in r and "min" in r[lyr]]
        if not vals:
            continue
        for label, key in (("min", "min"), ("max", "max"),
                           ("p2", "p2"), ("p98", "p98")):
            pass
        mn = np.array([v["min"] for v in vals])
        mx = np.array([v["max"] for v in vals])
        p2 = np.array([v["p2"] for v in vals])
        p98 = np.array([v["p98"] for v in vals])
        tl = np.array([v["tail_share"] for v in vals])
        print("%-9s %10.3f %10.3f %10.3f %10.3f %8.1f%%"
              % (lyr, mn.mean(), mx.mean(), p2.mean(), p98.mean(), 100 * tl.mean()))
        print("%-9s %10s %10s %10s %10s"
              % ("  spread", "%.3f" % np.ptp(mn), "%.3f" % np.ptp(mx),
                 "%.3f" % np.ptp(p2), "%.3f" % np.ptp(p98)))
        rng = mx - mn
        prng = p98 - p2
        print("%-9s min-max range varies %.1f%% of its mean; "
              "p2-p98 range varies %.1f%%"
              % ("", 100 * np.ptp(rng) / rng.mean() if rng.mean() else 0,
                 100 * np.ptp(prng) / prng.mean() if prng.mean() else 0))
    print("\nwrote %s" % OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
