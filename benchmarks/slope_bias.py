r"""Does Clear's cleanup demote ground for being on a hillside?

Clear measures each ground return against the 10th percentile of raw elevations
within 0.75 m. On a slope that percentile sits at the downhill edge of the
disc, so an upslope point is compared against ground that is genuinely lower
than it. The disc spans about 1.5*tan(theta) of elevation; the percentile sits
near its bottom; so a point at the upslope rim stands roughly 1.2*tan(theta)
above the reference on terrain alone. Against the 0.20 m threshold that gives
about 9.5 degrees, past which Clear starts removing ground for being on a hill.

That is arithmetic. This measures it.

Terrain slope comes from the ARCHIVE's own reference DEM, not from ours, so the
test does not use our own output to judge our own output. Slope is taken at the
scale of the cleanup patch -- the DEM smoothed over 3 m before the gradient --
because slope measured across a single 0.5 m cell is mostly noise and is not
what the rule responds to.

  --stage window     find a window with both flat and steep ground
  --stage classify   clip it, run SMRF, keep the ground cloud on disk
  --stage measure    demotion rate against slope, for each rule

Everything is read-only with respect to the pipeline: no configuration is
changed and no existing file is written.
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import SRC, TESTS
sys.path.insert(0, SRC)

OUT = os.path.join(TESTS, "slope_bias")
WINDOW_M = 500.0
BUFFER_M = 120.0          # SMRF's window runs to 70 m; give it room
SMOOTH_M = 3.0            # slope at the scale the cleanup patch sees


def reference_slope(ref_path, smooth_m=SMOOTH_M):
    """Slope in degrees from the archive DEM, smoothed to the patch scale."""
    import rasterio
    from scipy import ndimage
    with rasterio.open(ref_path) as r:
        z = r.read(1).astype(np.float32)
        nd, tf = r.nodata, r.transform
    if nd is not None:
        z[z == nd] = np.nan
    z[z == 0] = np.nan
    cell = abs(tf.a)
    ok = np.isfinite(z)
    # a boxcar over NaN would poison every cell near a gap, so smooth the sum
    # and the count separately and divide
    k = max(1, int(round(smooth_m / cell)))
    zf = np.where(ok, z, 0.0)
    num = ndimage.uniform_filter(zf, k, mode="nearest")
    den = ndimage.uniform_filter(ok.astype(np.float32), k, mode="nearest")
    sm = np.where(den > 0.3, num / np.maximum(den, 1e-6), np.nan)
    gy, gx = np.gradient(sm, cell)
    slope = np.degrees(np.arctan(np.hypot(gx, gy)))
    slope[~ok] = np.nan
    return slope, tf, cell


def sample(raster, tf, cell, x, y):
    """Nearest-cell lookup, NaN outside."""
    c = ((x - tf.c) / cell).astype(np.int64)
    r = ((tf.f - y) / cell).astype(np.int64)
    out = np.full(len(x), np.nan, np.float32)
    good = (r >= 0) & (r < raster.shape[0]) & (c >= 0) & (c < raster.shape[1])
    out[good] = raster[r[good], c[good]]
    return out


def stage_window(ref_path, n=6):
    """Rank candidate windows by how much steep ground they hold."""
    slope, tf, cell = reference_slope(ref_path)
    k = int(WINDOW_M / cell)
    h, w = slope.shape
    rows = []
    for r0 in range(0, h - k, k // 2):
        for c0 in range(0, w - k, k // 2):
            s = slope[r0:r0 + k, c0:c0 + k]
            v = s[np.isfinite(s)]
            if v.size < 0.6 * k * k:
                continue
            steep = float((v > 15).mean())
            flat = float((v < 5).mean())
            # want both, not one: the test needs a control population
            rows.append((min(steep, flat) * 2, steep, flat,
                         float(tf.c + (c0 + k / 2) * cell),
                         float(tf.f - (r0 + k / 2) * cell)))
    rows.sort(reverse=True)
    print("%4s %8s %8s %12s %12s" % ("rank", ">15deg", "<5deg", "easting", "northing"))
    for i, (_s, st, fl, x, y) in enumerate(rows[:n], 1):
        print("%4d %7.1f%% %7.1f%% %12.2f %12.2f" % (i, 100 * st, 100 * fl, x, y))
    return rows[:n]


def stage_classify(las_path, centre, tag):
    import pdal
    from replicalm import classify
    from replicalm.config import PRESETS
    os.makedirs(OUT, exist_ok=True)
    cfg = PRESETS["clear"]
    half = WINDOW_M / 2 + BUFFER_M
    x, y = centre
    b = (x - half, y - half, x + half, y + half)
    clip = os.path.join(OUT, "%s_clip.las" % tag)
    pl = pdal.Pipeline(json.dumps({"pipeline": [
        classify.las_reader(las_path),
        {"type": "filters.crop",
         "bounds": "([%f,%f],[%f,%f])" % (b[0], b[2], b[1], b[3])},
        {"type": "writers.las", "filename": clip, "forward": "all"}]}))
    npts = pl.execute()
    print("clipped %d points into %s" % (npts, clip), flush=True)
    ground = os.path.join(OUT, "%s_ground.laz" % tag)
    classify.classify_tile(clip, ground, cfg, verbose=True)
    print("wrote %s" % ground, flush=True)
    return ground


def stage_measure(ground_path, ref_path, centre, tag):
    from replicalm import classify, cleanup
    arr, _ = classify.read_points(ground_path)
    g = arr[arr["Classification"] == 2]
    x = g["X"].astype("f8"); y = g["Y"].astype("f8"); z = g["Z"].astype("f8")
    # drop the buffer: it exists to feed SMRF, not to be measured
    half = WINDOW_M / 2
    keep = ((np.abs(x - centre[0]) <= half) & (np.abs(y - centre[1]) <= half))
    print("%d ground returns, %d inside the %.0f m window"
          % (len(x), keep.sum(), WINDOW_M), flush=True)

    print("computing height above the 10th percentile floor ...", flush=True)
    above = cleanup.height_above_floor(x, y, z, patch=0.75, percentile=10.0)

    slope, tf, cell = reference_slope(ref_path)
    s = sample(slope, tf, cell, x, y)

    m = keep & np.isfinite(s)
    s, above = s[m], above[m]
    print("%d returns with a reference slope" % m.sum())

    edges = [0, 3, 6, 9, 12, 15, 20, 25, 30, 90]
    print("")
    print("%-12s %10s %10s %10s %10s"
          % ("slope band", "returns", "demoted %", "median h", "p90 h"))
    print("-" * 56)
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = (s >= lo) & (s < hi)
        if sel.sum() < 200:
            continue
        a = above[sel]
        dem = 100.0 * (a > 0.20).mean()
        rows.append({"lo": lo, "hi": hi, "n": int(sel.sum()),
                     "demoted_pct": round(float(dem), 3),
                     "median_h": round(float(np.median(a)), 4),
                     "p90_h": round(float(np.percentile(a, 90)), 4)})
        print("%3d - %-6d %10d %9.2f%% %10.4f %10.4f"
              % (lo, hi, sel.sum(), dem, np.median(a), np.percentile(a, 90)))
    with open(os.path.join(OUT, "%s_slope_bias.json" % tag), "w",
              encoding="utf-8") as fh:
        json.dump({"tag": tag, "window_m": WINDOW_M, "centre": list(centre),
                   "reference": ref_path, "bands": rows}, fh, indent=1)
    print("")
    print("wrote %s" % os.path.join(OUT, "%s_slope_bias.json" % tag))
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=("window", "classify", "measure"),
                    required=True)
    ap.add_argument("--tile", default="l0s444")
    ap.add_argument("--centre", nargs=2, type=float, default=None)
    ap.add_argument("--tag", default=None)
    a = ap.parse_args(argv)
    from clear_fulltile import TILES
    las, ref = TILES[a.tile]
    tag = a.tag or a.tile
    if a.stage == "window":
        stage_window(ref)
        return 0
    if a.centre is None:
        raise SystemExit("--centre EASTING NORTHING is required for this stage")
    if a.stage == "classify":
        stage_classify(las, tuple(a.centre), tag)
    else:
        stage_measure(os.path.join(OUT, "%s_ground.laz" % tag), ref,
                      tuple(a.centre), tag)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
