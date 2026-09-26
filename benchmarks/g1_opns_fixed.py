r"""EXPERIMENT, 2026-09-26: fixed corpus bounds on the openness layer only.

Nothing existing is modified. This is a separate renderer; reverting means
deleting this file. The four source layers are read from an rvt directory that
GLiHT_rvt.py already wrote, so no DEM is recomputed and the terrain under test
is byte-for-byte the terrain in the current images.

WHAT IT CHANGES, AND THE ONE REASON
-----------------------------------
Layer 3 is positive openness, stretched per tile between its own minimum and
maximum, then brightness -30 and contrast 55. On l0s395 that stretch runs
46.39 to 140.52 degrees, while 96% of the corpus's openness values lie between
81.5 and 90.7 -- nine degrees inside ninety-four. Ordinary ground therefore
lands near 0.39 to 0.43 of the stretch, the brightness subtracts 0.30, the
contrast multiplies the distance from 0.5 by 1.55, and the result goes
negative and clips. Measured on the window at E752163 N2099157: openness
transforms to exactly 0.000 on dark cells and on mid-tone cells alike. The
layer carries no information at all.

That is what produces the truncation. Overlay with a top layer of zero is flat
below base 0.5 -- output is zero whatever the base does -- and rises as
2*base - 1 above it. Every cell whose running composite falls below the hinge
receives the same contribution, so the transfer curve has no gradient there.

Fixing this layer's bounds to the pooled p2 and p98 puts ordinary ground near
the middle of the stretch, where brightness -30 and contrast 55 leave it
inside range, and the layer contributes again.

The other three layers keep per-tile min-max, unchanged, so the experiment has
one variable. --control renders all four per-tile for a side-by-side.
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Pooled from 60 archive tiles across six relief strata; benchmarks/layer_bounds.py
OPNS_P2, OPNS_P98 = 81.517, 90.697

# Table 3, Britton et al. 2025 -- unchanged
G1_LAYERS = [
    ("MultiHS",  0,  "normal",      10,  20, 2.0),
    ("Slope",    33, "difference",  10,  20, 3.0),
    ("OpnsPos",  33, "overlay",    -30,  55, 1.0),
    ("SVF",      33, "multiply",     1,  30, 0.5),
]


def minmax(arr):
    lo, hi = np.nanmin(arr), np.nanmax(arr)
    if hi == lo:
        return np.zeros_like(arr), float(lo), float(hi)
    return np.clip((arr - lo) / (hi - lo), 0.0, 1.0), float(lo), float(hi)


def bcg(arr, brightness, contrast, gamma):
    out = np.clip(arr + brightness / 100.0, 0.0, 1.0)
    if contrast != 0:
        out = np.clip((out - 0.5) * (1.0 + contrast / 100.0) + 0.5, 0.0, 1.0)
    if gamma != 1.0 and gamma > 0:
        out = np.power(np.clip(out, 1e-6, 1.0), 1.0 / gamma)
    return out.astype(np.float32)


def blend(base, top, mode, opacity):
    if mode == "normal":
        blended = top
    elif mode == "multiply":
        blended = base * top
    elif mode == "difference":
        blended = np.abs(base - top)
    elif mode == "overlay":
        blended = np.where(base < 0.5, 2.0 * base * top,
                           1.0 - 2.0 * (1.0 - base) * (1.0 - top))
    else:
        raise ValueError("unknown blend mode %r" % mode)
    return np.clip(base * (1.0 - opacity) + blended * opacity, 0.0, 1.0)


def find_layers(rvt_dir):
    found = {}
    for root, _dirs, files in os.walk(rvt_dir):
        for f in files:
            if not f.endswith(".tif"):
                continue
            for name in ("MultiHS", "Slope", "OpnsPos", "SVF"):
                if "_%s_" % name in f:
                    found[name] = os.path.join(root, f)
    missing = [n for n in ("MultiHS", "Slope", "OpnsPos", "SVF") if n not in found]
    if missing:
        raise SystemExit("missing layers under %s: %s" % (rvt_dir, missing))
    return found


def render(rvt_dir, out_png, opns="corpus", opns_file=None):
    import rasterio
    from PIL import Image
    found = find_layers(rvt_dir)
    if opns_file:
        # a positive-openness layer computed at another search radius
        found["OpnsPos"] = opns_file

    # The kernel layers come off disk a cell or two smaller than the slope
    # layer, trimmed at the border. GLiHT_rvt.py pads each one back into the
    # DEM's own shape, centred, with NaN; the same is done here against the
    # largest layer present, so nothing is resampled and the padding falls
    # outside coverage where it is written black.
    shapes = {}
    for name in ("MultiHS", "Slope", "OpnsPos", "SVF"):
        with rasterio.open(found[name]) as r:
            shapes[name] = (r.height, r.width)
    shape = (max(s[0] for s in shapes.values()), max(s[1] for s in shapes.values()))

    canvas, valid, slope_raw = None, None, None
    print("layer      stretch source            low       high    mean out")
    print("-" * 66)
    for name, transp, mode, b, c, g in G1_LAYERS:
        with rasterio.open(found[name]) as r:
            arr = r.read(1).astype(np.float32)
            nd = r.nodata
        if nd is not None:
            arr[arr == nd] = np.nan
        if arr.shape != shape:
            lh, lw = arr.shape
            full = np.full(shape, np.nan, np.float32)
            ro, co = (shape[0] - lh) // 2, (shape[1] - lw) // 2
            full[ro:ro + lh, co:co + lw] = arr
            arr = full
        here = np.isfinite(arr)
        valid = here if valid is None else (valid & here)

        if name == "OpnsPos" and opns != "minmax":
            if opns == "corpus":
                lo, hi, src = OPNS_P2, OPNS_P98, "corpus p2-p98"
            else:
                # "p2-98": this tile's own percentiles, so the recipe stays
                # per-tile and only the endpoint definition changes
                p_lo, p_hi = (float(t) for t in opns[1:].split("-"))
                v = arr[np.isfinite(arr)]
                lo, hi = (float(np.percentile(v, p_lo)),
                          float(np.percentile(v, p_hi)))
                src = "per-tile p%g-p%g" % (p_lo, p_hi)
            s = np.clip((arr - lo) / max(hi - lo, 1e-9), 0.0, 1.0)
        else:
            s, lo, hi = minmax(arr)
            src = "per-tile min-max"
        if name == "Slope":
            slope_raw = arr
        s = bcg(s, b, c, g)
        m = float(np.nanmean(s[here])) if here.any() else float("nan")
        print("%-9s  %-18s %9.3f %10.3f %10.4f" % (name, src, lo, hi, m))
        canvas = s if canvas is None else blend(canvas, s, mode,
                                                1.0 - transp / 100.0)

    rgb = (np.clip(np.nan_to_num(canvas, nan=0.0), 0, 1) * 255).astype(np.uint8)
    rgb[~valid] = 0
    Image.fromarray(rgb).save(out_png)
    v = rgb[valid]
    print("")
    print("%d x %d, %d data cells" % (shape[0], shape[1], int(valid.sum())))
    print("value mean %.1f   below 28: %.2f%%   below 10: %.3f%%"
          % (v.mean(), 100.0 * (v < 28).mean(), 100.0 * (v < 10).mean()))
    if slope_raw is not None:
        steep = valid & (slope_raw > 20.0)
        flat = valid & (slope_raw < 5.0)
        if steep.sum() > 1000 and flat.sum() > 1000:
            st = rgb[steep].astype(np.float64)
            fl = rgb[flat].astype(np.float64)
            print("steep >20deg: %d cells, mean %.1f, sd %.2f, iqr %.1f"
                  % (steep.sum(), st.mean(), st.std(),
                     np.percentile(st, 75) - np.percentile(st, 25)))
            print("flat  < 5deg: %d cells, mean %.1f   separation %.1f"
                  % (flat.sum(), fl.mean(), fl.mean() - st.mean()))
    print("wrote %s" % out_png)
    return out_png


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("rvt_dir")
    ap.add_argument("out_png")
    ap.add_argument("--opns-file", default=None,
                    help="use this positive-openness raster instead of the "
                         "one in rvt_dir (e.g. a wider search radius)")
    ap.add_argument("--opns", default="corpus",
                    help="openness stretch: minmax (the current build), "
                         "corpus (pooled p2-p98), or a per-tile percentile "
                         "pair such as p2-98 or p1-99")
    a = ap.parse_args(argv)
    os.makedirs(os.path.dirname(os.path.abspath(a.out_png)), exist_ok=True)
    render(a.rvt_dir, a.out_png, opns=a.opns, opns_file=a.opns_file)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
