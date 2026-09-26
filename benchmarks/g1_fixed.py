r"""The G1 composite with fixed bounds, so tiles render on one scale.

A copy of the G1 assembly in GLiHT_rvt.py, not a modification of it. That script
implements the published Table 3 recipe and produced the archive's own imagery;
changing it would alter a published product and break comparability with
everything already rendered. This sits beside it.

WHAT DIFFERS, AND WHY
---------------------
The original normalises each source layer to its own minimum and maximum within
each tile:

    lo, hi = np.nanmin(arr), np.nanmax(arr)

Measured across sixty tiles spanning six relief strata, that gives the outer 4%
of cells between 68% and 90% of the stretch, and moves the endpoints by 50 to
74% of their mean from one tile to the next. Real terrain is compressed into a
fraction of the output range -- which is the darkness -- and every tile lands on
a different scale, so brightness differences between tiles carry no information.

Here the bounds are fixed, pooled from that stratified sample at the 2nd and
98th percentiles of the corpus rather than of any tile. The blend that follows is
unchanged: the same four layers, the same Table 3 brightness, contrast, gamma,
opacity and blend modes, in the same order. Only what the layers are measured
against changes.

A consequence worth stating: flat tiles will now look flat. Under per-tile
normalisation every tile was stretched to full range, so a low-relief tile
rendered with the same apparent contrast as a rugged one. That was an artifact.

NODATA
------
The original carries NaN through to a three-channel uint8 with nodata unset, so
gaps and the trimmed edge take whatever the arithmetic produced. Here they are
written as black and declared, which is what makes the edge treatment visible
rather than accidental.
"""
import argparse
import json
import os
import sys

import numpy as np

# Pooled from 60 archive tiles across six relief strata, 3.7 to 216 m relief.
# benchmarks/layer_bounds.py rebuilds them; layer_histograms.npz holds the
# distributions if a different cut is ever wanted.
DEFAULT_BOUNDS = {
    "MultiHS": (0.525, 0.574),
    # Slope alone takes p99.5, not p98. Its upper tail climbs 39% between the
    # two where the other three have hard shoulders, and steep ground is
    # concentrated rather than spread: nine of the nineteen steepest tiles in
    # the sample have their own p98 above the pooled p98, some as high as 41
    # degrees, so a 23.8 cap flattens terrace risers and platform flanks on
    # exactly the tiles where they matter. The cost on ordinary ground is 29%
    # of the linear range but only about 11% after the gamma of 3.0 that the
    # recipe already applies to this layer.
    "Slope":   (0.484, 33.131),
    "OpnsPos": (81.517, 90.697),
    "SVF":     (0.815, 0.995),
}

# (source layer, transparency %, blend mode, brightness, contrast, gamma)
# Table 3, Britton et al. 2025 -- unchanged from the original recipe.
G1_LAYERS = [
    ("MultiHS",  0,  "normal",      10,  20, 2.0),
    ("Slope",    33, "difference",  10,  20, 3.0),
    ("OpnsPos",  33, "overlay",    -30,  55, 1.0),
    ("SVF",      33, "multiply",     1,  30, 0.5),
]


def stretch(arr, lo, hi):
    """Fixed bounds, not per-tile extremes. Values outside are clipped."""
    if hi <= lo:
        return np.zeros_like(arr)
    return np.clip((arr - lo) / (hi - lo), 0.0, 1.0)


def bcg(arr, brightness, contrast, gamma):
    """Brightness, contrast, gamma exactly as the original applies them."""
    out = np.clip(arr + brightness / 100.0, 0.0, 1.0)
    if contrast != 0:
        out = np.clip((out - 0.5) * (1.0 + contrast / 100.0) + 0.5, 0.0, 1.0)
    if gamma != 1.0 and gamma > 0:
        out = np.power(np.clip(out, 1e-6, 1.0), 1.0 / gamma)
    return out


def blend(base, top, mode, opacity):
    if mode == "normal":
        blended = top
    elif mode == "multiply":
        blended = base * top
    elif mode == "difference":
        blended = np.abs(base - top)
    elif mode == "overlay":
        blended = np.where(base < 0.5, 2 * base * top,
                           1 - 2 * (1 - base) * (1 - top))
    else:
        raise ValueError("unknown blend mode %r" % mode)
    return np.clip(base * (1.0 - opacity) + blended * opacity, 0.0, 1.0)


def build(rvt_root, bounds, out_path=None, suffix="_G1fixed"):
    """Rebuild G1 from layers already on disk under `rvt_root`."""
    import rasterio
    found = {}
    for root, _dirs, files in os.walk(rvt_root):
        for f in files:
            if not f.endswith(".tif"):
                continue
            for name, _lo, _hi in ((n, 0, 0) for n in DEFAULT_BOUNDS):
                if "_%s_" % name in f:
                    found[name] = os.path.join(root, f)
    missing = [n for n in DEFAULT_BOUNDS if n not in found]
    if missing:
        raise SystemExit("missing layers under %s: %s" % (rvt_root, missing))

    canvas, profile, valid = None, None, None
    for name, transp, mode, b, c, g in G1_LAYERS:
        with rasterio.open(found[name]) as r:
            arr = r.read(1).astype(np.float32)
            nd = r.nodata
            if profile is None:
                profile = r.profile.copy()
        if nd is not None:
            arr[arr == nd] = np.nan
        here = np.isfinite(arr)
        valid = here if valid is None else (valid & here)
        lo, hi = bounds[name]
        a = bcg(stretch(np.nan_to_num(arr, nan=lo), lo, hi), b, c, g)
        canvas = a if canvas is None else blend(canvas, a, mode,
                                                1.0 - transp / 100.0)

    rgb = (np.clip(canvas, 0, 1) * 255).astype(np.uint8)
    rgb[~valid] = 0                      # nodata is black, and declared below
    if out_path is None:
        sample = os.path.basename(found["MultiHS"])
        out_path = os.path.join(os.path.dirname(found["MultiHS"]),
                                sample.replace("_MultiHS_", suffix + "_"))
    profile.update(dtype="uint8", count=3, compress="lzw", nodata=0)
    with rasterio.open(out_path, "w", **profile) as dst:
        dst.write(np.stack([rgb, rgb, rgb], axis=0))
    return out_path, int(valid.sum()), int((~valid).sum())


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Rebuild the G1 composite with fixed, corpus-wide bounds.")
    ap.add_argument("rvt_roots", nargs="+",
                    help="directories holding the four source layer GeoTIFFs")
    ap.add_argument("--bounds", default=None,
                    help="layer_bounds.json; omit to use the pooled defaults")
    ap.add_argument("--percentiles", nargs=2, type=float, default=(2.0, 98.0),
                    metavar=("LOW", "HIGH"),
                    help="which pooled percentiles to read from --bounds")
    a = ap.parse_args(argv)

    bounds = dict(DEFAULT_BOUNDS)
    if a.bounds:
        with open(a.bounds, encoding="utf-8") as fh:
            doc = json.load(fh)
        pp = doc["pooled_percentiles"]
        lo_q, hi_q = ("%g" % a.percentiles[0], "%g" % a.percentiles[1])
        for name in bounds:
            q = pp[name]
            bounds[name] = (float(q[lo_q]), float(q[hi_q]))
    print("bounds in use:")
    for name, (lo, hi) in bounds.items():
        print("  %-9s %9.3f to %9.3f" % (name, lo, hi))

    for root in a.rvt_roots:
        try:
            out, nvalid, nnod = build(root, bounds)
        except SystemExit as exc:
            print("  %-50s %s" % (root[-50:], exc)); continue
        print("  wrote %s  (%d cells, %d nodata)"
              % (out, nvalid, nnod), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
