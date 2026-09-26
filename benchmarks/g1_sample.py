r"""A side-by-side sample: the published stretch against fixed bounds.

The point is to see how the two taste, not to compute how they should. Table 3's
brightness, contrast and gamma were tuned against per-tile min-max inputs; with
fixed bounds those same numbers act on a differently scaled input, and Slope's
gamma of 3.0 in particular was compensating for a squash that no longer happens
the same way. Whether the combination looks right is a judgement, so this makes
the judgement possible.

Both panels come from the same four layers computed once per tile, so nothing
differs between them except what the layers are measured against. Tiles are
chosen to span relief, because the argument for fixed bounds is that flat and
rugged ground should not both be stretched to full range -- and that claim is
only visible if both are on screen.

Written as PNG at a viewable size alongside the full-resolution GeoTIFFs.
"""
import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from g1_fixed import DEFAULT_BOUNDS, G1_LAYERS, bcg, blend, stretch
from halo import build_index, conform, strip_halo
from layer_bounds import compute, read_dem

DEM_ROOT = os.environ.get("REPLICALM_DEM_ROOT", r"D:\_DEMs")
HALO = 1
OUT = os.environ.get("REPLICALM_G1_SAMPLE",
                     r"C:\Replicalm\tests\g1_sample")


def per_tile_minmax(arr):
    """The published stretch: each layer to its own extremes, per tile."""
    v = arr[np.isfinite(arr)]
    if not v.size:
        return np.zeros_like(arr)
    lo, hi = float(v.min()), float(v.max())
    return np.zeros_like(arr) if hi <= lo else np.clip((arr - lo) / (hi - lo),
                                                       0.0, 1.0)


def pad_to(arr, h, w):
    """Centre a border-cropped layer on the full canvas.

    RVT trims a different border per algorithm -- slope loses a cell, sky-view
    factor loses more -- so the four layers come back at different sizes. The
    published script pads them to a common canvas before blending; not doing so
    is a broadcast error, or worse, a silent misalignment.
    """
    lh, lw = arr.shape
    if (lh, lw) == (h, w):
        return arr
    out = np.full((h, w), np.nan, dtype=np.float32)
    r0, c0 = (h - lh) // 2, (w - lw) // 2
    out[r0:r0 + lh, c0:c0 + lw] = arr
    return out


def composite(layers, bounds=None):
    """Assemble G1. Fixed bounds when given, per-tile extremes when not."""
    h = max(a.shape[0] for a in layers.values())
    w = max(a.shape[1] for a in layers.values())
    canvas, valid = None, None
    for name, transp, mode, b, c, g in G1_LAYERS:
        arr = pad_to(layers[name], h, w)
        here = np.isfinite(arr)
        valid = here if valid is None else (valid & here)
        if bounds is None:
            a = per_tile_minmax(np.where(here, arr, np.nan))
            a = np.nan_to_num(a, nan=0.0)
        else:
            lo, hi = bounds[name]
            a = stretch(np.nan_to_num(arr, nan=lo), lo, hi)
        a = bcg(a, b, c, g)
        canvas = a if canvas is None else blend(canvas, a, mode,
                                                1.0 - transp / 100.0)
    rgb = (np.clip(canvas, 0, 1) * 255).astype(np.uint8)
    rgb[~valid] = 0
    return rgb, valid


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--tiles", type=int, default=3)
    ap.add_argument("--max-cells", type=int, default=4_000_000)
    ap.add_argument("--stats", default=r"C:\Replicalm\tests\layer_bounds\layer_bounds.json")
    a = ap.parse_args(argv)
    os.makedirs(OUT, exist_ok=True)
    from PIL import Image
    import rasterio

    # pick tiles spanning relief, from the characterisation already done
    doc = json.load(open(a.stats, encoding="utf-8"))
    per = sorted(doc["per_tile"], key=lambda r: r["relief_m"])
    picks = [per[0], per[len(per) // 2], per[-1]][:a.tiles]
    print("bounds in use:")
    for k, (lo, hi) in DEFAULT_BOUNDS.items():
        print("  %-9s %9.3f to %9.3f" % (k, lo, hi))

    import glob
    all_dems = sorted(glob.glob(os.path.join(DEM_ROOT, "**", "*0p5m*.tif"),
                                recursive=True))
    print("indexing %d tiles for neighbour lookup" % len(all_dems), flush=True)
    index = build_index(all_dems)

    for rec in picks:
        name = rec["tile"]
        hits = glob.glob(os.path.join(DEM_ROOT, "**", name), recursive=True)
        if not hits:
            print("  %s not found" % name); continue
        t0 = time.time()
        # read with a halo so the kernels have context, compute, then strip
        # back to the tile's own grid; anything still off is conformed within
        # two cells and fatal beyond that
        dem, cell, rep = read_dem(hits[0], a.max_cells, index=index, halo=HALO)
        if np.isfinite(dem).sum() < 10000:
            print("  %s: no data in the sampled window" % name[:44]); continue
        th, tw = dem.shape[0] - 2 * HALO, dem.shape[1] - 2 * HALO
        layers = {}
        for n in DEFAULT_BOUNDS:
            arr = compute(n, dem, cell).astype(np.float32)
            # a kernel layer comes back short by its own border, so strip only
            # what the halo added and let conform settle the remainder
            off = (dem.shape[0] - arr.shape[0]) // 2
            arr = arr[off:arr.shape[0] + off - 0] if False else arr
            inner = arr[max(0, HALO - off):arr.shape[0] - max(0, HALO - off),
                        max(0, HALO - off):arr.shape[1] - max(0, HALO - off)]
            layers[n] = conform(inner, th, tw, name=n)
        print("     halo: %d cells from %d neighbours, %d reflected"
              % (rep["from_neighbours"], len(set(rep["neighbours_used"])),
                 rep["reflected"]), flush=True)
        old, _ = composite(layers, bounds=None)
        new, valid = composite(layers, bounds=DEFAULT_BOUNDS)

        stem = name.replace(".tif", "")
        for tag, img in (("published", old), ("fixed", new)):
            Image.fromarray(img).save(os.path.join(OUT, "%s_%s.png" % (stem, tag)))
        # the pair, side by side, with a divider
        h, w = old.shape
        pair = np.full((h, w * 2 + 8), 255, np.uint8)
        pair[:, :w] = old
        pair[:, w + 8:] = new
        Image.fromarray(pair).save(os.path.join(OUT, "%s_compare.png" % stem))
        print("  %-46s relief %6.1f m  mean level %3d -> %3d  %4.0f s"
              % (stem[:46], rec["relief_m"], int(old[valid].mean()),
                 int(new[valid].mean()), time.time() - t0), flush=True)

    print("\nwrote %s  (left panel published, right panel fixed bounds)" % OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
