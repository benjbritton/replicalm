r"""Is the black rim on steep flanks the cleanup, or the renderer?

Two mechanisms put dark in the same place and need opposite remedies.

  the cleanup   `height_above_floor` measures each return against the 10th
                percentile of a 0.75 m radius, which on a slope sits at the
                downhill side. From about 9 degrees the upslope fringe of every
                neighbourhood exceeds the 0.20 m threshold on terrain alone, so
                Clear removes ground that is not vegetation.

  the renderer  Sky-view factor collapses at the foot of a scarp, and layer 4
                squares it -- gamma 0.5, applied as a power of 1/gamma -- then
                multiplies. A correct DEM still renders black there.

So: build the same window from Freckles, which has no cleanup, and from Clear,
and put them side by side. If the rim is in both, the cleanup is innocent.

The two panels share one stretch, taken from the Clear window and applied to
both. Rendering each to its own min and max is what the published recipe does
per tile, and doing it here would give the two panels different scales and make
the comparison meaningless -- a trap this project has fallen into before.

Two stages, because PDAL and RVT live in different environments:
  --stage clip    in the interpreter that has PDAL: locate the feature, clip,
                  and process the window through Freckles and Clear
  --stage render  in the interpreter that has RVT: build both composites
"""
import argparse
import glob
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import TESTS

OUT = os.path.join(TESTS, "slope_artifact")
FULL = os.path.join(TESTS, "clear_fulltile")
WINDOW_M = 400.0        # side of the window put under the microscope
BUFFER_M = 120.0        # context for SMRF, whose window runs to 70 m


def find_dark_feature(g1_path):
    """Map coordinates of the largest saturated-dark cluster inside coverage."""
    from osgeo import gdal
    from scipy import ndimage
    gdal.UseExceptions()
    # the dataset is held while the band is read: chaining
    # gdal.Open(...).GetRasterBand(1) frees it underneath the band
    ds = gdal.Open(g1_path)
    band = ds.GetRasterBand(1).ReadAsArray()
    gt = ds.GetGeoTransform()
    ds = None
    data = band > 0                                  # 0 is nodata here
    # erode the coverage mask so the boundary's own darkness is excluded: the
    # question is about interior scarps, not about the edge
    inside = ndimage.binary_erosion(data, np.ones((41, 41), bool))
    dark = (band < 28) & inside
    lab, n = ndimage.label(dark)
    if not n:
        return None
    sizes = ndimage.sum(dark, lab, range(1, n + 1))
    big = int(np.argmax(sizes)) + 1
    r, c = ndimage.center_of_mass(dark, lab, big)
    x = gt[0] + (c + 0.5) * gt[1]
    y = gt[3] + (r + 0.5) * gt[5]
    # sizes is zero-based and `big` is already argmax + 1
    return {"x": float(x), "y": float(y), "cells": int(sizes[big - 1])}


def stage_clip(tiles, at=None):
    import pdal
    from dataclasses import replace
    from replicalm import pipeline
    from replicalm.config import PRESETS
    from replicalm.classify import las_reader

    # The source LAS comes from clear_fulltile's TILES, which is derived from
    # tests/clips/index.json. Not from clear_fulltile.json: that file is
    # rewritten by whichever run finishes last, and it currently names the wrong
    # l0s395 because a three-tile run overwrote a corrected single-tile one.
    from clear_fulltile import TILES
    with open(os.path.join(FULL, "clear_fulltile.json"), encoding="utf-8") as fh:
        done = {r["tile"]: r for r in json.load(fh)}
    os.makedirs(OUT, exist_ok=True)
    # merged, not replaced: rebuilding this from empty on every invocation is
    # how a single-tile run silently discards what earlier runs recorded
    wj = os.path.join(OUT, "windows.json")
    found = json.load(open(wj, encoding="utf-8")) if os.path.exists(wj) else {}
    for tile in tiles:
        g1 = glob.glob(os.path.join(FULL, tile, "rvt", "**", "*_G1_*.tif"),
                       recursive=True)
        if not g1:
            print("%-8s no G1 to locate a feature in" % tile); continue
        if at and tile in at:
            spot = {"x": at[tile][0], "y": at[tile][1], "cells": -1}
        else:
            spot = find_dark_feature(g1[0])
        if not spot:
            print("%-8s no dark cluster found" % tile); continue
        cell = done[tile]["cell_m"] if tile in done else 0.5
        half = WINDOW_M / 2 + BUFFER_M
        b = (spot["x"] - half, spot["y"] - half, spot["x"] + half, spot["y"] + half)
        las = TILES[tile][0]
        clip = os.path.join(OUT, "%s_E%d_N%d.las"
                            % (tile, round(spot["x"]), round(spot["y"])))
        print("%-8s feature at %.1f, %.1f (%d dark cells); clipping %.0f m box"
              % (tile, spot["x"], spot["y"], spot["cells"], 2 * half), flush=True)
        pl = pdal.Pipeline(json.dumps({"pipeline": [
            las_reader(las),
            {"type": "filters.crop",
             "bounds": "([%f,%f],[%f,%f])" % (b[0], b[2], b[1], b[3])},
            {"type": "writers.las", "filename": clip, "forward": "all"}]}))
        n = pl.execute()
        print("         %d points in the window" % n, flush=True)
        if n < 20000:
            print("         too few to process"); continue
        for profile in ("baseline", "clear"):
            d = os.path.join(OUT, "%s_E%d_N%d_%s"
                             % (tile, round(spot["x"]), round(spot["y"]), profile))
            os.makedirs(d, exist_ok=True)
            # same cell for both, so only the cleanup differs
            pipeline.process(clip, d, cfg=PRESETS[profile], cell_m=cell,
                             make_g1=False, keep_ground=False,
                             progress=lambda s, m, f=None: None)
            print("         %-8s DEM written" % profile, flush=True)
        key = "%s_E%d_N%d" % (tile, round(spot["x"]), round(spot["y"]))
        found[key] = {"tile": tile, "centre": spot, "cell_m": cell,
                      "window_m": WINDOW_M}
    with open(os.path.join(OUT, "windows.json"), "w", encoding="utf-8") as fh:
        json.dump(found, fh, indent=1)
    print("\nwrote %s" % os.path.join(OUT, "windows.json"))


def stage_render(tiles):
    import rasterio
    from PIL import Image
    from g1_fixed import G1_LAYERS, bcg, blend
    from layer_bounds import compute

    with open(os.path.join(OUT, "windows.json"), encoding="utf-8") as fh:
        windows = json.load(fh)

    def layers_for(dem_path, centre, side):
        with rasterio.open(dem_path) as r:
            arr = r.read(1).astype(np.float32)
            nd, tf = r.nodata, r.transform
        if nd is not None:
            arr[arr == nd] = np.nan
        arr[arr == 0] = np.nan
        cell = abs(tf.a)
        c = int((centre["x"] - tf.c) / cell)
        rr = int((tf.f - centre["y"]) / cell)
        k = int(side / cell / 2)
        sl = arr[max(0, rr - k):rr + k, max(0, c - k):c + k]
        return {n: compute(n, sl, cell).astype(np.float32) for n in
                ("MultiHS", "Slope", "OpnsPos", "SVF")}, sl

    for key, w in windows.items():
        tile = w.get("tile", key)
        panels, shared = {}, {}
        for profile in ("baseline", "clear"):
            dem = glob.glob(os.path.join(OUT, "%s_%s" % (key, profile),
                                         "*_DEM.tif"))
            if not dem:
                print("%-8s %s: no DEM" % (tile, profile)); continue
            panels[profile], _ = layers_for(dem[0], w["centre"], w["window_m"])
        if len(panels) != 2:
            continue
        # one stretch for both panels, taken from Clear
        for n in ("MultiHS", "Slope", "OpnsPos", "SVF"):
            v = panels["clear"][n]
            v = v[np.isfinite(v)]
            shared[n] = (float(v.min()), float(v.max())) if v.size else (0.0, 1.0)

        imgs = {}
        for profile, lay in panels.items():
            h = min(a.shape[0] for a in lay.values())
            wd = min(a.shape[1] for a in lay.values())
            canvas, valid = None, None
            for name, transp, mode, b, c, g in G1_LAYERS:
                a = lay[name][:h, :wd]
                here = np.isfinite(a)
                valid = here if valid is None else (valid & here)
                lo, hi = shared[name]
                s = np.clip((np.nan_to_num(a, nan=lo) - lo) / max(hi - lo, 1e-9),
                            0, 1)
                s = bcg(s, b, c, g)
                canvas = s if canvas is None else blend(canvas, s, mode,
                                                        1.0 - transp / 100.0)
            rgb = (np.clip(canvas, 0, 1) * 255).astype(np.uint8)
            rgb[~valid] = 0
            imgs[profile] = rgb
            Image.fromarray(rgb).save(os.path.join(
                OUT, "%s_%s.png" % (key, "freckles" if profile == "baseline" else "clear")))

        a_, b_ = imgs["baseline"], imgs["clear"]
        h = min(a_.shape[0], b_.shape[0]); wd = min(a_.shape[1], b_.shape[1])
        pair = np.full((h, wd * 2 + 8), 255, np.uint8)
        pair[:, :wd] = a_[:h, :wd]
        pair[:, wd + 8:] = b_[:h, :wd]
        Image.fromarray(pair).save(os.path.join(OUT, "%s_freckles_vs_clear.png" % key))
        dk = lambda im: 100.0 * (im[im > 0] < 28).mean()
        print("%-26s dark cells: freckles %5.2f%%   clear %5.2f%%"
              % (key, dk(a_), dk(b_)), flush=True)
    print("\nwrote %s  (left panel Freckles, right panel Clear)" % OUT)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=("clip", "render"), required=True)
    ap.add_argument("--tiles", nargs="*", default=["l0s395", "l0s444"])
    ap.add_argument("--at", nargs="*", default=None,
                    help="tile=easting,northing to centre a window explicitly")
    a = ap.parse_args(argv)
    if a.stage == "clip":
        sys.path.insert(0, os.path.join(os.path.dirname(
            os.path.dirname(os.path.abspath(__file__))), "src"))
        at = None
        if a.at:
            at = {}
            for item in a.at:
                t, xy = item.split("=")
                at[t] = tuple(float(v) for v in xy.split(","))
        stage_clip(a.tiles, at=at)
    else:
        stage_render(a.tiles)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
