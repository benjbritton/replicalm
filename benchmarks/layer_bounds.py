r"""Fixed display bounds for the four G1 source layers, from a stratified sample.

The twelve-tile measurement showed that per-tile min-max normalisation gives the
outer 4% of cells between 68% and 90% of the stretch, and that the endpoints move
by 50 to 74% of their mean from tile to tile. So the layers need bounds that are
fixed once and applied everywhere. This computes them.

Two passes. The first reads every archive DEM at 0.5 m and takes cheap statistics
-- relief and elevation scatter -- which cost a raster read each. Tiles are then
binned on relief and sampled evenly across the bins, so the sample spans flat
ground and rugged ground in proportion rather than by luck. The second pass
computes the four layers on each sampled tile and accumulates a pooled histogram
per layer.

Percentiles are taken from the pooled histogram, not averaged across per-tile
percentiles: the question is where the corpus's values lie, and a mean of
per-tile quantiles answers a different one.

Reads the archive's own DEMs, so bounds derived here can be applied to both the
archetype's renderings and ours, putting them on one scale for the first time.
Nothing is modified.
"""
import argparse
import glob
import json
import os
import sys
import time

import numpy as np

DEM_ROOT = os.environ.get("REPLICALM_DEM_ROOT", r"D:\_DEMs")
OUT_DIR = os.environ.get("REPLICALM_LAYER_STATS_DIR",
                         r"C:\Replicalm\tests\layer_bounds")
# (name, histogram range) -- ranges chosen from the layers' definitions, wide
# enough that nothing clips: hillshade and sky-view factor are [0, 1] by
# construction, slope cannot exceed 90 degrees, positive openness is an angle.
LAYERS = (("MultiHS", 0.0, 1.0), ("Slope", 0.0, 90.0),
          ("OpnsPos", 0.0, 180.0), ("SVF", 0.0, 1.0))
NBINS = 4000


def compute(vis, dem, cell):
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
    key = "svf" if vis == "SVF" else "opns"
    r = rvt.vis.sky_view_factor(dem=dem, resolution=cell,
                                compute_svf=(vis == "SVF"), compute_asvf=False,
                                compute_opns=(vis == "OpnsPos"),
                                svf_n_dir=16, svf_r_max=10, svf_noise=0,
                                no_data=np.nan)
    return r[key] if isinstance(r, dict) else r


def read_dem(path, max_cells, index=None, halo=1):
    """A DEM with kernel context, cropped to a window that actually has data.

    The halo comes from the neighbouring tile where one adjoins and shares the
    grid, and from symmetric reflection at the survey perimeter, so kernels
    compute on continuous terrain and nothing is lost at the border.

    The crop is placed on the densest part of the tile rather than at its
    geometric centre. Centring is what it was, and on one tile the centre fell
    entirely outside coverage: every layer came back empty and the failure
    surfaced three steps later as a NaN.
    """
    from halo import read_with_halo
    big, cell, transform, rep = read_with_halo(path, index=index, halo=halo)
    h, w = big.shape
    if h * w > max_cells:
        k = int(np.sqrt(max_cells))
        finite = np.isfinite(big)
        if finite.sum() < 1000:
            return big[:k, :k], cell, rep
        # centre the window on the centre of mass of the data, not of the grid
        rows, cols = np.nonzero(finite)
        r0 = int(np.clip(rows.mean() - k // 2, 0, max(0, h - k)))
        c0 = int(np.clip(cols.mean() - k // 2, 0, max(0, w - k)))
        big = big[r0:r0 + k, c0:c0 + k]
    return big, cell, rep


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=60, help="tiles in the sample")
    ap.add_argument("--bins", type=int, default=6, help="relief strata")
    ap.add_argument("--max-cells", type=int, default=9_000_000)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args(argv)
    os.makedirs(OUT_DIR, exist_ok=True)
    import rasterio

    tiles = sorted(glob.glob(os.path.join(DEM_ROOT, "**", "*0p5m*.tif"),
                             recursive=True))
    print("pass 1: relief of %d archive DEMs" % len(tiles), flush=True)
    traits, skipped, t0 = [], [], time.time()
    for i, p in enumerate(tiles, 1):
        try:
            with rasterio.open(p) as r:
                # decimated read: relief does not need every cell
                # out_shape is two-dimensional for a single-band read; a
                # three-tuple silently returns a 1-D array and every tile then
                # fails the data check without saying so.
                z = r.read(1, out_shape=(min(r.height, 512),
                                         min(r.width, 512))).astype(np.float32)
                nd = r.nodata
            if nd is not None:
                z[z == nd] = np.nan
            z[z == 0] = np.nan
            v = z[np.isfinite(z)]
            if v.size < 500:
                continue
            traits.append({"path": p,
                           "relief_m": float(np.percentile(v, 98) - np.percentile(v, 2)),
                           "sd_m": float(v.std())})
        except Exception as exc:
            # reported, not swallowed: a silent skip here once emptied the whole
            # sample and the failure only surfaced three steps later
            skipped.append((os.path.basename(p), "%s: %s"
                            % (type(exc).__name__, str(exc)[:70])))
            continue
        if i % 100 == 0:
            print("   %d/%d  %.0f s" % (i, len(tiles), time.time() - t0), flush=True)
    print("   %d tiles characterised, %d skipped, %.0f s"
          % (len(traits), len(skipped), time.time() - t0), flush=True)
    for name, why in skipped[:5]:
        print("     skipped %-44s %s" % (name[:44], why), flush=True)
    if not traits:
        raise SystemExit("pass 1 characterised nothing; see the skips above")

    # stratify on relief, sample evenly across strata
    traits.sort(key=lambda t: t["relief_m"])
    rng = np.random.default_rng(a.seed)
    strata = np.array_split(np.arange(len(traits)), a.bins)
    per = max(1, a.n // a.bins)
    picked = []
    for s in strata:
        take = min(per, len(s))
        picked += [traits[i] for i in rng.choice(s, size=take, replace=False)]
    print("\nsample: %d tiles across %d relief strata, %.1f to %.1f m relief"
          % (len(picked), a.bins, picked[0]["relief_m"],
             max(p["relief_m"] for p in picked)), flush=True)

    from halo import build_index
    print("\nindexing %d tiles for neighbour lookup" % len(traits), flush=True)
    index = build_index([t["path"] for t in traits])

    hists = {name: np.zeros(NBINS, dtype=np.int64) for name, _lo, _hi in LAYERS}
    edges = {name: np.linspace(lo, hi, NBINS + 1) for name, lo, hi in LAYERS}
    per_tile, t0 = [], time.time()
    for i, t in enumerate(picked, 1):
        try:
            dem, cell, rep = read_dem(t["path"], a.max_cells, index=index)
            if np.isfinite(dem).sum() < 10000:
                print("   skip %s: no data in the sampled window"
                      % os.path.basename(t["path"])[:44], flush=True)
                continue
            rec = {"tile": os.path.basename(t["path"]), "relief_m": t["relief_m"]}
            for name, _lo, _hi in LAYERS:
                arr = compute(name, dem, cell).astype(np.float32)
                v = arr[np.isfinite(arr)]
                if not v.size:
                    continue
                hists[name] += np.histogram(v, bins=edges[name])[0]
                rec[name] = {"min": float(v.min()), "max": float(v.max()),
                             "p2": float(np.percentile(v, 2)),
                             "p98": float(np.percentile(v, 98))}
            per_tile.append(rec)
        except Exception as exc:
            print("   skip %s: %s" % (os.path.basename(t["path"])[:40],
                                      str(exc)[:60]), flush=True)
            continue
        print("   [%2d/%d] %-44s relief %6.1f m  %4.0f s"
              % (i, len(picked), os.path.basename(t["path"])[:44],
                 t["relief_m"], time.time() - t0), flush=True)

    def pooled(name, q):
        h = hists[name]
        c = np.cumsum(h) / max(h.sum(), 1)
        j = int(np.searchsorted(c, q / 100.0))
        j = min(j, NBINS - 1)
        return float(0.5 * (edges[name][j] + edges[name][j + 1]))

    bounds = {}
    print("\npooled across %d tiles -- percentiles of the corpus, not of tiles"
          % len(per_tile))
    print("%-9s %9s %9s %9s %9s %9s %9s"
          % ("layer", "p0.5", "p1", "p2", "p98", "p99", "p99.5"))
    print("-" * 68)
    for name, _lo, _hi in LAYERS:
        qs = {q: pooled(name, q) for q in (0.5, 1, 2, 98, 99, 99.5)}
        bounds[name] = qs
        print("%-9s %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f"
              % (name, qs[0.5], qs[1], qs[2], qs[98], qs[99], qs[99.5]))

    out = {"dem_root": DEM_ROOT, "tiles_sampled": len(per_tile),
           "relief_strata": a.bins, "max_cells_per_tile": a.max_cells,
           "pooled_percentiles": bounds, "per_tile": per_tile}
    with open(os.path.join(OUT_DIR, "layer_bounds.json"), "w",
              encoding="utf-8") as fh:
        json.dump(out, fh, indent=1)
    np.savez(os.path.join(OUT_DIR, "layer_histograms.npz"),
             **{k: v for k, v in hists.items()},
             **{k + "_edges": v for k, v in edges.items()})
    print("\nwrote %s and layer_histograms.npz" % os.path.join(OUT_DIR, "layer_bounds.json"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
