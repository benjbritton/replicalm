r"""Move the crossing that makes the black bands, and see what it costs.

Layer 2 of the G1 composite is a DIFFERENCE blend: |hillshade - slope|. Both
inputs are stretched to 0..1 and then put through brightness, contrast and
gamma, and the slope layer's gamma of 3.0 is applied as a power of 1/gamma, so
it is a cube root -- which drives steep ground hard towards 1. The hillshade on
the same ground sits a little below 1. Where the two meet, the difference is
zero and the composite goes black regardless of what the terrain is doing.

Measured on l0s395 at a rim: hillshade 0.795, slope 0.879, difference 0.084.
On ordinary ground the same subtraction gives 0.530.

The crossing point is set by that gamma. Lower it and the slope layer stops
saturating, so it never climbs to meet the hillshade and the band does not form.
This sweeps it and reports the cost at each value, with the panels to look at.

It is a deviation from Table 3, not a replication of it. Nothing here changes
the pipeline -- it renders from a DEM already on disk and writes images.
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from g1_opns_fixed import bcg, blend

# Table 3, Britton et al. 2025. Only the slope layer's gamma is varied.
LAYERS = [("MultiHS", 0, "normal", 10, 20, 2.0),
          ("Slope", 33, "difference", 10, 20, 3.0),
          ("OpnsPos", 33, "overlay", -30, 55, 1.0),
          ("SVF", 33, "multiply", 1, 30, 0.5)]
GAMMAS = [3.0, 2.0, 1.5, 1.0, 0.7, 0.5]


def layers_from_dem(dem, cell):
    import rvt.vis
    out = {}
    hs = rvt.vis.multi_hillshade(dem=dem, resolution_x=cell, resolution_y=cell,
                                 nr_directions=8, sun_elevation=35,
                                 no_data=np.nan)
    out["MultiHS"] = np.nanmean(hs, axis=0) if np.ndim(hs) == 3 else hs
    sl = rvt.vis.slope_aspect(dem=dem, resolution_x=cell, resolution_y=cell,
                              output_units="degree", no_data=np.nan)
    out["Slope"] = sl["slope"] if isinstance(sl, dict) else sl
    for name, key, svf, opns in (("SVF", "svf", True, False),
                                 ("OpnsPos", "opns", False, True)):
        r = rvt.vis.sky_view_factor(dem=dem, resolution=cell, compute_svf=svf,
                                    compute_asvf=False, compute_opns=opns,
                                    svf_n_dir=16, svf_r_max=10, svf_noise=0,
                                    no_data=np.nan)
        out[name] = r[key] if isinstance(r, dict) else r
    # The kernel layers come back a cell smaller each side than the DEM; pad
    # them back into its shape, centred, with NaN. GLiHT_rvt.py does the same.
    shape = dem.shape
    padded = {}
    for k, v in out.items():
        v = np.asarray(v, np.float32)
        if v.shape != shape:
            full = np.full(shape, np.nan, np.float32)
            ro, co = (shape[0] - v.shape[0]) // 2, (shape[1] - v.shape[1]) // 2
            full[ro:ro + v.shape[0], co:co + v.shape[1]] = v
            v = full
        padded[k] = v
    return padded


def compose(lay, slope_gamma, bounds):
    canvas, valid = None, None
    for name, transp, mode, b, c, g in LAYERS:
        a = lay[name]
        here = np.isfinite(a)
        valid = here if valid is None else (valid & here)
        lo, hi = bounds[name]
        s = np.clip((np.nan_to_num(a, nan=lo) - lo) / max(hi - lo, 1e-9), 0, 1)
        s = bcg(s, b, c, g if name != "Slope" else slope_gamma)
        canvas = s if canvas is None else blend(canvas, s, mode,
                                                1.0 - transp / 100.0)
    rgb = (np.clip(np.nan_to_num(canvas, nan=0.0), 0, 1) * 255).astype(np.uint8)
    rgb[~valid] = 0
    return rgb, valid


def load_layers(rvt_dir):
    """The four layers as GLiHT_rvt.py wrote them, full tile."""
    import rasterio
    found = {}
    for root, _dirs, files in os.walk(rvt_dir):
        for f in files:
            if not f.endswith(".tif") or "R30" in f:
                continue
            for name in ("MultiHS", "Slope", "OpnsPos", "SVF"):
                if "_%s_" % name in f:
                    found[name] = os.path.join(root, f)
    missing = [n for n in ("MultiHS", "Slope", "OpnsPos", "SVF") if n not in found]
    if missing:
        raise SystemExit("missing layers under %s: %s" % (rvt_dir, missing))
    lay, tf = {}, None
    for name, path in found.items():
        with rasterio.open(path) as r:
            arr = r.read(1).astype(np.float32)
            if r.nodata is not None:
                arr[arr == r.nodata] = np.nan
            if tf is None or (r.height, r.width) > (tf[1], tf[2]):
                tf = (r.transform, r.height, r.width)
        lay[name] = arr
    shape = (tf[1], tf[2])
    for name, v in list(lay.items()):
        if v.shape != shape:
            full = np.full(shape, np.nan, np.float32)
            ro, co = (shape[0] - v.shape[0]) // 2, (shape[1] - v.shape[1]) // 2
            full[ro:ro + v.shape[0], co:co + v.shape[1]] = v
            lay[name] = full
    return lay, tf[0]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--dem", default=None)
    ap.add_argument("--rvt-dir", default=None,
                    help="use layers already written for this tile, so the "
                         "stretch is the real whole-tile one")
    ap.add_argument("--centre", nargs=2, type=float,
                    default=[316622.34, 1970507.42])
    ap.add_argument("--side", type=float, default=120.0)
    ap.add_argument("--out", default=None)
    ap.add_argument("--tag", default="")
    a = ap.parse_args(argv)

    import glob
    import rasterio
    from rasterio.windows import Window
    from PIL import Image, ImageDraw

    cx, cy = a.centre
    if a.rvt_dir:
        lay, tf = load_layers(a.rvt_dir)
        cell = abs(tf.a)
        print("layers from %s, %d x %d at %.3f m"
              % (a.rvt_dir, lay["Slope"].shape[1], lay["Slope"].shape[0], cell),
              flush=True)
    else:
        dem = a.dem or glob.glob(
            r"C:\NR_Lidar\outputs\**\NR_316500_1970500_DEM_0p25m_v1.tif",
            recursive=True)[0]
        with rasterio.open(dem) as r:
            tf = r.transform
            cell = abs(tf.a)
            z = r.read(1).astype(np.float32)
        z[z <= 0] = np.nan
        print("%s whole tile %d x %d at %.3f m"
              % (os.path.basename(dem), z.shape[1], z.shape[0], cell), flush=True)
        lay = layers_from_dem(z, cell)

    # The stretch comes from the WHOLE tile, because that is what the recipe
    # does and it is what makes the artifact. Stretching a small window to its
    # own extremes gives contrast the delivered product does not have, and the
    # black band never appears -- which would be a fix for a problem that was
    # not there.
    bounds = {}
    for name in lay:
        v = lay[name][np.isfinite(lay[name])]
        bounds[name] = (float(v.min()), float(v.max()))
        print("  %-9s tile range %10.3f to %10.3f" % (name, *bounds[name]))

    k = int(a.side / cell)
    col0 = int((cx - tf.c) / cell) - k // 2
    row0 = int((tf.f - cy) / cell) - k // 2
    h, w = lay["Slope"].shape
    r0 = max(0, min(row0, h - k)); c0 = max(0, min(col0, w - k))
    win = (slice(r0, r0 + k), slice(c0, c0 + k))
    cut = {n: v[win] for n, v in lay.items()}
    sw = cut["Slope"]
    steep = np.isfinite(sw) & (sw > 20.0)
    flat = np.isfinite(sw) & (sw < 5.0)
    print("  window %.0f m: %d cells, %d steep, %d flat"
          % (a.side, sw.size, steep.sum(), flat.sum()), flush=True)

    print("\n%7s %9s %9s %9s %9s %9s"
          % ("gamma", "dark %", "black %", "steep mu", "steep sd", "flat mu"))
    print("-" * 58)
    panels = []
    for gam in GAMMAS:
        rgb, valid = compose(cut, gam, bounds)
        v = rgb[valid]
        st = rgb[steep & valid].astype(float)
        fl = rgb[flat & valid].astype(float)
        print("%7.1f %8.1f%% %8.1f%% %9.1f %9.2f %9.1f"
              % (gam, 100 * (v < 28).mean(), 100 * (v < 10).mean(),
                 st.mean() if st.size else float("nan"),
                 st.std() if st.size else float("nan"),
                 fl.mean() if fl.size else float("nan")), flush=True)
        panels.append((gam, rgb))

    ph, pw = panels[0][1].shape
    cols = 3
    rows = (len(panels) + cols - 1) // cols
    sheet = np.full((rows * (ph + 26) + 6, cols * (pw + 6) + 6), 255, np.uint8)
    for i, (gam, im) in enumerate(panels):
        rr, cc = divmod(i, cols)
        sheet[rr * (ph + 26) + 26:rr * (ph + 26) + 26 + ph,
              cc * (pw + 6) + 6:cc * (pw + 6) + 6 + pw] = im
    pic = Image.fromarray(sheet).convert("RGB")
    dr = ImageDraw.Draw(pic)
    for i, (gam, _im) in enumerate(panels):
        rr, cc = divmod(i, cols)
        dr.text((cc * (pw + 6) + 8, rr * (ph + 26) + 8),
                "slope gamma %.1f%s" % (gam, "   (Table 3)" if gam == 3.0 else ""),
                fill=(0, 0, 0))
    out = a.out or os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "results", "nr_block",
        "slope_gamma%s_E%d_N%d.png" % (a.tag, round(cx), round(cy)))
    os.makedirs(os.path.dirname(out), exist_ok=True)
    pic.save(out)
    print("\nwrote %s" % out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
