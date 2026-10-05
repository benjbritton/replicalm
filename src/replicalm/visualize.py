"""The G1 composite, computed in process.

WHY THIS EXISTS
---------------
The image step used to shell out to `GLiHT_rvt.py`, found by looking beside the
interpreter, in REPLICALM_RVT_SCRIPT, or at a hardcoded `C:\\g1\\tools\\`. That
script is not in the installer, so on every user's machine the step returned
`{"skipped": "GLiHT_rvt.py not found"}`, left the DEM alone, and said nothing
an operator would notice. The checkbox did nothing.

The dependency it needed was never the problem: `rvt-py` is pinned in
`environment.yml` and ships in the installer already. Only the recipe was
missing. So the recipe lives here now, the external script is no longer
consulted, and the step works on a clean install with no new dependency and no
change to the installer's size.

WHAT IT PRODUCES
----------------
The G1 composite of Table 3, and the four visualizations it is blended from.
Nothing else: no VAT, no RVT product set. Replicalm is a ground-filtering tool
that renders one image to look at its own output, and the Relief Visualization
Toolbox remains the tool for relief visualization. The four layers are written
because they have to be computed anyway and are useful to inspect when a
composite looks wrong.

    <out>/G1/<stem>_G1.tif                  3-band 8-bit, the composite
    <out>/G1/layers/<stem>_MultiHS.tif      float32, as rvt computed them
    <out>/G1/layers/<stem>_Slope.tif
    <out>/G1/layers/<stem>_OpnsPos.tif
    <out>/G1/layers/<stem>_SVF.tif

THE RECIPE
----------
Table 3 of Britton et al. 2025, bottom layer first. Each source layer is
normalised to 0-1 over its own min and max, then brightness, contrast and gamma
are applied, then it is blended onto what is beneath it:

    MultiHS   normal      opacity 100%   B  10   C 20   G 2.0
    Slope     difference  opacity  67%   B  10   C 20   G 3.0
    OpnsPos   overlay     opacity  67%   B -30   C 55   G 1.0
    SVF       multiply    opacity  67%   B   1   C 30   G 0.5

The arithmetic is reproduced from the recipe rather than reinterpreted, because
a composite that is 18 levels darker than the archive is a difference nobody
can see the cause of, and the only defence is that the steps are the same
steps. Gamma is applied as `power(x, 1/gamma)`; min-max is per layer and per
tile, which is a known limitation -- a tile whose openness range is narrow gets
stretched the same as one whose range is wide -- and is what the published
recipe does.

RVT is Apache-2.0, ZRC SAZU and the University of Ljubljana (UL FGG). Its
functions are called as a library and are not modified. See
packaging/THIRD-PARTY-NOTICES.md.
"""
import os
import warnings

import numpy as np

# (source layer, transparency %, blend mode, brightness, contrast, gamma)
G1_LAYERS = [
    ("MultiHS",  0,  "normal",      10,  20, 2.0),
    ("Slope",    33, "difference",  10,  20, 3.0),
    ("OpnsPos",  33, "overlay",    -30,  55, 1.0),
    ("SVF",      33, "multiply",     1,  30, 0.5),
]
LAYER_NAMES = [spec[0] for spec in G1_LAYERS]

# rvt's horizon scan parameters, as the recipe sets them
SVF_DIRECTIONS = 16
SVF_RADIUS_CELLS = 10
HS_DIRECTIONS = 8
HS_SUN_ELEVATION = 35


class VisualizeError(RuntimeError):
    pass


def _rvt():
    try:
        import rvt.vis
    except ImportError as e:                       # pragma: no cover
        raise VisualizeError(
            "rvt-py is not importable, so the G1 image cannot be built. It is "
            "pinned in packaging/environment.yml; install it with "
            "`pip install rvt-py==2.2.1`. The DEM is unaffected."
        ) from e
    return rvt.vis


def normalise(arr):
    """Per-layer min-max to 0-1, ignoring NaN. Flat input returns zeros."""
    lo = np.nanmin(arr)
    hi = np.nanmax(arr)
    if not np.isfinite(lo) or not np.isfinite(hi) or hi == lo:
        return np.zeros_like(arr, dtype=np.float32)
    return np.clip((arr - lo) / (hi - lo), 0.0, 1.0).astype(np.float32)


def adjust(arr, brightness, contrast, gamma):
    """Brightness, then contrast about mid grey, then gamma as power(x, 1/g)."""
    out = np.clip(arr + brightness / 100.0, 0.0, 1.0)
    if contrast:
        out = np.clip((out - 0.5) * (1.0 + contrast / 100.0) + 0.5, 0.0, 1.0)
    if gamma != 1.0 and gamma > 0:
        out = np.power(np.clip(out, 1e-6, 1.0), 1.0 / gamma)
    return out.astype(np.float32)


def blend(base, top, mode, opacity):
    """One layer onto what is beneath it, at `opacity` in 0-1."""
    if mode == "normal":
        mixed = top
    elif mode == "multiply":
        mixed = base * top
    elif mode == "difference":
        mixed = np.abs(base - top)
    elif mode == "overlay":
        mixed = np.where(base < 0.5, 2.0 * base * top,
                         1.0 - 2.0 * (1.0 - base) * (1.0 - top))
    else:
        raise VisualizeError("unknown blend mode %r" % mode)
    return np.clip(base * (1.0 - opacity) + mixed * opacity, 0.0, 1.0)


def compute_layer(name, dem, cell):
    """One visualization, as the recipe parameterises it."""
    vis = _rvt()
    if name == "MultiHS":
        out = vis.multi_hillshade(dem=dem, resolution_x=cell, resolution_y=cell,
                                  nr_directions=HS_DIRECTIONS,
                                  sun_elevation=HS_SUN_ELEVATION,
                                  no_data=np.nan)
        # one image per direction; the recipe averages them. A column that is
        # nodata in every direction averages an empty slice, which is correct
        # and warns, so the warning is silenced rather than the NaN avoided.
        if np.ndim(out) == 3:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                out = np.nanmean(out, axis=0)
        return np.asarray(out, dtype=np.float32)
    if name == "Slope":
        out = vis.slope_aspect(dem=dem, resolution_x=cell, resolution_y=cell,
                               output_units="degree", ve_factor=1,
                               no_data=np.nan)
        return np.asarray(out["slope"], dtype=np.float32)
    if name in ("SVF", "OpnsPos"):
        want_opns = name == "OpnsPos"
        out = vis.sky_view_factor(
            dem=dem, resolution=cell,
            compute_svf=not want_opns, compute_asvf=False,
            compute_opns=want_opns,
            svf_n_dir=SVF_DIRECTIONS, svf_r_max=SVF_RADIUS_CELLS, svf_noise=0,
            no_data=np.nan)
        return np.asarray(out["opns" if want_opns else "svf"], dtype=np.float32)
    raise VisualizeError("unknown layer %r" % name)


def _centre_into(arr, shape):
    """rvt crops a border off some layers; put it back on the DEM's canvas."""
    h, w = arr.shape
    H, W = shape
    if (h, w) == (H, W):
        return arr
    if h > H or w > W:
        return arr[(h - H) // 2:(h - H) // 2 + H, (w - W) // 2:(w - W) // 2 + W]
    out = np.full(shape, np.nan, dtype=np.float32)
    r0, c0 = (H - h) // 2, (W - w) // 2
    out[r0:r0 + h, c0:c0 + w] = arr
    return out


def _read_dem(path):
    import rasterio
    with rasterio.open(str(path)) as r:
        dem = r.read(1).astype(np.float64)
        nodata = r.nodata
        profile = r.profile
        crs = r.crs
        cell = abs(r.transform.a)
    if nodata is not None:
        dem[dem == nodata] = np.nan
    dem[dem > 1e30] = np.nan
    return dem, cell, profile, crs


def _write_float(arr, path, profile, crs):
    import rasterio
    import rasterio.transform as rt
    os.makedirs(os.path.dirname(path), exist_ok=True)
    p = dict(profile)
    h, w = arr.shape
    if h != p["height"] or w != p["width"]:
        # a cropped layer keeps its own extent; move the origin to match
        t = p["transform"]
        p["transform"] = rt.from_origin(
            t.c + ((p["width"] - w) // 2) * t.a,
            t.f + ((p["height"] - h) // 2) * t.e, t.a, -t.e)
        p["height"], p["width"] = h, w
    p.update(dtype="float32", count=1, compress="lzw", crs=crs, nodata=None)
    with rasterio.open(str(path), "w", **p) as dst:
        dst.write(arr.astype(np.float32), 1)
    return str(path)


def _write_rgb(canvas, path, profile, crs):
    import rasterio
    os.makedirs(os.path.dirname(path), exist_ok=True)
    grey = (np.clip(np.nan_to_num(canvas, nan=0.0), 0, 1) * 255).astype(np.uint8)
    p = dict(profile)
    p.update(dtype="uint8", count=3, compress="lzw", crs=crs, nodata=None)
    with rasterio.open(str(path), "w", **p) as dst:
        dst.write(np.stack([grey, grey, grey], axis=0))
    return str(path)


def g1(dem_path, out_dir, keep_layers=True, progress=None, stem=None):
    """Build the G1 composite from a DEM. Returns what was written.

    `keep_layers` writes the four source visualizations beside the composite.
    They are computed either way, so turning it off saves disk and no time.
    """
    progress = progress or (lambda *a, **k: None)
    stem = stem or os.path.splitext(os.path.basename(str(dem_path)))[0]
    if stem.endswith("_DEM"):
        stem = stem[:-4]

    dem, cell, profile, crs = _read_dem(dem_path)
    g1_dir = os.path.join(out_dir, "G1")
    layer_dir = os.path.join(g1_dir, "layers")
    written = {"layers": {}}

    canvas = np.zeros(dem.shape, dtype=np.float32)
    for name, transparency, mode, brightness, contrast, gamma in G1_LAYERS:
        progress("visualize", "computing %s" % name)
        raw = compute_layer(name, dem, cell)
        if keep_layers:
            written["layers"][name] = _write_float(
                raw, os.path.join(layer_dir, "%s_%s.tif" % (stem, name)),
                profile, crs)
        layer = normalise(_centre_into(raw, dem.shape))
        layer = adjust(layer, brightness, contrast, gamma)
        canvas = blend(canvas, layer, mode, 1.0 - transparency / 100.0)

    written["g1"] = _write_rgb(canvas,
                               os.path.join(g1_dir, "%s_G1.tif" % stem),
                               profile, crs)
    progress("visualize", "G1 written: %s" % os.path.basename(written["g1"]))
    return written
