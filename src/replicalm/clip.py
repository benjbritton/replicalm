"""Cut matched windows from a cloud and its reference DEM.

WHY A WINDOW AND NOT THE TILE
-----------------------------
The G-LiHT clouds run from several hundred megabytes to nearly two gigabytes per
transect. A parameter sweep processes the same ground once per combination, so
sweeping a whole tile means reading two gigabytes twenty times to learn something
a few hundred metres of it would have told us. The window is cut once, both from
the cloud and from the reference, and every combination reads the small version.

THE WINDOW IS CHOSEN ON COVERAGE, NOT AT RANDOM
-----------------------------------------------
A window placed arbitrarily along a transect can land on a stretch with no
returns, and a sweep scored there measures nothing. `best_window` scans the
reference DEM for the densest covered area, which is where the comparison has
the most to say. The location is recorded so the choice is inspectable.

WHY COVERAGE ALONE WAS THE WRONG RULE
-------------------------------------
Maximum coverage is found on flat, open ground, because that is where returns
reach the floor everywhere. So the calibration windows were, without anyone
choosing it, the terrain least able to tell two ground filters apart.

Measured on South_GLAS_l0s395: below 10 degrees of slope, every configuration
tried agreed to within 0.08% of cells past half a metre. Above 30 degrees they
ranged from 2.25% to 40.38%. Ninety-eight per cent of all large error sat on
slopes above 10 degrees, which are 15.5% of the tile. A sweep scored on a
coverage-selected window is therefore scored almost entirely on the 85% of
ground where the answer does not matter, and the parameters it picked -- CSF
over SMRF, one ground pass instead of the source's two, noise filtering off --
were all wrong where it does.

`best_window` now scores a candidate on coverage AND relief together. Coverage
is still required, because a window of empty cells measures nothing; relief is
required too, because a window of flat cells measures nothing that matters.
"""
import json
import os

import numpy as np


class ClipError(RuntimeError):
    pass


def dem_coverage_profile(dem_path, win_m, step_m=None, min_coverage=0.80):
    """Coverage and relief of each candidate window along the raster.

    Returns tuples of (score, coverage, steep_fraction, row, col), sorted best
    first. `steep_fraction` is the share of a window's valid cells whose slope
    exceeds twenty degrees -- the band where ground filters begin to disagree.
    """
    from osgeo import gdal
    gdal.UseExceptions()
    d = gdal.Open(str(dem_path))
    b = d.GetRasterBand(1)
    a = b.ReadAsArray().astype("f4")
    nod = b.GetNoDataValue()
    valid = np.isfinite(a)
    if nod is not None:
        valid &= a != nod
    valid &= a != 0                      # the black-means-nodata convention
    gt = d.GetGeoTransform()
    cell = abs(gt[1])
    cy = abs(gt[5]) or cell

    # slope of the reference, so steepness is a property of the ground rather
    # than of whichever surface is being judged against it
    z = np.where(valid, a, np.nan).astype("f8")
    gyy, gxx = np.gradient(z, cy, cell)
    steep = np.isfinite(gyy) & (np.degrees(np.arctan(np.hypot(gxx, gyy))) >= 20.0)

    w = max(int(round(win_m / cell)), 8)
    step = max(int(round((step_m or win_m / 2) / cell)), 1)

    out = []
    H, W = valid.shape
    for r in range(0, max(H - w, 1), step):
        for c in range(0, max(W - w, 1), step):
            vv = valid[r:r + w, c:c + w]
            cov = float(vv.mean())
            if cov < min_coverage:
                continue
            n = int(vv.sum())
            st = float(steep[r:r + w, c:c + w].sum() / n) if n else 0.0
            # coverage is a gate, not a score: past the threshold, more of it
            # adds nothing, while more relief keeps adding discriminating power
            out.append((st, cov, st, r, c))
    if not out:
        raise ClipError("no window of %.0f m reaches %.0f%% coverage in %s"
                        % (win_m, 100 * min_coverage, dem_path))
    out.sort(reverse=True)
    return out, gt, cell, w


def _window_bounds(gt, cell, w, r, c):
    minx = gt[0] + c * cell
    maxy = gt[3] - r * cell
    return {"minx": minx, "maxy": maxy,
            "maxx": minx + w * cell, "miny": maxy - w * cell,
            "cell": cell, "window_m": w * cell}


def ground_density(las_path, bounds, tmp_dir=None, cell_m=1.0):
    """Ground returns per square metre inside a window, by a fast classification.

    Coarse on purpose: one SMRF pass on a 1 m working grid, which is too coarse
    to render with and ample to answer "does this window have a floor". The
    gate is a factor of three away from the values it has to separate, so it
    does not need the production settings.
    """
    import tempfile
    from . import classify
    from .config import ReplicalmConfig, GroundPass

    tmp_dir = tmp_dir or tempfile.gettempdir()
    stem = "dens_%d_%d" % (int(bounds["minx"]), int(bounds["maxy"]))
    win = os.path.join(tmp_dir, stem + ".las")
    gnd = os.path.join(tmp_dir, stem + "_g.las")
    try:
        if not os.path.exists(gnd):
            clip_las(las_path, win, bounds, verbose=False)
            cfg = ReplicalmConfig(
                passes=[GroundPass(algorithm="smrf", slope=0.1584,
                                   threshold_m=0.5, cell_m=cell_m)],
                remove_low_noise=False, remove_outliers=False)
            classify.classify_tile(win, gnd, cfg, verbose=False)
        arr, _ = classify.read_points(gnd)
        n = int((arr["Classification"] == 2).sum())
    finally:
        for p in (win,):
            if os.path.exists(p):
                try:
                    os.remove(p)
                except OSError:
                    pass
    area = (bounds["maxx"] - bounds["minx"]) * (bounds["maxy"] - bounds["miny"])
    return (n / area) if area > 0 else 0.0


def best_window(dem_path, win_m=400.0, min_coverage=0.80, prefer_relief=True,
                las_path=None, min_ground_density=0.0, min_steep_fraction=0.0,
                max_candidates=8, tmp_dir=None, verbose=False):
    """World coordinates of the window a sweep has most to learn from.

    With `prefer_relief` the choice is the steepest window clearing the coverage
    gate. Set it False for the old coverage-only behaviour, kept so earlier
    results can be regenerated rather than silently orphaned.

    With `las_path` and `min_ground_density`, candidates are walked in order and
    the first one whose ground returns clear the floor is taken. Relief with no
    returns beneath it measures the survey rather than the filter, which is how
    l8s431 entered a calibration at 0.92 ground points per square metre and
    failed on 93% of its steep cells under every configuration tried.
    """
    prof, gt, cell, w = dem_coverage_profile(dem_path, win_m,
                                             min_coverage=min_coverage)
    if not prefer_relief:
        prof = sorted(prof, key=lambda t: -t[1])
    tried = []
    for _, cov, steep, r, c in prof[:max(1, max_candidates)]:
        if steep < min_steep_fraction:
            continue
        b = _window_bounds(gt, cell, w, r, c)
        dens = None
        if las_path and min_ground_density > 0:
            try:
                dens = ground_density(las_path, b, tmp_dir=tmp_dir)
            except Exception as e:
                if verbose:
                    print("  density check failed at %d,%d: %s" % (c, r, str(e)[:40]))
                continue
            tried.append((dens, cov, steep))
            if dens < min_ground_density:
                if verbose:
                    print("  window at %d,%d: %.2f gnd/m2, below the %.2f floor"
                          % (c, r, dens, min_ground_density))
                continue
        b.update({"coverage": cov, "steep_fraction": steep,
                  "ground_density": dens,
                  "selected_for": "relief" if prefer_relief else "coverage"})
        return b
    raise ClipError(
        "no %.0f m window in %s clears coverage %.0f%%, steep %.1f%% and "
        "ground density %.2f/m2 (best densities tried: %s)"
        % (win_m, os.path.basename(str(dem_path)), 100 * min_coverage,
           100 * min_steep_fraction, min_ground_density,
           ", ".join("%.2f" % d for d, _, _ in tried) or "none"))


def clip_las(las_path, out_path, bounds, verbose=True):
    """Crop a point cloud to a bounding box, writing an uncompressed LAS."""
    import json as _json
    import pdal
    pipeline = {"pipeline": [
        str(las_path),
        {"type": "filters.crop",
         "bounds": "([%f,%f],[%f,%f])" % (bounds["minx"], bounds["maxx"],
                                          bounds["miny"], bounds["maxy"])},
        {"type": "writers.las", "filename": str(out_path),
         "minor_version": 2, "dataformat_id": 1, "compression": "false"}]}
    pl = pdal.Pipeline(_json.dumps(pipeline))
    n = pl.execute()
    if not n:
        raise ClipError("no points fall inside the window for %s" % las_path)
    if verbose:
        print("  clipped %-42s %8d points" % (os.path.basename(str(las_path)), n))
    return str(out_path), int(n)


def clip_raster(dem_path, out_path, bounds, verbose=True):
    """Crop a raster to the same box, without resampling."""
    from osgeo import gdal
    gdal.UseExceptions()
    gdal.Translate(str(out_path), str(dem_path),
                   projWin=[bounds["minx"], bounds["maxy"],
                            bounds["maxx"], bounds["miny"]],
                   projWinSRS=None, noData=None)
    if verbose:
        d = gdal.Open(str(out_path))
        print("  clipped %-42s %4d x %-4d cells"
              % (os.path.basename(str(dem_path)), d.RasterXSize, d.RasterYSize))
    return str(out_path)


def prepare(tile, out_dir, win_m=400.0, verbose=True):
    """Cut a matched window from one tile's cloud and reference DEM.

    `tile` is an entry from calibration_tiles.json. Returns the paths and the
    window, written beside the outputs so the selection is recoverable.
    """
    os.makedirs(out_dir, exist_ok=True)
    win = best_window(tile["reference_dem"], win_m)
    if verbose:
        print("%s %s: window %.0f m at %.0f, %.0f -- %.1f%% covered"
              % (tile["region"], tile["tile"], win["window_m"],
                 win["minx"], win["maxy"], 100 * win["coverage"]))
    las_out = os.path.join(out_dir, "%s_clip.las" % tile["tile"])
    dem_out = os.path.join(out_dir, "%s_reference.tif" % tile["tile"])
    _, n = clip_las(tile["las"], las_out, win, verbose)
    clip_raster(tile["reference_dem"], dem_out, win, verbose)
    rec = {"tile": tile["tile"], "region": tile["region"], "window": win,
           "source_las": tile["las"], "source_dem": tile["reference_dem"],
           "clipped_las": las_out, "clipped_reference": dem_out,
           "points": n}
    with open(os.path.join(out_dir, "%s_window.json" % tile["tile"]), "w",
              encoding="utf-8") as fh:
        json.dump(rec, fh, indent=1)
    return rec
