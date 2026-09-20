"""Tune the translation against DEMs the original workflow already produced.

WHY CALIBRATION IS NOT OPTIONAL HERE
------------------------------------
The classification parameters in config.py are a translation, not a port.
TerraScan's window, terrain angle, iteration angle and iteration distance have
no counterparts in SMRF, and the defaults were derived by reasoning about what
each term controls -- tan(9 deg) for a slope tolerance, and so on. Reasoning is
where a translation starts, not where it ends. This measures it.

The reference is the DEM the TerraScan and Surfer pipeline produced for the same
ground. Agreement with that is the only evidence available that the substitution
preserved the result.

THE SCORE, AND WHY NOT RMSE ALONE
---------------------------------
Elevation RMSE rewards smoothing. A ground filter that strips the mounds and
platforms out of the terrain produces a smoother surface that agrees better with
anything, and scores well, while destroying exactly the signal this data exists
to carry. Liao, Dong and He (2024, Remote Sensing 16:4563) address the same
problem for landslide micro-topography and score DEM quality on elevation RMSE
combined with a terrain-complexity term, so that losing detail is penalised
rather than rewarded.

The score here follows that structure:

    quality = w_rmse * normalised_rmse + w_complexity * complexity_loss

where complexity is measured as the standard deviation of local slope, and
complexity_loss is how much of the reference's complexity the candidate failed
to reproduce. Both terms are dimensionless and lower is better.

THE DENSITY FLOOR
-----------------
Liao et al. also report the ground point density at which DEM quality stops
improving: 2.08 pts/m2 for a 1 m product, 2.43 at 0.2 m, 1.84 at 2 m. Those
figures come from landslide terrain in Sichuan, not Maya karst, so they are used
here as a reported reference rather than as a law -- a tile falling below the
floor is flagged, not rejected. `density_sweep` runs their dilution experiment
on your own data, which is the way to find out whether the number transfers.
"""
import itertools
import json
import os
import time
from dataclasses import replace

import numpy as np


class CalibrationError(RuntimeError):
    pass


# Liao et al. 2024, Remote Sensing 16:4563, Table: optimal ground point density
# by DEM resolution. Reported, not assumed to transfer.
LIAO_OPTIMAL_DENSITY = {0.2: 2.43, 0.5: 2.08, 1.0: 2.08, 2.0: 1.84}


def read_dem(path):
    """A DEM and its grid, as arrays, with nodata masked."""
    from osgeo import gdal
    gdal.UseExceptions()
    d = gdal.Open(str(path))
    b = d.GetRasterBand(1)
    a = b.ReadAsArray().astype("f8")
    nod = b.GetNoDataValue()
    m = np.isfinite(a)
    if nod is not None:
        m &= a != nod
    return a, m, d.GetGeoTransform(), (d.RasterXSize, d.RasterYSize)


def terrain_complexity(dem, mask, cell):
    """Standard deviation of local slope, the detail a smoother would remove.

    Slope rather than elevation, because elevation variation is dominated by the
    regional surface and would swamp the metre-scale relief that matters.
    """
    z = np.where(mask, dem, np.nan)
    gy, gx = np.gradient(z, cell)
    slope = np.sqrt(gx ** 2 + gy ** 2)
    s = slope[np.isfinite(slope)]
    if s.size < 10:
        return float("nan")
    return float(np.nanstd(s))


def compare_to_reference(candidate_path, reference_path, w_rmse=0.5,
                         w_complexity=0.5):
    """Score one candidate DEM against the reference, on shared cells only.

    Resamples nothing. If the two grids do not overlap on a common lattice the
    comparison is refused, because a score computed across a half-cell offset
    measures the offset.
    """
    ca, cm, cgt, csz = read_dem(candidate_path)
    ra, rm, rgt, rsz = read_dem(reference_path)

    if abs(cgt[1] - rgt[1]) > 1e-6:
        raise CalibrationError(
            "cell sizes differ by %.2e m (%.6f vs %.6f). Build the candidate "
            "with grid.grid_from_raster(reference) so the comparison is cell "
            "for cell." % (abs(cgt[1] - rgt[1]), cgt[1], rgt[1]))
    if abs(abs(cgt[5]) - abs(rgt[5])) > 1e-6:
        raise CalibrationError("north-south cell sizes differ: %.6f vs %.6f"
                               % (abs(cgt[5]), abs(rgt[5])))
    cell = abs(cgt[1])
    dx = (cgt[0] - rgt[0]) / cell
    dy = (rgt[3] - cgt[3]) / abs(cgt[5])
    if abs(dx - round(dx)) > 1e-3 or abs(dy - round(dy)) > 1e-3:
        raise CalibrationError(
            "grids are offset by a fraction of a cell (%.3f, %.3f); a score "
            "would measure the offset rather than the surface" % (dx, dy))

    # align by integer cell offset
    ox, oy = int(round(dx)), int(round(dy))
    rx0, cx0 = max(0, oy), max(0, ox)
    ry0, cy0 = max(0, -oy), max(0, -ox)
    h = min(rsz[1] - rx0, csz[1] - ry0)
    w = min(rsz[0] - cx0, csz[0] - cy0)
    if h <= 0 or w <= 0:
        raise CalibrationError("candidate and reference do not overlap")

    R = ra[rx0:rx0 + h, cx0:cx0 + w]
    Rm = rm[rx0:rx0 + h, cx0:cx0 + w]
    C = ca[ry0:ry0 + h, cy0:cy0 + w]
    Cm = cm[ry0:ry0 + h, cy0:cy0 + w]
    both = Rm & Cm
    if both.sum() < 100:
        raise CalibrationError("only %d cells in common" % int(both.sum()))

    diff = C[both] - R[both]
    rmse = float(np.sqrt((diff ** 2).mean()))
    bias = float(diff.mean())
    mae = float(np.abs(diff).mean())

    rc = terrain_complexity(R, Rm, cell)
    cc = terrain_complexity(C, Cm, cell)
    # how much of the reference's detail the candidate failed to reproduce;
    # negative means it produced MORE structure, which is noise, so it is
    # penalised by magnitude either way
    complexity_loss = abs(rc - cc) / rc if rc and np.isfinite(rc) else float("nan")

    # normalise RMSE by the reference's own relief so the two terms combine
    relief = float(np.nanstd(R[Rm]))
    nrmse = rmse / relief if relief else float("nan")

    score = w_rmse * nrmse + w_complexity * complexity_loss
    return {"rmse_m": rmse, "bias_m": bias, "mae_m": mae,
            "nrmse": nrmse, "relief_m": relief,
            "complexity_reference": rc, "complexity_candidate": cc,
            "complexity_loss": complexity_loss,
            "cells_compared": int(both.sum()),
            "coverage_candidate": float(Cm.mean()),
            "coverage_reference": float(Rm.mean()),
            "score": score}


def check_density(n_ground, area_m2, cell):
    """Ground density against Liao's reported optimum for this resolution."""
    d = n_ground / area_m2 if area_m2 else 0.0
    want = LIAO_OPTIMAL_DENSITY.get(cell)
    out = {"density_pts_m2": d, "liao_optimum": want}
    if want:
        out["ratio"] = d / want
        out["below_optimum"] = bool(d < want)
    return out


def sweep(las_path, reference_dem, out_dir, base_cfg, param_grid,
          cell=1.0, radius=20.0, max_points=32, verbose=True):
    """Run the pipeline once per parameter combination and score each result.

    `param_grid` maps a GroundPass field to the values to try, for example
    {"slope": [0.1, 0.158, 0.213], "threshold_m": [1.0, 2.0, 3.0]}. Every
    combination is applied to every pass, since the translation sets both from
    the same reasoning.
    """
    from . import classify, grid as G, kriging as K, interpolate

    os.makedirs(out_dir, exist_ok=True)
    keys = sorted(param_grid)
    combos = list(itertools.product(*(param_grid[k] for k in keys)))
    if verbose:
        print("%d combinations over %s" % (len(combos), ", ".join(keys)))

    # Build on the reference's own lattice. Its cells are 0.500042 m on an
    # unsnapped origin, so a production grid would sit a fraction of a cell
    # away and the comparison would measure that offset rather than the
    # surface. Production output still uses the declared grid; calibration
    # borrows the reference's so the difference is cell for cell.
    g = G.grid_from_raster(reference_dem)
    wkt = interpolate.source_srs(las_path)
    area = (g.bounds[2] - g.bounds[0]) * (g.bounds[3] - g.bounds[1])
    rows = []

    for n, values in enumerate(combos, 1):
        tag = "_".join("%s%g" % (k[:4], v) for k, v in zip(keys, values))
        las_out = os.path.join(out_dir, "gnd_%s.las" % tag)
        dem_out = os.path.join(out_dir, "dem_%s.tif" % tag)
        cfg = replace(base_cfg,
                      passes=[replace(p, **dict(zip(keys, values)))
                              for p in base_cfg.passes])
        t0 = time.time()
        try:
            c = classify.classify_tile(las_path, las_out, cfg, verbose=False)
            arr, _ = classify.read_points(las_out)
            gnd = arr[arr["Classification"] == 2]
            if len(gnd) < cfg.min_ground_points:
                raise CalibrationError("only %d ground points" % len(gnd))
            dem, vg = K.krige_grid(gnd["X"], gnd["Y"], gnd["Z"], g,
                                   radius=radius, max_points=max_points,
                                   verbose=False)
            K.write_geotiff(dem, g, wkt, dem_out)
            sc = compare_to_reference(dem_out, reference_dem)
            dens = check_density(len(gnd), area, cell)
            row = {"combination": dict(zip(keys, values)), "tag": tag,
                   "ground_points": int(len(gnd)), "seconds": round(time.time() - t0, 1),
                   **dens, **sc}
        except Exception as e:
            row = {"combination": dict(zip(keys, values)), "tag": tag,
                   "error": "%s: %s" % (type(e).__name__, e)}
        rows.append(row)
        if verbose:
            if "error" in row:
                print("  %2d/%d  %-28s FAILED %s" % (n, len(combos), tag,
                                                     row["error"][:46]))
            else:
                print("  %2d/%d  %-28s score %.4f  rmse %.3f m  "
                      "complexity loss %.3f  %.1f pts/m2"
                      % (n, len(combos), tag, row["score"], row["rmse_m"],
                         row["complexity_loss"], row["density_pts_m2"]))

    ok = [r for r in rows if "error" not in r]
    ok.sort(key=lambda r: r["score"])
    result = {"las": str(las_path), "reference": str(reference_dem),
              "cell_m": cell, "radius_m": radius, "max_points": max_points,
              "ranked": ok, "failed": [r for r in rows if "error" in r]}
    with open(os.path.join(out_dir, "sweep.json"), "w", encoding="utf-8") as fh:
        json.dump(result, fh, indent=1, default=float)
    if verbose and ok:
        b = ok[0]
        print("\nbest: %s  score %.4f  (rmse %.3f m, complexity loss %.3f)"
              % (b["tag"], b["score"], b["rmse_m"], b["complexity_loss"]))
    return result


def density_sweep(las_path, out_dir, fractions=(1.0, 0.75, 0.5, 0.35, 0.25,
                                                0.15, 0.1, 0.05),
                  cell=1.0, radius=20.0, seed=0, verbose=True):
    """Liao's dilution experiment, on your own terrain.

    Thins the ground returns to each fraction, builds a DEM from each, and
    scores every one against the DEM built from the full set. The knee in that
    curve is the density below which quality falls away -- Liao's 2.08 pts/m2
    for a 1 m product is the equivalent figure for landslide terrain, and this
    is how to find out whether Maya karst agrees.
    """
    from . import classify, grid as G, kriging as K, interpolate

    os.makedirs(out_dir, exist_ok=True)
    arr, _ = classify.read_points(las_path)
    gnd = arr[arr["Classification"] == 2]
    if not len(gnd):
        raise CalibrationError("no ground-classified points in %s" % las_path)
    g = G.grid_for_las(las_path, cell=cell)
    wkt = interpolate.source_srs(las_path)
    area = (g.bounds[2] - g.bounds[0]) * (g.bounds[3] - g.bounds[1])
    rng = np.random.default_rng(seed)

    full_path = os.path.join(out_dir, "dilute_1.00.tif")
    dem, _ = K.krige_grid(gnd["X"], gnd["Y"], gnd["Z"], g, radius=radius,
                          verbose=False)
    K.write_geotiff(dem, g, wkt, full_path)
    if verbose:
        print("reference from all %d ground points (%.2f pts/m2)"
              % (len(gnd), len(gnd) / area))
        print("%8s %10s %10s %9s %9s" % ("fraction", "points", "pts/m2",
                                         "rmse", "cplx loss"))

    rows = []
    for f in fractions:
        if f >= 1.0:
            continue
        k = max(int(len(gnd) * f), 10)
        idx = rng.choice(len(gnd), k, replace=False)
        sub = gnd[idx]
        p = os.path.join(out_dir, "dilute_%.2f.tif" % f)
        try:
            d2, _ = K.krige_grid(sub["X"], sub["Y"], sub["Z"], g,
                                 radius=radius, verbose=False)
            K.write_geotiff(d2, g, wkt, p)
            sc = compare_to_reference(p, full_path)
            row = {"fraction": f, "points": int(k),
                   "density_pts_m2": k / area, **sc}
        except Exception as e:
            row = {"fraction": f, "points": int(k),
                   "error": "%s: %s" % (type(e).__name__, e)}
        rows.append(row)
        if verbose:
            if "error" in row:
                print("%8.2f %10d  %s" % (f, k, row["error"][:44]))
            else:
                print("%8.2f %10d %10.2f %9.3f %9.3f"
                      % (f, k, row["density_pts_m2"], row["rmse_m"],
                         row["complexity_loss"]))

    out = {"las": str(las_path), "cell_m": cell,
           "full_density_pts_m2": len(gnd) / area,
           "liao_optimum_for_cell": LIAO_OPTIMAL_DENSITY.get(cell),
           "levels": rows}
    with open(os.path.join(out_dir, "density_sweep.json"), "w",
              encoding="utf-8") as fh:
        json.dump(out, fh, indent=1, default=float)
    return out


def sweep_interpolation(las_path, reference_dem, out_dir, cfg,
                        radii=(5.0, 10.0, 20.0), max_points=(8, 16, 32, 64),
                        verbose=True):
    """Vary the kriging while holding the classification fixed.

    Classification is the expensive stage and does not depend on these
    parameters, so it runs once and every combination reuses its ground points.
    The previous sweep repeated it identically for each combination, which cost
    most of the nine minutes and told us nothing.

    Radius and neighbour count are what control smoothing, and smoothing is what
    the complexity term measures. A candidate can match the reference to seven
    centimetres while reproducing under half its terrain detail -- the elevation
    is right and the relief that carries the archaeology is not.
    """
    from . import classify, grid as G, kriging as K, interpolate

    os.makedirs(out_dir, exist_ok=True)
    g = G.grid_from_raster(reference_dem)
    wkt = interpolate.source_srs(las_path)

    gnd_path = os.path.join(out_dir, "ground.las")
    c = classify.classify_tile(las_path, gnd_path, cfg, verbose=False)
    arr, _ = classify.read_points(gnd_path)
    gnd = arr[arr["Classification"] == 2]
    area = (g.bounds[2] - g.bounds[0]) * (g.bounds[3] - g.bounds[1])
    if verbose:
        print("classified once: %d ground points, %.2f pts/m2"
              % (len(gnd), len(gnd) / area))
        print("%8s %7s %9s %9s %9s %9s %9s"
              % ("radius", "points", "score", "rmse m", "cplx loss",
                 "cplx cand", "filled%"))

    ref_c = None
    rows = []
    for radius in radii:
        vg = None
        for k in max_points:
            tag = "r%g_k%d" % (radius, k)
            dem_out = os.path.join(out_dir, "dem_%s.tif" % tag)
            t0 = time.time()
            try:
                if vg is None:
                    vg = K.fit_variogram(gnd["X"], gnd["Y"], gnd["Z"],
                                         max_lag=radius)
                dem, _ = K.krige_grid(gnd["X"], gnd["Y"], gnd["Z"], g,
                                      radius=radius, max_points=k,
                                      variogram=vg, verbose=False)
                K.write_geotiff(dem, g, wkt, dem_out)
                sc = compare_to_reference(dem_out, reference_dem)
                ref_c = sc["complexity_reference"]
                row = {"radius_m": radius, "max_points": k, "tag": tag,
                       "seconds": round(time.time() - t0, 1), **sc}
            except Exception as e:
                row = {"radius_m": radius, "max_points": k, "tag": tag,
                       "error": "%s: %s" % (type(e).__name__, e)}
            rows.append(row)
            if verbose:
                if "error" in row:
                    print("%8.0f %7d  %s" % (radius, k, row["error"][:52]))
                else:
                    print("%8.0f %7d %9.4f %9.3f %9.3f %9.4f %8.1f%%"
                          % (radius, k, row["score"], row["rmse_m"],
                             row["complexity_loss"], row["complexity_candidate"],
                             100 * row["coverage_candidate"]))

    ok = [r for r in rows if "error" not in r]
    ok.sort(key=lambda r: r["score"])
    out = {"las": str(las_path), "reference": str(reference_dem),
           "ground_points": int(len(gnd)),
           "complexity_reference": ref_c, "ranked": ok,
           "failed": [r for r in rows if "error" in r]}
    with open(os.path.join(out_dir, "interp_sweep.json"), "w",
              encoding="utf-8") as fh:
        json.dump(out, fh, indent=1, default=float)
    if verbose and ok:
        b = ok[0]
        print("" + chr(10) + "reference complexity %.4f" % ref_c)
        print("best: radius %.0f m, %d nearest -- score %.4f, rmse %.3f m, "
              "complexity %.4f vs %.4f"
              % (b["radius_m"], b["max_points"], b["score"], b["rmse_m"],
                 b["complexity_candidate"], ref_c))
    return out
