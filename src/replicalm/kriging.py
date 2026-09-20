"""Local ordinary kriging onto a declared grid.

WHAT THE SOURCE METHOD SPECIFIES
--------------------------------
"interpolated into a Digital Elevation Model raster using the Kriging routine in
Golden Software Surfer considering a search radius of 20 m". That is local
ordinary kriging: for each cell, take the points within the radius, solve a small
system, produce one estimate. Not global kriging, which on a million-cell tile
would be an intractable dense solve and is not what Surfer does either.

WHAT THE SOURCE METHOD DOES NOT SPECIFY
---------------------------------------
Two things it needed and did not state, both of which are set here and recorded:

  the variogram model and its parameters
      Spherical, fitted per tile over a lag equal to the kriging radius. That
      combination was chosen by measurement, not preference. On the first test
      tile an unbounded linear model put 1% of cells outside the range of the
      data, the worst of them six metres below the lowest return: a near-flat
      semivariance makes the kriging system ill-conditioned and the weights
      diverge. Spherical and exponential fitted to 20 m both gave zero cells
      outside the data range. Fitting to 60 m instead was as bad as linear,
      because the sill it finds sits at 87 m and is nearly flat across a 20 m
      neighbourhood. So the fitting lag tracks the search radius.

  the maximum number of points per solve
      A 20 m radius at these densities holds a few thousand returns. Kriging
      that many per cell is both intractable and pointless, since distant
      points carry almost no weight. Surfer caps this; the supplement does not
      say where. The cap here is explicit, defaults to 32 nearest, and is one
      of the parameters the calibration harness varies.

WHY NOT SCIKIT-LEARN OR PYKRIGE
-------------------------------
Neither is installed, both are heavy for the job, and the job is small: a KD-tree
query and a dense solve of order 32. Writing it keeps the dependency list at
NumPy and SciPy, which the environment already has, and keeps the search radius
exactly the documented figure rather than approximately whatever a library's
defaults produce.
"""
import math

import numpy as np


class KrigingError(RuntimeError):
    pass


# --- variogram models -------------------------------------------------------
# Each returns semivariance at lag h. Ordinary kriging solves in terms of these.

def _linear(h, nugget, slope, _range=None):
    return nugget + slope * h


def _spherical(h, nugget, sill, rng):
    out = np.full_like(h, nugget + sill, dtype=float)
    inside = h < rng
    hh = h[inside] / rng
    out[inside] = nugget + sill * (1.5 * hh - 0.5 * hh ** 3)
    out[h == 0] = 0.0
    return out


def _exponential(h, nugget, sill, rng):
    out = nugget + sill * (1.0 - np.exp(-3.0 * h / rng))
    out[h == 0] = 0.0
    return out


MODELS = {"linear": _linear, "spherical": _spherical,
          "exponential": _exponential}


def fit_variogram(x, y, z, model="spherical", n_lags=12, max_lag=None,
                  sample=20000, seed=0):
    """Fit an empirical variogram, so the model is measured rather than guessed.

    Subsamples before computing pairwise distances: the empirical variogram is a
    summary and does not improve past a few thousand points, whereas the pair
    count grows with the square.
    """
    rng_ = np.random.default_rng(seed)
    n = len(z)
    if n > sample:
        idx = rng_.choice(n, sample, replace=False)
        x, y, z = x[idx], y[idx], z[idx]
    from scipy.spatial import cKDTree
    if max_lag is None:
        max_lag = 20.0
    tree = cKDTree(np.c_[x, y])
    pairs = tree.query_pairs(max_lag, output_type="ndarray")
    if len(pairs) < 100:
        raise KrigingError("too few point pairs within %.1f m to fit a "
                           "variogram (%d)" % (max_lag, len(pairs)))
    d = np.hypot(x[pairs[:, 0]] - x[pairs[:, 1]],
                 y[pairs[:, 0]] - y[pairs[:, 1]])
    g = 0.5 * (z[pairs[:, 0]] - z[pairs[:, 1]]) ** 2

    edges = np.linspace(0, max_lag, n_lags + 1)
    which = np.clip(np.digitize(d, edges) - 1, 0, n_lags - 1)
    lag, semi, cnt = [], [], []
    for i in range(n_lags):
        m = which == i
        if m.sum() >= 30:
            lag.append(d[m].mean()); semi.append(g[m].mean()); cnt.append(int(m.sum()))
    lag, semi = np.array(lag), np.array(semi)
    if len(lag) < 3:
        raise KrigingError("variogram has too few populated lags to fit")

    if model == "linear":
        A = np.c_[np.ones_like(lag), lag]
        coef, *_ = np.linalg.lstsq(A, semi, rcond=None)
        nugget, slope = max(coef[0], 0.0), max(coef[1], 1e-12)
        params = (nugget, slope, None)
    else:
        sill0 = float(semi.max())
        rng0 = float(lag[-1])
        best, params = np.inf, (0.0, sill0, rng0)
        for rng_try in np.linspace(lag[1], lag[-1] * 1.5, 24):
            for nug in (0.0, 0.05 * sill0, 0.1 * sill0):
                pred = MODELS[model](lag, nug, sill0 - nug, rng_try)
                err = float(((pred - semi) ** 2).sum())
                if err < best:
                    best, params = err, (nug, sill0 - nug, rng_try)
    return {"model": model, "params": params, "lags": lag.tolist(),
            "semivariance": semi.tolist(), "counts": cnt,
            "pairs_used": int(len(pairs))}


def krige_grid(x, y, z, grid, radius=20.0, max_points=32, min_points=3,
               model="spherical", variogram=None, nodata=-9999.0,
               verbose=True):
    """Ordinary kriging of scattered points onto a declared grid.

    Cells with fewer than `min_points` neighbours inside the radius are left as
    nodata rather than extrapolated. A kriged value from two points twenty
    metres away is not a measurement, and filling it would hide the hole from
    everything downstream.
    """
    from scipy.spatial import cKDTree

    x = np.asarray(x, float); y = np.asarray(y, float); z = np.asarray(z, float)
    if len(z) < min_points:
        raise KrigingError("only %d points supplied" % len(z))

    if variogram is None:
        variogram = fit_variogram(x, y, z, model=model, max_lag=radius)
    mdl = MODELS[variogram["model"]]
    p = variogram["params"]

    tree = cKDTree(np.c_[x, y])
    gx, gy = grid.cell_centres()
    flat = np.c_[gx.ravel(), gy.ravel()]
    out = np.full(len(flat), nodata, float)

    # one query for every cell, capped at max_points, bounded by the radius
    dist, idx = tree.query(flat, k=min(max_points, len(z)),
                           distance_upper_bound=radius)
    if dist.ndim == 1:
        dist, idx = dist[:, None], idx[:, None]

    var = np.full(len(flat), nodata, float)
    filled = 0
    n_fallback = 0
    for c in range(len(flat)):
        d = dist[c]
        good = np.isfinite(d)
        k = int(good.sum())
        if k < min_points:
            continue
        ii = idx[c][good]
        px, py, pz = x[ii], y[ii], z[ii]

        # pairwise lags among the neighbours, and to the cell centre
        dxy = np.hypot(px[:, None] - px[None, :], py[:, None] - py[None, :])
        G = mdl(dxy, *p)
        g0 = mdl(d[good], *p)

        # ordinary kriging system with the unbiasedness constraint
        A = np.empty((k + 1, k + 1))
        A[:k, :k] = G
        A[:k, k] = 1.0
        A[k, :k] = 1.0
        A[k, k] = 0.0
        b = np.empty(k + 1)
        b[:k] = g0
        b[k] = 1.0
        # An ordinary kriging estimate interpolates its neighbours and cannot
        # legitimately fall outside their range. When it does, the system was
        # degenerate: this happens where the search radius is large relative to
        # the correlation range and the neighbours are packed tightly, so every
        # pairwise semivariance is nearly identical and the matrix loses rank.
        # On two of three calibration tiles -- ranges of 1.5 and 1.8 m against a
        # 3 m floor, one of them at 10 points per square metre -- the lstsq
        # fallback returned weights that produced values of order 1e15.
        #
        # The fallback is inverse-distance over the same neighbours: bounded by
        # construction, and a defensible estimate for a cell whose neighbours
        # carry no distinguishable spatial structure.
        try:
            w = np.linalg.solve(A, b)
            est = float(w[:k] @ pz)
            degenerate = not np.isfinite(est) or est < pz.min() or est > pz.max()
        except np.linalg.LinAlgError:
            degenerate = True
        if degenerate:
            dd = np.maximum(d[good], 1e-9)
            ww = 1.0 / dd ** 2
            est = float((ww @ pz) / ww.sum())
            n_fallback += 1
            # inverse distance carries no covariance model, so there is no
            # variance to report. Left as nodata rather than filled with a
            # number that would look like an uncertainty estimate.
        else:
            # ordinary kriging variance: the weights against the cell-to-point
            # semivariances, plus the Lagrange multiplier. This is the quantity
            # that distinguishes a modelled estimate from a weighted average,
            # and it is what a fallback cell does not have.
            var[c] = float(w[:k] @ g0 + w[k])
        out[c] = est
        filled += 1

    dem = out.reshape(grid.height, grid.width)
    variance = var.reshape(grid.height, grid.width)
    frac = n_fallback / filled if filled else 0.0
    if verbose:
        print("  kriged %d of %d cells (%.1f%%), radius %.1f m, "
              "%d nearest, %s variogram%s"
              % (filled, len(flat), 100 * filled / len(flat), radius,
                 max_points, variogram["model"],
                 "" if not n_fallback
                 else ", %d cells (%.1f%%) fell back to inverse distance"
                      % (n_fallback, 100 * frac)))
    if frac > 0.25:
        print("  WARNING: %.0f%% of cells were degenerate. The search radius "
              "(%.1f m) is far larger than the correlation range (%.2f m); "
              "kriging is not adding anything here."
              % (100 * frac, radius,
                 variogram["params"][2] if variogram["params"][2] else float("nan")))
    variogram = dict(variogram)
    variogram["fallback_cells"] = n_fallback
    variogram["fallback_fraction"] = frac
    variogram["variance"] = variance
    return dem, variogram


def write_geotiff(dem, grid, wkt, path, nodata=-9999.0):
    """Write the array with the grid's own transform. Nothing is recomputed."""
    from osgeo import gdal, osr
    gdal.UseExceptions()
    drv = gdal.GetDriverByName("GTiff")
    ds = drv.Create(str(path), grid.width, grid.height, 1, gdal.GDT_Float32,
                    options=["COMPRESS=DEFLATE", "TILED=YES", "PREDICTOR=3"])
    ds.SetGeoTransform(grid.geotransform)
    if wkt:
        srs = osr.SpatialReference()
        srs.ImportFromWkt(wkt)
        ds.SetProjection(srs.ExportToWkt())
    b = ds.GetRasterBand(1)
    b.SetNoDataValue(nodata)
    b.WriteArray(dem.astype("f4"))
    b.FlushCache()
    ds = None
    return str(path)


def radius_for_density(density, max_points=16, factor=4.0, floor=3.0,
                       ceiling=20.0):
    """A search radius scaled to how far away the neighbours actually are.

    WHY NOT A FIXED RADIUS
    ----------------------
    A fixed radius means different things at different point densities. At 3
    points per square metre a 20 m circle holds about 3,800 returns and the
    `max_points` cap picks the nearest sixteen, which sit within roughly 1.3 m;
    the radius does nothing and 5 m and 20 m give identical answers. At 8 points
    per square metre the same circle holds 10,000, the nearest sixteen sit
    within 0.8 m, every pairwise semivariance is nearly identical, and the
    kriging matrix loses rank: on l0s444 that put 77 to 88% of cells onto the
    inverse-distance fallback.

    The quantity that matters is the distance at which `max_points` neighbours
    are found, which is sqrt(k / (pi * density)). Multiplying by `factor` leaves
    room for the clustering that real scan lines produce, so the cap still binds
    on well-covered ground while the radius stays in proportion to the spacing.

    Bounded at both ends: below the floor too few points are found to krige at
    all, and above the ceiling the one-sided fringe at the coverage boundary
    grows faster than the estimate improves, since that fringe is one radius
    wide by construction.
    """
    d = float(density)
    if not np.isfinite(d) or d <= 0:
        return float(ceiling)
    r = float(np.sqrt(max_points / (np.pi * d)) * factor)
    return float(min(max(r, floor), ceiling))


def radius_from_variogram(variogram, floor=3.0, ceiling=50.0, factor=1.0):
    """Derive a search radius from the fitted correlation range.

    Beyond the range, two points carry no shared information, so including them
    dilutes the ones that do. A 2x2 on South_GLAS_l0s395 showed this directly:
    terrain whose fitted range was 6.89 m gave 0.041 m RMSE with a 5 m search
    and 0.080 m with the 20 m search the source method specifies -- while the
    same 20 m search was the better of the two when the variogram's range was
    28.77 m. The factors interact, and the radius belongs with the range.

    The source method's fixed 20 m is therefore not wrong so much as terrain-
    specific. It suits ground with a correlation range near 20-30 m, which this
    tile does not have.

    Bounded because a degenerate fit should not produce a degenerate search:
    below the floor too few points are found to krige, above the ceiling the
    solve grows without improving.
    """
    model = variogram.get("model")
    p = variogram.get("params") or ()
    if model in ("spherical", "exponential") and len(p) >= 3 and p[2]:
        r = float(p[2]) * factor
    else:
        # a linear model has no range; fall back to the documented radius
        r = 20.0
    return float(min(max(r, floor), ceiling))


def thin_declustered(x, y, target_density, seed=0, keep="centre"):
    """Thin to a target density by taking one point per cell of a coarse grid.

    Random thinning would preserve the clustering that is the problem: at ten
    points per square metre the sixteen nearest neighbours of a cell sit within
    about 0.7 m of each other, every pairwise semivariance is then nearly
    identical, and the kriging matrix loses rank. Removing points at random
    thins the clusters without spreading them.

    Grid thinning spreads them by construction. One point per cell of side
    1/sqrt(target) gives a set whose spacing is bounded below, which is what
    the solve needs.

    Whether this costs anything is the open question. Liao, Dong and He (2024)
    put the optimal ground density for a metre-scale DEM at 2.08 points per
    square metre, above which quality stops improving -- but that was landslide
    terrain in Sichuan, and Maya karst carries plazas, platform edges and
    retaining walls at scales their study had no reason to consider. Applying
    their figure here tests it rather than assumes it.

    `keep` chooses which point represents a cell: "centre" takes the one
    nearest the cell centre, which spreads the result most evenly; "random"
    takes any, which preserves the elevation distribution better where cells
    hold a mix of surfaces.
    """
    x = np.asarray(x, float); y = np.asarray(y, float)
    if target_density <= 0:
        raise KrigingError("target density must be positive")
    cell = 1.0 / np.sqrt(target_density)
    ix = np.floor((x - x.min()) / cell).astype(np.int64)
    iy = np.floor((y - y.min()) / cell).astype(np.int64)
    key = ix * (iy.max() + 1) + iy

    order = np.argsort(key, kind="stable")
    ks = key[order]
    starts = np.flatnonzero(np.r_[True, ks[1:] != ks[:-1]])
    ends = np.r_[starts[1:], len(ks)]

    rng = np.random.default_rng(seed)
    pick = np.empty(len(starts), dtype=np.int64)
    if keep == "centre":
        cxs = (ix[order] + 0.5) * cell + x.min()
        cys = (iy[order] + 0.5) * cell + y.min()
        d2 = (x[order] - cxs) ** 2 + (y[order] - cys) ** 2
        for i, (a, b) in enumerate(zip(starts, ends)):
            pick[i] = order[a + int(np.argmin(d2[a:b]))]
    else:
        for i, (a, b) in enumerate(zip(starts, ends)):
            pick[i] = order[a + int(rng.integers(b - a))]
    return np.sort(pick)
