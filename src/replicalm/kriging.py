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
import os

import numpy as np


def _query_workers():
    """Cores for the neighbour query: all of them unless told otherwise."""
    try:
        return int(os.environ.get("REPLICALM_QUERY_WORKERS", "-1"))
    except ValueError:
        return -1


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
                  sample=20000, seed=0, min_nugget_fraction=0.0):
    """Fit an empirical variogram, so the model is measured rather than guessed.

    Subsamples before computing pairwise distances: the empirical variogram is a
    summary and does not improve past a few thousand points, whereas the pair
    count grows with the square.

    `min_nugget_fraction` floors the nugget at a fraction of the total sill.

    WHY A FLOOR IS AVAILABLE AT ALL
    -------------------------------
    The fit chooses among nuggets of 0, 5% and 10% of the sill, and on this data
    it picks exactly zero more often than not: 34 of 60 fits across ten tiles. A
    variogram with no nugget makes ordinary kriging an exact interpolator -- it
    honors every observation, so per-return ranging noise becomes single-cell
    spikes in the raster. At four points per square meter with sixteen neighbors
    inside about 1.2 m, each cell is decided by one or two returns and there is
    nothing to average against.

    Roughness on flat ground, against the reference surface from the original
    workflow, measured on the l0s395 trench window:

        0%  1.26x    2%  1.09x    4%  0.99x    6%  0.92x   20%  0.72x
        1%  1.16x    3%  1.04x    5%  0.96x   10%  0.83x

    Four per cent is where the output carries the same micro-relief the original
    workflow did; below it the surface is sharper than the reference, above it
    smoother.

    The default is 0.0, so nothing changes unless a caller asks. This is not in
    the locked baseline: it was measured on one window of one tile, and the
    visual half of that comparison is not yet valid, because G1 normalizes each
    raster by its own extremes and one outlier cell restretches the whole image.
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
        # The floor is applied after fitting rather than inside the search, so
        # the reported fit stays the fit and the floor stays a stated choice.
        if min_nugget_fraction > 0:
            nug, psill, rng_ = params
            total = nug + psill
            floor = min_nugget_fraction * total
            if nug < floor:
                params = (floor, total - floor, rng_)
    return {"model": model, "params": params, "lags": lag.tolist(),
            "semivariance": semi.tolist(), "counts": cnt,
            "pairs_used": int(len(pairs)),
            "min_nugget_fraction": float(min_nugget_fraction)}


def krige_grid(x, y, z, grid, radius=20.0, max_points=32, min_points=3,
               model="spherical", variogram=None, nodata=-9999.0,
               verbose=True, chunk_cells=1000000, solve_batch=50000,
               ridge_fraction=1e-6):
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

    # DUPLICATE LOCATIONS MAKE THE SYSTEM EXACTLY SINGULAR
    #
    # Two returns at one XY give the kriging matrix two identical rows, which no
    # nugget, ridge or neighbour count can rescue -- the matrix has no inverse,
    # and every cell whose neighbourhood contains such a pair falls out of the
    # solve. The delivered G-LiHT tiles carry them: 13.05% of the returns on one
    # tile are exact XYZ repeats, same return number, same GPS time, same source,
    # while another tile has none. It is a property of the archive, not of the
    # processing, and it accounted for 89% of that tile's cells abandoning
    # kriging for inverse distance.
    #
    # Exact repeats carry no information, so removing them costs nothing. Where
    # the same XY holds different elevations -- 0.1% of groups -- the lowest is
    # kept, because a ground surface cannot pass through two heights at one
    # place and the lower return is the ground candidate.
    order = np.lexsort((z, y, x))
    x, y, z = x[order], y[order], z[order]
    first = np.ones(len(z), bool)
    if len(z) > 1:
        first[1:] = (x[1:] != x[:-1]) | (y[1:] != y[:-1])
    n_dropped = int(len(z) - first.sum())
    if n_dropped:
        x, y, z = x[first], y[first], z[first]
    if len(z) < min_points:
        raise KrigingError("only %d distinct locations after de-duplication"
                           % len(z))

    if variogram is None:
        variogram = fit_variogram(x, y, z, model=model, max_lag=radius)
    mdl = MODELS[variogram["model"]]
    p = variogram["params"]

    # THE SYSTEM IS BUILT FROM COVARIANCES, NOT SEMIVARIANCES
    #
    # The two are algebraically equivalent -- C(h) = C(0) - gamma(h) gives the
    # same weights -- but not numerically. With semivariances the diagonal is
    # gamma(0) = 0 while the off-diagonals carry the nugget, so a nugget makes
    # the matrix more nearly constant and conditions it worse. With covariances
    # the diagonal is the total sill and exceeds the off-diagonals by exactly
    # the nugget, which is the diagonal dominance that conditions it. This is
    # where a nugget belongs, and why one had no effect while the system was
    # written the other way.
    c0_sill = float(p[0] + p[1]) if len(p) > 1 else 1.0
    ridge = ridge_fraction * c0_sill

    tree = cKDTree(np.c_[x, y])
    n_cells = grid.height * grid.width
    out = np.full(n_cells, nodata, float)
    var = np.full(n_cells, nodata, float)

    # THE NEIGHBOUR QUERY IS WHAT BOUNDS MEMORY, SO IT IS CHUNKED
    #
    # Querying every cell at once allocates two arrays of cells x max_points.
    # On a full G-LiHT tile -- 46 million cells at 0.5 m, sixteen neighbours --
    # that is about 12 GB of distances and indices before a single cell is
    # solved, which is the reason tile processing was blocked at 1 km in the
    # first place. Chunking by rows bounds it at roughly `chunk_cells` x
    # max_points regardless of how large the grid is, so the block size becomes
    # a choice about the method rather than a limit of the machine.
    #
    # Results are identical: the same cells, the same neighbours, the same
    # solves, in the same order.
    k_use = min(max_points, len(z))
    rows_per_chunk = max(1, int(chunk_cells // max(grid.width, 1)))
    xs = grid.origin_x + (np.arange(grid.width) + 0.5) * grid.cell

    filled = 0
    n_fallback = 0
    for r0 in range(0, grid.height, rows_per_chunk):
        r1 = min(r0 + rows_per_chunk, grid.height)
        ys = grid.origin_y - (np.arange(r0, r1) + 0.5) * grid.dy
        cxs, cys = np.meshgrid(xs, ys)
        chunk = np.c_[cxs.ravel(), cys.ravel()]
        # The query threads across cores. That is right for one tile and wrong
        # inside a batch run, where every worker process would claim every core
        # and they would spend their time contending rather than working. The
        # batch driver sets this to 1 and parallelises over tiles instead.
        dist, idx = tree.query(chunk, k=k_use, distance_upper_bound=radius,
                               workers=_query_workers())
        if dist.ndim == 1:
            dist, idx = dist[:, None], idx[:, None]
        base = r0 * grid.width

        # THE SOLVES ARE BATCHED BY NEIGHBOUR COUNT
        #
        # Every cell builds the same shape of system, so cells that found the
        # same number of neighbours can be solved together: numpy's solve takes
        # a stack of matrices and hands it to LAPACK in one call, where the
        # per-cell loop paid Python's overhead on every one of tens of millions
        # of 17x17 systems. The arithmetic is identical -- same neighbours, same
        # semivariances, same system -- only the dispatch changes.
        #
        # cKDTree returns distances ascending with infinity padding, so a cell
        # with k finite neighbours has them in the first k columns and `[:, :k]`
        # selects exactly what `d[good]` selected before.
        n_good = np.isfinite(dist).sum(1)
        usable = n_good >= min_points
        for kval in np.unique(n_good[usable]):
            sel = np.flatnonzero(usable & (n_good == kval))
            kval = int(kval)
            # Sub-batched because the stacked systems are the memory high water
            # mark: k+1 squared doubles per cell, ~2.3 kB at sixteen neighbours.
            for s0 in range(0, len(sel), solve_batch):
                cells = sel[s0:s0 + solve_batch]
                ii = idx[cells, :kval]
                dd = dist[cells, :kval]
                px, py, pz = x[ii], y[ii], z[ii]
                m = len(cells)

                dxy = np.hypot(px[:, :, None] - px[:, None, :],
                               py[:, :, None] - py[:, None, :])
                # covariances, and the cell-to-point right-hand side with them
                c_0 = c0_sill - mdl(dd, *p)
                A = np.empty((m, kval + 1, kval + 1))
                A[:, :kval, :kval] = c0_sill - mdl(dxy, *p)
                # A small ridge on the diagonal, a fraction of the sill, guards
                # whatever ill-conditioning survives de-duplication. Inert where
                # the system already solves.
                if ridge:
                    dg = np.arange(kval)
                    A[:, dg, dg] += ridge
                A[:, :kval, kval] = 1.0
                A[:, kval, :kval] = 1.0
                A[:, kval, kval] = 0.0
                # Trailing axis kept explicit: numpy 2 reads a two-dimensional
                # right-hand side as one matrix rather than a stack of vectors,
                # so (m, k+1) would be misread as a single (m, k+1) system.
                b = np.empty((m, kval + 1, 1))
                b[:, :kval, 0] = c_0
                b[:, kval, 0] = 1.0

                try:
                    w = np.linalg.solve(A, b)[:, :, 0]
                except np.linalg.LinAlgError:
                    # One singular system in the stack fails the whole call, so
                    # the batch is redone one at a time and the singular ones
                    # are marked rather than losing the batch that contained
                    # them. Rare: degeneracy here is usually caught by the
                    # range test below, not by a raise.
                    w = np.empty((m, kval + 1))
                    singular = np.zeros(m, bool)
                    for j in range(m):
                        try:
                            w[j] = np.linalg.solve(A[j], b[j])[:, 0]
                        except np.linalg.LinAlgError:
                            w[j] = 0.0
                            singular[j] = True
                else:
                    singular = np.zeros(m, bool)

                est = np.einsum("ij,ij->i", w[:, :kval], pz)
                # An ordinary kriging estimate interpolates its neighbours and
                # cannot legitimately fall outside their range. When it does the
                # system was degenerate -- the search radius large against the
                # correlation range, neighbours packed tightly, every pairwise
                # semivariance nearly identical and the matrix short of rank.
                # The fallback is inverse distance over the same neighbours:
                # bounded by construction, and defensible for a cell whose
                # neighbours carry no distinguishable spatial structure.
                degenerate = (singular | ~np.isfinite(est)
                              | (est < pz.min(1)) | (est > pz.max(1)))
                if degenerate.any():
                    ddg = np.maximum(dd[degenerate], 1e-9)
                    ww = 1.0 / ddg ** 2
                    est[degenerate] = ((ww * pz[degenerate]).sum(1)
                                       / ww.sum(1))
                    n_fallback += int(degenerate.sum())

                good_cells = base + cells
                out[good_cells] = est
                # inverse distance carries no covariance model, so a fallback
                # cell has no variance to report. Left as nodata rather than
                # filled with a number that would look like an uncertainty.
                keep = ~degenerate
                if keep.any():
                    # ordinary kriging variance in covariance form:
                    # C(0) - sum(w_i C_i0) - mu
                    var[good_cells[keep]] = (
                        c0_sill
                        - np.einsum("ij,ij->i", w[keep, :kval], c_0[keep])
                        - w[keep, kval])
                filled += m

    dem = out.reshape(grid.height, grid.width)
    variance = var.reshape(grid.height, grid.width)
    frac = n_fallback / filled if filled else 0.0
    if verbose:
        print("  kriged %d of %d cells (%.1f%%), radius %.1f m, "
              "%d nearest, %s variogram%s"
              % (filled, n_cells, 100 * filled / max(n_cells, 1), radius,
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
