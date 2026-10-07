"""Minimum curvature with tension, as an alternative surface engine.

WHY A SECOND ENGINE
-------------------
Ordinary kriging estimates each cell as a weighted average of its neighbours
under an assumption of second-order stationarity. Where ground returns arrive in
clumps separated by voids -- which is what a canopy produces -- the weights
concentrate on whichever clump is nearest and the surface relaxes toward the
local mean between them. On sloped ground that relaxation reads as curvilinear
grooves following the clump geometry, and no setting of the neighbourhood,
nugget or output cell removes it, because the mechanism is in the estimator
rather than in its parameters.

This functional has no local mean to relax toward. The surface is the one of
minimum tensioned curvature passing through the observations, and between them
it is governed by a smoothness condition rather than by an average.

Measured on Centro GLAS l0s441, interior cells only, identical accepted ground
and identical grid, against a TerraScan/Surfer reference of 21.7 degrees at the
99th percentile of local slope:

    tension   0.01   0.03   0.05   0.10   0.25   1.00
    p99       20.9   20.1   21.0   22.8   24.0   29.3
    ordinary kriging on the same points and cells: 24.1

Taut surfaces are much worse than kriging and freely bending ones match the
reference, which is the opposite of what a smoothing argument would predict and
is the signature of the void-bridging problem: the surface needs curvature
continuity across a gap, not a membrane stretched over it.

0.03 is the operating point. Ringing was looked for directly rather than
inferred from the slope statistic, on a 60 m window over the hillocks, by three
measures a halo shows in and a smooth surface does not: the elevation range, the
departure of each cell from the median of its 5 m neighbourhood, and the tail of
|curvature|. Against the reference on the same window:

    surface                p99    elevation range    excursion p99.9 / max
    TerraScan / Surfer    22.8     97.23 - 101.78        0.382 / 0.600
    tension 0.03          20.9     97.28 - 101.78        0.346 / 0.443
    tension 0.02          21.1     97.28 - 101.77        0.354 / 0.457

Neither bends outside the reference's elevation range, and both depart less from
their own neighbourhoods than the reference does -- the opposite of overshoot.
What the lower tension does is stop paying: 0.02 scores slightly worse than 0.03
on every column, so the response has converged and there is nothing below it to
reach for. The earlier reading that this engine rings below 0.03 came from the
unregularised form, before `smoothing` was added; with the Tikhonov term in
place the regularised surface does not.

THE FORMULATION
---------------
Regularised spline with tension, after Mitasova and Mitas (1993). In two
dimensions the Green's function of the tensioned functional is

    R(r) = -[ E1(rho) + ln(rho) + Ce ],      rho = (phi r / 2)^2

with E1 the exponential integral and Ce the Euler-Mascheroni constant. The
surface is

    z(s) = a + sum_j lambda_j R(|s - s_j|)      subject to  sum_j lambda_j = 0

solved over the k nearest observations to each cell. Locality keeps the system
small and keeps the surface responsive to real terrain; a global solve on
millions of returns is a dense system of that order and is not available.

CONDITIONING IS NOT OPTIONAL HERE
---------------------------------
R(0) is zero and R(r) goes as -rho for small r, so two observations a few
centimetres apart give near-identical rows and the coefficients diverge --
elevation ranges in the thousands of metres, from terrain spanning nine. Exact
coincidences are removed upstream by `preprocess.dedupe_xy`, which on one clip
took out 13.56% of the accepted returns, but clustering at centimetre scale
survives that and still has to be handled. `smoothing` adds a Tikhonov term
scaled to the kernel block, relaxing exact interpolation by a stated fraction.

That term is the same statement the kriging nugget makes: returns closer
together than the ranging precision are not independent observations of
different places.
"""
import numpy as np

EULER_MASCHERONI = 0.5772156649015329

DEFAULT_TENSION = 0.03        # inverse metres; see the tables above
DEFAULT_SMOOTHING = 1e-2      # fraction of the kernel block's magnitude
DEFAULT_NEIGHBOURS = 48


class SplineError(RuntimeError):
    pass


def kernel(r2, tension):
    """R(r) for the tensioned functional, taking SQUARED distance."""
    rho = (tension * tension / 4.0) * np.asarray(r2, dtype=float)
    # R(0) is zero in the limit: E1(rho) -> -Ce - ln(rho) as rho -> 0, so the
    # bracket vanishes. Clamping keeps the logarithm finite at the diagonal.
    rho = np.maximum(rho, 1e-12)
    from scipy.special import exp1
    return -(exp1(rho) + np.log(rho) + EULER_MASCHERONI)


def spline_grid(x, y, z, grid, tension=DEFAULT_TENSION,
                neighbours=DEFAULT_NEIGHBOURS, smoothing=DEFAULT_SMOOTHING,
                min_points=1, nodata=-9999.0, chunk_cells=1000000,
                solve_batch=8000, verbose=True, progress=None):
    """Minimum-curvature-with-tension surface on a declared grid.

    Signature mirrors `kriging.krige_grid` so the two are interchangeable at the
    call site. Returns (dem, info) the same way.

    Points must already be distinct in XY; pass them through
    `preprocess.dedupe_xy` first. The caller does that once, upstream, so that
    every engine sees the same conditioned set.
    """
    from scipy.spatial import cKDTree

    x = np.asarray(x, float); y = np.asarray(y, float); z = np.asarray(z, float)
    if len(z) < min_points:
        raise SplineError("only %d points supplied" % len(z))
    k = int(min(neighbours, len(z)))
    if k < 3:
        raise SplineError("a surface needs at least three observations; got %d"
                          % len(z))

    tree = cKDTree(np.column_stack([x, y]))
    dem = np.full((grid.height, grid.width), np.nan)

    # Rows per chunk, so memory is bounded by the neighbour block rather than by
    # the raster, exactly as the kriging path bounds it.
    per_row = max(grid.width, 1)
    rows_per_chunk = max(1, int(chunk_cells // per_row))
    singular = 0
    solved = 0

    for r0 in range(0, grid.height, rows_per_chunk):
        r1 = min(r0 + rows_per_chunk, grid.height)
        rows = np.arange(r0, r1)
        cols = np.arange(grid.width)
        cx = grid.origin_x + (cols + 0.5) * grid.cell
        cy = grid.origin_y - (rows + 0.5) * grid.cell
        QX, QY = np.meshgrid(cx, cy)
        q = np.column_stack([QX.ravel(), QY.ravel()])

        est = np.full(len(q), np.nan)
        for s0 in range(0, len(q), solve_batch):
            s1 = min(s0 + solve_batch, len(q))
            qs = q[s0:s1]
            _d, idx = tree.query(qs, k=k, workers=-1)
            idx = np.atleast_2d(idx)
            sx, sy, sz = x[idx], y[idx], z[idx]
            m = len(qs)

            dx = sx[:, :, None] - sx[:, None, :]
            dy = sy[:, :, None] - sy[:, None, :]
            A = np.empty((m, k + 1, k + 1))
            A[:, :k, :k] = kernel(dx * dx + dy * dy, tension)
            A[:, :k, k] = 1.0
            A[:, k, :k] = 1.0
            A[:, k, k] = 0.0
            if smoothing:
                scale = np.mean(np.abs(A[:, :k, :k]), axis=(1, 2))
                d = np.arange(k)
                A[:, d, d] += (smoothing * scale)[:, None]

            # The right-hand side carries a trailing axis so this is a stack of
            # systems under both numpy 1 and numpy 2. Numpy 2 reads a 2-D `b`
            # as one matrix rather than as a stack of vectors, so the bare form
            # works on the pinned runtime and raises on a newer one.
            b = np.zeros((m, k + 1, 1))
            b[:, :k, 0] = sz
            try:
                sol = np.linalg.solve(A, b)[:, :, 0]
            except np.linalg.LinAlgError:
                # one bad neighbourhood should not lose a chunk
                sol = np.full((m, k + 1), np.nan)
                for j in range(m):
                    try:
                        sol[j] = np.linalg.solve(A[j], b[j, :, 0])
                    except np.linalg.LinAlgError:
                        singular += 1
            lam, a0 = sol[:, :k], sol[:, k]
            ddx = qs[:, 0][:, None] - sx
            ddy = qs[:, 1][:, None] - sy
            est[s0:s1] = a0 + np.einsum(
                "ij,ij->i", lam, kernel(ddx * ddx + ddy * ddy, tension))
            solved += m

        dem[r0:r1, :] = est.reshape(r1 - r0, grid.width)
        if progress is not None:
            progress("interpolate", "spline rows %d of %d" % (r1, grid.height),
                     0.35 + 0.35 * r1 / grid.height)

    info = {"engine": "spline", "tension": float(tension),
            "neighbours": int(k), "smoothing": float(smoothing),
            "singular_cells": int(singular),
            # the kriging path reports this and callers read it; a spline has no
            # inverse-distance fallback, so it is zero rather than absent
            "fallback_fraction": 0.0,
            "params": None, "model": "tensioned_spline"}
    if verbose:
        print("  spline: tension %.3f, %d neighbours, %d cells, %d singular"
              % (tension, k, solved, singular))
    return dem, info
