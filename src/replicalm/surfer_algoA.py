"""algoA — slope-adaptive constraint relaxation (saved variant, 2026-10-07).

Preserved verbatim as the slope-focused arm of the A/B comparison. The active
engine is `surfer_spline`; this module is kept importable so algoA can be
re-run without reconstructing it. Measured on the l0s441 benchmark: p95 19.3,
p99 25.6, max 48.2 deg against the archive 18.5 / 22.8 / 32.1, anisotropy
0.076, RMSE 0.136, MAE 0.076, gentle-ground high-pass rms held to within
-0.93 percent of the exact-constraint control.

Surfer paradigm: one biharmonic surface, solved globally on the raster.

WHY THIS IS A DIFFERENT PARADIGM
--------------------------------
The shipped engines solve one small system per output cell. This solves one
large system for the whole raster at once. Every cell appears as an unknown in a
single linear system, and the coupling between them is the differential operator
rather than a neighbour list, so there is no radius, no neighbour count, and no
stiffness term to choose. A change to any observation propagates everywhere,
which is what "global" means here and is the property the local engines lack.

THE FORMULATION
---------------
Minimum curvature, after Briggs (1974), which is the classical basis of Surfer's
minimum-curvature gridder. The surface minimises total squared curvature, whose
Euler-Lagrange equation is the biharmonic

    laplacian^2 z = 0

subject to passing through the data. Discretised on the regular output grid with
the thirteen-point stencil obtained by applying the five-point Laplacian twice:

           1
       2  -8   2
    1 -8  20  -8   1
       2  -8   2
           1

Cells not holding observations get the stencil above with right-hand side zero.
Cells holding observations get the stencil too, plus a data term pulling them
toward the mean of the returns inside them. The result is one square sparse
system in width*height unknowns, solved once.

WHY THE DATA TERM IS WEIGHTED RATHER THAN EXACT
-----------------------------------------------
Pinning each observed cell by an identity row makes the surface pass through
every cell mean exactly, which sounds right and is not. Adjacent accepted
returns disagree: across this benchmark the step between neighbouring
constrained cells reaches 0.520 m over 0.4128 m of ground, and an identity row
admits no residual, so the solver honours it and the surface reads 52.6 degrees
where the reference reads 32.1.

Those steps are not noise and not a sampling artefact. Single-return cells are
enriched only 1.19 times among the steepest pairs over the window's own rate,
the median such cell holds two returns, and evaluating the constraint from a
plane at the cell centre rather than the raw mean -- which removes the centroid
bias entirely -- left the maximum step at 0.520 m and the maximum slope at 52.6
degrees. The returns are real and they disagree.

What they disagree by is set by the ground. On flat terrain two returns a cell
apart should agree to the ranging precision, and any disagreement there is
signal worth keeping -- that is where micro-relief lives. On a slope they differ
by the terrain's own rise across the cell as well, and that difference carries
no information the surface needs to honour cell by cell.

So the tolerance is piecewise, with a knee: nothing at all below it, growing
above it.

    tolerance(theta) = 0                                     theta <= knee
                     = beta * (tan theta - tan knee)^gamma    theta >  knee

Below the knee the cell keeps an exact identity row and the surface passes
through it, so gentle ground is untouched and its fine texture is retained in
full. Above it the cell instead carries its stencil plus a data row

    w_i * (z_i - d_i) = 0,   w_i = tightness * cell / tolerance_i

so the surface may depart from the datum by about the tolerance. With the
defaults the tolerance is zero to 5 degrees, 0.07 m at 8, 0.12 m at 10, 0.25 m
at 15 and 0.28 m at 16 -- enough to absorb the 0.2 to 0.35 m steps that produce
the slope tail, while nothing below the knee moves at all.

theta_i is measured from the engine's own first pass: the exact-constraint
surface is solved, smoothed over five metres so the spikes themselves cannot
inflate the gradient, and its slope is read per cell. Two solves, no external
input.

Setting `tightness` to None restores exact identity rows everywhere and the
behaviour this replaced.

BOUNDARY CONDITIONS
-------------------
The thirteen-point stencil needs two cells of margin, which the raster edge does
not have. Those cells take the natural condition of a free thin plate,
laplacian z = 0, built from whichever four-neighbours exist. That is the
standard treatment and introduces no parameter; it does mean the outermost two
rings are governed by a weaker condition than the interior, so interior
statistics are the ones to compare.

WHAT BINNING THIS IMPLIES
-------------------------
A finite-difference formulation has one unknown per cell, so observations are
assigned to cells and averaged where several share one. That is inherent to the
paradigm rather than a choice: sub-cell structure cannot be represented by a
grid of this spacing. Points should still arrive distinct in XY from
`preprocess.dedupe_xy`, so the averaging is over genuinely separate returns.

THE SOLVE DOMAIN IS NOT THE OUTPUT DOMAIN
-----------------------------------------
A grid has no cells outside itself, so a solve run on exactly the requested
extent simply discards every point beyond it and then applies a free-plate
condition along a boundary that has real data sitting just across it. On the
l0s441 benchmark that threw away 12,562 of 17,860 returns and left the output
edge governed by the boundary condition instead of by the observations.

So the system is assembled over whatever extent the supplied points occupy, and
the requested window is cut out of the solution afterwards. The caller's choice
of how much margin to pass is what sets the buffer; there is no separate
parameter, and the surface inside the requested window is the one the full point
set implies.
"""
import numpy as np

NODATA = -9999.0

DEFAULT_KNEE_DEG = 5.0       # below this the constraint stays exact
DEFAULT_BETA = 1.4           # metres of tolerance per unit of tangent above it
DEFAULT_GAMMA = 1.0          # curvature of the tolerance ramp
DEFAULT_TIGHTNESS = 10.0     # data weight, in units of cell / tolerance
SLOPE_SMOOTH_M = 5.0         # scale at which the local gradient is read


class SurferError(RuntimeError):
    pass


def surfer_grid(x, y, z, grid, knee_deg=DEFAULT_KNEE_DEG, beta=DEFAULT_BETA,
                gamma=DEFAULT_GAMMA, tightness=DEFAULT_TIGHTNESS,
                pool_sigma_m=0.0, nodata=NODATA, verbose=True,
                progress=None):
    """Global biharmonic minimum-curvature surface on a declared grid.

    Returns (dem, info), matching the other engines' call signature.

    The data term is weighted by local terrain gradient; see the module
    docstring. `tightness=None` gives exact identity constraints instead.
    """
    import scipy.sparse as sp
    import scipy.sparse.linalg as spla

    x = np.asarray(x, float)
    y = np.asarray(y, float)
    z = np.asarray(z, float)
    if int(grid.height) * int(grid.width) == 0:
        raise SurferError("empty grid")
    if len(z) == 0:
        raise SurferError("no observations supplied")

    # ---- extend the solve domain to cover every point supplied ----------
    cell = grid.cell
    pad_w = int(np.ceil(max(0.0, (grid.origin_x - x.min()) / cell)))
    pad_e = int(np.ceil(max(0.0, (x.max() - (grid.origin_x
                                             + grid.width * cell)) / cell)))
    pad_n = int(np.ceil(max(0.0, (y.max() - grid.origin_y) / cell)))
    pad_s = int(np.ceil(max(0.0, ((grid.origin_y - grid.height * cell)
                                  - y.min()) / cell)))
    ox = grid.origin_x - pad_w * cell
    oy = grid.origin_y + pad_n * cell
    w = int(grid.width) + pad_w + pad_e
    h = int(grid.height) + pad_n + pad_s
    n = h * w
    if verbose and (pad_w or pad_e or pad_n or pad_s):
        print("  surfer: solving on %d x %d to cover the supplied points, "
              "cutting %d x %d out of it (buffer W%d E%d N%d S%d cells)"
              % (h, w, grid.height, grid.width, pad_w, pad_e, pad_n, pad_s))

    # ---- observations onto the finite-difference grid -------------------
    col = np.floor((x - ox) / cell).astype(np.int64)
    row = np.floor((oy - y) / cell).astype(np.int64)
    inside = (col >= 0) & (col < w) & (row >= 0) & (row < h)
    if not inside.any():
        raise SurferError("no observation falls inside the grid")
    flat = row[inside] * w + col[inside]
    tot = np.bincount(flat, weights=z[inside], minlength=n)
    cnt = np.bincount(flat, minlength=n)
    has = cnt > 0
    fixed = np.zeros(n)
    fixed[has] = tot[has] / cnt[has]

    # Distance-weighted pooling of the returns behind each constraint.
    #
    # At this density a 0.5 m cell holds about half a return, so most
    # constrained cells carry a single one and inherit that return's ranging
    # error in full, independently of their neighbours. That is white noise at
    # the cell scale. Pooling over a Gaussian support averages it down while
    # leaving structure coarser than the support alone, because the noise is
    # uncorrelated between cells and real form is not.
    #
    # Blurring the summed heights and the counts separately and dividing is
    # exactly a Gaussian-weighted local mean of the returns, at grid cost.
    if pool_sigma_m and pool_sigma_m > 0:
        from scipy.ndimage import gaussian_filter as _gf
        s = pool_sigma_m / cell
        sg = _gf(tot.reshape(h, w), s, mode="constant").ravel()
        ng = _gf(cnt.astype(float).reshape(h, w), s, mode="constant").ravel()
        ok_pool = has & (ng > 1e-9)
        fixed[ok_pool] = sg[ok_pool] / ng[ok_pool]
        if verbose:
            print("  surfer: constraints pooled over a %.2f m Gaussian "
                  "support (%d cells)" % (pool_sigma_m, int(ok_pool.sum())))
    n_fixed = int(has.sum())
    if verbose:
        print("  surfer: %d of %d cells carry observations (%.2f%%), %d points "
              "dropped outside the grid"
              % (n_fixed, n, 100.0 * n_fixed / n, int((~inside).sum())))
    if progress is not None:
        progress("interpolate", "%d of %d cells constrained" % (n_fixed, n),
                 0.40)

    rr = np.arange(n).reshape(h, w)
    hasg = has.reshape(h, w)

    # cells far enough from the edge for the thirteen-point stencil
    deep = np.zeros((h, w), bool)
    if h > 4 and w > 4:
        deep[2:h - 2, 2:w - 2] = True
    edge = ~deep
    STENCIL = [(0, 0, 20.0),
               (-1, 0, -8.0), (1, 0, -8.0), (0, -1, -8.0), (0, 1, -8.0),
               (-1, -1, 2.0), (-1, 1, 2.0), (1, -1, 2.0), (1, 1, 2.0),
               (-2, 0, 1.0), (2, 0, 1.0), (0, -2, 1.0), (0, 2, 1.0)]

    def assemble(weights, hard):
        """Operator everywhere except pinned cells; soft cells add a data term.

        `hard` marks observed cells kept on an exact identity row. `weights`
        supplies w_i for the observed cells that are not hard; pass None for
        the all-exact case.
        """
        hardg = hard.reshape(h, w)
        softg = hasg & ~hardg
        rows_i, cols_i, vals = [], [], []
        # a pinned cell carries no operator at all; everything else does
        op_deep, op_edge = deep & ~hardg, edge & ~hardg
        fi = rr[hardg]
        if len(fi):
            rows_i.append(fi); cols_i.append(fi)
            vals.append(np.ones(len(fi)))
        si = rr[softg]
        if len(si):
            rows_i.append(si); cols_i.append(si)
            vals.append(weights.reshape(h, w)[softg])

        di = rr[op_deep]
        if len(di):
            dr, dc = np.nonzero(op_deep)
            for dy, dx, v in STENCIL:
                rows_i.append(di)
                cols_i.append(rr[dr + dy, dc + dx])
                vals.append(np.full(len(di), v))

        ei = rr[op_edge]
        if len(ei):
            er, ec = np.nonzero(op_edge)
            deg = np.zeros(len(ei))
            for dy, dx in ((-1, 0), (1, 0), (0, -1), (0, 1)):
                ny, nx = er + dy, ec + dx
                ok = (ny >= 0) & (ny < h) & (nx >= 0) & (nx < w)
                if ok.any():
                    rows_i.append(ei[ok])
                    cols_i.append(rr[ny[ok], nx[ok]])
                    vals.append(np.full(int(ok.sum()), -1.0))
                    deg += ok
            rows_i.append(ei); cols_i.append(ei); vals.append(deg)

        A = sp.coo_matrix((np.concatenate(vals),
                           (np.concatenate(rows_i), np.concatenate(cols_i))),
                          shape=(n, n)).tocsr()
        b = np.zeros(n)
        b[hard] = fixed[hard]
        soft = has & ~hard
        if soft.any():
            b[soft] = fixed[soft] * weights[soft]
        return A, b

    def solve(A, b):
        with np.errstate(all="ignore"):
            sol = spla.spsolve(A, b)
        if np.all(np.isfinite(sol)):
            return sol, "spsolve"
        # a free region with no constraint anywhere in it leaves the operator
        # rank deficient; least squares resolves it rather than returning nan
        sol = spla.lsqr(A, b, atol=1e-10, btol=1e-10, iter_lim=20000)[0]
        if not np.all(np.isfinite(sol)):
            raise SurferError("the global system did not resolve")
        return sol, "lsqr"

    if progress is not None:
        progress("interpolate", "solving %d unknowns" % n, 0.50)

    A, b = assemble(None, has)
    sol, method = solve(A, b)
    passes = 1
    weights = None
    tol_stats = None
    n_soft = 0

    if tightness is not None:
        # local gradient from the first pass, smoothed so that the spikes the
        # weighting exists to remove cannot set their own tolerance
        from scipy.ndimage import gaussian_filter
        sm = gaussian_filter(sol.reshape(h, w), SLOPE_SMOOTH_M / cell / 2.0)
        gy, gx = np.gradient(sm, cell)
        tan_theta = np.hypot(gx, gy)
        over = np.maximum(tan_theta - np.tan(np.radians(knee_deg)), 0.0)
        tol = (beta * over ** gamma).ravel()
        hard = has & (tol <= 0.0)
        soft = has & ~hard
        n_soft = int(soft.sum())
        weights = np.zeros(n)
        weights[soft] = tightness * cell / tol[soft]
        ts = tol[soft]
        tol_stats = ((float(np.median(ts)), float(np.percentile(ts, 99)),
                      float(ts.max())) if n_soft else None)
        if verbose:
            print("  surfer: %d observed cells pinned exactly (theta <= %.0f "
                  "deg), %d relaxed" % (int(hard.sum()), knee_deg, n_soft))
            if n_soft:
                print("  surfer: tolerance %.3f to %.3f m (median %.3f); data "
                      "weight %.0f to %.0f"
                      % (ts.min(), ts.max(), np.median(ts),
                         weights[soft].min(), weights[soft].max()))
        A, b = assemble(weights, hard)
        sol, method = solve(A, b)
        passes = 2
    else:
        hard = has

    full = sol.reshape(h, w)
    dem = full[pad_n:pad_n + int(grid.height), pad_w:pad_w + int(grid.width)]

    misfit = None
    if tightness is not None and n_soft:
        d = np.abs(full.ravel()[has & ~hard] - fixed[has & ~hard])
        misfit = (float(np.median(d)), float(np.percentile(d, 99)),
                  float(d.max()))

    info = {"engine": "surfer_spline",
            "model": ("global_biharmonic_minimum_curvature" if tightness is None
                      else "global_biharmonic_slope_adaptive_constraints"),
            "constraint": ("exact_identity" if tightness is None
                           else "slope_adaptive_piecewise"),
            "knee_deg": None if tightness is None else float(knee_deg),
            "beta": None if tightness is None else float(beta),
            "gamma": None if tightness is None else float(gamma),
            "tightness": None if tightness is None else float(tightness),
            "cells_pinned": int(hard.sum()),
            "cells_relaxed": n_soft,
            "tolerance_m": tol_stats,
            "constraint_misfit_m": misfit,
            "passes": passes,
            "points": int(len(z)),
            "points_used": int(inside.sum()),
            "buffer_cells": [int(pad_w), int(pad_e), int(pad_n), int(pad_s)],
            "solve_shape": [int(h), int(w)],
            "unknowns": int(n),
            "nonzeros": int(A.nnz),
            "constrained_cells": n_fixed,
            "constrained_fraction": float(n_fixed) / n,
            "interior_biharmonic_cells": int(deep.sum()),
            "edge_laplacian_cells": int(edge.sum()),
            "solver": method,
            "fallback_fraction": 0.0,
            "singular_cells": 0,
            "params": None}
    if verbose:
        print("  surfer: %d unknowns, %d nonzeros, %s, %.2f to %.2f m"
              % (n, A.nnz, method, np.nanmin(dem), np.nanmax(dem)))
    return dem, info
