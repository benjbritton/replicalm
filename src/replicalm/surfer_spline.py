"""Surfer paradigm: one biharmonic surface, solved globally on the raster.

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
Cells holding observations get the stencil too, plus a weighted data row pulling
them toward the mean of the returns inside them. The result is one square sparse
system in width*height unknowns, solved once.

WHY THE DATA ROW IS WEIGHTED, AND WHY SLOPE ALONE CANNOT SET THE WEIGHT
-----------------------------------------------------------------------
Pinning each observed cell exactly makes the surface pass through every cell
mean, which sounds right and is not. Adjacent accepted returns disagree: on this
benchmark the step between neighbouring constrained cells reaches 0.520 m over
0.4128 m of ground, and an identity row admits no residual, so the solver
honours it and the surface reads 52.6 degrees where the reference reads 32.1.
Those steps are not sparse acceptances and not centroid bias -- both were tested
and ruled out -- so the returns are real and they disagree.

The obvious response is to relax the constraint wherever the ground is steep.
That is wrong for this landscape. Terrace faces, platform walls, escarpments and
other engineered earthworks occupy exactly the steep part of the surface, and a
purely slope-driven rule would relax the constraint hardest on precisely the
features the survey exists to record.

What separates them is curvature, not slope. A hillside is steep and smooth: a
large gradient with a near-zero second derivative. A terrace edge or platform
wall is a break in slope -- the gradient changes over a metre or two, which is a
large second derivative. So the weight is set by both.

    protection = logistic(slope)  OR  gate(curvature)

    logistic(theta) = 1 / (1 + exp(k (theta - theta0)))

centred at theta0 = 5.5 degrees and spanning roughly 3 to 8, so there is no
threshold anywhere and no seam in the solution. It falls from about 0.95 at 3
degrees to about 0.05 at 8.

    gate(kappa) = 1 / (1 + exp(-kc (|kappa| / kappa0 - 1)))

with kappa the Laplacian of the surface read over a 1.5 m window and kappa0 the
curvature of a break of `break_relief_m` across that window. The gate rises
toward one wherever the ground bends sharply, and the two combine as a smooth
union, protection = L + G - L*G, so a cell is held firmly if it is either gentle
or sharply curved, and released only where high slope meets low curvature --
broad monotonic hillsides, which is where the step artefacts live.

The weight interpolates geometrically between the two regimes,

    w_i = w_relaxed * (w_max / w_relaxed) ** protection_i

so it is smooth in the exponent across four orders of magnitude. w_max is capped
at 1e4 rather than taken to infinity: against a stencil diagonal of 20 that
tracks the data to within a fraction of a millimetre while keeping the matrix
far from ill-conditioned.

Curvature is read from the engine's own first pass, smoothed over 1.5 m before
the Laplacian is taken. That smoothing is what separates the two populations: an
isolated single-cell spike, 0.41 m across, is heavily attenuated by it, while a
terrace edge running for metres survives. Two solves, no external input.

Setting `tightness` to None restores exact identity rows everywhere.

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

DEFAULT_THETA0 = 5.5          # degrees; centre of the slope logistic
DEFAULT_K_SLOPE = 1.2         # per degree; spans roughly 3 to 8 degrees
DEFAULT_W_MAX = 1.0e4         # capped, to keep the matrix well conditioned
DEFAULT_W_RELAXED = 4.0       # weight on a broad, smooth, steep hillside
DEFAULT_CURV_WINDOW_M = 1.5   # window the Laplacian is read over
DEFAULT_BREAK_RELIEF_M = 0.10 # relief across that window counting as a break
DEFAULT_K_CURV = 6.0          # sharpness of the curvature gate


class SurferError(RuntimeError):
    pass


def surfer_grid(x, y, z, grid, theta0=DEFAULT_THETA0,
                k_slope=DEFAULT_K_SLOPE, w_max=DEFAULT_W_MAX,
                w_relaxed=DEFAULT_W_RELAXED,
                curv_window_m=DEFAULT_CURV_WINDOW_M,
                break_relief_m=DEFAULT_BREAK_RELIEF_M,
                k_curv=DEFAULT_K_CURV, tightness=1.0, pool_sigma_m=0.0, nodata=NODATA,
                verbose=True, progress=None):
    """Global biharmonic surface with curvature-aware constraint relaxation.

    Returns (dem, info), matching the other engines' call signature.

    See the module docstring for how the observation weight is set.
    `tightness=None` gives exact identity constraints instead.
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

    deep = np.zeros((h, w), bool)
    if h > 4 and w > 4:
        deep[2:h - 2, 2:w - 2] = True
    edge = ~deep
    STENCIL = [(0, 0, 20.0),
               (-1, 0, -8.0), (1, 0, -8.0), (0, -1, -8.0), (0, 1, -8.0),
               (-1, -1, 2.0), (-1, 1, 2.0), (1, -1, 2.0), (1, 1, 2.0),
               (-2, 0, 1.0), (2, 0, 1.0), (0, -2, 1.0), (0, 2, 1.0)]

    def assemble(weights):
        """Operator on every free cell; observed cells add their data row.

        `weights` None means observed cells are pinned by identity rows and
        carry no operator, which is the exact-constraint behaviour.
        """
        rows_i, cols_i, vals = [], [], []
        if weights is None:
            op_deep, op_edge = deep & ~hasg, edge & ~hasg
            fi = rr[hasg]
            rows_i.append(fi); cols_i.append(fi)
            vals.append(np.ones(len(fi)))
        else:
            op_deep, op_edge = deep, edge
            fi = rr[hasg]
            rows_i.append(fi); cols_i.append(fi)
            vals.append(weights.reshape(h, w)[hasg])

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
        b[has] = fixed[has] * (1.0 if weights is None else weights[has])
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

    A, b = assemble(None)
    sol, method = solve(A, b)
    passes = 1
    diag = None

    if tightness is not None:
        from scipy.ndimage import gaussian_filter, laplace
        first = sol.reshape(h, w)

        # slope, read over 5 m so the spikes cannot inflate it
        sm5 = gaussian_filter(first, 5.0 / cell / 2.0)
        gy, gx = np.gradient(sm5, cell)
        theta = np.degrees(np.arctan(np.hypot(gx, gy)))

        # curvature, read over 1.5 m: wide enough to damp a single-cell spike,
        # narrow enough to keep a terrace edge
        sm_c = gaussian_filter(first, curv_window_m / cell / 2.0)
        kappa = np.abs(laplace(sm_c) / (cell * cell))
        kappa0 = 2.0 * break_relief_m / (curv_window_m ** 2)

        L = 1.0 / (1.0 + np.exp(np.clip(k_slope * (theta - theta0), -60, 60)))
        G = 1.0 / (1.0 + np.exp(np.clip(-k_curv * (kappa / kappa0 - 1.0),
                                        -60, 60)))
        prot = L + G - L * G                      # smooth union
        weights = (w_relaxed * (w_max / w_relaxed) ** prot).ravel()

        hs_ = has
        diag = dict(theta=theta, kappa=kappa, kappa0=float(kappa0), prot=prot,
                    weights=weights)
        if verbose:
            p = prot.ravel()[hs_]
            print("  surfer: slope %.1f to %.1f deg; curvature threshold %.4f "
                  "1/m, observed-cell |kappa| median %.4f"
                  % (theta.min(), theta.max(), kappa0,
                     float(np.median(kappa.ravel()[hs_]))))
            print("  surfer: protection %.2f median; %d observed cells held "
                  "(>0.5), %d relaxed; weight %.1f to %.0f"
                  % (np.median(p), int((p > 0.5).sum()), int((p <= 0.5).sum()),
                     weights[hs_].min(), weights[hs_].max()))
            held_by_curv = int(((L.ravel()[hs_] <= 0.5)
                                & (G.ravel()[hs_] > 0.5)).sum())
            print("  surfer: %d steep cells held by the curvature gate alone"
                  % held_by_curv)
        A, b = assemble(weights)
        sol, method = solve(A, b)
        passes = 2

    full = sol.reshape(h, w)
    dem = full[pad_n:pad_n + int(grid.height), pad_w:pad_w + int(grid.width)]

    misfit = None
    held = relaxed = None
    if diag is not None:
        d = np.abs(full.ravel()[has] - fixed[has])
        misfit = (float(np.median(d)), float(np.percentile(d, 99)),
                  float(d.max()))
        p = diag["prot"].ravel()[has]
        held, relaxed = int((p > 0.5).sum()), int((p <= 0.5).sum())

    info = {"engine": "surfer_spline",
            "model": ("global_biharmonic_minimum_curvature" if diag is None
                      else "global_biharmonic_curvature_aware_logistic"),
            "constraint": ("exact_identity" if diag is None
                           else "logistic_slope_with_curvature_gate"),
            "theta0_deg": None if diag is None else float(theta0),
            "k_slope": None if diag is None else float(k_slope),
            "w_max": None if diag is None else float(w_max),
            "w_relaxed": None if diag is None else float(w_relaxed),
            "curv_window_m": None if diag is None else float(curv_window_m),
            "break_relief_m": None if diag is None else float(break_relief_m),
            "kappa_threshold_per_m": None if diag is None else diag["kappa0"],
            "cells_held": held,
            "cells_relaxed": relaxed,
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
