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

Cells holding observations are fixed, by an identity row, to the mean of the
observations falling inside them. Cells not holding observations get the stencil
above, with right-hand side zero. The result is one square sparse system in
width*height unknowns, solved once.

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


class SurferError(RuntimeError):
    pass


def surfer_grid(x, y, z, grid, nodata=NODATA, verbose=True, progress=None):
    """Global biharmonic minimum-curvature surface on a declared grid.

    Returns (dem, info), matching the other engines' call signature.

    There are no tuning parameters. The surface is determined by the grid, the
    operator and the data.
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
    n_fixed = int(has.sum())
    if verbose:
        print("  surfer: %d of %d cells carry observations (%.2f%%), %d points "
              "dropped outside the grid"
              % (n_fixed, n, 100.0 * n_fixed / n, int((~inside).sum())))
    if progress is not None:
        progress("interpolate", "%d of %d cells constrained" % (n_fixed, n),
                 0.40)

    rr = np.arange(n).reshape(h, w)
    free = (~has).reshape(h, w)

    # interior free cells take the full thirteen-point biharmonic
    deep = np.zeros((h, w), bool)
    if h > 4 and w > 4:
        deep[2:h - 2, 2:w - 2] = True
    deep &= free
    edge = free & ~deep                      # free, but too close to an edge

    rows_i, cols_i, vals = [], [], []

    # ---- fixed cells: identity -------------------------------------------
    fi = rr.ravel()[has]
    rows_i.append(fi); cols_i.append(fi); vals.append(np.ones(len(fi)))

    # ---- interior: biharmonic ---------------------------------------------
    di = rr[deep]
    if len(di):
        STENCIL = [(0, 0, 20.0),
                   (-1, 0, -8.0), (1, 0, -8.0), (0, -1, -8.0), (0, 1, -8.0),
                   (-1, -1, 2.0), (-1, 1, 2.0), (1, -1, 2.0), (1, 1, 2.0),
                   (-2, 0, 1.0), (2, 0, 1.0), (0, -2, 1.0), (0, 2, 1.0)]
        dr, dc = np.nonzero(deep)
        for dy, dx, v in STENCIL:
            rows_i.append(di)
            cols_i.append(rr[dr + dy, dc + dx])
            vals.append(np.full(len(di), v))

    # ---- near the edge: laplacian z = 0 over the neighbours that exist ----
    ei = rr[edge]
    if len(ei):
        er, ec = np.nonzero(edge)
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
    b[has] = fixed[has]

    if progress is not None:
        progress("interpolate", "solving %d unknowns, %d nonzeros"
                 % (n, A.nnz), 0.50)

    method = "spsolve"
    with np.errstate(all="ignore"):
        sol = spla.spsolve(A, b)
    if not np.all(np.isfinite(sol)):
        # a free region with no constraint anywhere in it leaves the operator
        # rank deficient; least squares resolves it rather than returning nan
        method = "lsqr"
        sol = spla.lsqr(A, b, atol=1e-10, btol=1e-10, iter_lim=20000)[0]
        if not np.all(np.isfinite(sol)):
            raise SurferError("the global system did not resolve")

    full = sol.reshape(h, w)
    dem = full[pad_n:pad_n + int(grid.height), pad_w:pad_w + int(grid.width)]

    info = {"engine": "surfer_spline",
            "model": "global_biharmonic_minimum_curvature",
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
