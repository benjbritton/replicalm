"""TerraScan paradigm: a triangulated surface, sampled directly.

WHY THIS IS A DIFFERENT PARADIGM
--------------------------------
Both shipped engines -- ordinary kriging and the regularised spline with tension
-- estimate each output cell from a neighbourhood gathered for that cell. Two
adjacent cells can hold different neighbour sets, so the surface is only as
continuous as the overlap between them, and the parameters that govern the
gathering (radius, count, stiffness) have to be chosen.

This engine has no such parameters. The ground returns define one Delaunay
triangulation, that triangulation defines one surface, and every output cell is
a sample of it. Continuity is a property of the construction rather than
something a neighbourhood size has to deliver: the Clough-Tocher element is a
cubic on each triangle, C1 across every shared edge, so the surface and its
first derivatives are continuous everywhere inside the hull.

THE FORMULATION
---------------
Clough and Tocher (1965). Each triangle is split at its centroid into three
sub-triangles, and a cubic polynomial is fitted on each so that value and
gradient match across all interior and shared edges. Vertex gradients are
estimated by minimising the curvature of the result, which is global over the
triangulation rather than per-cell.

WHAT IT CANNOT DO
-----------------
A triangulation interpolates; it does not extrapolate. Cells outside the convex
hull of the ground returns have no triangle and are returned as nodata. Inside a
canopy void the surface is whatever the triangle spanning that void says, which
is a flat-to-cubic bridge between its three corners -- the engine has no
mechanism to prefer any other shape there, and that is the honest behaviour of
the paradigm rather than a defect to patch.

Points must be distinct in XY. `preprocess.dedupe_xy` does that upstream for
every engine; coincident points make the triangulation degenerate.
"""
import numpy as np

NODATA = -9999.0


class TinError(RuntimeError):
    pass


def tin_grid(x, y, z, grid, nodata=NODATA, verbose=True, progress=None):
    """Delaunay plus Clough-Tocher, sampled on a declared grid.

    Returns (dem, info), matching `kriging.krige_grid` and `spline.spline_grid`
    so the three are interchangeable at the call site.

    There are no tuning parameters. Nothing here is chosen; the triangulation
    is determined by the points and the surface by the triangulation.
    """
    from scipy.spatial import Delaunay
    from scipy.interpolate import CloughTocher2DInterpolator

    x = np.asarray(x, float)
    y = np.asarray(y, float)
    z = np.asarray(z, float)
    if len(z) < 3:
        raise TinError("a triangulation needs at least three points; got %d"
                       % len(z))

    pts = np.column_stack([x, y])
    try:
        tri = Delaunay(pts)
    except Exception as exc:                        # QhullError and friends
        raise TinError("triangulation failed: %s" % exc)

    if progress is not None:
        progress("interpolate", "triangulated %d points into %d triangles"
                 % (len(z), len(tri.simplices)), 0.40)

    interp = CloughTocher2DInterpolator(tri, z)

    cols = np.arange(grid.width)
    rows = np.arange(grid.height)
    cx = grid.origin_x + (cols + 0.5) * grid.cell
    cy = grid.origin_y - (rows + 0.5) * grid.cell
    QX, QY = np.meshgrid(cx, cy)

    dem = interp(QX, QY)
    outside = ~np.isfinite(dem)

    if progress is not None:
        progress("interpolate", "sampled %d cells" % dem.size, 0.70)

    info = {"engine": "terrascan_tin",
            "model": "delaunay_clough_tocher",
            "points": int(len(z)),
            "triangles": int(len(tri.simplices)),
            "outside_hull": int(outside.sum()),
            "outside_fraction": float(outside.mean()),
            # the other engines report these and callers read them; a
            # triangulation has no fallback estimator and no local solves
            "fallback_fraction": 0.0,
            "singular_cells": 0,
            "params": None}
    if verbose:
        print("  tin: %d points, %d triangles, %d cells outside the hull "
              "(%.2f%%)" % (len(z), len(tri.simplices), outside.sum(),
                            100.0 * outside.mean()))
    return dem, info
