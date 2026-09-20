"""Universal kriging: model the trend instead of hiding from it.

WHY THIS EXISTS
---------------
Ordinary kriging assumes the field is stationary -- that its mean does not change
across the neighbourhood being solved. The three calibration tiles say this
terrain is not. Fitting variograms at windows of 5, 10, 20, 40 and 80 m gave
ranges that climbed with every window and sills that never plateaued:

    l8s431   range  1.76 -> 59.99 m    sill 0.0381 -> 0.0712
    l0s444   range  1.46 -> 55.35 m    sill 0.0219 -> 0.0361
    l0s395   range  6.94 -> 92.52 m    sill 0.0142 -> 0.3880

Semivariance still growing at eighty metres is the signature of a regional
gradient: points further apart are further down the slope, not less correlated.
There is no correlation range to find, and a bounded model fitted to such data
reports whichever range the window implies.

Ordinary kriging copes by staying local. A small neighbourhood on a trending
surface is approximately stationary, which is why a 5 m radius beat a 20 m one
on l0s395 -- not because it preserved detail, but because it stayed inside the
region where the assumption holds.

Universal kriging removes the need for that bargain. It fits a low-order
polynomial drift within each neighbourhood and krige the residual, so the
neighbourhood can widen without the assumption failing.

WHAT IT COSTS
-------------
Three extra unknowns per solve for a linear drift, six for quadratic, and the
system can be ill-conditioned where the neighbours are nearly collinear -- along
a scan line, for instance, which is exactly how these returns are distributed.
Whether it is worth that is a measurement, not an argument, which is what
`compare_methods` is for.
"""
import numpy as np

from .kriging import MODELS, fit_variogram, KrigingError


def _drift_basis(px, py, cx, cy, order):
    """Polynomial terms in local coordinates, centred on the cell.

    Local coordinates keep the design matrix well scaled: UTM eastings near
    seven hundred thousand squared would swamp the variogram terms and make the
    system numerically hopeless.
    """
    dx, dy = px - cx, py - cy
    if order == 0:
        return np.ones((len(px), 1))
    if order == 1:
        return np.c_[np.ones(len(px)), dx, dy]
    if order == 2:
        return np.c_[np.ones(len(px)), dx, dy, dx * dx, dx * dy, dy * dy]
    raise KrigingError("drift order %r not supported" % order)


def krige_universal(x, y, z, grid, radius=20.0, max_points=32, min_points=6,
                    model="spherical", variogram=None, drift_order=1,
                    nodata=-9999.0, verbose=True):
    """Kriging with a polynomial drift, solved per cell.

    With `drift_order=0` this reduces to ordinary kriging and the two paths can
    be compared on identical code. Order 1 fits a tilted plane within each
    neighbourhood, which is what a regional gradient looks like locally.
    """
    from scipy.spatial import cKDTree

    x = np.asarray(x, float); y = np.asarray(y, float); z = np.asarray(z, float)
    n_terms = {0: 1, 1: 3, 2: 6}[drift_order]
    need = max(min_points, n_terms + 1)
    if len(z) < need:
        raise KrigingError("only %d points supplied" % len(z))

    if variogram is None:
        variogram = fit_variogram(x, y, z, model=model, max_lag=radius)
    mdl = MODELS[variogram["model"]]
    p = variogram["params"]

    tree = cKDTree(np.c_[x, y])
    gx, gy = grid.cell_centres()
    flat = np.c_[gx.ravel(), gy.ravel()]
    out = np.full(len(flat), nodata, float)

    dist, idx = tree.query(flat, k=min(max_points, len(z)),
                           distance_upper_bound=radius)
    if dist.ndim == 1:
        dist, idx = dist[:, None], idx[:, None]

    filled = n_fallback = 0
    for c in range(len(flat)):
        d = dist[c]
        good = np.isfinite(d)
        k = int(good.sum())
        if k < need:
            continue
        ii = idx[c][good]
        px, py, pz = x[ii], y[ii], z[ii]
        cx, cy = flat[c]

        F = _drift_basis(px, py, cx, cy, drift_order)      # k x n_terms
        dxy = np.hypot(px[:, None] - px[None, :], py[:, None] - py[None, :])
        G = mdl(dxy, *p)
        g0 = mdl(d[good], *p)

        m = n_terms
        A = np.zeros((k + m, k + m))
        A[:k, :k] = G
        A[:k, k:] = F
        A[k:, :k] = F.T
        b = np.zeros(k + m)
        b[:k] = g0
        # the drift is evaluated at the cell, where local coordinates are zero,
        # so only the constant term is one
        b[k] = 1.0

        # An estimate that leaves the range of its own neighbours means the
        # system was degenerate, which happens where the neighbours are nearly
        # collinear -- common along a scan line. Inverse distance is bounded by
        # construction and is the fallback.
        try:
            w = np.linalg.solve(A, b)
            est = float(w[:k] @ pz)
            bad = not np.isfinite(est) or est < pz.min() or est > pz.max()
        except np.linalg.LinAlgError:
            bad = True
        if bad:
            ww = 1.0 / np.maximum(d[good], 1e-9) ** 2
            est = float((ww @ pz) / ww.sum())
            n_fallback += 1
        out[c] = est
        filled += 1

    dem = out.reshape(grid.height, grid.width)
    frac = n_fallback / filled if filled else 0.0
    if verbose:
        print("  universal kriging, drift order %d: %d of %d cells (%.1f%%), "
              "radius %.1f m, %d nearest%s"
              % (drift_order, filled, len(flat), 100 * filled / len(flat),
                 radius, max_points,
                 "" if not n_fallback else ", %.1f%% fell back" % (100 * frac)))
    info = dict(variogram)
    info["drift_order"] = drift_order
    info["fallback_fraction"] = frac
    return dem, info
