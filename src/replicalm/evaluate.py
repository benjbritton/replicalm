"""Judge a surface against the point cloud, with no reference DEM involved.

WHY
---
Every measure used in this project so far -- roughness ratio, RMSE, speck count
-- is defined against NCALM's output. They answer "how closely does this
resemble the delivered product", which is the right question for a replication
and the wrong one for asking whether a surface is accurate. If a surface were
better than the delivered product, all three would get worse.

These three do not mention a reference. They ask whether the surface is a good
account of the returns that were actually measured.

    fidelity      buffered leave-one-out. Withhold a return and everything
                  within `buffer` of it, predict at its location from what is
                  left, and measure the residual. How well does the surface
                  predict real measurements it has never seen?

    penetration   returns lying below the surface. A pulse can be stopped early
                  by vegetation but cannot arrive below the soil, so a return
                  beneath the modelled ground is evidence that the surface is
                  floating on clutter. This is the measure with the sharpest
                  teeth, and the one an over-smoothed or clutter-contaminated
                  surface fails.

    sharpness     elevation change over a fixed one-metre baseline, sampled
                  bilinearly. Measured in metres rather than cells so that grids
                  of different resolution are comparable -- a finer grid does
                  not score higher merely for having smaller cells.

THE TENSION IS THE POINT
------------------------
Fidelity alone rewards a surface that honours every return including
vegetation, which is the exact-interpolation failure this project started from.
Penetration alone rewards a surface pushed down to the lowest return anywhere.
Sharpness alone rewards noise. A surface has to do well on all three, and
reading them together is what makes the judgement possible without a reference.

Fidelity must therefore be scored against returns believed to be ground. The
belief cannot come from TerraScan's labels without reintroducing the benchmark
this module exists to avoid, so it comes from the physical asymmetry instead:
the lower envelope of returns is the evidence, and a return far above its
neighbours is not a ground measurement to be predicted.
"""
import numpy as np


class EvaluateError(RuntimeError):
    pass


def buffered_loo(x, y, z, variogram, radius=20.0, max_points=16, min_points=3,
                 buffer=0.75, n_samples=2000, seed=0, trust=None):
    """Predict withheld returns from returns that are not near them.

    A plain leave-one-out is meaningless on this data: returns cluster along
    scan lines, so the nearest neighbour of any withheld point is centimetres
    away and the surface reproduces it almost exactly. Withholding the point's
    whole neighbourhood is what makes the question non-trivial.

    `trust` optionally restricts which returns are used as test targets, without
    restricting which are used to predict them.
    """
    from scipy.spatial import cKDTree
    from .kriging import MODELS

    x = np.asarray(x, float); y = np.asarray(y, float); z = np.asarray(z, float)
    n = len(z)
    if n < 100:
        raise EvaluateError("only %d returns" % n)
    rng = np.random.default_rng(seed)
    pool = np.flatnonzero(trust) if trust is not None else np.arange(n)
    if len(pool) == 0:
        raise EvaluateError("no trusted returns to test against")
    test = rng.choice(pool, min(n_samples, len(pool)), replace=False)

    tree = cKDTree(np.c_[x, y])
    mdl = MODELS[variogram["model"]]
    p = variogram["params"]

    res, used = [], 0
    for i in test:
        near = tree.query_ball_point([x[i], y[i]], buffer)
        excluded = np.zeros(n, bool)
        excluded[near] = True
        cand = tree.query_ball_point([x[i], y[i]], radius)
        cand = np.array([c for c in cand if not excluded[c]], dtype=int)
        if len(cand) < min_points:
            continue
        d = np.hypot(x[cand] - x[i], y[cand] - y[i])
        if len(cand) > max_points:
            keep = np.argsort(d)[:max_points]
            cand, d = cand[keep], d[keep]
        k = len(cand)
        px, py, pz = x[cand], y[cand], z[cand]
        dxy = np.hypot(px[:, None] - px[None, :], py[:, None] - py[None, :])
        A = np.empty((k + 1, k + 1))
        A[:k, :k] = mdl(dxy, *p)
        A[:k, k] = 1.0; A[k, :k] = 1.0; A[k, k] = 0.0
        b = np.empty(k + 1)
        b[:k] = mdl(d, *p); b[k] = 1.0
        try:
            w = np.linalg.solve(A, b)
            est = float(w[:k] @ pz)
            if not np.isfinite(est) or est < pz.min() - (pz.max() - pz.min()) \
                    or est > pz.max() + (pz.max() - pz.min()):
                raise np.linalg.LinAlgError
        except np.linalg.LinAlgError:
            ww = 1.0 / np.maximum(d, 1e-9) ** 2
            est = float((ww @ pz) / ww.sum())
        res.append(abs(est - z[i])); used += 1
    if not res:
        raise EvaluateError("no test point had enough support outside its buffer")
    res = np.array(res)
    return {"tested": used, "buffer_m": buffer,
            "median": float(np.median(res)), "mae": float(res.mean()),
            "rmse": float(np.sqrt((res ** 2).mean())),
            "p90": float(np.percentile(res, 90))}


def returns_below(dem, grid, x, y, z, nodata=-9999.0, tolerance=0.05):
    """Returns lying beneath the modelled surface, which none should.

    Reported as a share of the returns that fall on valid cells, plus how far
    the offenders sit below. A surface resting on vegetation shows this
    immediately: the true ground returns beneath the clutter end up under it.
    """
    dem = np.asarray(dem)
    x = np.asarray(x, float); y = np.asarray(y, float); z = np.asarray(z, float)
    col = np.floor((x - grid.origin_x) / grid.cell).astype(int)
    row = np.floor((grid.origin_y - y) / grid.dy).astype(int)
    ok = ((col >= 0) & (col < grid.width) & (row >= 0) & (row < grid.height))
    if not ok.any():
        raise EvaluateError("no returns fall on the grid")
    surf = dem[row[ok], col[ok]]
    good = np.isfinite(surf) & (surf != nodata) & (surf != 0.0)
    if not good.any():
        raise EvaluateError("no returns fall on a filled cell")
    depth = surf[good] - z[ok][good]          # positive where a return is below
    under = depth > tolerance
    return {"considered": int(good.sum()), "below": int(under.sum()),
            "below_fraction": float(under.mean()),
            "median_depth_m": float(np.median(depth[under])) if under.any() else 0.0,
            "p95_depth_m": float(np.percentile(depth[under], 95)) if under.any() else 0.0}


def sharpness(dem, grid, wkt=None, baseline_m=1.0, nodata=-9999.0,
              n_samples=200000, seed=0):
    """Elevation change over a fixed baseline in metres, not in cells.

    Sampled bilinearly at a physical offset so that a 0.25 m grid and a 1.0 m
    grid are asked the same question. Reported as percentiles of the gradient
    magnitude: a surface that resolves real breaklines carries a heavier upper
    tail, and one that is merely noisy carries a heavier middle.
    """
    dem = np.asarray(dem, float)
    valid = np.isfinite(dem) & (dem != nodata) & (dem != 0.0)
    if valid.sum() < 1000:
        raise EvaluateError("only %d valid cells" % int(valid.sum()))
    rng = np.random.default_rng(seed)
    rows, cols = np.nonzero(valid)
    pick = rng.choice(len(rows), min(n_samples, len(rows)), replace=False)
    rows, cols = rows[pick], cols[pick]
    dx = baseline_m / grid.cell
    dy = baseline_m / grid.dy

    def at(r, c):
        r0 = np.clip(np.floor(r).astype(int), 0, grid.height - 1)
        c0 = np.clip(np.floor(c).astype(int), 0, grid.width - 1)
        r1 = np.clip(r0 + 1, 0, grid.height - 1)
        c1 = np.clip(c0 + 1, 0, grid.width - 1)
        fr, fc = r - r0, c - c0
        v = (dem[r0, c0] * (1 - fr) * (1 - fc) + dem[r0, c1] * (1 - fr) * fc +
             dem[r1, c0] * fr * (1 - fc) + dem[r1, c1] * fr * fc)
        m = (valid[r0, c0] & valid[r0, c1] & valid[r1, c0] & valid[r1, c1])
        return v, m

    east, me = at(rows.astype(float), cols + dx)
    west, mw = at(rows.astype(float), cols - dx)
    north, mn = at(rows - dy, cols.astype(float))
    south, ms = at(rows + dy, cols.astype(float))
    good = me & mw & mn & ms
    if good.sum() < 100:
        raise EvaluateError("too few sampleable positions")
    gx = (east[good] - west[good]) / (2 * baseline_m)
    gy = (south[good] - north[good]) / (2 * baseline_m)
    g = np.hypot(gx, gy)
    return {"sampled": int(good.sum()), "baseline_m": baseline_m,
            "median": float(np.median(g)), "p90": float(np.percentile(g, 90)),
            "p99": float(np.percentile(g, 99))}
