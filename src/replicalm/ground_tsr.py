"""Ground classification that reproduces the NCALM delivery from raw returns.

The target is the ground population a TerraScan/Surfer chain produces -- the
Class 2 and Class 8 union of an NCALM delivery, treated as one class, since the
delivery carries no attribute by which the two could be separated. Those labels
are used to calibrate and validate this routine. They are never read while it
runs, and nothing here requires them to be present.

WHY NOT A SINGLE FILTER

A morphological opening of a minimum surface is conservative: run at the
Replicalm baseline it keeps almost no vegetation but discards real ground on
closed landforms, worst on mounds, and the resulting surface sits high over
depression floors. A multiscale curvature test is the opposite: it recovers
most of that ground and admits enough low vegetation to make the surface worse
overall. Neither error is reachable by tuning the other filter, because they
are different mistakes -- one drops true ground, the other keeps false ground.

So this routine separates the two decisions. A conservative filter supplies
points that are almost certainly ground; a permissive one supplies everything
that might be; and a third stage adjudicates between them against a surface
built only from the confident set, with a tolerance scaled to how rough and how
densely sampled the neighbourhood actually is. The adjudication is what the
single filters cannot do, because it is the only stage that sees both answers.

STAGES

  1  density       local returns per square metre on a coarse grid, which sets
                   the working cell and the floor under every tolerance
  2  confident     morphological opening, left conservative on purpose
  3  candidates    multiscale curvature tolerance, permissive on purpose
  4  reference     robust surface from the confident set alone: a low
                   quantile per cell, holes filled by normalised convolution
  5  adjudicate    accept a candidate whose residual to that surface falls
                   inside a locally scaled, asymmetric gate
  6  refine        rebuild the surface from the accepted set and adjudicate
                   once more

Stage 5 is asymmetric because the two sides are not symmetric in nature.
Anything standing above the ground is vegetation, structure or noise; anything
below is either real microtopography or a dropout that the noise filter should
already have taken. The gate is therefore tighter above than below, by a
factor of about five at the calibrated settings.

The gate adds to the confident set and never subtracts from it. Letting it
reject confident points as well was implemented and swept, and it fails badly:
recall falls to between 0.39 and 0.74 while precision goes to 1.000, because
the confident points scatter about the quantile surface by more than a gate
tuned for candidates allows. Every pruning configuration scored worse than
every additive one. `prune_confident` is kept, defaulting off, so the result
is reproducible rather than merely asserted.

What the gate is worth, on eight development blocks: the opening alone scores
0.0470 m against the reference surface, the curvature test alone 0.0518 m, and
the adjudicated combination 0.0434 m. The gain is about 8 percent over the
better of the two filters, and it is consistent across every gate setting in
the additive family, which is what distinguishes it from a lucky parameter.

OPERATIONAL LIMITS

Calibrated on G-LiHT Riegl VQ-480 coverage over tropical karst with closed
canopy, at a measured 19 to 24 returns per square metre on the calibration
blocks -- not the 4 to 11 per square metre of ground returns, which is a
different quantity and the one earlier Replicalm work quotes. The density scaling in
stage 1 is what carries it to other acquisitions, and it is interpolation, not
extrapolation: outside the calibrated density band the working cell is clamped
rather than projected. Buildings are not modelled -- this is a bare-earth
routine for terrain, and a roof large enough to defeat the opening window will
survive into the ground class. Water is not detected.
"""
from dataclasses import dataclass, field

import numpy as np

try:
    from scipy import ndimage
except ImportError as exc:                          # pragma: no cover
    raise ImportError("ground_tsr needs scipy") from exc


class GroundTSRError(RuntimeError):
    pass


@dataclass
class GroundTSRParams:
    """Defaults are the calibrated configuration. See `scale_for_density`.

    Every length is in metres and every tolerance in metres of elevation.
    """

    # stage 1, density
    density_cell_m: float = 10.0
    density_ref: float = 20.0            # returns/m2 the defaults were set at
    cell_min_m: float = 0.5
    cell_max_m: float = 2.0

    # stage 2, the confident set: a deliberately conservative opening
    smrf_window_m: float = 32.0
    smrf_slope: float = 0.158
    smrf_threshold_m: float = 0.25
    smrf_scalar: float = 0.0
    smrf_cell_m: float = 1.0

    # stage 3, the candidate set: a deliberately permissive curvature test
    mcc_scales_m: tuple = (0.5, 1.0, 1.5)
    mcc_tolerance_m: float = 0.30
    mcc_iterations: int = 4

    # stage 4, the reference surface built from the confident set
    surface_cell_m: float = 2.0
    surface_quantile: float = 0.10
    surface_smooth_m: float = 3.0

    # stage 5, the gate
    gate_scale_m: float = 8.0            # window the local spread is read over
    gate_k_above: float = 0.75           # multiples of local spread, upward
    gate_k_below: float = 4.0            # and downward
    gate_floor_above_m: float = 0.05
    gate_floor_below_m: float = 0.15
    gate_ceiling_m: float = 1.00
    slope_allowance: float = 0.30        # metres of gate per unit of slope

    # stage 6
    refine_iterations: int = 1
    prune_confident: bool = False        # measured and rejected: see below
    min_ground_fraction: float = 0.02    # safety net under pruning

    # noise handling, applied before anything else
    outlier_neighbours: int = 8
    outlier_multiplier: float = 2.5

    def scale_for_density(self, density):
        """Return a copy with lengths adapted to measured return density.

        Point spacing goes as the inverse square root of density, so the
        working cell is scaled that way and then clamped. Clamping rather than
        extrapolating is deliberate: beyond the calibrated band the right
        behaviour is the nearest configuration that was actually measured, not
        a projection of a trend that was never tested there.
        """
        if not np.isfinite(density) or density <= 0:
            return self
        f = float(np.sqrt(self.density_ref / density))
        out = GroundTSRParams(**{k: getattr(self, k)
                                 for k in self.__dataclass_fields__})
        out.smrf_cell_m = float(np.clip(self.smrf_cell_m * f,
                                        self.cell_min_m, self.cell_max_m))
        out.surface_cell_m = float(np.clip(self.surface_cell_m * f,
                                           2 * self.cell_min_m,
                                           2 * self.cell_max_m))
        out.mcc_scales_m = tuple(float(np.clip(s * f, 0.25, 4.0))
                                 for s in self.mcc_scales_m)
        return out


# --------------------------------------------------------------------------
# stage 1


def return_density(x, y, cell_m=10.0):
    """Median returns per square metre over occupied cells of `cell_m`.

    The median rather than the mean, so a strip of overlap between flight
    lines does not set the working cell for the whole block.
    """
    if len(x) == 0:
        return 0.0
    ci = np.floor((x - x.min()) / cell_m).astype(np.int64)
    cj = np.floor((y - y.min()) / cell_m).astype(np.int64)
    k = ci * 1000003 + cj
    cnt = np.bincount(np.unique(k, return_inverse=True)[1])
    return float(np.median(cnt) / (cell_m * cell_m))


# --------------------------------------------------------------------------
# stages 2 and 3


def _pdal():
    from . import classify
    return classify._pdal()


def denoise(x, y, z, p):
    """Statistical outlier removal, the same step the pipeline applies."""
    import json
    pdal = _pdal()
    a = np.zeros(len(z), dtype=[("X", "f8"), ("Y", "f8"), ("Z", "f8"),
                                ("Classification", "u1")])
    a["X"], a["Y"], a["Z"] = x, y, z
    stg = [{"type": "filters.assign", "value": "Classification = 0"},
           {"type": "filters.outlier", "method": "statistical",
            "mean_k": p.outlier_neighbours,
            "multiplier": p.outlier_multiplier, "class": 7},
           {"type": "filters.expression", "expression": "Classification != 7"}]
    pl = pdal.Pipeline(json.dumps({"pipeline": stg}), arrays=[a])
    pl.execute()
    g = pl.arrays[0]
    return (g["X"].astype("f8"), g["Y"].astype("f8"), g["Z"].astype("f8"))


def confident_ground(x, y, z, p):
    """Morphological opening, kept conservative. Returns a boolean mask."""
    import json
    pdal = _pdal()
    a = np.zeros(len(z), dtype=[("X", "f8"), ("Y", "f8"), ("Z", "f8"),
                                ("Classification", "u1"),
                                ("PointSourceId", "u2")])
    a["X"], a["Y"], a["Z"] = x, y, z
    a["PointSourceId"] = np.arange(len(z), dtype="u2") if len(z) < 65536 \
        else 0
    stg = [{"type": "filters.assign", "value": "Classification = 0"},
           {"type": "filters.smrf", "window": p.smrf_window_m,
            "slope": p.smrf_slope, "threshold": p.smrf_threshold_m,
            "scalar": p.smrf_scalar, "cell": p.smrf_cell_m},
           {"type": "filters.expression", "expression": "Classification == 2"}]
    pl = pdal.Pipeline(json.dumps({"pipeline": stg}), arrays=[a])
    pl.execute()
    g = pl.arrays[0]
    # match back on coordinates; SMRF neither reorders nor edits geometry, but
    # it does drop points, so an index is rebuilt rather than assumed
    return _match_mask(x, y, z, g["X"], g["Y"], g["Z"])


def _match_mask(x, y, z, gx, gy, gz):
    key = _coord_key(x, y, z)
    gk = _coord_key(gx.astype("f8"), gy.astype("f8"), gz.astype("f8"))
    o = np.argsort(key, kind="stable")
    ks = key[o]
    pos = np.searchsorted(ks, gk)
    pos = np.clip(pos, 0, len(ks) - 1)
    hit = ks[pos] == gk
    out = np.zeros(len(z), bool)
    out[o[pos[hit]]] = True
    return out


def _coord_key(x, y, z):
    xi = np.round(x * 1000).astype(np.int64)
    yi = np.round(y * 1000).astype(np.int64)
    zi = np.round(z * 1000).astype(np.int64)
    return (xi * 1000003 + yi) * 1000003 + zi


def candidate_ground(x, y, z, p):
    """Multiscale curvature tolerance, kept permissive.

    At each scale a mean surface is formed from the points still in play,
    holes are closed by nearest-neighbour fill, the surface is smoothed, and
    points standing more than the tolerance above it are dropped. Coarse
    scales remove canopy, fine scales remove what is left in the understorey.
    """
    keep = np.ones(len(z), bool)
    ox, oy = x.min(), y.max()
    for lam in p.mcc_scales_m:
        for _ in range(p.mcc_iterations):
            w = int((x.max() - ox) / lam) + 2
            h = int((oy - y.min()) / lam) + 2
            c = np.clip(((x - ox) / lam).astype(np.int64), 0, w - 1)
            r = np.clip(((oy - y) / lam).astype(np.int64), 0, h - 1)
            S = np.zeros((h, w))
            N = np.zeros((h, w))
            np.add.at(S, (r[keep], c[keep]), z[keep])
            np.add.at(N, (r[keep], c[keep]), 1.0)
            empty = N < 1
            if empty.all():
                break
            with np.errstate(invalid="ignore"):
                M = np.where(empty, np.nan, S / np.maximum(N, 1))
            idx = ndimage.distance_transform_edt(
                empty, return_distances=False, return_indices=True)
            M = np.where(empty, M[tuple(idx)], M)
            M = ndimage.gaussian_filter(np.nan_to_num(
                M, nan=float(np.nanmean(M))), 1.0)
            drop = keep & (z - M[r, c] > p.mcc_tolerance_m)
            if not drop.any():
                break
            keep &= ~drop
    return keep


# --------------------------------------------------------------------------
# stage 4


def _grid_shape(x, y, cell):
    ox, oy = x.min(), y.max()
    w = int((x.max() - ox) / cell) + 2
    h = int((oy - y.min()) / cell) + 2
    return ox, oy, w, h


def _fill(M, ok, smooth_cells):
    """Normalised convolution: smooth the surface and its mask separately.

    Filling holes with a scene constant before smoothing would put a step at
    every hole edge and feed it into everything computed downstream, which is
    a mistake this code has made before.
    """
    num = ndimage.gaussian_filter(np.where(ok, M, 0.0), smooth_cells)
    den = ndimage.gaussian_filter(ok.astype("f8"), smooth_cells)
    sm = np.where(den > 1e-6, num / np.maximum(den, 1e-6), np.nan)
    if not np.isfinite(sm).all():
        bad = ~np.isfinite(sm)
        if bad.all():
            return np.zeros_like(sm)
        idx = ndimage.distance_transform_edt(
            bad, return_distances=False, return_indices=True)
        sm = np.where(bad, sm[tuple(idx)], sm)
    return sm


def reference_surface(x, y, z, mask, p):
    """A low quantile of the confident set per cell, smoothed and gap-filled.

    A quantile rather than a minimum: the minimum of a cell follows the lowest
    noise spike that survived denoising, while a low quantile tracks the soil
    and is stable at the densities this runs at.
    """
    cell = p.surface_cell_m
    ox, oy, w, h = _grid_shape(x, y, cell)
    c = np.clip(((x[mask] - ox) / cell).astype(np.int64), 0, w - 1)
    r = np.clip(((oy - y[mask]) / cell).astype(np.int64), 0, h - 1)
    zz = z[mask]
    flat = r * w + c
    order = np.argsort(flat, kind="stable")
    fs, zs = flat[order], zz[order]
    uniq, start = np.unique(fs, return_index=True)
    end = np.append(start[1:], len(fs))
    M = np.full(h * w, np.nan)
    for u, a_, b_ in zip(uniq, start, end):
        M[u] = np.quantile(zs[a_:b_], p.surface_quantile)
    M = M.reshape(h, w)
    ok = np.isfinite(M)
    if not ok.any():
        raise GroundTSRError("no confident ground to build a surface from")
    sm = _fill(M, ok, max(p.surface_smooth_m / cell / 2.0, 0.5))
    return sm, (ox, oy, cell, w, h)


def _sample(S, geom, x, y):
    ox, oy, cell, w, h = geom
    c = np.clip(((x - ox) / cell).astype(np.int64), 0, w - 1)
    r = np.clip(((oy - y) / cell).astype(np.int64), 0, h - 1)
    return S[r, c]


def _local_spread(x, y, resid, mask, geom, p):
    """Median absolute residual of the confident set, per neighbourhood.

    This is what makes the gate adaptive. Where the confident points scatter
    about the surface -- rough ground, sparse coverage, steep faces -- the
    gate opens; on smooth well-sampled ground it closes to its floor.
    """
    ox, oy, cell, w, h = geom
    c = np.clip(((x[mask] - ox) / cell).astype(np.int64), 0, w - 1)
    r = np.clip(((oy - y[mask]) / cell).astype(np.int64), 0, h - 1)
    a = np.abs(resid[mask])
    S = np.zeros((h, w))
    N = np.zeros((h, w))
    np.add.at(S, (r, c), a)
    np.add.at(N, (r, c), 1.0)
    sig = max(p.gate_scale_m / cell / 2.0, 0.5)
    num = ndimage.gaussian_filter(S, sig)
    den = ndimage.gaussian_filter(N, sig)
    mean_abs = np.where(den > 1e-6, num / np.maximum(den, 1e-6), np.nan)
    ok = np.isfinite(mean_abs)
    if not ok.any():
        return np.zeros((h, w))
    # 1.4826 takes a mean absolute deviation to a Gaussian-equivalent sigma
    return _fill(1.4826 * mean_abs, ok, 1.0)


def _surface_slope(S, geom):
    _ox, _oy, cell, _w, _h = geom
    gy, gx = np.gradient(S, cell)
    return np.hypot(gx, gy)


# --------------------------------------------------------------------------
# stages 5 and 6


def classify_ground(x, y, z, params=None, adapt_density=True, verbose=False,
                    keep_masks=False):
    """Ground mask over the input returns. Raw input only.

    Returns (mask, report). `mask` is boolean over x/y/z as given. `report`
    carries the measured density, the effective parameters and the per-stage
    counts, so a run can be audited without being repeated.

    `keep_masks` additionally puts the stage 2 and stage 3 masks into the
    report under "confident_mask" and "candidate_mask". It is off by default
    because the arrays are the size of the input; it exists so the stages can
    be scored against each other on identical data in one pass, rather than by
    classifying the same block two or three times and hoping nothing else
    varied between runs.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    if not (len(x) == len(y) == len(z)):
        raise GroundTSRError("x, y and z must be the same length")
    if len(z) < 32:
        raise GroundTSRError("need at least 32 returns; got %d" % len(z))

    p = params or GroundTSRParams()
    dens = return_density(x, y, p.density_cell_m)
    if adapt_density:
        p = p.scale_for_density(dens)
    rep = {"returns": int(len(z)), "density_per_m2": dens,
           "smrf_cell_m": p.smrf_cell_m, "surface_cell_m": p.surface_cell_m,
           "mcc_scales_m": list(p.mcc_scales_m)}
    if verbose:
        print("  ground_tsr: %d returns, %.2f returns/m2, cell %.2f m"
              % (len(z), dens, p.smrf_cell_m), flush=True)

    dx, dy, dz = denoise(x, y, z, p)
    alive = _match_mask(x, y, z, dx, dy, dz)
    rep["after_denoise"] = int(alive.sum())

    conf = np.zeros(len(z), bool)
    conf[alive] = confident_ground(x[alive], y[alive], z[alive], p)
    rep["confident"] = int(conf.sum())
    if conf.sum() < 16:
        raise GroundTSRError("the conservative filter returned %d points; "
                             "the block is too small or too sparse"
                             % int(conf.sum()))

    cand = np.zeros(len(z), bool)
    cand[alive] = candidate_ground(x[alive], y[alive], z[alive], p)
    cand |= conf
    rep["candidates"] = int(cand.sum())

    accepted = conf.copy()
    basis = conf
    for it in range(1 + max(p.refine_iterations, 0)):
        S, geom = reference_surface(x, y, z, basis, p)
        resid = z - _sample(S, geom, x, y)
        spread = _local_spread(x, y, resid, basis, geom, p)
        sl = _surface_slope(S, geom)
        sp = _sample(spread, geom, x, y)
        slp = _sample(sl, geom, x, y)
        extra = p.slope_allowance * slp
        hi = np.clip(p.gate_k_above * sp + extra, p.gate_floor_above_m,
                     p.gate_ceiling_m)
        lo = np.clip(p.gate_k_below * sp + extra, p.gate_floor_below_m,
                     p.gate_ceiling_m)
        accepted = cand & (resid <= hi) & (resid >= -lo)
        if not p.prune_confident:
            accepted |= conf
        elif accepted.sum() < p.min_ground_fraction * len(z):
            # the gate has rejected almost everything, which means the
            # surface it was measured against is wrong rather than that the
            # ground is absent. Fall back rather than emit a ruined block.
            accepted = accepted | conf
        rep["accepted_iter%d" % it] = int(accepted.sum())
        if verbose:
            print("    iteration %d: gate above %.3f-%.3f m, below "
                  "%.3f-%.3f m, accepted %d"
                  % (it, float(np.nanmin(hi)), float(np.nanmax(hi)),
                     float(np.nanmin(lo)), float(np.nanmax(lo)),
                     int(accepted.sum())), flush=True)
        basis = accepted

    rep["ground"] = int(accepted.sum())
    rep["ground_fraction"] = float(accepted.mean())
    if keep_masks:
        rep["confident_mask"] = conf
        rep["candidate_mask"] = cand
    return accepted, rep
