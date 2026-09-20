"""Remove residual low vegetation from the ground class.

THE CLEAR RULE
--------------
`clear()` is the working rule and the one the "Clear" product uses: demote a
ground return standing more than 0.20 m above the tenth percentile of ground
returns within 0.75 m. Measured on the Pixoyal window of South_GLAS_l0s395,
against labels taken from the delivered cloud's own TerraScan classification:

    recall of what TerraScan rejected     92.5%
    flat-ground roughness    1.26x the archive -> 0.99x
    specks                             3,901 -> 426
    RMSE                              0.0825 -> 0.0451
    error above 20 degrees              1.86% -> 0.29%

Validated unchanged on three windows it was not fitted to -- another part of
l0s395, and the relief windows of l0s417 and l2s444 -- improving RMSE and
steep-ground error on all three.

WHAT IT COSTS
-------------
It removes 17.75% of ground returns where an oracle using TerraScan's own
labels removes 3.88% for a better result (roughness 1.00x, specks 253, RMSE
0.0140). The difference is collateral: genuine returns, including some on
platform edges, which is visible as slightly softer architecture. Closing that
gap is what the "Deep" work is for.

NOT IN THE LOCKED BASELINE. Adopting it changes the baseline and needs the
usual re-verification.

THE RULES THAT DID NOT WORK, kept so they are not re-derived: The module is kept because the
negative results are worth not re-deriving, and because the diagnosis they
failed to act on is sound.

    height above the patch floor   demotes 72% of synthetic brush and 100% of a
                                   synthetic rock of the same footprint. A
                                   feature narrower than the patch has
                                   neighbours on the flat ground around it, so
                                   the floor stays low and the whole feature
                                   reads as high.
    spread, slope-compensated      54% of the brush, 92% of the rock, and 0% of
                                   brush sitting on a 20 degree slope. A rock's
                                   own edge produces spread.
    source disagreement            on the real window: demotes 3.1% of ground
                                   returns, touches 3.4% of speck cells and
                                   10.2% of steep cells. It damages architecture
                                   three times more than it cleans clutter,
                                   because on a slope two source chunks sample
                                   different parts of the same real surface and
                                   are recorded as disagreeing.

WHY THIS IS HARD
----------------
At four ground returns per square metre, a 0.4 m bush and a 0.4 m rock are each
described by about five points, and the geometry is the same. The discriminator
that does exist statistically -- source chunks disagreeing by 0.47 m at specks
against 0.11 m on control ground -- is available at only 27.4% of speck cells,
and slope confounds it everywhere else.

WHAT HAS NOT BEEN TRIED
-----------------------
Every test here looks only at points already classified as ground. The full
return stack has not been examined: a pulse that stopped on brush should have
been preceded by nothing and followed by nothing, whereas genuine ground under
scrub often shows a return above it within a short vertical distance. That
information is discarded before any of these rules see the data.

THE DIAGNOSIS THESE RULES WERE BUILT ON
---------------------------------------


WHAT THIS IS FOR
----------------
Measured on the l0s395 trench window, the cells that render as dark specks on
otherwise flat ground carry a median of five classified ground returns each.
Those returns sit 5.8 cm above the surrounding 3 m median, 76.2% of them
positive -- a one-sided bias, so not ranging noise, which would be symmetric.
But they scatter among themselves by 8.8 cm, more than the feature they
describe, and where two source chunks both contribute they disagree by 0.47 m
against 0.11 m on control ground. Mean scan angle is identical to control
(11.5 against 11.2 degrees), so it is not incidence geometry, and the returns
are 11% darker than control ground.

Something at those locations intercepts pulses at variable heights, regardless
of viewing angle, and reflects less than soil or limestone. Low vegetation --
scrub, brush, root mass -- does exactly that.

WHY HEIGHT ALONE CANNOT BE THE TEST
-----------------------------------
The obvious rule -- demote a return standing above the floor of its
neighbourhood -- destroys real terrain. Tested on synthetic ground carrying a
0.4 m brush clump and a 0.4 m rock of the same footprint, it demoted 72% of the
brush and 100% of the rock, because a feature narrower than the search radius
has neighbours on the flat ground around it, so the local floor stays low and
the whole feature reads as high.

Height above the neighbourhood is not what distinguishes clutter from terrain.
Agreement is. A rock raises every return that lands on it to the same height;
brush returns whatever height each pulse happened to stop at. That is precisely
what was measured -- 0.47 m of disagreement at specks against 0.11 m on control
ground -- and it is the property this module tests.

THE RULE
--------
Within each patch, compare the spread of ground-return heights against what the
local slope can account for. Where the spread is larger than terrain explains,
the patch is inconsistent, and the returns standing above its floor are demoted;
the lowest survive, because a pulse can be stopped early by vegetation but
cannot arrive below the surface. Where the spread is consistent -- flat ground
or a rock, both -- nothing is touched at any height.

Slope is estimated from the floor surface rather than from all returns, so that
clutter does not inflate the slope estimate and protect itself.
"""
import numpy as np


class CleanupError(RuntimeError):
    pass


def _cell_stats(x, y, z, cell, percentile):
    """Per-cell floor, spread and index, on a regular lattice."""
    ix = np.floor((x - x.min()) / cell).astype(np.int64)
    iy = np.floor((y - y.min()) / cell).astype(np.int64)
    nx = int(ix.max()) + 1
    flat_idx = iy * nx + ix
    order = np.argsort(flat_idx, kind="stable")
    sorted_idx = flat_idx[order]
    starts = np.flatnonzero(np.r_[True, sorted_idx[1:] != sorted_idx[:-1]])
    ends = np.r_[starts[1:], len(sorted_idx)]
    cell_id = sorted_idx[starts]

    floor = np.empty(len(starts))
    spread = np.empty(len(starts))
    count = (ends - starts).astype(np.int64)
    for k, (a, b) in enumerate(zip(starts, ends)):
        zz = z[order[a:b]]
        if len(zz) >= 3:
            lo, hi = np.percentile(zz, [percentile, 100.0 - percentile])
        else:
            lo, hi = zz.min(), zz.max()
        floor[k] = lo
        spread[k] = hi - lo
    lookup = {int(c): k for k, c in enumerate(cell_id)}
    which = np.array([lookup[int(c)] for c in flat_idx], dtype=np.int64)
    return which, floor, spread, count, cell_id, nx


def inconsistent_mask(x, y, z, cell=0.6, tolerance=0.15, spread_limit=0.20,
                      percentile=10.0, min_points=4, slope_factor=1.5):
    """True where a ground return stands above an inconsistent patch.

    `cell`          patch size, metres. 0.6 m is the scale at which the
                    disagreement was measured.
    `spread_limit`  how much height disagreement a patch may carry, beyond what
                    its slope explains, before it is judged inconsistent.
    `tolerance`     how far above the patch floor a return may sit and be kept,
                    once the patch is judged inconsistent.
    `percentile`    which return defines the floor. The minimum is too sensitive
                    to a single low outlier.
    `slope_factor`  multiplier on the slope-explained spread, for the fact that
                    a patch is wider than one cell across its diagonal.
    """
    x = np.asarray(x, float); y = np.asarray(y, float); z = np.asarray(z, float)
    if not (len(x) == len(y) == len(z)):
        raise CleanupError("x, y and z differ in length")
    if len(z) == 0:
        return np.zeros(0, bool)

    which, floor, spread, count, cell_id, nx = _cell_stats(
        x, y, z, cell, percentile)

    # local slope from the floor surface: clutter sits above the floor, so it
    # cannot inflate the slope estimate and thereby excuse itself
    ny = int(cell_id.max()) // nx + 1
    grid = np.full(nx * ny, np.nan)
    grid[cell_id] = floor
    grid = grid.reshape(ny, nx)
    gy, gx = np.gradient(np.where(np.isfinite(grid), grid,
                                  np.nanmean(grid)), cell, cell)
    grad = np.hypot(gx, gy).ravel()[cell_id]
    explained = slope_factor * grad * cell

    inconsistent = (spread > explained + spread_limit) & (count >= min_points)
    above = z - floor[which]
    return inconsistent[which] & (above > tolerance)


def demote(arr, cell=0.6, tolerance=0.15, spread_limit=0.20, percentile=10.0,
           min_points=4, ground_class=2, to_class=1, verbose=True):
    """Reclassify inconsistent high ground returns. Returns (array, report).

    Demoted rather than deleted: the points stay in the file under `to_class`,
    so the decision is inspectable and reversible rather than silently
    destructive.
    """
    arr = arr.copy()
    is_ground = arr["Classification"] == ground_class
    n_ground = int(is_ground.sum())
    if n_ground == 0:
        return arr, {"ground": 0, "demoted": 0, "fraction": 0.0}

    idx = np.flatnonzero(is_ground)
    bad = inconsistent_mask(arr["X"][idx], arr["Y"][idx], arr["Z"][idx],
                            cell=cell, tolerance=tolerance,
                            spread_limit=spread_limit, percentile=percentile,
                            min_points=min_points)
    arr["Classification"][idx[bad]] = to_class
    report = {"ground": n_ground, "demoted": int(bad.sum()),
              "fraction": float(bad.sum() / n_ground), "cell_m": cell,
              "tolerance_m": tolerance, "spread_limit_m": spread_limit,
              "percentile": percentile}
    if verbose:
        print("  demoted %d of %d ground returns (%.2f%%): patches "
              "disagreeing by more than %.2f m beyond their slope, returns "
              "more than %.2f m above the patch floor"
              % (report["demoted"], n_ground, 100 * report["fraction"],
                 spread_limit, tolerance))
    return arr, report

def source_disagreement_mask(x, y, z, source, cell=0.6, tolerance=0.12,
                             disagreement=0.20, min_sources=2, min_points=4):
    """True where independent source chunks disagree about the ground height.

    THE ONE DISCRIMINATOR THAT SURVIVED TESTING
    -------------------------------------------
    Height above the local floor and height spread within a patch both fail:
    tested on synthetic ground carrying a 0.4 m brush clump and a 0.4 m rock of
    the same footprint, they demoted 72% and 54% of the brush but also 100% and
    92% of the rock. At this scale a small rock and a small bush are the same
    shape, and a rock's own edge produces spread.

    What a rock cannot do is give two different answers. Every pulse that lands
    on it returns the same height. Vegetation returns whatever height it
    happened to stop the pulse at, which is why specks showed 0.47 m of
    disagreement between source chunks against 0.11 m on control ground.

    So the patch is judged by whether its sources agree, and only the returns
    standing above the lowest-reading source are demoted -- a pulse can be
    stopped early but cannot arrive below the surface.

    LIMIT, STATED PLAINLY
    ---------------------
    This can only judge a patch where two or more sources actually overlap. On
    the measured window that was 27.4% of speck cells. The remaining 72.6% are
    not addressed by this rule and are left alone, which is the correct
    behaviour for a test that cannot see them, but it means this is a partial
    fix rather than a complete one.
    """
    x = np.asarray(x, float); y = np.asarray(y, float); z = np.asarray(z, float)
    source = np.asarray(source)
    if len(z) == 0:
        return np.zeros(0, bool)

    ix = np.floor((x - x.min()) / cell).astype(np.int64)
    iy = np.floor((y - y.min()) / cell).astype(np.int64)
    nx = int(ix.max()) + 1
    key = iy * nx + ix
    order = np.argsort(key, kind="stable")
    ks = key[order]
    starts = np.flatnonzero(np.r_[True, ks[1:] != ks[:-1]])
    ends = np.r_[starts[1:], len(ks)]

    out = np.zeros(len(z), bool)
    for a, b in zip(starts, ends):
        members = order[a:b]
        if len(members) < min_points:
            continue
        src = source[members]
        uniq = np.unique(src)
        if len(uniq) < min_sources:
            continue
        means = np.array([z[members[src == u]].mean() for u in uniq])
        if means.max() - means.min() <= disagreement:
            continue                      # the sources agree: real surface
        floor = means.min()
        out[members] = z[members] - floor > tolerance
    return out

def height_above_floor(x, y, z, patch=0.75, percentile=10.0):
    """Height of each return above the low percentile of its neighbourhood."""
    from scipy.spatial import cKDTree
    x = np.asarray(x, float); y = np.asarray(y, float); z = np.asarray(z, float)
    if len(z) == 0:
        return np.zeros(0, float)
    tree = cKDTree(np.c_[x, y])
    nb = tree.query_ball_point(np.c_[x, y], patch)
    floor = np.array([np.percentile(z[i], percentile) if len(i) >= 3 else z[j]
                      for j, i in enumerate(nb)])
    return z - floor


def clear(arr, patch=0.75, height=0.20, percentile=10.0, ground_class=2,
          to_class=1, verbose=True):
    """The Clear rule. Returns (array, report).

    Demoted rather than deleted, so the decision stays inspectable.

    `height` is the one number worth tuning per survey. Lower removes more
    clutter and more genuine ground with it; the measured trade on the Pixoyal
    window was 0.20 m at 92.5% recall and 17.75% removed, 0.30 m at 73.2% and
    7.66%, 0.40 m at 51.8% and 3.60%.
    """
    arr = arr.copy()
    is_ground = arr["Classification"] == ground_class
    n = int(is_ground.sum())
    if n == 0:
        return arr, {"ground": 0, "demoted": 0, "fraction": 0.0}
    idx = np.flatnonzero(is_ground)
    above = height_above_floor(arr["X"][idx], arr["Y"][idx], arr["Z"][idx],
                               patch=patch, percentile=percentile)
    bad = above > height
    arr["Classification"][idx[bad]] = to_class
    report = {"ground": n, "demoted": int(bad.sum()),
              "fraction": float(bad.sum() / n), "patch_m": patch,
              "height_m": height, "percentile": percentile, "rule": "clear"}
    if verbose:
        print("  Clear: demoted %d of %d ground returns (%.2f%%) standing more "
              "than %.2f m above the %gth percentile within %.2f m"
              % (report["demoted"], n, 100 * report["fraction"], height,
                 percentile, patch))
    return arr, report
