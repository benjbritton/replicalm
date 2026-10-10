"""Conditioning that every surface solver needs, in one place.

WHY THIS MODULE EXISTS
----------------------
Coincident XY locations make a linear system exactly singular. Two returns at
one place give two identical rows, which no nugget, ridge term or neighbour
count can rescue -- the matrix has no inverse. The delivered G-LiHT clouds carry
them in quantity: 13.05% of the returns on one tile are exact XYZ repeats, same
return number, same GPS time, same source, while another tile has none. It is a
property of the archive rather than of the processing, and on that tile it put
89% of cells onto the inverse-distance fallback before it was found.

Ordinary kriging carried its own copy of the fix. Then a thin-plate spline was
tried and failed with `LinAlgError: Singular matrix` on the same data, for the
same reason, because the fix lived inside the kriging solver rather than ahead
of it. A radial basis solver, a minimum-curvature solver and a Delaunay
triangulation all share this requirement and none of them should have to know
about it.

So the conditioning happens once, upstream, and every solver takes points that
are already distinct.

WHAT IT DOES NOT DO
-------------------
It does not smooth, thin, filter or reweight. Removing an exact repeat loses no
information. Where one XY carries different elevations -- about 0.1% of groups
-- the lowest is kept, because a ground surface cannot pass through two heights
at one place and the lower return is the ground candidate. That is the only
judgement in here, it is stated, and it is the same judgement the kriging path
has always made.
"""
import numpy as np


class PreprocessError(RuntimeError):
    pass


def dedupe_xy(x, y, z, keep="lowest", round_to=None):
    """One point per XY location. Returns (x, y, z, report).

    `keep` chooses which elevation survives a collision: "lowest" (the ground
    candidate, and the default), "highest", or "first" in input order.

    `round_to` snaps coordinates to a grid before comparing, in the same units
    as the coordinates. Leave it None for exact matching, which is what the
    delivered clouds need -- their repeats are bit-identical. Set it where a
    solver is sensitive to near-coincidence rather than to exact coincidence;
    a spline with a small smoothing term still conditions badly on two points
    a micron apart.

    The report says what happened, so a caller can log it rather than discover
    the thinning later as a surprise.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    if not (len(x) == len(y) == len(z)):
        raise PreprocessError("x, y and z must be the same length")
    n_in = len(z)
    if n_in == 0:
        return x, y, z, {"points_in": 0, "points_out": 0, "dropped": 0,
                         "collision_groups": 0, "keep": keep,
                         "round_to": round_to}

    kx = x if round_to is None else np.round(x / round_to) * round_to
    ky = y if round_to is None else np.round(y / round_to) * round_to

    # Sort by location, then by elevation in the direction that puts the
    # survivor first, so taking the first of each group is the whole selection.
    if keep == "lowest":
        order = np.lexsort((z, ky, kx))
    elif keep == "highest":
        order = np.lexsort((-z, ky, kx))
    elif keep == "first":
        order = np.lexsort((np.arange(n_in), ky, kx))
    elif keep == "mean":
        # For a set that is already ground, "lowest" is the wrong survivor.
        # A delivered ground population can hold two returns at one XY that
        # straddle the surface -- in the NCALM deliveries class 2 sits a
        # median 16 mm below it and class 8 12 mm above -- and taking the
        # minimum silently keeps the lower class on every collision. That is
        # a systematic downward pull on an eighth of the points, not a
        # tie-break. Averaging coincident returns has no such preference.
        order = np.lexsort((ky, kx))
    else:
        raise PreprocessError(
            "keep must be 'lowest', 'highest', 'first' or 'mean'")

    kxs, kys = kx[order], ky[order]
    first = np.ones(n_in, dtype=bool)
    if n_in > 1:
        first[1:] = (kxs[1:] != kxs[:-1]) | (kys[1:] != kys[:-1])

    sel = order[first]
    if keep == "mean":
        starts = np.flatnonzero(first)
        counts = np.diff(np.append(starts, n_in)).astype(float)
        zbar = np.add.reduceat(z[order], starts) / counts
        back = np.argsort(sel)
        xs, ys, zs_ = x[sel][back], y[sel][back], zbar[back]
        dropped = n_in - len(sel)
        return xs, ys, zs_, {
            "points_in": n_in, "points_out": int(len(sel)),
            "dropped": int(dropped),
            "collision_groups": int((counts > 1).sum()),
            "keep": keep, "round_to": round_to}
    sel.sort()                      # back to input order; solvers need no order
    dropped = n_in - len(sel)
    groups = int(np.count_nonzero(~first))      # points that lost a collision
    report = {"points_in": n_in, "points_out": int(len(sel)),
              "dropped": int(dropped),
              "collision_groups": int(n_in - first.sum()) and groups,
              "keep": keep, "round_to": round_to}
    return x[sel], y[sel], z[sel], report


def describe(report):
    """One line for a log, or None when nothing was removed."""
    if not report or not report.get("dropped"):
        return None
    return ("de-duplicated %d of %d returns sharing an XY location (%.2f%%); "
            "kept the %s of each" % (report["dropped"], report["points_in"],
                                     100.0 * report["dropped"] / report["points_in"],
                                     report["keep"]))
