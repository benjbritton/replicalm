"""Is the difference between two DEMs noise, or is it smoothing?

THE QUESTION
------------
Every interpolation setting tried so far reproduces the reference's elevations to
a few centimetres while carrying 1.6 to 2.8 times its terrain complexity. Two
explanations fit that equally well from the magnitude alone:

  we are adding noise      low returns taken as ground, rendered as pits
  they removed terrain     the original smoothed, and we are keeping detail

The magnitude cannot separate them. The spatial structure can. Noise is
uncorrelated: scattered single cells, no relationship between neighbours.
Smoothing is correlated over the width of whatever kernel did it.

WHAT THE SHAPE SAYS
-------------------
The residual's variogram answers it directly, and measures the answer:

  pure nugget            semivariance flat from the first bin, nugget/sill near
                         1, no range. The residual is noise.
  rises to a sill        the range is the correlation length of the difference,
                         which for smoothing is the kernel width. A range of
                         fifteen metres is a statement about the original
                         workflow, not an inference about it.

COMPUTING IT ON A LATTICE
-------------------------
These are regular grids, so pairwise lags need no distance computation. By the
Wiener-Khinchin theorem the whole lag structure comes from three Fourier
transforms in O(N log N), against O(N^2) for explicit pairs -- 640,000 cells
would otherwise be 2x10^11 pairs.

The complication is nodata. FFT autocovariance assumes a complete field, and a
third of these cells are empty. Expanding the semivariance definition solves it,
because every term is then a correlation:

    2 gamma(h) N(h) = sum z_i^2 + sum z_j^2 - 2 sum z_i z_j

with z set to zero outside the mask and the mask carried alongside. N(h), the
count of valid pairs at each lag, is the mask correlated with itself. Nothing is
infilled and nothing is approximated; empty cells simply contribute nothing to
any term, which is what they should do.
"""
import numpy as np


class ResidualError(RuntimeError):
    pass


def _corr_fft(a, b, shape):
    """Cross-correlation of two padded arrays, by FFT."""
    fa = np.fft.rfft2(a, s=shape)
    fb = np.fft.rfft2(b, s=shape)
    return np.fft.irfft2(fa * np.conj(fb), s=shape)


def gridded_variogram(field, mask, cell=1.0, max_lag_cells=None, n_bins=24,
                      min_pairs=100):
    """Empirical variogram of a masked raster, exactly, via FFT.

    Returns lag distances in metres, semivariance, and the pair count behind
    each bin, so a bin resting on twenty pairs can be recognised as such.
    """
    z = np.where(mask, field, 0.0).astype(np.float64)
    m = mask.astype(np.float64)
    h, w = z.shape
    if max_lag_cells is None:
        max_lag_cells = min(h, w) // 3
    shape = (h + max_lag_cells + 1, w + max_lag_cells + 1)

    # three correlations give every term of the semivariance
    zz = _corr_fft(z, z, shape)            # sum z_i z_j
    z2m = _corr_fft(z * z, m, shape)       # sum z_i^2 over valid partners
    mz2 = _corr_fft(m, z * z, shape)       # sum z_j^2
    nn = _corr_fft(m, m, shape)            # pair count

    # lag vectors: rows 0..max, columns -max..max, taken from the wrapped result
    ky = np.arange(0, max_lag_cells + 1)
    kx = np.concatenate([np.arange(0, max_lag_cells + 1),
                         np.arange(-max_lag_cells, 0)])
    sub = np.ix_(ky % shape[0], kx % shape[1])
    num = z2m[sub] + mz2[sub] - 2.0 * zz[sub]
    cnt = nn[sub]

    dy = ky[:, None] * np.ones_like(kx)[None, :]
    dx = np.ones_like(ky)[:, None] * kx[None, :]
    dist = np.hypot(dy, dx) * cell

    good = (cnt > 0.5) & (dist > 0)
    d = dist[good]
    g = np.maximum(num[good], 0.0) / (2.0 * cnt[good])
    n = cnt[good]

    edges = np.linspace(0, d.max(), n_bins + 1)
    which = np.clip(np.digitize(d, edges) - 1, 0, n_bins - 1)
    lags, semis, counts = [], [], []
    for i in range(n_bins):
        sel = which == i
        if not sel.any():
            continue
        tot = n[sel].sum()
        if tot < min_pairs:
            continue
        # weight by pair count: a lag bin is the pooled estimate, not a mean of
        # means, since bins hold wildly different numbers of pairs
        lags.append(float((d[sel] * n[sel]).sum() / tot))
        semis.append(float((g[sel] * n[sel]).sum() / tot))
        counts.append(int(tot))
    if len(lags) < 3:
        raise ResidualError("only %d usable lag bins" % len(lags))
    return np.array(lags), np.array(semis), np.array(counts)


def morans_i(field, mask):
    """Moran's I at lag one, using the four rook neighbours.

    One number: near zero for an uncorrelated field, near one for a structured
    one. The variogram says more, but this says it in a single figure.
    """
    z = np.where(mask, field, np.nan)
    v = z[mask]
    if v.size < 10:
        raise ResidualError("too few valid cells")
    mu = v.mean()
    d = np.where(mask, z - mu, 0.0)
    num = w = 0.0
    for a, b in (((slice(1, None), slice(None)), (slice(None, -1), slice(None))),
                 ((slice(None), slice(1, None)), (slice(None), slice(None, -1)))):
        pair = mask[a] & mask[b]
        num += float((d[a][pair] * d[b][pair]).sum()) * 2
        w += float(pair.sum()) * 2
    denom = float((d[mask] ** 2).sum())
    if denom == 0 or w == 0:
        return float("nan")
    return (v.size / w) * (num / denom)


def characterise(lags, semis, counts):
    """Nugget, sill and range, read off the empirical curve.

    Read rather than fitted: fitting a model presumes the shape that is in
    question. The nugget is the intercept of a line through the first three
    bins, the sill is the plateau of the far bins, and the range is where the
    curve first reaches 95% of it.
    """
    k = min(3, len(lags))
    A = np.c_[np.ones(k), lags[:k]]
    nugget = float(max(np.linalg.lstsq(A, semis[:k], rcond=None)[0][0], 0.0))
    tail = max(len(semis) // 3, 1)
    sill = float(np.mean(semis[-tail:]))
    ratio = nugget / sill if sill > 0 else float("nan")
    rng, resolved = None, False
    if sill > 0:
        hit = np.where(semis >= 0.95 * sill)[0]
        if hit.size:
            rng = float(lags[hit[0]])
            # A range equal to the first bin is not a range. It means the curve
            # was already at its sill by the shortest lag measured, so whatever
            # structure exists lies below the resolution of these bins. Saying
            # "5.4 m" there reports the instrument, not the terrain.
            resolved = hit[0] > 0
    if ratio >= 0.9:
        verdict = "pure nugget: the residual is spatially uncorrelated noise"
    elif not resolved:
        verdict = ("sill reached at the first lag (%.1f m): any structure is "
                   "finer than these bins resolve. Re-run with a shorter "
                   "max_lag to see below it." % (rng if rng else float("nan")))
    elif ratio <= 0.35:
        verdict = ("structured: correlated over about %.1f m, the scale of "
                   "whatever produced it" % rng)
    else:
        verdict = ("mixed: a nugget over structure with a range of about "
                   "%.1f m" % rng)
    return {"nugget": nugget, "sill": sill, "nugget_to_sill": ratio,
            "range_m": rng, "range_resolved": bool(resolved),
            "lag_resolution_m": float(lags[0]), "verdict": verdict}


def compare(candidate_path, reference_path, max_lag_m=60.0):
    """Full diagnosis of the difference between two aligned DEMs."""
    from osgeo import gdal
    gdal.UseExceptions()

    def read(p):
        d = gdal.Open(str(p))
        b = d.GetRasterBand(1)
        a = b.ReadAsArray().astype("f8")
        nod = b.GetNoDataValue()
        m = np.isfinite(a)
        if nod is not None:
            m &= a != nod
        m &= a != 0.0
        return a, m, d.GetGeoTransform(), (d.RasterXSize, d.RasterYSize)

    ca, cm, cgt, csz = read(candidate_path)
    ra, rm, rgt, rsz = read(reference_path)
    if csz != rsz or abs(cgt[0] - rgt[0]) > 1e-6 or abs(cgt[3] - rgt[3]) > 1e-6:
        raise ResidualError("rasters are not on the same grid; align them first")

    cell = abs(cgt[1])
    both = cm & rm
    resid = np.where(both, ca - ra, 0.0)
    if both.sum() < 1000:
        raise ResidualError("only %d cells in common" % int(both.sum()))

    lags, semis, counts = gridded_variogram(
        resid, both, cell=cell, max_lag_cells=int(max_lag_m / cell))
    out = characterise(lags, semis, counts)
    out.update({"morans_i": morans_i(resid, both),
                "residual_mean_m": float(resid[both].mean()),
                "residual_sd_m": float(resid[both].std()),
                "cells": int(both.sum()),
                "lags_m": lags.tolist(), "semivariance": semis.tolist(),
                "pair_counts": counts.tolist()})
    return out
