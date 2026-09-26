r"""Kernel-safe tile ingestion: real neighbour data first, reflection second.

WHY
---
Kernel-based layers lose their border. Multi-directional hillshade runs a 3x3
gradient, so a 2000 x 2000 DEM returns 1998 x 1998, while slope, sky-view factor
and positive openness come back full size. Reconciling that afterwards leaves a
choice between a nodata ring, a shifted geotransform, or resampling -- and
resampling a layer that is about to be differenced and multiplied against three
others smears it by a sub-pixel amount for no gain.

Better to not lose the border in the first place. Read the tile with a halo, run
the kernel over the enlarged array, then strip the halo. The output is exactly
the target grid and the geotransform never moves.

WHERE THE HALO COMES FROM, IN ORDER
-----------------------------------
1. The neighbouring tile, when one exists. The kernel then computes on real,
   continuous terrain across the seam, and cells at the tile edge are as valid as
   cells in the middle. G-LiHT transects overlap at their joins, so a neighbour is
   usually there.
2. Reflection, at the survey perimeter where nothing adjoins. Symmetric reflection
   continues the local gradient instead of inventing a wall -- zero-padding or
   nodata would give the kernel a cliff at the boundary and stamp a bright or dark
   rim along the edge of every outer tile.

Neighbour data is only used when the two rasters share a pixel grid to within a
hundredth of a cell. Tiles at the same resolution but on unaligned origins would
need resampling to combine, which is the thing this exists to avoid, so those
edges fall back to reflection and say so.

GEOTRANSFORMS
-------------
Never modified. The halo is a read-time device; everything returned to the caller
is on the tile's own grid with the tile's own transform.
"""
import math
import os
import warnings

import numpy as np


class DimensionMismatch(ValueError):
    """A layer came back too far from the target grid to correct safely."""


def build_index(paths):
    """Path -> grid description, from headers only. Cheap enough for a corpus."""
    import rasterio
    index = []
    for p in paths:
        try:
            with rasterio.open(p) as r:
                index.append({"path": p, "transform": r.transform,
                              "width": r.width, "height": r.height,
                              "bounds": r.bounds, "nodata": r.nodata,
                              "cell": abs(r.transform.a)})
        except Exception as exc:
            warnings.warn("cannot index %s: %s" % (os.path.basename(p), exc))
    return index


def _aligned(a, b, tol=0.01):
    """Do two transforms share a pixel grid, to within a hundredth of a cell?"""
    if abs(abs(a.a) - abs(b.a)) > tol * abs(a.a):
        return False
    if abs(abs(a.e) - abs(b.e)) > tol * abs(a.e):
        return False
    dx = (b.c - a.c) / a.a
    dy = (b.f - a.f) / a.e
    return (abs(dx - round(dx)) < tol) and (abs(dy - round(dy)) < tol)


def read_with_halo(path, index=None, halo=1, clean_nodata=True):
    """A DEM plus `halo` cells of context on every side.

    Returns (array, cell, transform, report). The array is
    (height + 2*halo, width + 2*halo); `transform` is the tile's own, unchanged.
    """
    import rasterio
    with rasterio.open(path) as r:
        dem = r.read(1).astype(np.float32)
        nd, transform = r.nodata, r.transform
        h, w, cell = r.height, r.width, abs(r.transform.a)
        bounds = r.bounds
    if clean_nodata:
        if nd is not None:
            dem[dem == nd] = np.nan
        dem[dem == 0] = np.nan          # the archive's other nodata convention

    big = np.full((h + 2 * halo, w + 2 * halo), np.nan, dtype=np.float32)
    big[halo:halo + h, halo:halo + w] = dem
    report = {"halo": halo, "from_neighbours": 0, "reflected": 0,
              "neighbours_used": [], "unaligned_skipped": []}

    # 1. real neighbour data, where a tile adjoins and shares the grid
    if index:
        want = (bounds.left - halo * cell, bounds.bottom - halo * cell,
                bounds.right + halo * cell, bounds.top + halo * cell)
        for nb in index:
            if os.path.abspath(nb["path"]) == os.path.abspath(path):
                continue
            b = nb["bounds"]
            if (b.left >= want[2] or b.right <= want[0]
                    or b.bottom >= want[3] or b.top <= want[1]):
                continue
            if not _aligned(transform, nb["transform"]):
                report["unaligned_skipped"].append(os.path.basename(nb["path"]))
                continue
            # where the neighbour sits in the enlarged grid
            col0 = int(round((b.left - (bounds.left - halo * cell)) / cell))
            row0 = int(round(((bounds.top + halo * cell) - b.top) / cell))
            r0, c0 = max(0, row0), max(0, col0)
            r1 = min(big.shape[0], row0 + nb["height"])
            c1 = min(big.shape[1], col0 + nb["width"])
            if r1 <= r0 or c1 <= c0:
                continue
            with rasterio.open(nb["path"]) as rr:
                win = rasterio.windows.Window(c0 - col0, r0 - row0,
                                              c1 - c0, r1 - r0)
                patch = rr.read(1, window=win).astype(np.float32)
                pnd = rr.nodata
            if clean_nodata:
                if pnd is not None:
                    patch[patch == pnd] = np.nan
                patch[patch == 0] = np.nan
            target = big[r0:r1, c0:c1]
            fill = np.isnan(target) & np.isfinite(patch)
            if fill.any():
                target[fill] = patch[fill]
                report["from_neighbours"] += int(fill.sum())
                report["neighbours_used"].append(os.path.basename(nb["path"]))

    # 2. reflection for whatever the neighbours did not supply
    if halo:
        ring = np.zeros_like(big, dtype=bool)
        ring[:halo, :] = ring[-halo:, :] = True
        ring[:, :halo] = ring[:, -halo:] = True
        todo = ring & ~np.isfinite(big)
        if todo.any():
            # reflect the interior outward; 'symmetric' mirrors the edge cell
            # itself, which continues the local gradient rather than inventing
            # a wall at the boundary
            filled = np.pad(dem, halo, mode="symmetric")
            take = todo & np.isfinite(filled)
            big[take] = filled[take]
            report["reflected"] = int(take.sum())

    return big, cell, transform, report


def strip_halo(arr, halo):
    """Back to the tile's own grid after the kernel has run."""
    if not halo:
        return arr
    return arr[halo:-halo, halo:-halo]


def conform(arr, target_h, target_w, name="layer", tol=2):
    """Tiered dimension validator against the reference grid.

    Within `tol` cells the difference is corrected deterministically -- symmetric
    padding when under-sized, an inner crop when over-sized -- and a warning is
    emitted, because a silent correction is how a one-cell georeferencing offset
    survives to publication. Beyond `tol` it is fatal: a layer that far from the
    target grid is not a border effect, and padding it would be guesswork
    dressed as arithmetic.
    """
    h, w = arr.shape
    dh, dw = target_h - h, target_w - w
    if dh == 0 and dw == 0:
        return arr
    if abs(dh) > tol or abs(dw) > tol:
        raise DimensionMismatch(
            "%s is %dx%d against a target of %dx%d, off by (%d, %d) which "
            "exceeds the %d-cell tolerance. This is not a kernel border; "
            "refusing to pad or crop it."
            % (name, h, w, target_h, target_w, dh, dw, tol))
    warnings.warn("%s is %dx%d against a target of %dx%d; correcting by (%d, %d)"
                  % (name, h, w, target_h, target_w, dh, dw), RuntimeWarning)
    if dh > 0 or dw > 0:
        top, left = max(dh, 0) // 2, max(dw, 0) // 2
        out = np.full((max(h, target_h), max(w, target_w)), np.nan, np.float32)
        out[top:top + h, left:left + w] = arr
        arr, h, w = out, out.shape[0], out.shape[1]
    if h > target_h or w > target_w:
        top, left = (h - target_h) // 2, (w - target_w) // 2
        arr = arr[top:top + target_h, left:left + target_w]
    return arr
