"""Turn a kriged surface into a deliverable raster: holes filled, edges trimmed,
black meaning nothing but no data.

THREE DEFECTS THIS FIXES, ALL MEASURED ON South_GLAS_l0s395
-----------------------------------------------------------
interior holes
    13,561 cells in 7 holes sat empty inside the covered area, every one of
    them carrying data in the reference. They are not voids in the survey; they
    are cells where fewer than `min_points` ground returns fell inside the
    search radius. Left as nodata they render as transparent gaps, and a reader
    cannot tell them from the edge of the flight line.

the edge fringe
    Cells within one search radius of the coverage boundary are supported from
    one side only, so they extrapolate. Measured against the reference, median
    error ran 3.2 to 6.9 m inside 20 cells of the edge with over 90% of cells
    past half a metre, then fell to 0.003 m beyond it. That fringe was about a
    third of every large error on the tile.

    `erode_cells: 1` cannot express this. The damage extends one search radius,
    so the erosion depth belongs in metres and is converted to cells here.

the black convention
    Elevation nodata written as -9999 is invisible to anything that lost the
    nodata tag, and zero is worse because it is also a legal elevation. The
    convention is the one the edge-fixed products already use: outside the mask
    every pixel is zero, inside it none is, and the mask travels in band 2.

WHAT FILLING IS AND IS NOT
--------------------------
Interior holes are filled by GDAL's inverse-distance fill, the same operation
the original workflow's kriging performed implicitly when it put values there.
A filled cell is interpolated terrain, not measured terrain, so the mask band
marks the covered area and `fill_report` records how many cells were invented
and how large the largest hole was. A hole wider than `max_fill_m` is left
alone: past some span, filling is drawing.
"""
import numpy as np


class FinaliseError(RuntimeError):
    pass


def hole_mask(valid):
    """Cells that are empty but enclosed by data: the ones worth filling."""
    from scipy import ndimage as ndi
    return ndi.binary_fill_holes(valid) & ~valid


def erosion_cells(radius_m, cell_m, factor=1.0, minimum=1):
    """How deep to trim, given the radius that did the extrapolating.

    The fringe is one search radius wide because that is how far a cell can
    reach for support that exists on one side only. Expressed in cells so the
    caller does not have to know the grid.
    """
    return max(minimum, int(np.ceil(factor * radius_m / cell_m)))


def fill_holes(dem, valid, nodata=-9999.0, max_fill_m=60.0, cell_m=0.5,
               smoothing=2, verbose=True):
    """Fill enclosed gaps by inverse distance, leaving large ones alone.

    Returns the filled array, the mask of cells that were filled, and a report.
    """
    from osgeo import gdal
    from scipy import ndimage as ndi
    gdal.UseExceptions()

    holes = hole_mask(valid)
    report = {"hole_cells": int(holes.sum()), "filled_cells": 0,
              "holes": 0, "skipped_holes": 0, "largest_hole_cells": 0}
    if not holes.any():
        return np.array(dem, copy=True), np.zeros_like(holes), report

    lab, n = ndi.label(holes)
    report["holes"] = int(n)
    sizes = ndi.sum(holes, lab, range(1, n + 1))
    report["largest_hole_cells"] = int(sizes.max()) if n else 0

    # a hole wider than the limit is left as nodata: filling it would be an
    # invention spanning further than the surrounding data can justify
    max_cells = (max_fill_m / cell_m) ** 2
    keep = np.zeros(n + 1, bool)
    for i, s in enumerate(sizes, 1):
        if s <= max_cells:
            keep[i] = True
        else:
            report["skipped_holes"] += 1
    fillable = keep[lab]

    work = np.where(valid, dem, np.nan).astype("f8")
    src = np.where(np.isfinite(work), work, 0.0)
    drv = gdal.GetDriverByName("MEM")
    ds = drv.Create("", dem.shape[1], dem.shape[0], 1, gdal.GDT_Float64)
    band = ds.GetRasterBand(1)
    band.WriteArray(src)
    mem_mask = drv.Create("", dem.shape[1], dem.shape[0], 1, gdal.GDT_Byte)
    mem_mask.GetRasterBand(1).WriteArray(valid.astype("u1"))
    gdal.FillNodata(band, mem_mask.GetRasterBand(1),
                    maxSearchDist=float(max_fill_m / cell_m),
                    smoothingIterations=int(smoothing))
    got = band.ReadAsArray()

    out = np.array(dem, copy=True)
    did = fillable & np.isfinite(got) & (~valid)
    out[did] = got[did]
    report["filled_cells"] = int(did.sum())
    if verbose:
        print("  filled %d of %d enclosed cells in %d holes (largest %d cells)%s"
              % (report["filled_cells"], report["hole_cells"], report["holes"],
                 report["largest_hole_cells"],
                 "" if not report["skipped_holes"]
                 else ", %d holes too large and left open"
                      % report["skipped_holes"]))
    return out, did, report


def infer_transect_axis(grid, name=None):
    """'x' strips left and right, 'y' strips top and bottom.

    G-LiHT names carry the flight direction: a GLAS transect runs roughly
    north-south, so its stitching seams are on the left and right. The
    fallback is the aspect ratio, which says the same thing for an unnamed
    raster.
    """
    if name:
        u = str(name).upper()
        if "_GLAS_" in u:
            return "x"
        if "_NFI_" in u:
            return "y"
    return "x" if grid.height >= grid.width else "y"


def transect_edge_mask(valid, cell_m, strip_m, axis="x",
                       protect_lo=False, protect_hi=False):
    """Drop `strip_m` inward from the outermost data in each row or column.

    This is the edge treatment the reference products use, from the G-LiHT
    reprocessing work (`clean_las_edges.py`, `mask_dem_edges.py`). It exists
    for a specific defect: where adjacent flight passes overlap and are merged,
    the ground classification at the overlap differs from the interior and
    leaves a seam about 21 m wide down both sides of every transect.

    Two properties matter and neither is shared by a morphological erosion:

    it is directional
        only the across-track sides are stripped. Eroding a ring takes the
        along-track ends too, where there is no seam and nothing wrong.

    it cannot open an interior hole
        only the first and last runs of each line are touched, so a void
        enclosed by data is never reached. An erosion grows every interior
        void by its depth, which is how a canopy gap becomes a black island.
    """
    v = np.asarray(valid, bool)
    out = v.copy()
    strip = int(round(float(strip_m) / float(cell_m)))
    if strip <= 0:
        return out
    if axis == "y":
        v = v.T
        out = out.T
    # `protect_lo` and `protect_hi` say that the lattice's own border is a
    # clip rather than the end of the survey, because the point cloud carries
    # on past it. Stripping there would take good data along a perfectly
    # straight line for a seam that is not present -- on Centro GLAS l0s454,
    # coverage reached the west border on 36.6 percent of rows and the strip
    # cut a 260 m vertical face out of it.
    #
    # They default to False, and that matters: when the grid is derived from
    # the cloud, the data touches the frame on every side by construction, and
    # protecting it unconditionally would disable the strip entirely. Only a
    # caller that knows points exist beyond the border should set them.
    last = v.shape[1] - 1
    lines = np.flatnonzero(v.any(axis=1))
    for r in lines:
        cols = np.flatnonzero(v[r])
        lo, hi = int(cols[0]), int(cols[-1])
        l_end = lo
        if not (protect_lo and lo == 0):
            l_end = min(lo + strip, hi + 1)
            out[r, lo:l_end] = False
        if not (protect_hi and hi == last):
            r_start = max(hi - strip + 1, l_end)
            out[r, r_start:hi + 1] = False
    return out.T if axis == "y" else out


def support_mask(x, y, grid, radius_m):
    """Cells within `radius_m` of an accepted ground return.

    Erosion trims inward from whatever the interpolator filled, which is the
    wrong reference when the interpolator extrapolates: a global solver run in
    blocks will fringe out to each block's own extent, and no erosion depth
    expressed in the coverage mask reaches that. Distance to the data does.

    Interior voids are filled before returning. Ground returns are sparse under
    canopy, so a plain distance threshold punches holes through the middle of
    perfectly good terrain; those holes are not edges and this mask has no
    business removing them.
    """
    from scipy import ndimage as ndi

    col = np.floor((np.asarray(x) - grid.origin_x) / grid.cell).astype(np.int64)
    row = np.floor((grid.origin_y - np.asarray(y)) / grid.cell).astype(np.int64)
    ok = ((col >= 0) & (col < grid.width) & (row >= 0) & (row < grid.height))
    occ = np.zeros((grid.height, grid.width), bool)
    if ok.any():
        occ[row[ok], col[ok]] = True
    if not occ.any():
        return np.ones_like(occ)
    dist = ndi.distance_transform_edt(~occ, sampling=grid.cell)
    near = ndi.binary_closing(dist <= float(radius_m),
                              structure=np.ones((5, 5)))
    return ndi.binary_fill_holes(near)


def _distance_to_data(x, y, grid):
    """Metres from each cell to the nearest return. Shared by the two below."""
    from scipy import ndimage as ndi

    col = np.floor((np.asarray(x) - grid.origin_x) / grid.cell).astype(np.int64)
    row = np.floor((grid.origin_y - np.asarray(y)) / grid.cell).astype(np.int64)
    ok = ((col >= 0) & (col < grid.width) & (row >= 0) & (row < grid.height))
    occ = np.zeros((grid.height, grid.width), bool)
    if not ok.any():
        return None
    occ[row[ok], col[ok]] = True
    return ndi.distance_transform_edt(~occ, sampling=grid.cell)


def auto_support_radius(x, y, grid, lo_m=1.0, hi_m=60.0, step_m=1.0,
                        verbose=False):
    """How far from a return a cell may sit and still be measured, per tile.

    A fixed radius cannot be right for every survey: it has to bridge the
    gaps between ground returns, and under closed canopy those are far wider
    than the mean spacing. On the calibration tile the mean ground spacing is
    0.43 m and the radius that reproduced the reference footprint was 20 m,
    forty-six times larger, because what it is really bridging is canopy
    occlusion rather than sampling interval.

    MEASURED AND REJECTED. The knee of the coverage-versus-radius curve was
    the obvious automatic rule and it does not reproduce the reference: on
    Centro GLAS l0s454 it chose 7.0 m, giving a footprint 7.4 percent short
    of the archive's and broken into four pieces. It is also sensitive to the
    range searched, which makes it a parameter wearing a disguise.

    Use `search_radius_m` instead. The reference footprint is not a property
    of the gap structure at all -- it is the interpolator's own reach. Surfer
    ran kriging at a 20 m search radius with a minimum of one point, so a
    cell received a value exactly when a return lay within 20 m of it, and
    the archive's maximum distance to a ground return is 19.73 m. A cell
    beyond the search radius was never interpolated; it was extrapolated.

    Kept because the negative result is worth having in one place.
    """
    dist = _distance_to_data(x, y, grid)
    if dist is None:
        return float(hi_m)
    radii = np.arange(lo_m, hi_m + step_m, step_m)
    area = np.array([float((dist <= r).sum()) for r in radii])
    if area[-1] <= area[0]:
        return float(lo_m)
    rn = (radii - radii[0]) / max(radii[-1] - radii[0], 1e-9)
    an = (area - area[0]) / max(area[-1] - area[0], 1e-9)
    # distance from each point to the straight line between the two ends
    knee = int(np.argmax(np.abs(an - rn)))
    r = float(radii[knee])
    if verbose:
        print("  support radius %.1f m chosen at the knee of the coverage "
              "curve (%.1f to %.1f m tested)" % (r, lo_m, hi_m))
    return r


def data_footprint(x, y, grid, radius_m, axis=None, strip_m=None,
                   ends_m=None):
    """The active area: near a return, closed, optionally stripped.

    Returns (mask, voids, report). `voids` marks cells the closing added --
    inside the survey but carrying no measurement -- which the caller fills
    and records rather than leaving as holes.

    Closing runs per connected area, so a detached island becomes its own
    closed polygon rather than being bridged to the main body.
    """
    from scipy import ndimage as ndi

    dist = _distance_to_data(x, y, grid)
    if dist is None:
        m = np.ones((grid.height, grid.width), bool)
        return m, np.zeros_like(m), {"radius_m": radius_m, "areas": 0}
    near = dist <= float(radius_m)
    lab, n = ndi.label(near)
    closed = np.zeros_like(near)
    for k in range(1, n + 1):
        closed |= ndi.binary_fill_holes(lab == k)
    voids = closed & ~near
    # Does the cloud carry on past the lattice? If it does, that border is a
    # clip and not a survey edge; if the grid was sized to the cloud, it is
    # the edge and must be stripped like any other.
    w = float(grid.origin_x)
    e = w + grid.width * grid.cell
    nth = float(grid.origin_y)
    sth = nth - grid.height * grid.cell
    past = {"w": float(np.min(x)) < w - 0.5 * grid.cell,
            "e": float(np.max(x)) > e + 0.5 * grid.cell,
            "n": float(np.max(y)) > nth + 0.5 * grid.cell,
            "s": float(np.min(y)) < sth - 0.5 * grid.cell}
    if strip_m:
        ax = axis or infer_transect_axis(grid)
        lo_key, hi_key = ("w", "e") if ax == "x" else ("n", "s")
        closed = transect_edge_mask(closed, grid.cell, strip_m, ax,
                                    protect_lo=past[lo_key],
                                    protect_hi=past[hi_key])
        # The ends as well as the sides. The sides carry the flight-pass
        # stitching seam; the ends carry the end-of-pass taper, where the
        # scan narrows and the solver has support on one side only. Two
        # different defects, and the reference treats both: measured on
        # Centro GLAS l0s454 the archive is inset exactly 42 rows, 21.0 m,
        # at the north and south ends as well as east and west.
        if ends_m:
            other = "y" if ax == "x" else "x"
            lo2, hi2 = ("n", "s") if other == "y" else ("w", "e")
            closed = transect_edge_mask(closed, grid.cell, ends_m, other,
                                        protect_lo=past[lo2],
                                        protect_hi=past[hi2])
        voids &= closed
    _l2, n2 = ndi.label(closed)
    return closed, voids, {"radius_m": float(radius_m),
                           "strip_m": float(strip_m or 0.0),
                           "ends_m": float(ends_m or 0.0),
                           "cloud_past_border": past,
                           "areas_before_strip": int(n),
                           "areas": int(n2),
                           "cells": int(closed.sum()),
                           "void_cells": int(voids.sum())}


def fill_voids(dem, have_value, want):
    """Give every wanted cell a continuous value. Returns (dem, filled).

    A void inside the survey must not be written as a sentinel elevation. A
    cell set to 1 sits a hundred metres below its neighbours, so every
    derived layer sees a cliff and the artifact spreads about five metres
    past the hole in each direction. The surface is kept continuous here and
    the void is marked in the image and recorded in the vector sidecar
    instead, which is where "no measurement" belongs.
    """
    from scipy import ndimage as ndi

    out = np.array(dem, dtype=float, copy=True)
    gap = want & ~have_value
    if not gap.any():
        return out, gap
    idx = ndi.distance_transform_edt(~have_value, return_distances=False,
                                     return_indices=True)
    out[gap] = out[tuple(i[gap] for i in idx)]
    # soften the nearest-neighbour patches so the fill carries no facets
    src = np.where(have_value | gap, out, 0.0)
    wgt = (have_value | gap).astype(float)
    num = ndi.gaussian_filter(src, 1.0)
    den = ndi.gaussian_filter(wgt, 1.0)
    smooth = np.where(den > 1e-6, num / np.maximum(den, 1e-6), out)
    out[gap] = smooth[gap]
    return out, gap


def return_footprint(x, y, grid, axis=None):
    """The span between the first and last ground return, row by row.

    This replaces `support_mask`, which thresholded a distance transform at a
    few metres and so drew the coverage boundary as a union of discs. Where
    returns thin toward a swath margin that union is literally a row of
    overlapping circles with outliers as detached islands, which is the
    scalloped boundary the production renders showed.

    The reference implementation -- mask_dem_edges.py in the G1 reprocessing
    work -- takes the first and last valid column of each row and strips a
    fixed distance inward from each. It can do that because Surfer bounded its
    own output; this interpolator has no distance limit and fills the whole
    lattice, so the footprint has to come from the returns instead. Taking one
    span per row reproduces the reference's geometry: a single contiguous run
    per row, so no islands and no scallops are representable at all.

    `axis` follows the transect convention -- "x" spans each row for a GLAS
    strip, "y" spans each column for an NFI strip.
    """
    col = np.floor((np.asarray(x) - grid.origin_x) / grid.cell).astype(np.int64)
    row = np.floor((grid.origin_y - np.asarray(y)) / grid.cell).astype(np.int64)
    ok = ((col >= 0) & (col < grid.width) & (row >= 0) & (row < grid.height))
    occ = np.zeros((grid.height, grid.width), bool)
    if not ok.any():
        return np.ones_like(occ)
    occ[row[ok], col[ok]] = True
    ax = axis or infer_transect_axis(grid)
    span = np.zeros_like(occ)
    if ax == "x":
        for r in np.flatnonzero(occ.any(axis=1)):
            c = np.flatnonzero(occ[r])
            span[r, c[0]:c[-1] + 1] = True
    else:
        for c in np.flatnonzero(occ.any(axis=0)):
            r = np.flatnonzero(occ[:, c])
            span[r[0]:r[-1] + 1, c] = True
    return span


def _affine(grid):
    """rasterio Affine from a Grid. `geotransform` is GDAL's ordering."""
    from rasterio.transform import Affine

    gt = grid.geotransform
    return Affine(gt[1], gt[2], gt[0], gt[4], gt[5], gt[3])


def write_footprint_vectors(mask, grid, wkt, path, voids=None):
    """The active area as closed polygons, one per piece, exterior ring only.

    GeoJSON because it needs no driver beyond what rasterio already carries;
    fiona is not in the shipped environment, so a GeoPackage cannot be
    written without adding a dependency for a sidecar.

    Interior rings are dropped deliberately: a void inside the survey is
    recorded as its own feature with kind "no_measurement" rather than as a
    hole in the active polygon, so a consumer reading only the outline gets
    the survey's extent and a consumer reading the attributes gets the
    caveat.
    """
    import json

    from rasterio import features

    feats = []

    def rings(m, kind):
        for geom, v in features.shapes(np.asarray(m, "uint8"), mask=m,
                                       transform=_affine(grid),
                                       connectivity=4):
            if v != 1:
                continue
            ring = geom["coordinates"][0]
            xs = [q[0] for q in ring]
            ys = [q[1] for q in ring]
            area = abs(sum(xs[i] * ys[i + 1] - xs[i + 1] * ys[i]
                           for i in range(len(xs) - 1))) / 2.0
            feats.append({
                "type": "Feature",
                "properties": {"poly_id": len(feats) + 1, "kind": kind,
                               "area_m2": round(area, 2),
                               "area_ha": round(area / 10000.0, 4)},
                "geometry": {"type": "Polygon", "coordinates": [ring]}})

    rings(np.asarray(mask, bool), "active")
    if voids is not None and np.asarray(voids, bool).any():
        rings(np.asarray(voids, bool), "no_measurement")
    with open(path, "w") as fh:
        json.dump({"type": "FeatureCollection", "features": feats}, fh)
    return {"path": path, "features": len(feats),
            "active": sum(1 for f in feats
                          if f["properties"]["kind"] == "active"),
            "voids": sum(1 for f in feats
                         if f["properties"]["kind"] == "no_measurement")}


def finalise(dem, grid, wkt, path, radius_m, nodata_in=-9999.0,
             erode_factor=1.0, max_fill_m=60.0, smoothing=2, support=None,
             strip_m=None, axis=None, verbose=True):
    """Fill, trim and write one DEM under the black-means-nodata convention.

    `support` is an optional boolean mask of cells close enough to a ground
    return to be written; see `support_mask`. It is applied alongside the
    erosion rather than instead of it, so neither relaxes the other.
    """
    from . import mask as M

    dem = np.asarray(dem)
    valid = M.coverage_mask(dem, nodata_in)
    filled, did, report = fill_holes(dem, valid, nodata_in,
                                     max_fill_m=max_fill_m, cell_m=grid.cell,
                                     smoothing=smoothing, verbose=verbose)
    covered = valid | did

    from scipy import ndimage as ndi

    if strip_m:
        # the reference convention: across-track only, interior untouched
        ax = axis or infer_transect_axis(grid)
        trimmed = transect_edge_mask(covered, grid.cell, strip_m, ax)
        cells = int(round(float(strip_m) / grid.cell))
    elif radius_m and radius_m > 0:
        cells = erosion_cells(radius_m, grid.cell, erode_factor)
        trimmed = M.erode(covered, cells)
    else:
        # Edge fixing off: keep every cell the interpolator produced. The
        # fringe where support is one-sided is left in deliberately, so the
        # engine can be judged on all of its active-data area rather than on
        # whatever a trim chose to leave behind.
        cells = 0
        trimmed = covered
    n_eroded = int(trimmed.sum())
    if support is not None:
        trimmed = trimmed & np.asarray(support, bool)

    # Nothing enclosed by data is ever dropped. A hole left here would render
    # as black inside the survey, and black means no data; an interior void is
    # not the edge of anything.
    holes_restored = int((ndi.binary_fill_holes(trimmed) & covered
                          & ~trimmed).sum())
    trimmed = ndi.binary_fill_holes(trimmed) & covered
    if verbose and holes_restored:
        print("  restored %d interior cells the trim had enclosed"
              % holes_restored)
    if verbose:
        print("  trimmed %d cells (%.1f m) from the coverage edge: "
              "%d of %d cells remain (%.2f%%)"
              % (cells, cells * grid.cell, trimmed.sum(), covered.sum(),
                 100 * trimmed.sum() / max(covered.sum(), 1)))
        if support is not None:
            print("  support mask removed a further %d cells beyond the "
                  "returns" % (n_eroded - int(trimmed.sum())))
    if not trimmed.any():
        raise FinaliseError("erosion of %d cells removed the whole raster"
                            % cells)

    # Single band. The mask band restated what band 1 already says -- zero is
    # no data, and no valid cell is exactly zero -- while failing to record the
    # one thing it could usefully have carried, which cells were interpolated.
    # That provenance is in the report below and belongs in the run's sidecar
    # JSON, not in a band every consumer has to know to ignore.
    M.write_masked_geotiff(filled, grid, wkt, path, valid=trimmed,
                           nodata_value=0.0, erode_cells=0,
                           write_mask_band=False, verbose=verbose)
    report.update({"erode_cells": cells, "erode_m": cells * grid.cell,
                   "support_mask": support is not None,
                   "edge_strip_m": strip_m,
                   "interior_cells_restored": holes_restored,
                   "cells_beyond_support": (None if support is None
                                            else n_eroded
                                            - int(trimmed.sum())),
                   "cells_before_trim": int(covered.sum()),
                   "cells_after_trim": int(trimmed.sum()),
                   "filled_fraction": (report["filled_cells"] /
                                       max(int(trimmed.sum()), 1)),
                   "bands": 1,
                   "path": str(path)})
    return report
