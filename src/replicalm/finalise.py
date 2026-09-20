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


def finalise(dem, grid, wkt, path, radius_m, nodata_in=-9999.0,
             erode_factor=1.0, max_fill_m=60.0, smoothing=2, verbose=True):
    """Fill, trim and write one DEM under the black-means-nodata convention."""
    from . import mask as M

    dem = np.asarray(dem)
    valid = M.coverage_mask(dem, nodata_in)
    filled, did, report = fill_holes(dem, valid, nodata_in,
                                     max_fill_m=max_fill_m, cell_m=grid.cell,
                                     smoothing=smoothing, verbose=verbose)
    covered = valid | did

    cells = erosion_cells(radius_m, grid.cell, erode_factor)
    trimmed = M.erode(covered, cells)
    if verbose:
        print("  trimmed %d cells (%.1f m) from the coverage edge: "
              "%d of %d cells remain (%.2f%%)"
              % (cells, cells * grid.cell, trimmed.sum(), covered.sum(),
                 100 * trimmed.sum() / max(covered.sum(), 1)))
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
                   "cells_before_trim": int(covered.sum()),
                   "cells_after_trim": int(trimmed.sum()),
                   "filled_fraction": (report["filled_cells"] /
                                       max(int(trimmed.sum()), 1)),
                   "bands": 1,
                   "path": str(path)})
    return report
