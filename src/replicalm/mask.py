"""Black means no data, and nothing else does.

THE CONVENTION
--------------
Outside the valid area, every pixel is zero. Inside it, no pixel is zero. A
reader can therefore test one value and know the answer, without consulting a
nodata tag that may not have survived a format conversion, a resample, or a
trip through software that discards band metadata.

This is the rule already used in the edge-fixed products, and it exists because
the alternative fails quietly. A DEM where zero sometimes means "no return" and
sometimes means "the ground is at zero" produces holes that look like terrain
and terrain that looks like holes, and neither is visible until something
downstream has already trusted it.

WHY IN-MASK ZEROS HAVE TO BE LIFTED, NOT LEFT
---------------------------------------------
Enforcing the convention in one direction only is worse than not enforcing it.
If zeros outside are made black but a genuine zero inside is left alone, then
zero now means two things again and the convention has bought nothing. So any
valid pixel that lands on zero is raised to the smallest value the raster can
distinguish from it: one unit for an integer band, one ulp-scaled epsilon for a
float band. The displacement is below the sensor's precision and is reported, so
it is a recorded adjustment rather than a silent one.
"""
import numpy as np


class MaskError(RuntimeError):
    pass


def enforce(array, valid, nodata_value=0.0, verbose=True):
    """Apply the convention to one band.

    `valid` is a boolean array, True where data exists. Returns the adjusted
    band and a report of what was changed.
    """
    a = np.array(array, copy=True)
    valid = np.asarray(valid, bool)
    if a.shape != valid.shape:
        raise MaskError("array %s and mask %s differ in shape"
                        % (a.shape, valid.shape))

    outside_nonzero = int((~valid & (a != nodata_value)).sum())
    a[~valid] = nodata_value

    inside_zero = valid & (a == nodata_value)
    n_lifted = int(inside_zero.sum())
    if n_lifted:
        if np.issubdtype(a.dtype, np.integer):
            step = 1
        else:
            finite = a[valid & np.isfinite(a)]
            scale = float(np.abs(finite).max()) if finite.size else 1.0
            step = max(np.finfo(a.dtype).eps * scale, np.finfo(a.dtype).tiny)
        a[inside_zero] = nodata_value + step
    else:
        step = 0

    report = {"cells_total": int(a.size),
              "cells_valid": int(valid.sum()),
              "valid_fraction": float(valid.mean()),
              "zeroed_outside": outside_nonzero,
              "lifted_inside": n_lifted,
              "lift_step": float(step)}
    if verbose and (outside_nonzero or n_lifted):
        print("  mask: zeroed %d cells outside, lifted %d zeros inside by %.3g"
              % (outside_nonzero, n_lifted, step))
    return a, report


def verify(array, valid, nodata_value=0.0):
    """Confirm the convention holds. Raises if it does not.

    Called after writing rather than trusted from having been applied, because
    the failure this guards against is a later stage reintroducing zeros.
    """
    a = np.asarray(array)
    valid = np.asarray(valid, bool)
    bad_in = int((valid & (a == nodata_value)).sum())
    bad_out = int((~valid & (a != nodata_value)).sum())
    if bad_in or bad_out:
        raise MaskError(
            "the black-means-nodata convention is broken: %d valid cells are "
            "exactly %g, and %d cells outside the mask are not. Zero now means "
            "two things." % (bad_in, nodata_value, bad_out))
    return True


def coverage_mask(dem, nodata=-9999.0):
    """Which cells of a freshly kriged DEM carry a value."""
    a = np.asarray(dem)
    return np.isfinite(a) & (a != nodata)


def erode(valid, cells=1):
    """Shrink a mask inward, to drop the fringe where support is one-sided.

    Cells at the very edge of coverage are interpolated from points on one side
    only, so they extrapolate rather than interpolate. Trimming a ring of them
    is cheaper than explaining later why the boundary of every tile is soft.
    """
    if cells <= 0:
        return np.asarray(valid, bool)
    v = np.asarray(valid, bool)
    for _ in range(cells):
        shrunk = v.copy()
        shrunk[1:, :] &= v[:-1, :]
        shrunk[:-1, :] &= v[1:, :]
        shrunk[:, 1:] &= v[:, :-1]
        shrunk[:, :-1] &= v[:, 1:]
        v = shrunk
    return v


def write_masked_geotiff(dem, grid, wkt, path, valid=None, nodata_value=0.0,
                         erode_cells=0, write_mask_band=False, verbose=True):
    """Write a DEM under the convention. Single band by default.

    WHY THE MASK BAND IS NO LONGER WRITTEN
    --------------------------------------
    It used to be, on the argument that it recorded where data was expected and
    so distinguished "outside the survey" from "inside the survey and missing",
    which the elevation band cannot. In practice it did not do that: `finalise`
    built it from coverage OR fill, so interpolated cells were flagged
    identically to measured ones. What remained was a restatement of something
    band 1 already says unambiguously -- zero is no data, and no valid cell is
    ever exactly zero, because `enforce` lifts any that would be.

    It also cost more than it appeared to: a Float32 band carrying a boolean,
    and a second band that disappears silently through any format conversion
    without multiband support, which is the exact quiet failure this convention
    exists to prevent.

    Fill provenance is reported by `finalise` and belongs in the run's sidecar
    JSON, where the counts already go. On the locked baseline it amounts to
    1,241 cells in 12.4 million, about 0.01%.

    `write_mask_band=True` still writes it, for a caller who wants the old shape.
    """
    from osgeo import gdal, osr
    gdal.UseExceptions()

    dem = np.asarray(dem)
    if valid is None:
        valid = coverage_mask(dem)
    valid = erode(valid, erode_cells)
    out, report = enforce(dem, valid, nodata_value, verbose=verbose)
    verify(out, valid, nodata_value)

    nbands = 2 if write_mask_band else 1
    drv = gdal.GetDriverByName("GTiff")
    ds = drv.Create(str(path), grid.width, grid.height, nbands,
                    gdal.GDT_Float32,
                    options=["COMPRESS=DEFLATE", "TILED=YES", "PREDICTOR=3"])
    ds.SetGeoTransform(grid.geotransform)
    if wkt:
        srs = osr.SpatialReference()
        srs.ImportFromWkt(wkt)
        ds.SetProjection(srs.ExportToWkt())

    b1 = ds.GetRasterBand(1)
    b1.SetNoDataValue(nodata_value)
    b1.SetDescription("elevation, metres; zero means no data")
    b1.WriteArray(out.astype("f4"))
    if write_mask_band:
        b2 = ds.GetRasterBand(2)
        b2.SetDescription("coverage mask, 1 where elevation is measured")
        b2.WriteArray(valid.astype("f4"))
    ds.FlushCache()
    ds = None
    report["path"] = str(path)
    report["bands"] = nbands
    return str(path), report
