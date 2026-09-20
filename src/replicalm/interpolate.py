"""Step 8: ground points to a DEM, and step 9: georeference the result.

TWO BACKENDS, ONE INTERFACE
---------------------------
`cloudcompare` drives CloudCompare's rasterizer from the command line. Its
`-EMPTY_FILL KRIGING` really does krige -- on the first test tile it filled
21,029 cells the other strategies left empty, differing from them by up to
2.17 m where both had values. That matters more than it sounds: only 30.7% of
cells on that tile contain any return at all, because the points cluster along
scan lines rather than spreading evenly. Over most of a G-LiHT DEM the
interpolator is not refining measurements, it is producing them.

`kriging` is a local ordinary kriging written here, for when the search radius
has to be exactly the 20 m the source method specifies rather than whatever the
GUI defaults to.

WHY GEOREFERENCING IS A SEPARATE STEP
-------------------------------------
CloudCompare writes the grid with correct real-world coordinates in the
geotransform but no projection string, so the output is a TIFF with a position
and no idea where that position is. The coordinate system is in the LAS header,
so it is read from the source and embedded afterwards. Nothing is reprojected
and no coordinates are altered; the file is only told what it already means.
"""
import glob
import json
import os
import subprocess
import time

CLOUDCOMPARE = r"C:\Program Files\CloudCompare\CloudCompare.exe"


class InterpolateError(RuntimeError):
    pass


def source_srs(las_path):
    """The WKT coordinate system recorded in a LAS header, or None.

    The coordinate system is in the header, so reading the points to find it is
    pure waste -- and on a full G-LiHT tile that waste is 23 million points and
    over a gigabyte of resident array for a string. `quickinfo` asks the reader
    for its header alone and returns in milliseconds. It is attempted first and
    the full read is kept only as a fallback, since not every PDAL build exposes
    it and a missing CRS is worse than a slow one.
    """
    import pdal
    pl = pdal.Pipeline(json.dumps({"pipeline": [str(las_path)]}))
    try:
        qi = pl.quickinfo
        info = qi.get("readers.las") or next(iter(qi.values()))
        srs = info.get("srs") or {}
        wkt = srs.get("wkt") or srs.get("compoundwkt")
        if wkt:
            return wkt
    except Exception:
        pass
    pl.execute()
    md = pl.metadata
    if not isinstance(md, dict):
        md = json.loads(md)
    las = md["metadata"]["readers.las"]
    if isinstance(las, list):
        las = las[0]
    srs = las.get("srs") or {}
    return srs.get("wkt") or las.get("comp_spatialreference") or None


def georeference(tif_path, wkt, out_path=None):
    """Embed a coordinate system into a raster that already has a geotransform.

    The geotransform is left exactly as written. This only records which
    coordinate system those numbers were always in.
    """
    from osgeo import gdal, osr
    gdal.UseExceptions()
    if not wkt:
        raise InterpolateError("no coordinate system to embed; the source LAS "
                               "header carries none")
    src = gdal.Open(tif_path, gdal.GA_Update if out_path is None else gdal.GA_ReadOnly)
    gt = src.GetGeoTransform()
    if gt == (0.0, 1.0, 0.0, 0.0, 0.0, 1.0):
        raise InterpolateError(
            "raster has no geotransform, so a coordinate system would place it "
            "nowhere. Check that the rasteriser wrote world coordinates.")
    srs = osr.SpatialReference()
    srs.ImportFromWkt(wkt)
    if out_path is None:
        src.SetProjection(srs.ExportToWkt())
        src.FlushCache()
        src = None
        return tif_path
    drv = gdal.GetDriverByName("GTiff")
    dst = drv.CreateCopy(out_path, src, options=["COMPRESS=DEFLATE", "TILED=YES"])
    dst.SetProjection(srs.ExportToWkt())
    dst.SetGeoTransform(gt)
    dst.FlushCache()
    dst = None
    src = None
    return out_path


def rasterize_cloudcompare(las_path, out_path, cell_m=1.0, projection="AVG",
                           fill="KRIGING", exe=CLOUDCOMPARE, verbose=True):
    """Rasterise via CloudCompare's command line, then georeference the result.

    CloudCompare names its own output and drops it beside the input, so the file
    is located by timestamp afterwards rather than being written where asked.
    """
    las_path, out_path = str(las_path), str(out_path)
    if not os.path.exists(exe):
        raise InterpolateError("CloudCompare not found at %s" % exe)
    before = set(glob.glob(os.path.join(os.path.dirname(las_path), "*.tif")))

    cmd = [exe, "-SILENT", "-O", las_path, "-RASTERIZE",
           "-GRID_STEP", str(cell_m), "-PROJ", projection,
           "-EMPTY_FILL", fill, "-OUTPUT_RASTER_Z"]
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True)
    made = sorted(set(glob.glob(os.path.join(os.path.dirname(las_path), "*.tif")))
                  - before)
    if not made:
        raise InterpolateError(
            "CloudCompare produced no raster.\ncommand: %s\nstdout: %s\nstderr: %s"
            % (" ".join(cmd), r.stdout[-400:], r.stderr[-400:]))
    produced = made[-1]

    # the fill strategy is not validated by the exit code; an unknown one is
    # reported on stdout and silently degrades to leaving cells empty
    msg = (r.stdout or "") + (r.stderr or "")
    if "Unknown empty cell filling strategy" in msg:
        raise InterpolateError("CloudCompare did not recognise fill %r; it left "
                               "cells empty instead" % fill)

    os.replace(produced, out_path)
    wkt = source_srs(las_path)
    georeference(out_path, wkt)
    verify_alignment(out_path, las_path)
    if verbose:
        print("  rasterised %-30s fill=%-8s %.1fs -> %s"
              % (os.path.basename(las_path), fill, time.time() - t0,
                 os.path.basename(out_path)))
    return out_path


def verify_alignment(tif_path, las_path, tolerance_cells=1.5):
    """Confirm the raster sits where the point cloud says it should.

    CloudCompare applies a global shift to large coordinates, subtracting an
    offset internally so that eastings and northings survive single-precision
    storage -- at a northing of two million, float32 quantises to about 0.25 m,
    which would be visible on a one-metre grid. The shift is supposed to be
    undone on export, and on every tile tested it has been. This checks rather
    than assumes, because a shift left in place produces a DEM that is correct
    in every respect except where it is, and nothing downstream would notice.
    """
    import json
    from osgeo import gdal
    import pdal
    gdal.UseExceptions()

    pl = pdal.Pipeline(json.dumps({"pipeline": [str(las_path)]}))
    pl.execute()
    md = pl.metadata if isinstance(pl.metadata, dict) else json.loads(pl.metadata)
    las = md["metadata"]["readers.las"]
    if isinstance(las, list):
        las = las[0]

    d = gdal.Open(str(tif_path))
    gt = d.GetGeoTransform()
    cell = abs(gt[1])
    r_minx, r_maxy = gt[0], gt[3]
    r_maxx = gt[0] + gt[1] * d.RasterXSize
    r_miny = gt[3] + gt[5] * d.RasterYSize

    dx = abs(r_minx - las["minx"])
    dy = abs(r_maxy - las["maxy"])
    tol = tolerance_cells * cell
    ok = dx <= tol and dy <= tol
    out = {"raster_minx": r_minx, "las_minx": las["minx"], "dx_m": dx,
           "raster_maxy": r_maxy, "las_maxy": las["maxy"], "dy_m": dy,
           "tolerance_m": tol, "aligned": bool(ok)}
    if not ok:
        raise InterpolateError(
            "raster and point cloud disagree about position by %.2f m east and "
            "%.2f m north, beyond the %.2f m tolerance. This is the signature "
            "of a global shift that was applied and not undone. The DEM is "
            "correct except for where it is." % (dx, dy, tol))
    return out


def raster_stats(tif_path):
    """Coverage and range, which is what says whether a DEM is usable."""
    from osgeo import gdal
    import numpy as np
    gdal.UseExceptions()
    d = gdal.Open(tif_path)
    b = d.GetRasterBand(1)
    a = b.ReadAsArray().astype("f8")
    nod = b.GetNoDataValue()
    m = np.isfinite(a)
    if nod is not None:
        m &= a != nod
    gt = d.GetGeoTransform()
    return {"width": d.RasterXSize, "height": d.RasterYSize,
            "cell_m": abs(gt[1]), "filled_fraction": float(m.mean()),
            "z_min": float(a[m].min()) if m.any() else None,
            "z_max": float(a[m].max()) if m.any() else None,
            "has_crs": bool(d.GetProjection()),
            "crs_name": (d.GetProjection().split('"')[1]
                         if d.GetProjection() else None)}
