"""Steps 3 to 7: reset, classify ground in two passes, height-filter, export.

Everything here runs through PDAL, which ships inside ArcGIS Pro's Python
environment, so no second environment is needed and no data crosses a process
boundary between classification and interpolation.

WHY THE PASSES ARE SEPARATE PIPELINES
-------------------------------------
The source method runs Classify Ground twice, the second pass reusing the first
pass's ground and searching again with a steeper tolerance. PDAL's ground filters
do not reuse prior classification; each one starts from the points it is given.
The second pass is therefore run on the FIRST pass's ground points only, which
is what "reuse previously classified ground" amounts to: refine the ground set,
do not re-examine what was already rejected. Points demoted by the second pass
return to unclassified rather than vanishing, so nothing is lost silently.
"""
import json
import os

import numpy as np


class ClassifyError(RuntimeError):
    pass


def _pdal():
    try:
        import pdal
    except ImportError as e:
        raise ClassifyError(
            "PDAL not importable. Run this with ArcGIS Pro's Python, which "
            "ships it: "
            r'"C:\Program Files\ArcGIS\Pro\bin\Python\envs\arcgispro-py3\python.exe"'
        ) from e
    return pdal


def _ground_stage(p, noise_class=None):
    """One GroundPass -> the PDAL stage dict that implements it."""
    if p.algorithm == "smrf":
        st = {"type": "filters.smrf", "window": p.window_m, "slope": p.slope,
              "threshold": p.threshold_m, "scalar": p.scalar, "cell": p.cell_m}
        if noise_class is not None:
            # without this the filter reads the points ELM just condemned
            st["ignore"] = "Classification[%d:%d]" % (noise_class, noise_class)
        return st
    if p.algorithm == "pmf":
        return {"type": "filters.pmf", "max_window_size": p.window_m,
                "slope": p.slope, "max_distance": p.threshold_m}
    if p.algorithm == "csf":
        return {"type": "filters.csf", "resolution": p.csf_resolution_m,
                "rigidness": p.csf_rigidness, "threshold": p.csf_threshold_m}
    raise ClassifyError("unknown ground algorithm %r" % p.algorithm)


# LAS 1.4 with point data record format 6 or higher must record its coordinate
# system as WKT, and must set a bit in the global encoding field saying so. Many
# files in circulation set the format and not the bit, and PDAL refuses them:
# "Global encoding WKT flag not set for point format 6 - 10". The file is
# malformed, but it is the file the user has, and refusing it means refusing
# much of the modern LAS a survey will hand over.
#
# `nosrs` tells the reader to skip the coordinate system, which gets the points
# in. It cannot invent the projection, so a run that falls back this way has no
# CRS unless the caller supplies one -- which is why it reports rather than
# quietly continuing.
WKT_FLAG_ERROR = "Global encoding WKT flag not set"


def las_reader(path, srs=None, nosrs=False):
    """The reader stage, with the options awkward files need."""
    stage = {"type": "readers.las", "filename": str(path)}
    if srs:
        stage["override_srs"] = str(srs)
    if nosrs:
        stage["nosrs"] = True
    return stage


def run_pipeline(stages, verbose=True):
    """Execute a PDAL pipeline, retrying without the SRS if the header is bad.

    Returns (pipeline, point_count, fell_back).
    """
    pdal = _pdal()
    try:
        pl = pdal.Pipeline(json.dumps({"pipeline": stages}))
        return pl, pl.execute(), False
    except RuntimeError as e:
        if WKT_FLAG_ERROR not in str(e):
            raise
    retry = list(stages)
    first = retry[0]
    if isinstance(first, dict) and first.get("type") == "readers.las":
        first = dict(first); first["nosrs"] = True
    else:
        first = las_reader(first, nosrs=True)
    retry[0] = first
    if verbose:
        print("  this file declares LAS 1.4 point format 6-10 without the WKT "
              "flag its own header requires; re-reading without the coordinate "
              "system. Supply one with srs= or the output will be unprojected.")
    pl = pdal.Pipeline(json.dumps({"pipeline": retry}))
    return pl, pl.execute(), True


def read_points(path, srs=None):
    """Read a LAS/LAZ into a structured array, with its header metadata."""
    pl, n, _fell_back = run_pipeline([las_reader(path, srs=srs)], verbose=False)
    if not n:
        raise ClassifyError("no points read from %s" % path)
    return pl.arrays[0], pl.metadata


def classify_tile(in_path, out_path, cfg, verbose=True, srs=None):
    """Run steps 3 to 7 on one tile. Returns a dict of counts for the log.

    The counts matter as much as the output: a tile whose ground fraction is far
    from its neighbours' is the signature of a classification that has gone
    wrong, and that is visible here before anything is interpolated.
    """
    pdal = _pdal()
    stages = [las_reader(in_path, srs=srs)]

    # step 3: every point back to unclassified
    if cfg.reset_classification:
        stages.append({"type": "filters.assign",
                       "value": "Classification = 0"})

    # noise removal, ahead of the ground passes.
    #
    # ELM marks low outliers so SMRF does not anchor the surface to them;
    # the statistical outlier filter removes isolated returns in all
    # directions. Both write to the noise class rather than deleting, so the
    # ground filter can be told to ignore them and the decision stays visible
    # in the output rather than being silently applied.
    if getattr(cfg, "remove_low_noise", False):
        stages.append({"type": "filters.elm", "cell": cfg.elm_cell_m,
                       "threshold": cfg.elm_threshold_m,
                       "class": cfg.noise_class})
    if getattr(cfg, "remove_outliers", False):
        stages.append({"type": "filters.outlier", "method": "statistical",
                       "mean_k": cfg.outlier_neighbours,
                       "multiplier": cfg.outlier_multiplier,
                       "class": cfg.noise_class})

    # steps 4 and 5: ground passes.
    #
    # The source runs Classify Ground twice, the second "using previously
    # classified ground returns and a new minimum elevation search". SMRF has
    # no such refinement mode: each invocation classifies the cloud it is given,
    # from scratch. An earlier version here approximated the reuse by filtering
    # to ground before the second pass, which did restrict pass 2 to pass 1's
    # ground -- and silently destroyed every other return, so step 6 had nothing
    # left to find and class 8 came out empty on a real tile.
    #
    # Both passes therefore run over the whole cloud, the second superseding the
    # first with a steeper tolerance. That is a real divergence from the source
    # and is recorded rather than papered over: with SMRF the two-pass structure
    # collapses to the second pass's parameters. If the calibration shows the
    # refinement mattered, the way to get it is a ground-seeded filter such as
    # CSF with a prior, not a pipeline trick.
    for p in cfg.passes:
        stages.append(_ground_stage(
            p, cfg.noise_class if getattr(cfg, 'remove_low_noise', False)
            or getattr(cfg, 'remove_outliers', False) else None))

    # step 6: height above the ground MODEL, not above the nearest ground return.
    #
    # The source method computes "the height (hag) of each return above the
    # ground model", then classifies returns between -0.2 and 0.2 m into class
    # 8. hag_nn measures to the nearest ground point instead, which is a
    # different quantity: on a slope, or anywhere ground returns are sparse,
    # the nearest one sits at a different elevation than the surface beneath
    # the point. The first run on a real tile showed the consequence -- a single
    # point landed in the class 8 band out of 59,738, because every ground
    # return measures zero against itself and little else fell in between.
    #
    # hag_delaunay triangulates the ground returns and measures to that
    # surface, which is what "above the ground model" means.
    stages.append({"type": "filters.hag_delaunay"})
    stages.append({"type": "filters.expression",
                   "expression": "HeightAboveGround >= %g && "
                                 "HeightAboveGround <= %g"
                                 % (cfg.hag_low_cut_m, cfg.hag_high_cut_m)})
    lo, hi = cfg.near_ground_band_m
    stages.append({"type": "filters.assign",
                   "value": "Classification = %d WHERE "
                            "HeightAboveGround > %g && HeightAboveGround < %g "
                            "&& Classification != 2"
                            % (cfg.near_ground_class, lo, hi)})

    # step 7: ground and near-ground only, LAS 1.2
    keep = " || ".join("Classification == %d" % c for c in cfg.keep_classes)
    stages.append({"type": "filters.expression", "expression": keep})
    # LAS 1.4 point format 6, written as LAZ. `forward` is held to scale and
    # offset so nothing else carries over from the input header -- clouds
    # delivered from TerraScan bring proprietary Terrasolid records that inflate
    # the file and mean nothing elsewhere. The coordinate system is written
    # explicitly from the source rather than forwarded, and PDAL sets the WKT
    # flag that 1.4 requires, so the output is the well-formed version of the
    # format that so often arrives malformed.
    writer = {"type": "writers.las", "filename": str(out_path),
              "minor_version": int(str(cfg.las_version).split(".")[1]),
              "dataformat_id": int(getattr(cfg, "las_point_format", 6)),
              "compression": ("laszip" if getattr(cfg, "las_compression", True)
                              else "false"),
              "forward": getattr(cfg, "las_forward", "scale,offset"),
              "software_id": "Replicalm"}
    if srs:
        writer["a_srs"] = str(srs)
    else:
        from .interpolate import source_srs
        found = source_srs(in_path)
        if found:
            writer["a_srs"] = found
    stages.append(writer)

    pl, kept, fell_back = run_pipeline(stages, verbose=verbose)

    arr = pl.arrays[0] if pl.arrays else np.zeros(0)
    cls = arr["Classification"] if len(arr) else np.zeros(0, "u1")
    out = {"input": str(in_path), "output": str(out_path),
           "points_out": int(kept),
           "ground": int((cls == 2).sum()),
           "near_ground": int((cls == cfg.near_ground_class).sum())}
    if verbose:
        print("  %-38s %8d kept  %8d ground  %7d class-%d"
              % (os.path.basename(str(in_path)), out["points_out"],
                 out["ground"], out["near_ground"], cfg.near_ground_class))
    if out["ground"] < cfg.min_ground_points:
        out["warning"] = ("only %d ground points; too few to rasterise honestly"
                          % out["ground"])
    return out
