"""One tile, end to end: cloud in, DEM out, and optionally the G1 image.

This is what the desktop launcher calls and what a command line drives. It
exists so that the GUI contains no processing logic -- a front end that also
decides things is a front end whose behaviour cannot be reproduced from a
script.

WHAT IT DOES
------------
    classify      the locked baseline ground filter
    interpolate   kriging onto a declared grid, radius scaled to point density
    finalise      fill enclosed holes, trim the one-sided fringe, write with
                  zero meaning no data
    visualize     optional: the G1 composite and the four layers it blends

THE VISUALIZATION USED TO BE A SHELL-OUT AND DID NOT WORK
---------------------------------------------------------
G1 is the four-layer blend published as Table 3 of Britton et al. 2025. This
module used to call `GLiHT_rvt.py` for it, found beside the interpreter or at a
path on the machine the script was written on, so that one definition of the
recipe existed and this package could not drift from it. That was the right
instinct and the wrong mechanism: the script is not in the installer, so on
every user's machine the step returned "GLiHT_rvt.py not found", left the DEM
alone, and reported it only in a log line nobody reads. The feature was absent
rather than optional.

The recipe is now in `replicalm.visualize`, which calls `rvt.vis` as a library.
`rvt-py` was already pinned in packaging/environment.yml and already shipped,
so this costs no new dependency and nothing in installer size. The import stays
deferred: a run that only wants elevation should not need the visualization
stack present, and its absence is reported rather than raised.

What it writes is G1 and its four source layers under `<out>/G1/`, and nothing
else. RVT remains the tool for relief visualization; this renders one image so
the operator can look at what the ground filter produced.
"""
import os
import time

import numpy as np


class PipelineError(RuntimeError):
    pass


def _noop(stage, message, fraction=None):
    pass


def process(las_path, out_dir, cfg=None, cell_m=None, make_g1=False,
            progress=None, keep_ground=True,
            derive_cell=True, grid=None, on_grid=None,
            reuse_ground=False, ground_stage="after"):
    """Process one tile. Returns a dict describing everything written.

    `progress(stage, message, fraction)` is called as work proceeds, so a GUI
    can show what is happening without this module knowing what a GUI is.

    `grid` overrides the lattice entirely; `on_grid` is a reference raster whose
    lattice is adopted, which is the usual way to ask for one. Either makes the
    output directly comparable to that raster, cell for cell.

    This matters because production derives its cell from measured density and
    snaps cell edges to whole multiples, so neighbouring tiles mosaic without
    resampling. Reference DEMs from the original workflow have neither property
    -- kriging to a data-derived extent in Surfer produced cells of 0.500042 m
    on origins that are multiples of nothing -- so a production raster and an
    archive raster never share a lattice, and `calibrate.compare_to_reference`
    refuses the comparison rather than score a half-cell offset. Building a
    second candidate on the archive's grid is the answer; resampling the first
    one afterwards is not, because that measures the resampler.

    `reuse_ground` takes the classified ground beside the output if it is
    already there, which is what makes building that second candidate cheap:
    classification is the expensive half and does not depend on the lattice.
    It reuses whatever is on disk without checking it came from this config, so
    it is for a second grid over the same run, not for resuming a changed one.

    `keep_ground` saves the ground returns as LAZ. `ground_stage` chooses when:

      "after"   (default) the ground Clear kept, written to
                `<stem>_ground_cleaned.laz`. This is the surface the DEM was
                actually built from.
      "before"  the ground as the filter classified it, residual low vegetation
                included, written to `<stem>_ground.laz`. Re-deriving the
                cleanup at another threshold starts here and skips
                classification, which is the expensive half.

    The two go to different filenames deliberately. `reuse_ground` looks for
    `<stem>_ground.laz`, so a cleaned file saved under that name would be
    cleaned a second time on the next run, quietly removing more ground each
    pass. Keeping the names apart makes that impossible rather than merely
    unlikely.
    """
    from . import classify, finalise, grid as G, interpolate, kriging as K
    from .config import PRESETS, verify_baseline

    progress = progress or _noop
    cfg = cfg or PRESETS["clear"]
    # None means derive it from measured density; cfg.dem_cell_m is the
    # source's fixed figure and is used only if explicitly asked for.
    if cell_m is None and getattr(cfg, "dem_cell_m", None) and derive_cell is False:
        cell_m = cfg.dem_cell_m
    os.makedirs(out_dir, exist_ok=True)
    stem = os.path.splitext(os.path.basename(str(las_path)))[0]
    t0 = time.time()

    # A run should say which configuration produced it, and stop if that
    # configuration is not the one the measurements were taken against.
    try:
        verify_baseline(cfg)
        baseline = True
    except ValueError as e:
        baseline = False
        progress("config", "NOT the locked baseline: %s" % str(e).splitlines()[-1].strip())
    cfg.to_json(os.path.join(out_dir, stem + "_config.json"))

    ground_las = os.path.join(out_dir, stem + "_ground.laz")
    if reuse_ground and os.path.exists(ground_las):
        progress("classify", "reusing the ground beside the output", 0.05)
    else:
        progress("classify", "classifying ground returns", 0.05)
        classify.classify_tile(las_path, ground_las, cfg, verbose=False)

    arr, _ = classify.read_points(ground_las)
    gnd = arr[arr["Classification"] == 2]
    cleaned = False

    # Residual low vegetation the ground filter accepted. Removing it is the
    # delivered method; set clean_vegetation False for the unfiltered surface,
    # which is the pure translation of the source and what the published
    # comparison figures are measured against.
    if getattr(cfg, "clean_vegetation", False) and len(gnd):
        from . import cleanup
        gx = gnd["X"].astype("f8"); gy = gnd["Y"].astype("f8")
        gz = gnd["Z"].astype("f8")
        if getattr(cfg, "clean_slope_aware", False):
            fit_r = (getattr(cfg, "clean_plane_patch_m", None)
                     or getattr(cfg, "clean_patch_m", 0.75))
            progress("classify", "cleanup measured against a plane fitted over "
                                 "%.2f m, so the rule does not fire on slope"
                     % fit_r)
            above = cleanup.height_above_plane(
                gx, gy, gz, patch=getattr(cfg, "clean_patch_m", 0.75),
                plane_patch=fit_r)
        else:
            above = cleanup.height_above_floor(
                gx, gy, gz,
                patch=getattr(cfg, "clean_patch_m", 0.75),
                percentile=getattr(cfg, "clean_percentile", 10.0))
        keep = above <= getattr(cfg, "clean_height_m", 0.20)
        progress("classify", "removed %d of %d ground returns (%.1f%%) standing "
                             "above the ground around them"
                 % ((~keep).sum(), len(keep), 100 * (~keep).mean()))
        gnd = gnd[keep]
        cleaned = True

    if len(gnd) < cfg.min_ground_points:
        raise PipelineError(
            "only %d ground points; below min_ground_points (%d). Rasterising "
            "this would not be meaningful." % (len(gnd), cfg.min_ground_points))
    x = gnd["X"].astype("f8"); y = gnd["Y"].astype("f8"); z = gnd["Z"].astype("f8")

    # A cell size can be given, or derived from the density actually measured.
    # 1/sqrt(density) is where rasterised block cross-validation found the
    # residual curve turning; see grid.cell_for_density.
    #
    # Density is measured over the ground the survey covers, not over the
    # bounding box. A flight strip crosses its box on the diagonal, so the box
    # is mostly empty: a full G-LiHT tile reads 1.24 returns per square metre
    # that way against about 4.2 where the returns are, and a cell derived from
    # the first figure comes out nearly twice too coarse.
    density, covered = G.covered_density(x, y)
    if on_grid is not None and grid is None:
        grid = G.grid_from_raster(on_grid)
    if grid is not None:
        g = grid
        cell_m = g.cell
        progress("interpolate",
                 "cell size %.4f m and origin adopted from a reference "
                 "lattice, so the result is comparable to it cell for cell"
                 % cell_m)
    else:
        if cell_m is None:
            cell_m = G.cell_for_density(density,
                                        factor=getattr(cfg, 'cell_factor', 1.0))
            progress("interpolate",
                     "cell size %.2f m derived from %.2f returns per m2 over "
                     "%.2f km2 of covered ground"
                     % (cell_m, density, covered / 1e6))
        g = G.grid_for_las(las_path, cell=cell_m)
    radius = (cfg.search_radius_m if cfg.search_radius_mode == "fixed"
              and cfg.search_radius_m else
              K.radius_for_density(density, cfg.max_points,
                                   cfg.search_radius_factor,
                                   cfg.search_radius_floor_m,
                                   cfg.search_radius_ceiling_m))
    progress("interpolate",
             "%d ground points, %.2f per m2, search radius %.2f m"
             % (len(z), density, radius), 0.35)

    # Coincident XY make every solver's system singular, so the points are
    # conditioned once here rather than inside whichever engine runs. See
    # preprocess.dedupe_xy: the kriging path used to carry its own copy, and a
    # spline tried on the same data failed for the same reason because the fix
    # was not reachable from outside kriging.
    from .preprocess import dedupe_xy, describe
    # "mean", not "lowest": by this point the set is already ground,
    # and a delivered or classified ground population can hold two
    # returns at one XY straddling the surface. Taking the minimum
    # biases every collision downward; averaging has no preference.
    x, y, z, dedupe_report = dedupe_xy(x, y, z, keep="mean")
    note = describe(dedupe_report)
    if note:
        progress("interpolate", note)

    engine = getattr(cfg, "interpolator", "kriging")
    if engine == "spline":
        from . import spline as S
        progress("interpolate",
                 "minimum curvature with tension %.3f, %d neighbours"
                 % (getattr(cfg, "spline_tension", S.DEFAULT_TENSION),
                    getattr(cfg, "spline_neighbours", S.DEFAULT_NEIGHBOURS)))
        dem, info = S.spline_grid(
            x, y, z, g,
            tension=getattr(cfg, "spline_tension", S.DEFAULT_TENSION),
            neighbours=getattr(cfg, "spline_neighbours", S.DEFAULT_NEIGHBOURS),
            smoothing=getattr(cfg, "spline_smoothing", S.DEFAULT_SMOOTHING),
            min_points=cfg.min_points, verbose=False,
            chunk_cells=getattr(cfg, "chunk_cells", 1000000))
    elif engine == "kriging":
        v = K.fit_variogram(x, y, z, model=cfg.variogram_model,
                            max_lag=cfg.variogram_fit_lag_m * 2)
        dem, info = K.krige_grid(x, y, z, g, radius=radius,
                                 max_points=cfg.max_points,
                                 min_points=cfg.min_points, variogram=v,
                                 verbose=False,
                                 chunk_cells=getattr(cfg, "chunk_cells",
                                                     1000000))
    elif engine in ("biharmonic", "surfer"):
        # The Surfer half of the reference chain is a global minimum-curvature
        # solve, and this is the engine that reproduces it. It was measured
        # against the delivered archive at 0.0524 m with a roughness ratio of
        # 0.939 while the local spline above loses half the coherent detail
        # below 1 m -- and until now it was not reachable from here at all,
        # only from analysis scripts, so no released build has ever run it.
        from . import surfer_algoA as SA
        progress("interpolate",
                 "global biharmonic minimum curvature, %.0f m cores with a "
                 "%.0f m buffer, pooling %.2f m"
                 % (getattr(cfg, "biharmonic_core_m", SA.DEFAULT_CORE_M),
                    getattr(cfg, "biharmonic_buffer_m", SA.DEFAULT_BUFFER_M),
                    getattr(cfg, "biharmonic_pool_sigma_m",
                            SA.DEFAULT_POOL_SIGMA_M)))
        dem, info = SA.surfer_grid_blocked(
            x, y, z, g,
            core_m=getattr(cfg, "biharmonic_core_m", SA.DEFAULT_CORE_M),
            buffer_m=getattr(cfg, "biharmonic_buffer_m", SA.DEFAULT_BUFFER_M),
            pool_sigma_m=getattr(cfg, "biharmonic_pool_sigma_m",
                                 SA.DEFAULT_POOL_SIGMA_M),
            verbose=False, progress=progress)
    else:
        raise PipelineError("unknown interpolator %r; expected 'kriging', "
                            "'spline' or 'biharmonic'" % engine)
    if info.get("fallback_fraction", 0.0) > 0.25:
        progress("interpolate",
                 "%.0f%% of cells fell back to inverse distance: the "
                 "neighbourhood is too tightly packed for the covariance model"
                 % (100 * info["fallback_fraction"]))

    progress("finalise", "filling holes, trimming the coverage edge", 0.75)
    wkt = interpolate.source_srs(las_path)

    # Save the ground, if asked, at the stage asked for. "after" is the surface
    # the DEM was built from; "before" is what the filter produced and is the
    # one to keep if the cleanup threshold may be revisited.
    kept_ground = None
    if keep_ground:
        if ground_stage == "before":
            kept_ground = ground_las
        else:
            cleaned_las = os.path.join(out_dir, stem + "_ground_cleaned.laz")
            try:
                classify.write_points(gnd, cleaned_las, cfg, srs=wkt)
                kept_ground = cleaned_las
                progress("classify", "ground after cleanup written: %s"
                         % os.path.basename(cleaned_las))
            except Exception as exc:
                progress("classify", "could not write the cleaned ground (%s); "
                         "keeping the classified ground instead"
                         % str(exc)[:60])
                kept_ground = ground_las
            if kept_ground != ground_las and not cleaned:
                # nothing was removed, so the two files would be identical
                progress("classify", "cleanup did not run for this profile; "
                         "the saved ground is the classified ground")

    dem_path = os.path.join(out_dir, stem + "_DEM.tif")
    # Two different trims, because they answer two different questions. The
    # erosion depth is the distance a cell can reach for support that exists on
    # one side only, which is the search radius and is a kriging quantity; the
    # spline and the global solver have no such radius, so for them the depth
    # is the support radius itself rather than a number borrowed from an engine
    # that is not running. The support mask is measured against the data and
    # applies either way: see config.edge_fix_m for how 4.0 m was arrived at.
    edge_fix = getattr(cfg, "edge_fix_m", None)
    strip_m = getattr(cfg, "edge_strip_m", None)
    axis = getattr(cfg, "edge_axis", None) or finalise.infer_transect_axis(
        g, os.path.basename(las_path))

    # The footprint is one span per row between the first and last ground
    # return, which is the geometry mask_dem_edges.py produces in the G1
    # reprocessing work. It replaces support_mask, whose thresholded distance
    # transform drew the boundary as a union of discs and rendered as a
    # scalloped margin with detached islands.
    #
    # The footprint cannot simply be dropped: this interpolator has no
    # distance limit and fills the whole lattice, so without it the across-
    # track strip would leave the entire bounding rectangle less its margins.
    # Surfer bounded its own output, which is why the reference script needed
    # nothing of the kind.
    support = finalise.return_footprint(x, y, g, axis)
    progress("finalise", "footprint: one span per row between the first and "
             "last ground return, %d of %d cells"
             % (int(support.sum()), support.size))
    trim_radius = radius if engine == "kriging" else (edge_fix or radius)
    if strip_m:
        progress("finalise", "edge strip: %.1f m across-track (%s), the "
                 "reference convention" % (strip_m, axis))
    report = finalise.finalise(dem, g, wkt, dem_path, radius_m=trim_radius,
                               nodata_in=-9999.0,
                               erode_factor=cfg.erode_factor,
                               support=support, strip_m=strip_m, axis=axis,
                               verbose=False)

    out = {"dem": dem_path, "ground_las": kept_ground,
           "ground_stage": ground_stage if kept_ground else None,
           "ground_points": int(len(z)), "density": density,
           "search_radius_m": radius, "cell_m": cell_m,
           "edge_fix_m": edge_fix, "edge_strip_m": strip_m,
           "edge_axis": axis,
           "fallback_fraction": info["fallback_fraction"],
           # The correlation range is what explains the fallback: when the
           # nearest max_points neighbours all sit inside it, their pairwise
           # semivariances are nearly identical, the system loses rank, and the
           # estimate is inverse distance rather than kriging. Recorded so the
           # two can be read together across a survey.
           "variogram_range_m": (float(info["params"][2])
                                 if info.get("params") is not None
                                 and len(info["params"]) > 2 else None),
           "variogram_model": info.get("model"),
           "locked_baseline": baseline, "seconds": round(time.time() - t0, 1)}
    out.update({k: report[k] for k in
                ("filled_cells", "holes", "erode_cells", "erode_m",
                 "cells_after_trim", "filled_fraction",
                 "interior_cells_restored")})
    if kept_ground != ground_las and os.path.exists(ground_las):
        os.remove(ground_las)

    if make_g1:
        progress("visualize", "building the G1 composite", 0.85)
        out["g1"] = visualize(dem_path, out_dir,
                              progress=progress)
    progress("done", "finished in %.0f s" % (time.time() - t0), 1.0)
    return out


def visualize(dem_path, out_dir, progress=None, keep_layers=True, **_legacy):
    """The G1 composite and its four source layers, under `<out_dir>/G1/`.

    Delegates to `replicalm.visualize`, which calls rvt as a library. Returns
    what was written, or a dict saying why nothing was -- the DEM is the
    product and a missing visualization stack must not fail a run that has
    already produced one.

    `**_legacy` swallows `rvt_script=` and `res=`, which the shell-out needed
    and this does not. Callers that still pass them keep working.
    """
    progress = progress or _noop
    try:
        from . import visualize as _vis
    except ImportError as exc:                       # pragma: no cover
        progress("visualize", "the image step is unavailable: %s" % exc)
        return {"skipped": str(exc)}
    try:
        return _vis.g1(dem_path, out_dir, keep_layers=keep_layers,
                       progress=progress)
    except Exception as exc:
        progress("visualize", "the image step failed: %s" % exc)
        return {"failed": str(exc)}


def main(argv=None):
    """Command line: the same work the launcher does, scriptable."""
    import argparse
    from dataclasses import replace
    from .config import PRESETS

    p = argparse.ArgumentParser(
        prog="replicalm",
        description="Bare-earth DEM from a lidar point cloud, at the locked "
                    "baseline. Optionally the G1 image as well.")
    p.add_argument("las", help="input LAS or LAZ")
    p.add_argument("-o", "--out", default=".", help="output directory")
    p.add_argument("--cell", type=float, default=None,
                   help="output cell size in metres; default is derived from "
                        "measured ground-return density")
    # The same control in the units people ask for it in. The relationship is
    # quadratic -- 4 px/m2 is a 0.50 m cell, not 0.25 -- so this exists to stop
    # that conversion being done wrongly by hand. Squares only, because those
    # are the values that land on round cell sizes.
    p.add_argument("--pixels-per-m2", type=int, default=None, metavar="N",
                   choices=(1, 4, 9, 16, 25),
                   help="output resolution as pixels per square metre: 1 is "
                        "1.00 m cells, 4 is 0.50 m, 9 is 0.33 m, 16 is 0.25 m, "
                        "25 is 0.20 m. An alternative to --cell")
    p.add_argument("--g1", action="store_true", help="also build the G1 image")
    p.add_argument("--keep-ground", action="store_true",
                   help="save the ground returns as LAZ beside the DEM")
    p.add_argument("--slope-aware-cleanup", action="store_true",
                   help="measure the cleanup height against a plane fitted to "
                        "each neighbourhood instead of the 10th percentile of "
                        "its raw elevations, so the rule stops removing ground "
                        "for being on a slope")
    p.add_argument("--ground-stage", choices=("after", "before"),
                   default="after",
                   help="which ground to save: 'after' the cleanup, the surface "
                        "the DEM was built from (default), or 'before' it, "
                        "which is where re-deriving the cleanup starts")
    # The interpolation neighbourhood. Offered as a pair because the two move
    # together: 64 points only condition well inside a radius wide enough to
    # hold them, and 16 inside a fixed 20 m lose matrix rank on dense data.
    # "published" is the source's own figure and measured better or neutral on
    # both test tiles; "fast" is 7 to 9 times quicker to interpolate.
    p.add_argument("--interpolation", choices=("published", "fast"),
                   default=None,
                   help="published: 64 points within a fixed 20 m, the source's "
                        "neighbourhood and the profile default. fast: 16 points "
                        "with the radius from density, 7-9x quicker")
    p.add_argument("--profile", "--preset", dest="preset",
                   default="clear", choices=sorted(PRESETS),
                   help="baseline | clear (default) | deep")
    # The ground filter's elevation tolerance. It decides how far a return may
    # stand above the provisional surface and still be called ground, so it is
    # the parameter that most often needs changing on a difficult tile: too
    # loose and low vegetation is accepted, too tight and real relief is cut.
    p.add_argument("--threshold", type=float, default=None, metavar="M",
                   help="SMRF elevation threshold in metres (profile default: "
                        "0.25 for clear and deep, 0.50 for baseline)")
    # The search neighbourhood. The kriging estimate weights neighbours by the
    # fitted variogram, so beyond the correlation range extra neighbours carry
    # almost no weight and the count stops mattering -- which is why 16 and 64
    # agree on surveys whose range is short against the search radius. Where the
    # range is long relative to point spacing they diverge, so the count is a
    # property of the survey rather than a constant, and is exposed.
    p.add_argument("--max-points", type=int, default=None, metavar="N",
                   help="most neighbours per estimate (profile default: 64 for "
                        "baseline, matching the source; 16 otherwise)")
    p.add_argument("--min-points", type=int, default=None, metavar="N",
                   help="fewest neighbours before a cell is left empty "
                        "(profile default: 1 for baseline, 3 otherwise)")
    p.add_argument("--reason", default=None, metavar="TEXT",
                   help="why this run departs from the method's defaults. "
                        "Recorded in the settings file beside the output")
    a = p.parse_args(argv)

    cfg = PRESETS[a.preset]
    if a.max_points is not None or a.min_points is not None:
        mx = a.max_points if a.max_points is not None else cfg.max_points
        mn = a.min_points if a.min_points is not None else cfg.min_points
        if mx < 1 or mn < 1:
            p.error("--max-points and --min-points must be at least 1")
        if mn > mx:
            p.error("--min-points (%d) cannot exceed --max-points (%d)" % (mn, mx))
        cfg = replace(cfg, max_points=mx, min_points=mn)

    # A departure from the defaults is a claim about the survey. The run is not
    # blocked for want of a reason -- that would only teach people to type
    # anything -- but an unexplained departure is called out here and left
    # blank in the settings file, where its absence is as visible as its
    # presence would have been.
    overrides = [n for n, v in (("--threshold", a.threshold),
                                ("--max-points", a.max_points),
                                ("--min-points", a.min_points)) if v is not None]
    if a.reason:
        cfg = replace(cfg, note=a.reason)
    elif overrides:
        print("NOTE: %s set without --reason. The defaults are what the source "
              "specifies and what the published measurements were taken at; a "
              "departure should say why, and the settings file will record that "
              "it did not." % ", ".join(overrides))
    if a.threshold is not None:
        if a.threshold <= 0:
            p.error("--threshold must be positive")
        cfg = replace(cfg, passes=[replace(cfg.passes[0],
                                           threshold_m=a.threshold)]
                                  + list(cfg.passes[1:]))

    def show(stage, message, fraction=None):
        print("[%-11s] %s" % (stage, message))

    if a.slope_aware_cleanup:
        cfg = replace(cfg, clean_slope_aware=True)
    if a.interpolation == "fast":
        cfg = replace(cfg, max_points=16, min_points=3,
                      search_radius_mode="density", search_radius_m=None)
    elif a.interpolation == "published":
        cfg = replace(cfg, max_points=64, min_points=1,
                      search_radius_mode="fixed", search_radius_m=20.0)

    cell = a.cell
    if a.pixels_per_m2 is not None:
        if a.cell is not None:
            p.error("give --cell or --pixels-per-m2, not both")
        cell = (1.0 / a.pixels_per_m2) ** 0.5
        print("resolution %d px/m2 is a %.4f m cell"
              % (a.pixels_per_m2, cell))

    out = process(a.las, a.out, cfg=cfg, cell_m=cell,
                  make_g1=a.g1, progress=show,
                  keep_ground=a.keep_ground, ground_stage=a.ground_stage)
    print("\nDEM: %s" % out["dem"])
    g1 = out.get("g1")
    if isinstance(g1, dict) and g1.get("g1"):
        print("G1: %s" % g1["g1"])
        if g1.get("layers"):
            print("G1 source layers: %s"
                  % os.path.dirname(next(iter(g1["layers"].values()))))
    elif isinstance(g1, dict):
        print("image step not completed: %s"
              % (g1.get("skipped") or g1.get("failed")))
    if out.get("ground_las"):
        print("Ground (%s cleanup): %s" % (out["ground_stage"],
                                           out["ground_las"]))
    print("%d ground points, %.2f per m2, radius %.2f m, %d cells kept"
          % (out["ground_points"], out["density"], out["search_radius_m"],
             out["cells_after_trim"]))
    if not out["locked_baseline"]:
        print("NOTE: this run did not use the locked baseline configuration.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
