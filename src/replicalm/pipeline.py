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
    visualize     optional, and a shell-out rather than a reimplementation

THE VISUALIZATION IS SOMEONE ELSE'S RECIPE
------------------------------------------
G1 is the four-layer blend published as Table 3 of Britton et al. 2025, and the
script that produces it (`GLiHT_rvt.py`) belongs to that work rather than to
this package. It is called, not copied, so that one definition of the recipe
exists and this package cannot drift from it.

It is also optional. It needs `rvt-py`, which is Apache 2.0 from ZRC SAZU and
the University of Ljubljana, and which pulls in GDAL and matplotlib. A run that
only wants elevation should not have to install any of that, so the import and
the call are deferred and their absence is reported rather than raised.

Point the shell-out at the script with REPLICALM_RVT_SCRIPT, or pass
`rvt_script=`. Without it, `process` produces the DEM and says plainly that the
image step was skipped.
"""
import os
import subprocess
import sys
import time

import numpy as np


class PipelineError(RuntimeError):
    pass


def _noop(stage, message, fraction=None):
    pass


def process(las_path, out_dir, cfg=None, cell_m=None, make_g1=False,
            rvt_script=None, progress=None, keep_ground=True,
            derive_cell=True):
    """Process one tile. Returns a dict describing everything written.

    `progress(stage, message, fraction)` is called as work proceeds, so a GUI
    can show what is happening without this module knowing what a GUI is.
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

    progress("classify", "classifying ground returns", 0.05)
    ground_las = os.path.join(out_dir, stem + "_ground.laz")
    classify.classify_tile(las_path, ground_las, cfg, verbose=False)

    arr, _ = classify.read_points(ground_las)
    gnd = arr[arr["Classification"] == 2]

    # Residual low vegetation the ground filter accepted. Removing it is the
    # delivered method; set clean_vegetation False for the unfiltered surface,
    # which is the pure translation of the source and what the published
    # comparison figures are measured against.
    if getattr(cfg, "clean_vegetation", False) and len(gnd):
        from . import cleanup
        above = cleanup.height_above_floor(
            gnd["X"].astype("f8"), gnd["Y"].astype("f8"),
            gnd["Z"].astype("f8"),
            patch=getattr(cfg, "clean_patch_m", 0.75),
            percentile=getattr(cfg, "clean_percentile", 10.0))
        keep = above <= getattr(cfg, "clean_height_m", 0.20)
        progress("classify", "removed %d of %d ground returns (%.1f%%) standing "
                             "above the ground around them"
                 % ((~keep).sum(), len(keep), 100 * (~keep).mean()))
        gnd = gnd[keep]

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

    v = K.fit_variogram(x, y, z, model=cfg.variogram_model,
                        max_lag=cfg.variogram_fit_lag_m * 2)
    dem, info = K.krige_grid(x, y, z, g, radius=radius,
                             max_points=cfg.max_points,
                             min_points=cfg.min_points, variogram=v,
                             verbose=False,
                             chunk_cells=getattr(cfg, "chunk_cells", 1000000))
    if info["fallback_fraction"] > 0.25:
        progress("interpolate",
                 "%.0f%% of cells fell back to inverse distance: the "
                 "neighbourhood is too tightly packed for the covariance model"
                 % (100 * info["fallback_fraction"]))

    progress("finalise", "filling holes, trimming the coverage edge", 0.75)
    wkt = interpolate.source_srs(las_path)
    dem_path = os.path.join(out_dir, stem + "_DEM.tif")
    report = finalise.finalise(dem, g, wkt, dem_path, radius_m=radius,
                               nodata_in=-9999.0,
                               erode_factor=cfg.erode_factor, verbose=False)

    out = {"dem": dem_path, "ground_las": ground_las if keep_ground else None,
           "ground_points": int(len(z)), "density": density,
           "search_radius_m": radius, "cell_m": cell_m,
           "fallback_fraction": info["fallback_fraction"],
           "locked_baseline": baseline, "seconds": round(time.time() - t0, 1)}
    out.update({k: report[k] for k in
                ("filled_cells", "holes", "erode_cells", "erode_m",
                 "cells_after_trim", "filled_fraction")})
    if not keep_ground and os.path.exists(ground_las):
        os.remove(ground_las)

    if make_g1:
        progress("visualize", "building the G1 composite", 0.85)
        out["g1"] = visualize(dem_path, out_dir, rvt_script=rvt_script,
                              progress=progress)
    progress("done", "finished in %.0f s" % (time.time() - t0), 1.0)
    return out


def find_rvt_script(explicit=None):
    """Locate GLiHT_rvt.py, or return None. Never raises."""
    for candidate in (explicit, os.environ.get("REPLICALM_RVT_SCRIPT"),
                      os.path.join(os.path.dirname(sys.executable),
                                   "GLiHT_rvt.py"),
                      r"C:\g1\tools\GLiHT_rvt.py"):
        if candidate and os.path.isfile(candidate):
            return candidate
    return None


def visualize(dem_path, out_dir, rvt_script=None, res="0p5m", progress=None):
    """Produce the G1 composite by calling the published recipe.

    Returns the output directory, or a dict explaining why nothing was made.
    The recipe expects DEMs named `*_DEM_<res>_v1.tif` under a root it walks,
    so the DEM is linked into that shape rather than renamed in place.
    """
    progress = progress or _noop
    script = find_rvt_script(rvt_script)
    if not script:
        progress("visualize", "GLiHT_rvt.py not found; set REPLICALM_RVT_SCRIPT "
                              "to enable the image step. DEM is unaffected.")
        return {"skipped": "GLiHT_rvt.py not found"}

    stem = os.path.splitext(os.path.basename(dem_path))[0]
    staged_root = os.path.join(out_dir, "_rvt_in")
    staged_dir = os.path.join(staged_root, "tile")
    os.makedirs(staged_dir, exist_ok=True)
    staged = os.path.join(staged_dir, "%s_DEM_%s_v1.tif" % (stem, res))
    if not os.path.exists(staged):
        try:
            os.link(dem_path, staged)          # cheap; same volume
        except OSError:
            import shutil
            shutil.copy2(dem_path, staged)

    rvt_root = os.path.join(out_dir, "rvt")
    cmd = [sys.executable, script, "--res", res, "--no-metadata",
           "--workers", "1", "--dem-root", staged_root, "--rvt-root", rvt_root]
    progress("visualize", "running %s" % os.path.basename(script))
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        tail = (proc.stderr or proc.stdout or "").strip().splitlines()[-3:]
        progress("visualize", "the image step failed: %s" % " / ".join(tail))
        return {"failed": proc.returncode, "stderr": "\n".join(tail)}
    return {"rvt_root": rvt_root}


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
                   help="output cell size in metres (default 1.0)")
    p.add_argument("--g1", action="store_true", help="also build the G1 image")
    p.add_argument("--rvt-script", default=None,
                   help="path to GLiHT_rvt.py for the image step")
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

    out = process(a.las, a.out, cfg=cfg, cell_m=a.cell,
                  make_g1=a.g1, rvt_script=a.rvt_script, progress=show)
    print("\nDEM: %s" % out["dem"])
    print("%d ground points, %.2f per m2, radius %.2f m, %d cells kept"
          % (out["ground_points"], out["density"], out["search_radius_m"],
             out["cells_after_trim"]))
    if not out["locked_baseline"]:
        print("NOTE: this run did not use the locked baseline configuration.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
