r"""One tile, one configuration, on the archive's grid, with its G1.

Built to settle which of Clear's two departures from Freckles costs the detail
the archive keeps: the kriging neighbourhood (16/3 against 64/1) or the
vegetation cleanup. Those are the only two differences between the presets, so
varying them independently gives a 2x2 and the answer is in the corners rather
than in an argument.

The classified ground from an earlier run is reused, so a variant costs the
interpolation only -- about 45 minutes rather than three hours. That cloud is
the filter's output, before any cleanup, which is the right starting point for
every cell of the factorial.
"""
import argparse
import os
import shutil
import subprocess
import sys
import time
from dataclasses import replace

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import DEM_ROOT, LAS_ROOT, SRC, TESTS
from reference import PairingError, reference_for
sys.path.insert(0, SRC)

from replicalm import pipeline
from replicalm.config import PRESETS

# Defaults are the tile the cleanup work was developed on; --tile points the
# same comparison at another one, which is the only way to tell a parameter that
# generalises from one that was fitted to a single window.
DEFAULT_TILE = "AMIGACarb_Yuc_South_GLAS_Apr2013_l0s395"


def locate(tile):
    """LAS, archive reference and working directory for a tile stem.

    The reference is resolved on region, sensor and tile id together. It used to
    be globbed on the tile id alone -- `*l0s441*0p5m*.tif` -- and tile ids repeat
    across regions, so that was ambiguous for 223 of the 458 corpus tiles and
    returned whichever the filesystem listed first. Concretely: the Chiapas-1
    tile l0s422 has no reference of its own, and the old glob matched
    `OutYuc_NFI_l0s422` instead and would have scored it against a survey from
    a different part of the peninsula, silently and with a plausible number.

    The working directory still uses the bare id, which is safe because a run
    names one tile, but it means two same-id tiles from different regions share
    a directory. Pass REPLICALM_WORK to separate them if that ever matters.
    """
    import glob as _glob
    las = _glob.glob(os.path.join(LAS_ROOT, "**", tile + ".la[sz]"),
                     recursive=True)
    if not las:
        raise SystemExit("no LAS found for %s" % tile)
    try:
        ref = reference_for(tile, root=DEM_ROOT)
    except PairingError as exc:
        raise SystemExit("%s" % exc)
    sid = tile.split("_")[-1]
    base = os.path.join(TESTS, "clear_fulltile", sid)
    return las[0], ref, base

G1_PY = os.environ.get("REPLICALM_G1_PYTHON",
                       r"C:\Users\benja\anaconda3\envs\g1\python.exe")
RVT = os.environ.get("REPLICALM_RVT_SCRIPT", r"C:\g1\tools\GLiHT_rvt.py")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True)
    ap.add_argument("--tile", default=DEFAULT_TILE)
    ap.add_argument("--neighbours", type=int, default=16)
    ap.add_argument("--min-points", type=int, default=3)
    ap.add_argument("--clean", dest="clean", action="store_true")
    ap.add_argument("--no-clean", dest="clean", action="store_false")
    ap.set_defaults(clean=True)
    ap.add_argument("--slope-aware", action="store_true",
                    help="measure the cleanup against a fitted plane")
    ap.add_argument("--plane-patch", type=float, default=None,
                    help="radius the cleanup plane is fitted over")
    ap.add_argument("--skip-rvt", action="store_true")
    a = ap.parse_args(argv)

    TILE = a.tile
    LAS, REF, BASE = locate(TILE)
    GROUND = os.path.join(BASE, "_ground", TILE + "_ground.laz")
    out = os.path.join(BASE, a.tag)
    os.makedirs(out, exist_ok=True)
    os.makedirs(os.path.dirname(GROUND), exist_ok=True)
    cfg = replace(PRESETS["clear"], max_points=a.neighbours,
                  min_points=a.min_points, clean_vegetation=a.clean,
                  clean_slope_aware=a.slope_aware,
                  clean_plane_patch_m=a.plane_patch)
    print("%s: %d neighbours, min %d, cleanup %s%s, threshold %.2f m"
          % (a.tag, cfg.max_points, cfg.min_points,
             "ON" if cfg.clean_vegetation else "OFF",
             " (slope-aware)" if cfg.clean_slope_aware else "",
             cfg.passes[0].threshold_m), flush=True)

    # The first run on a tile classifies and banks the ground; every run after
    # reuses it, so a sweep over cleanup parameters costs interpolation only.
    dst = os.path.join(out, TILE + "_ground.laz")
    first = not os.path.exists(GROUND)
    if first:
        print("no banked ground for this tile: classifying once, then reusing",
              flush=True)
    elif not os.path.exists(dst):
        try:
            os.link(GROUND, dst)
        except OSError:
            shutil.copy2(GROUND, dst)

    def show(stage, message, fraction=None):
        print("  [%-11s] %s" % (stage, message), flush=True)

    t0 = time.time()
    res = pipeline.process(LAS, out, cfg=cfg, make_g1=False, keep_ground=True,
                           ground_stage="before", reuse_ground=not first,
                           on_grid=REF, progress=show)
    if first and os.path.exists(dst) and not os.path.exists(GROUND):
        try:
            os.link(dst, GROUND)
        except OSError:
            shutil.copy2(dst, GROUND)
        print("banked the classified ground at %s" % GROUND, flush=True)
    print("DEM in %.0f s: %s" % (time.time() - t0, res["dem"]), flush=True)
    if a.skip_rvt:
        return 0

    staged = os.path.join(out, "_rvt_in", "tile")
    os.makedirs(staged, exist_ok=True)
    link = os.path.join(staged, "%s_%s_DEM_0p5m_v1.tif"
                        % (TILE.split("_")[-1], a.tag))
    if not os.path.exists(link):
        shutil.copy2(res["dem"], link)
    env = {k: v for k, v in os.environ.items()
           if k not in ("PYTHONPATH", "PYTHONHOME", "GDAL_DATA", "PROJ_LIB",
                        "GDAL_DRIVER_PATH", "PROJ_DATA")}
    home = os.path.dirname(G1_PY)
    env["PATH"] = os.pathsep.join(
        [home, os.path.join(home, "Library", "bin"),
         os.path.join(home, "Library", "usr", "bin"),
         os.path.join(home, "Scripts"), env.get("PATH", "")])
    t0 = time.time()
    p = subprocess.run([G1_PY, RVT, "--res", "0p5m", "--no-metadata",
                        "--workers", "1",
                        "--dem-root", os.path.join(out, "_rvt_in"),
                        "--rvt-root", os.path.join(out, "rvt")],
                       env=env, text=True)
    print("RVT exit %d in %.0f s; products under %s"
          % (p.returncode, time.time() - t0, os.path.join(out, "rvt")),
          flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
