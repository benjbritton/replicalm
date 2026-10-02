r"""Is Clear's cleanup what costs the shadow detail?

Measured on South GLAS l0s395: our composite renders the same ground 18 levels
darker than the archive's and carries 2.7x less texture in its darkest band,
while the DEMs agree to a median of 2 mm and their micro-relief is identical.
The difference traces to the positive-openness layer, whose tile minimum is 39
in ours against 18 in the archive -- enough that Table 3's brightness of -30
clips our layer to zero where the archive's survives.

This asks whether the cleanup is responsible, by building the same tile with it
switched off and nothing else changed.

ONE VARIABLE. Clear's own configuration with `clean_vegetation` set False --
not the Freckles preset, which also carries 64 neighbours against Clear's 16
and would make this a two-variable test.

The classified cloud from the on-grid run is reused, so this pays for the
interpolation only. That cloud is the filter's output, before the cleanup, which
is exactly the input wanted here.
"""
import argparse
import os
import shutil
import subprocess
import sys
import time
from dataclasses import replace

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import SRC, TESTS
sys.path.insert(0, SRC)

from replicalm import pipeline
from replicalm.config import PRESETS

TILE = "AMIGACarb_Yuc_South_GLAS_Apr2013_l0s395"
LAS = (r"D:\GLiHT_LAS_orig\Yuc_South\South_Glas\South_Glas\%s.las" % TILE)
REF = (r"D:\_Archive_EdgeFixed\Yuc_South\South_Glas\South_Glas"
       r"\South_GLAS_l0s395_DEM_0p5m_v1.tif")
SRC_RUN = os.path.join(TESTS, "clear_fulltile", "l0s395", "on_archive_grid")
OUT = os.path.join(TESTS, "clear_fulltile", "l0s395", "no_cleanup")

G1_PY = os.environ.get("REPLICALM_G1_PYTHON",
                       r"C:\Users\benja\anaconda3\envs\g1\python.exe")
RVT = os.environ.get("REPLICALM_RVT_SCRIPT", r"C:\g1\tools\GLiHT_rvt.py")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-dem", action="store_true")
    a = ap.parse_args(argv)
    os.makedirs(OUT, exist_ok=True)

    cfg = replace(PRESETS["clear"], clean_vegetation=False)
    print("Clear, cleanup OFF: %d neighbours, min %d, threshold %.2f m"
          % (cfg.max_points, cfg.min_points, cfg.passes[0].threshold_m),
          flush=True)

    dem = os.path.join(OUT, TILE + "_DEM.tif")
    if not (a.skip_dem and os.path.exists(dem)):
        src_g = os.path.join(SRC_RUN, TILE + "_ground.laz")
        dst_g = os.path.join(OUT, TILE + "_ground.laz")
        if not os.path.exists(dst_g):
            if not os.path.exists(src_g):
                raise SystemExit("no classified ground to reuse at %s" % src_g)
            try:
                os.link(src_g, dst_g)
            except OSError:
                shutil.copy2(src_g, dst_g)
            print("reusing the classified ground from the on-grid run",
                  flush=True)

        def show(stage, message, fraction=None):
            print("  [%-11s] %s" % (stage, message), flush=True)

        t0 = time.time()
        out = pipeline.process(LAS, OUT, cfg=cfg, make_g1=False,
                               keep_ground=True, ground_stage="before",
                               reuse_ground=True, on_grid=REF, progress=show)
        print("DEM in %.0f s: %s" % (time.time() - t0, out["dem"]), flush=True)
        dem = out["dem"]

    # RVT, in the interpreter that has it and without this one's environment
    staged = os.path.join(OUT, "_rvt_in", "tile")
    os.makedirs(staged, exist_ok=True)
    link = os.path.join(staged, "l0s395_noclean_DEM_0p5m_v1.tif")
    if not os.path.exists(link):
        shutil.copy2(dem, link)
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
                        "--dem-root", os.path.join(OUT, "_rvt_in"),
                        "--rvt-root", os.path.join(OUT, "rvt")],
                       env=env, text=True)
    print("RVT exit %d in %.0f s" % (p.returncode, time.time() - t0), flush=True)
    print("products under %s" % os.path.join(OUT, "rvt"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
