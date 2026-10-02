r"""G1 and its source layers for a chosen handful of finished DEMs.

The survey run produced bare-earth rasters only. The relief products are built
from the DEM, not from the point cloud, so they can be made afterwards without
re-classifying or re-kriging anything.

TWO THINGS THIS HAS TO GET RIGHT
--------------------------------
GLiHT_rvt.py discovers its input by globbing `*_DEM_{res}_v1.tif`, so the files
have to be staged under a name of that shape. The token is only a label: the
script reads the real cell size from the raster's own transform (line 125). Our
tiles derive their cell from measured point density and no two of these six
share one -- 0.423 to 0.615 m -- so staging under a single token is fine for the
run but would leave every output file claiming a resolution it does not have.
The outputs are therefore renamed afterwards to carry each tile's true cell.

The second thing is the environment. GLiHT_rvt.py runs under the interpreter
that has RVT, and it must not inherit this one's: PYTHONPATH points at
Replicalm's src, GDAL_DATA and PROJ_LIB point into ArcGIS, and PATH leads with
the ArcGIS binaries, which makes the child's GDAL fail with "No module named
'_gdal'".
"""
import argparse
import glob
import os
import re
import shutil
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import TESTS

RUN = os.path.join(TESTS, "full_clear")
STAGE = os.path.join(RUN, "_viz_stage", "tile")
RVT_OUT = os.path.join(RUN, "_viz_rvt")
VIZ = os.path.join(RUN, "viz")
TOKEN = "0p5m"                      # staging label only; see the note above

G1_PY = os.environ.get("REPLICALM_G1_PYTHON",
                       r"C:\Users\benja\anaconda3\envs\g1\python.exe")
RVT = os.environ.get("REPLICALM_RVT_SCRIPT", r"C:\g1\tools\GLiHT_rvt.py")

# The six asked for: two named, four drawn from Centro with a fixed seed so the
# selection is reproducible rather than merely arbitrary.
NAMED = ["AMIGACarb_Yuc_South_GLAS_Apr2013_l0s395_DEM.tif",
         "AMIGACarb_Yuc_Centro_NFI_Apr2013_l5s410_DEM.tif"]
SEED = 395410
N_RANDOM = 4
WANTED = ("G1", "Slope", "MultiHS", "OpnsPos", "SVF")


def cell_token(path):
    import rasterio
    with rasterio.open(path) as r:
        c = abs(r.transform.a)
    return ("%.3f" % c).replace(".", "p").rstrip("0").rstrip("p") + "m", c


def select():
    import random
    named = [os.path.join(RUN, n) for n in NAMED]
    for p in named:
        if not os.path.exists(p):
            raise SystemExit("missing requested DEM: %s" % p)
    centro = sorted(glob.glob(os.path.join(RUN, "*Centro*_DEM.tif")))
    pool = [p for p in centro if p not in named]
    random.seed(SEED)
    return named + random.sample(pool, N_RANDOM)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--keep-stage", action="store_true")
    a = ap.parse_args(argv)

    picks = select()
    for d in (STAGE, RVT_OUT, VIZ):
        os.makedirs(d, exist_ok=True)

    print("staging %d DEMs" % len(picks), flush=True)
    cells = {}
    for p in picks:
        stem = os.path.basename(p)[:-len("_DEM.tif")]
        tok, c = cell_token(p)
        cells[stem] = (tok, c)
        dst = os.path.join(STAGE, "%s_DEM_%s_v1.tif" % (stem, TOKEN))
        if not os.path.exists(dst):
            try:
                os.link(p, dst)          # same volume: costs nothing
            except OSError:
                shutil.copy2(p, dst)
        print("  %-54s cell %.3f m" % (stem[:54], c), flush=True)

    # the child must not inherit this interpreter's GDAL or PATH
    env = {k: v for k, v in os.environ.items()
           if k not in ("PYTHONPATH", "PYTHONHOME", "GDAL_DATA", "PROJ_LIB",
                        "GDAL_DRIVER_PATH", "PROJ_DATA")}
    home = os.path.dirname(G1_PY)
    env["PATH"] = os.pathsep.join(
        [home, os.path.join(home, "Library", "bin"),
         os.path.join(home, "Library", "usr", "bin"),
         os.path.join(home, "Scripts"), env.get("PATH", "")])

    t0 = time.time()
    print("\nrunning RVT on %d tiles, %d workers" % (len(picks), a.workers),
          flush=True)
    proc = subprocess.run(
        [G1_PY, RVT, "--res", TOKEN, "--no-metadata",
         "--workers", str(a.workers),
         "--dem-root", os.path.dirname(STAGE),
         "--rvt-root", RVT_OUT],
        env=env, text=True)
    print("RVT exit %d after %.0f s" % (proc.returncode, time.time() - t0),
          flush=True)

    # collect, renaming the staging token to each tile's real cell size
    made = []
    for root, _dirs, files in os.walk(RVT_OUT):
        for f in files:
            if not f.endswith(".tif"):
                continue
            m = re.match(r"(.+?)_(%s)_%s_v1\.tif$" % ("|".join(WANTED), TOKEN), f)
            if not m:
                continue
            stem, vis = m.group(1), m.group(2)
            tok = cells.get(stem, (TOKEN, 0))[0]
            out = os.path.join(VIZ, "%s_%s_%s_v1.tif" % (stem, vis, tok))
            shutil.copy2(os.path.join(root, f), out)
            made.append(out)

    print("\n%d products written to %s" % (len(made), VIZ))
    by = {}
    for p in made:
        by.setdefault(re.sub(r"_(%s)_.*$" % "|".join(WANTED), "",
                             os.path.basename(p)), []).append(
            re.search(r"_(%s)_" % "|".join(WANTED), os.path.basename(p)).group(1))
    for k in sorted(by):
        print("  %-54s %s" % (k[:54], " ".join(sorted(by[k]))))
    if not a.keep_stage:
        shutil.rmtree(os.path.dirname(STAGE), ignore_errors=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
