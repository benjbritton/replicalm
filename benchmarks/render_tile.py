r"""Process one whole G-LiHT tile with this pipeline, for visual comparison.

WHY A BLOCKED DRIVER
--------------------
Everything measured so far ran on 400 m clips of 640,000 cells. South_GLAS_l0s395
is 2890 x 15950 cells at 0.5 m -- 46.1 million, seventy-two times larger -- and
neither stage survives being handed the whole thing at once:

  classification  CSF simulates a cloth over the entire cloud it is given. The
                  full tile is roughly 23 million returns.
  interpolation   krige_grid queries its KD-tree for every cell in one call. At
                  46 million cells and 16 neighbours that is about 12 GB of
                  distances and indices before a single cell is solved.

So the tile is processed in 1 km blocks with a 10 m buffer, which is what the
source method prescribes at steps 1 and 2 rather than a workaround invented
here. The buffer is widened to the search radius for the interpolation, since a
cell at a block edge needs neighbours from across it; points outside the block
proper contribute to its solve and are then discarded, so no cell is solved from
a truncated neighbourhood and no point is counted twice.

THE GRID IS THE REFERENCE'S
---------------------------
Output lands on the reference DEM's exact lattice, irregular cells and all, via
grid_from_raster. Production would use snap_outward; this run is for comparison
against the archive G1, and a comparison is only cell-for-cell if both rasters
sit on the same cells.

THE RADIUS IS A CONSTANT, AND IS NOT DERIVED FROM THE TERRAIN
-------------------------------------------------------------
The refit across ten tiles found no correlation length on any of them: fitted
range tracks the fitting window at 0.30 to 1.44 times its width, and the sill
grows with the window by 3x to 140x, so there is no plateau to read a range off.
A variogram is still fitted here, once, on a subsample, and reused for every
block -- fitting per block would let the radius drift between blocks and print
the block boundaries into the output. It is recorded as what it is: a constant
chosen by calibration, not a quantity measured from this tile.

CLASSIFIER
----------
CSF with rigidness 2 at a 0.5 m threshold, which is the setting calibration
chose for this tile specifically (0.029 m RMSE, 1.01x complexity). The pilot's
broader default, smrf at a 0.25 m threshold, lost here.
"""
import json, os, sys, time
from dataclasses import replace
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import ROOT, SRC, DEM_ROOT, LAS_ROOT, TESTS
sys.path.insert(0, SRC)
import numpy as np
from replicalm import classify, clip, grid as G, kriging as K, interpolate
from replicalm.config import PRESETS, GroundPass

TILE = "l0s395"
LAS = (LAS_ROOT + r"\\Yuc_South\South_Glas\South_Glas"
       r"\AMIGACarb_Yuc_South_GLAS_Apr2013_l0s395.las")
REF = (DEM_ROOT + r"\\Yuc_South\South_Glas\South_Glas"
       r"\South_GLAS_l0s395_DEM_0p5m_v1.tif")
SUFFIX = sys.argv[2] if len(sys.argv) > 2 else ""
OUT = ROOT + r"\\render" + SUFFIX
WORK = os.path.join(OUT, "work")
# the published G1s live under <root>\<Region>\<Sub>\<Sub>\, and GLiHT_rvt.py
# mirrors whatever structure it finds, so the DEM is written into the same shape
DEM_OUT = os.path.join(OUT, "dem", "Yuc_South", "South_Glas", "South_Glas",
                       "South_GLAS_l0s395_DEM_0p5m_v1.tif")

BLOCK_M = 1000.0        # source method step 1
BUFFER_M = 10.0         # source method step 2
MAX_LAG = 10.0          # only the variogram SHAPE is taken from this
# The search radius is set here, not derived. The refit across ten tiles found
# no correlation length to derive it from, and a radius ladder on this tile
# measured what the derived value cost: at 14.4 m the surface chords across
# slope breaks and extrapolates past the data boundary, and window RMSE falls
# from 0.551 to 0.065 m as the radius comes down to 1.5 m. Below that there are
# not 16 neighbours to find at 4 points per square metre, and coverage fails
# instead of accuracy.
RADIUS_M = float(sys.argv[1]) if len(sys.argv) > 1 else 1.5
# Nugget floor as a fraction of sill, argv[3]. Zero reproduces the locked
# baseline. Four per cent is where flat-ground roughness matched the reference
# on the l0s395 trench window; see kriging.fit_variogram.
NUGGET_FRACTION = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0
MIN_GROUND = 500
for d in (WORK, os.path.dirname(DEM_OUT)):
    os.makedirs(d, exist_ok=True)

g = G.grid_from_raster(REF)
wkt = interpolate.source_srs(LAS)
print("tile grid: %s" % g.describe())
print("%d cells total\n" % (g.width * g.height))

# cells per block, from the reference's own cell size
bw = max(1, int(round(BLOCK_M / g.cell)))
bh = max(1, int(round(BLOCK_M / g.dy)))
# The source's structure, restored. Calibration on flat 400 m clips had led to
# CSF, to one pass instead of two, and to the noise filters switched off,
# because none of those choices cost anything on flat ground. On the trench
# window at 1500,5250 they cost a great deal: two SMRF passes at the source's
# 9 and 12 degrees, with noise removal on, scored 0.083 m RMSE and 0.25%
# terracing against the single pass's 0.161 m and 2.44%.
#
# `cell` is the one parameter the source cannot specify, since TerraScan has no
# SMRF working grid. Left at 1 m it quantizes steep faces into steps while the
# output is rendered at 0.5 m; tied to the output cell the terracing halves.
import math
_CELL = 0.5
cfg = replace(PRESETS["ncalm"], passes=[
    GroundPass(algorithm="smrf", slope=round(math.tan(math.radians(9.0)), 4),
               threshold_m=0.5, cell_m=_CELL, window_m=70.0,
               stands_in_for="TerraScan pass 1: 70 m, 88 deg, 9 deg, 3 m"),
    GroundPass(algorithm="smrf", slope=round(math.tan(math.radians(12.0)), 4),
               threshold_m=0.5, cell_m=_CELL, window_m=70.0, reuse_ground=True,
               stands_in_for="TerraScan pass 2: 70 m, 88 deg, 12 deg, 3 m"),
])

from osgeo import gdal
gdal.UseExceptions()
ref_ds = gdal.Open(REF)
ref_band = ref_ds.GetRasterBand(1)
ref_nod = ref_band.GetNoDataValue()

dem = np.full((g.height, g.width), -9999.0, dtype="f4")
variogram = None
radius = None
stats = []
t_start = time.time()
blocks = [(i, j) for j in range(0, g.height, bh) for i in range(0, g.width, bw)]
print("%d blocks of up to %d x %d cells\n" % (len(blocks), bw, bh))

for n, (i, j) in enumerate(blocks, 1):
    w = min(bw, g.width - i)
    h = min(bh, g.height - j)
    sub = G.Grid(g.origin_x + i * g.cell, g.origin_y - j * g.dy,
                 g.cell, w, h, cell_y=g.cell_y)

    # a block the reference never covered has nothing to compare against and
    # nothing to interpolate; skipping it is most of the saving on this tile,
    # which is 73% empty bounding box
    ref_blk = ref_band.ReadAsArray(i, j, w, h).astype("f8")
    valid = np.isfinite(ref_blk) & (ref_blk != 0.0)
    if ref_nod is not None:
        valid &= ref_blk != ref_nod
    if valid.sum() < 100:
        print("[%2d/%d] %5d,%-5d empty in the reference, skipped" % (n, len(blocks), i, j))
        continue

    mnx, mny, mxx, mxy = sub.bounds
    pad = BUFFER_M + (radius if radius else 25.0)
    blk_las = os.path.join(WORK, "%s_%05d_%05d.las" % (TILE, i, j))
    gnd_las = os.path.join(WORK, "%s_%05d_%05d_gnd.las" % (TILE, i, j))
    t0 = time.time()
    try:
        if not os.path.exists(gnd_las):
            clip.clip_las(LAS, blk_las,
                          {"minx": mnx - pad, "miny": mny - pad,
                           "maxx": mxx + pad, "maxy": mxy + pad},
                          verbose=False)
            classify.classify_tile(blk_las, gnd_las, cfg, verbose=False)
        arr, _ = classify.read_points(gnd_las)
        gnd = arr[arr["Classification"] == 2]
        if len(gnd) < MIN_GROUND:
            print("[%2d/%d] %5d,%-5d only %d ground points, left as nodata"
                  % (n, len(blocks), i, j, len(gnd)))
            continue
        x, y, z = (gnd["X"].astype("f8"), gnd["Y"].astype("f8"),
                   gnd["Z"].astype("f8"))

        # one variogram for the whole tile, fitted on the first block that has
        # enough ground to fit on, then reused everywhere
        if variogram is None and len(z) >= 20000:
            sel = np.random.default_rng(20260919).choice(len(z), 20000,
                                                         replace=False)
            variogram = K.fit_variogram(x[sel], y[sel], z[sel],
                                        model="spherical", max_lag=MAX_LAG,
                                        min_nugget_fraction=NUGGET_FRACTION)
            radius = RADIUS_M
            print("    variogram fitted on block %d,%d: range %.2f m; "
                  "search radius set to %.2f m (not derived), every block"
                  % (i, j, variogram["params"][2], radius))
            if NUGGET_FRACTION:
                print("    nugget floored at %.0f%% of sill -> %.5f (kriging "
                      "smooths instead of interpolating exactly)"
                      % (100 * NUGGET_FRACTION, variogram["params"][0]))
        if variogram is None:
            print("[%2d/%d] %5d,%-5d too few points to fit on yet, deferred"
                  % (n, len(blocks), i, j))
            continue

        blk, info = K.krige_grid(x, y, z, sub, radius=radius, max_points=16,
                                 variogram=variogram, verbose=False)
        dem[j:j + h, i:i + w] = blk.astype("f4")
        filled = float((blk > -9000).mean())
        stats.append({"i": i, "j": j, "ground": int(len(gnd)),
                      "filled": filled, "fallback": info["fallback_fraction"],
                      "seconds": round(time.time() - t0, 1)})
        print("[%2d/%d] %5d,%-5d %8d ground  %5.1f%% filled  %4.1f%% fallback  "
              "%5.0f s  (%.0f min elapsed)"
              % (n, len(blocks), i, j, len(gnd), 100 * filled,
                 100 * info["fallback_fraction"], time.time() - t0,
                 (time.time() - t_start) / 60))
    except Exception as e:
        print("[%2d/%d] %5d,%-5d FAILED %s" % (n, len(blocks), i, j, str(e)[:60]))
    finally:
        if os.path.exists(blk_las):
            try:
                os.remove(blk_las)
            except OSError:
                pass

K.write_geotiff(dem, g, wkt, DEM_OUT)
with open(os.path.join(OUT, "render_%s.json" % TILE), "w", encoding="utf-8") as fh:
    json.dump({"tile": TILE, "grid": g.describe(), "radius_m": radius,
               "variogram": {k: v for k, v in (variogram or {}).items()
                             if k in ("model", "params")},
               "blocks": stats}, fh, indent=1, default=float)

covered = float((dem > -9000).mean())
print("\nwrote %s" % DEM_OUT)
print("%.1f%% of the tile grid carries a value, %.0f minutes total"
      % (100 * covered, (time.time() - t_start) / 60))
