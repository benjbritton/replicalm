r"""Re-run the classifier calibration on ground that can tell the arms apart.

WHY THIS EXISTS
---------------
Every earlier calibration window was chosen for maximum reference coverage,
which finds flat, open ground. Measured after the fact, those windows contained
between 0.00% and 0.20% cells steeper than twenty degrees. The sweeps that chose
CSF over SMRF, one ground pass over the source's two, and noise filtering off
were therefore scored on terrain where all of those choices are equivalent.

On a window with relief, they are not equivalent at all: the same three choices
cost 40.38% of cells above thirty degrees past half a metre, against 2.25% for
the source's structure, and produced the bead and terrace artifacts visible in
the G1.

WHAT IS SCORED
--------------
RMSE is reported but is not the criterion, because flat ground dominates it and
flat ground is where the arms agree. The criteria are:

  steep_support   share of cells above 20 degrees with at least three ground
                  returns within 1.5 m. This is the classification question:
                  did the filter keep the flank at all. CSF scored 67.8% here
                  against SMRF's 100%.
  over30          share of cells above 30 degrees whose elevation is more than
                  half a metre from the reference.
  terrace         share of steep cells whose own gradient has collapsed below a
                  quarter of the reference's -- a step where the reference has
                  a face. The archive scores 0.00%; a 1 m SMRF working grid
                  scored 5.06%.

THE ARMS
--------
  ncalm_2pass   the corrected default: the source's two passes at 9 and 12
                degrees, 0.5 m threshold, working grid tied to the output cell,
                noise filters on.
  smrf_1pass    one pass, otherwise the same. Isolates the second pass.
  csf_r2        what the flat calibration chose. Included so the size of that
                error is on the record per tile rather than asserted.
  ncalm_asdoc   the source's parameters verbatim, 3 m threshold on a 1 m grid.
                Shows which parts of the translation are load-bearing.
"""
import json, math, os, sys, time
from dataclasses import replace
sys.path.insert(0, r"C:\Replicalm\src")
import numpy as np
from scipy.spatial import cKDTree
from osgeo import gdal
from replicalm import calibrate, classify, clip, grid as G, kriging as K, interpolate
from replicalm.config import PRESETS, GroundPass
gdal.UseExceptions()

OUT = r"C:\Replicalm\tests\recal"
CLIPS = os.path.join(OUT, "clips")
os.makedirs(CLIPS, exist_ok=True)
WIN_M = 400.0
MIN_COVERAGE = 0.70
RADIUS = 20.0


def two_pass(cell=0.5, thr=0.5):
    return [GroundPass(algorithm="smrf", slope=round(math.tan(math.radians(9.0)), 4),
                       threshold_m=thr, cell_m=cell, window_m=70.0),
            GroundPass(algorithm="smrf", slope=round(math.tan(math.radians(12.0)), 4),
                       threshold_m=thr, cell_m=cell, window_m=70.0,
                       reuse_ground=True)]


ARMS = [
    ("ncalm_2pass", replace(PRESETS["ncalm"], passes=two_pass())),
    ("smrf_1pass", replace(PRESETS["ncalm"], passes=two_pass()[:1])),
    ("csf_r2", replace(PRESETS["ncalm"],
                       passes=[GroundPass(algorithm="csf", csf_rigidness=2,
                                          csf_threshold_m=0.5)],
                       remove_low_noise=False, remove_outliers=False)),
    ("ncalm_asdoc", replace(PRESETS["ncalm"], passes=two_pass(cell=1.0, thr=3.0))),
]

tiles = json.load(open(r"C:\Replicalm\tests\calibration_tiles.json"))["tiles"]
tiles += json.load(open(r"C:\Replicalm\tests\pilot_tiles.json"))["tiles"]

results = []
for t in tiles:
    tile = t["tile"]
    try:
        win = clip.best_window(t["reference_dem"], WIN_M,
                               min_coverage=MIN_COVERAGE, prefer_relief=True)
    except Exception as e:
        print("\n%s: %s" % (tile, str(e)[:70]))
        continue
    las_p = os.path.join(CLIPS, "%s_clip.las" % tile)
    ref_p = os.path.join(CLIPS, "%s_reference.tif" % tile)
    print("\n%s %s" % ("=" * 78, tile))
    print("  window %.1f%% covered, %.2f%% of cells above 20 deg"
          % (100 * win["coverage"], 100 * win["steep_fraction"]))
    try:
        if not os.path.exists(las_p):
            clip.clip_las(t["las"], las_p, win, verbose=False)
        if not os.path.exists(ref_p):
            clip.clip_raster(t["reference_dem"], ref_p, win, verbose=False)
    except Exception as e:
        print("  clip failed: %s" % str(e)[:60])
        continue

    g = G.grid_from_raster(ref_p)
    wkt = interpolate.source_srs(las_p)
    rds = gdal.Open(ref_p); rb = rds.GetRasterBand(1)
    ref = rb.ReadAsArray().astype("f8"); nod = rb.GetNoDataValue()
    rm = np.isfinite(ref) & (ref != 0.0)
    if nod is not None:
        rm &= ref != nod
    gyy, gxx = np.gradient(np.where(rm, ref, np.nan), g.dy, g.cell)
    slope = np.degrees(np.arctan(np.hypot(gxx, gyy)))
    steep_ref = rm & (slope >= 20)
    if steep_ref.sum() < 500:
        print("  only %d steep cells; this tile cannot discriminate" % steep_ref.sum())
    gx_, gy_ = g.cell_centres(); flat = np.c_[gx_.ravel(), gy_.ravel()]

    def terrace(z, m):
        gy, gx = np.gradient(np.where(m, z, np.nan), g.dy, g.cell)
        s = np.degrees(np.arctan(np.hypot(gx, gy)))
        k = m & (slope >= 20) & np.isfinite(s)
        return float((s[k] < 0.25 * slope[k]).mean()) if k.sum() > 200 else float("nan")

    print("  %-13s %9s %9s %9s %9s %9s" %
          ("arm", "gnd/m2", "support", "over30", "terrace", "rmse"))
    for name, cfg in ARMS:
        t0 = time.time()
        gp = os.path.join(OUT, "%s_%s.las" % (tile, name))
        try:
            if not os.path.exists(gp):
                classify.classify_tile(las_p, gp, cfg, verbose=False)
            arr, _ = classify.read_points(gp)
            gnd = arr[arr["Classification"] == 2]
            if len(gnd) < 500:
                print("  %-13s only %d ground points" % (name, len(gnd)))
                continue
            x, y, z = (gnd["X"].astype("f8"), gnd["Y"].astype("f8"),
                       gnd["Z"].astype("f8"))
            area = (g.bounds[2] - g.bounds[0]) * (g.bounds[3] - g.bounds[1])
            ng = np.array(cKDTree(np.c_[x, y]).query_ball_point(
                flat, 1.5, return_length=True)).reshape(g.height, g.width)
            support = (float((ng[steep_ref] >= 3).mean())
                       if steep_ref.sum() else float("nan"))
            v = K.fit_variogram(x, y, z, model="spherical", max_lag=10.0)
            dem, _ = K.krige_grid(x, y, z, g, radius=RADIUS, max_points=16,
                                  variogram=v, verbose=False)
            K.write_geotiff(dem, g, wkt, os.path.join(OUT, "%s_%s.tif" % (tile, name)))
            m = rm & (dem > -9000) & np.isfinite(slope)
            e = np.abs(dem - ref)
            b30 = m & (slope >= 30)
            row = {"tile": tile, "arm": name,
                   "steep_fraction": win["steep_fraction"],
                   "density": len(z) / area, "steep_support": support,
                   "over30": float((e[b30] > 0.5).mean()) if b30.sum() > 200 else None,
                   "terrace": terrace(dem, m),
                   "rmse": float(np.sqrt((e[m] ** 2).mean())),
                   "coverage": float(m.sum() / rm.sum()),
                   "seconds": round(time.time() - t0, 1)}
            results.append(row)
            print("  %-13s %9.2f %8.1f%% %8s %8.2f%% %9.4f"
                  % (name, row["density"], 100 * support,
                     "n/a" if row["over30"] is None else "%.2f%%" % (100 * row["over30"]),
                     100 * row["terrace"], row["rmse"]))
        except Exception as ex:
            print("  %-13s FAILED %s" % (name, str(ex)[:50]))

with open(os.path.join(OUT, "recal.json"), "w", encoding="utf-8") as fh:
    json.dump(results, fh, indent=1, default=float)

print("\n\n%s\nAGGREGATE across tiles with steep ground\n%s" % ("=" * 78, "=" * 78))
print("%-13s %8s %10s %10s %10s %10s" %
      ("arm", "tiles", "support", "over30", "terrace", "rmse"))
print("-" * 68)
for name, _ in ARMS:
    rs = [r for r in results if r["arm"] == name and r["steep_fraction"] > 0.01]
    if not rs:
        continue
    o30 = [r["over30"] for r in rs if r["over30"] is not None]
    print("%-13s %8d %9.1f%% %9s %9.2f%% %10.4f"
          % (name, len(rs), 100 * np.nanmean([r["steep_support"] for r in rs]),
             "n/a" if not o30 else "%.2f%%" % (100 * np.mean(o30)),
             100 * np.nanmean([r["terrace"] for r in rs]),
             float(np.mean([r["rmse"] for r in rs]))))
print("\nwrote %s" % os.path.join(OUT, "recal.json"))
