r"""Do the noise filters do anything on steep ground, with pass count held fixed?

THE UNRESOLVED CLAIM
--------------------
Two measurements point in opposite directions and both cannot be right.

The trench comparison showed two SMRF passes with the noise filters on scoring
0.083 m RMSE and 0.25% terracing, against one pass with the filters off at
0.161 m and 2.44%. Two variables moved at once. The recalibration then held the
filters on and varied only the pass count across seven tiles, and found the two
indistinguishable -- steep support 98.1% either way, error above 30 degrees
10.34% against 10.18%. By elimination that left the filters as the cause.

But the flat micro-sweep found ELM completely inert at every threshold tried,
removing nothing, and the statistical outlier filter marginally harmful. If the
filters do nothing, they cannot be what the trench comparison measured.

WHAT THIS RUN DOES
------------------
Varies only the filters, with one SMRF pass everywhere, on the three windows
that actually carry relief: l0s395 at 33.2% of cells above twenty degrees,
l0s417 at 13.1% and l2s444 at 11.3%. Four arms, one variable.

Either the filters earn their place on steep ground -- in which case the flat
sweep simply showed that flat ground does not need them, which is unsurprising
-- or they do not, and something else explains the trench result, in which case
the default should be off and the trench comparison needs re-running arm by arm.
"""
import json, math, os, sys, time
from dataclasses import replace
sys.path.insert(0, r"C:\Replicalm\src")
import numpy as np
from scipy.spatial import cKDTree
from osgeo import gdal
from replicalm import classify, grid as G, kriging as K, interpolate
from replicalm.config import PRESETS, GroundPass
gdal.UseExceptions()

CLIPS = r"C:\Replicalm\tests\recal\clips"
OUT = r"C:\Replicalm\tests\filt"
os.makedirs(OUT, exist_ok=True)
TILES = ["l0s395", "l0s417", "l2s444"]      # 33.2%, 13.1%, 11.3% steep

ONE_PASS = [GroundPass(algorithm="smrf",
                       slope=round(math.tan(math.radians(9.0)), 4),
                       threshold_m=0.5, cell_m=0.5, window_m=70.0)]

ARMS = [
    ("both_off",     dict(remove_low_noise=False, remove_outliers=False)),
    ("elm_only",     dict(remove_low_noise=True,  remove_outliers=False)),
    ("outlier_only", dict(remove_low_noise=False, remove_outliers=True)),
    ("both_on",      dict(remove_low_noise=True,  remove_outliers=True)),
]

results = []
for tile in TILES:
    las = os.path.join(CLIPS, "%s_clip.las" % tile)
    ref_p = os.path.join(CLIPS, "%s_reference.tif" % tile)
    if not (os.path.exists(las) and os.path.exists(ref_p)):
        print("%s: clip missing" % tile)
        continue
    g = G.grid_from_raster(ref_p)
    wkt = interpolate.source_srs(las)
    area = (g.bounds[2] - g.bounds[0]) * (g.bounds[3] - g.bounds[1])
    rds = gdal.Open(ref_p); rb = rds.GetRasterBand(1)
    ref = rb.ReadAsArray().astype("f8"); nod = rb.GetNoDataValue()
    rm = np.isfinite(ref) & (ref != 0.0)
    if nod is not None:
        rm &= ref != nod
    gyy, gxx = np.gradient(np.where(rm, ref, np.nan), g.dy, g.cell)
    slope = np.degrees(np.arctan(np.hypot(gxx, gyy)))
    steep_ref = rm & (slope >= 20)
    gx_, gy_ = g.cell_centres(); flat = np.c_[gx_.ravel(), gy_.ravel()]

    def terrace(z, m):
        gy, gx = np.gradient(np.where(m, z, np.nan), g.dy, g.cell)
        s = np.degrees(np.arctan(np.hypot(gx, gy)))
        k = m & (slope >= 20) & np.isfinite(s)
        return float((s[k] < 0.25 * slope[k]).mean()) if k.sum() > 200 else float("nan")

    print("\n%s %s   %.1f%% of cells above 20 deg"
          % ("=" * 60, tile, 100 * steep_ref.sum() / rm.sum()))
    print("  %-14s %8s %8s %9s %9s %9s %9s"
          % ("arm", "gnd/m2", "removed", "support", "over30", "terrace", "rmse"))
    base_density = None
    for name, params in ARMS:
        t0 = time.time()
        cfg = replace(PRESETS["ncalm"], passes=list(ONE_PASS), **params)
        gp = os.path.join(OUT, "%s_%s.las" % (tile, name))
        try:
            if not os.path.exists(gp):
                classify.classify_tile(las, gp, cfg, verbose=False)
            arr, _ = classify.read_points(gp)
            gnd = arr[arr["Classification"] == 2]
            x, y, z = (gnd["X"].astype("f8"), gnd["Y"].astype("f8"),
                       gnd["Z"].astype("f8"))
            dens = len(z) / area
            if base_density is None:
                base_density = dens
            removed = 100.0 * (1.0 - dens / base_density)
            ng = np.array(cKDTree(np.c_[x, y]).query_ball_point(
                flat, 1.5, return_length=True)).reshape(g.height, g.width)
            support = (float((ng[steep_ref] >= 3).mean())
                       if steep_ref.sum() else float("nan"))
            radius = K.radius_for_density(dens, 16, 4.0, 3.0, 20.0)
            v = K.fit_variogram(x, y, z, model="spherical", max_lag=10.0)
            dem, _ = K.krige_grid(x, y, z, g, radius=radius, max_points=16,
                                  variogram=v, verbose=False)
            K.write_geotiff(dem, g, wkt, os.path.join(OUT, "%s_%s.tif" % (tile, name)))
            m = rm & (dem > -9000) & np.isfinite(slope)
            e = np.abs(dem - ref)
            b30 = m & (slope >= 30)
            row = {"tile": tile, "arm": name, "density": dens,
                   "removed_pct": removed, "steep_support": support,
                   "over30": float((e[b30] > 0.5).mean()) if b30.sum() > 200 else None,
                   "terrace": terrace(dem, m),
                   "rmse": float(np.sqrt((e[m] ** 2).mean())),
                   "seconds": round(time.time() - t0, 1)}
            results.append(row)
            print("  %-14s %8.2f %7.2f%% %8.1f%% %8s %8.2f%% %9.4f"
                  % (name, dens, removed, 100 * support,
                     "n/a" if row["over30"] is None else "%.2f%%" % (100 * row["over30"]),
                     100 * row["terrace"], row["rmse"]))
        except Exception as ex:
            print("  %-14s FAILED %s" % (name, str(ex)[:46]))

with open(os.path.join(OUT, "filters.json"), "w", encoding="utf-8") as fh:
    json.dump(results, fh, indent=1, default=float)

print("\n\n%s\nAGGREGATE over the three relief windows\n%s" % ("=" * 70, "=" * 70))
print("%-14s %9s %9s %10s %10s %10s"
      % ("arm", "removed", "support", "over30", "terrace", "rmse"))
print("-" * 66)
for name, _ in ARMS:
    rs = [r for r in results if r["arm"] == name]
    if not rs:
        continue
    o = [r["over30"] for r in rs if r["over30"] is not None]
    print("%-14s %8.2f%% %8.1f%% %9s %9.2f%% %10.4f"
          % (name, float(np.mean([r["removed_pct"] for r in rs])),
             100 * np.nanmean([r["steep_support"] for r in rs]),
             "n/a" if not o else "%.2f%%" % (100 * np.mean(o)),
             100 * np.nanmean([r["terrace"] for r in rs]),
             float(np.mean([r["rmse"] for r in rs]))))
print("\n'removed' is the share of ground points the arm dropped relative to")
print("both_off on the same window. A filter removing 0.00%% is inert.")
