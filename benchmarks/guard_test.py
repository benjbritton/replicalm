r"""Is the bead artifact on steep slopes caused by the fallback guard?

THE CLAIM IN THE CODE
---------------------
krige_grid rejects a kriging estimate that falls outside the range of its own
neighbours, on the stated grounds that "an ordinary kriging estimate interpolates
its neighbours and cannot legitimately fall outside their range".

That is not a property of ordinary kriging. The weights sum to one but are not
constrained to be positive, and negative weights are ordinary -- the screen
effect, where a near point shields a farther one, produces them routinely. A
cell near the top of a slope, surrounded by neighbours that are mostly downhill,
SHOULD estimate above all of them. Ordinary kriging is not a convex combination
and extrapolation within a neighbourhood is not evidence of a degenerate system.

WHY THAT WOULD PRODUCE EXACTLY THIS ARTIFACT
--------------------------------------------
Inverse distance IS a convex combination: its estimate is always within the
range of the points it averages. So every cell the guard diverts to the fallback
gets clamped to its neighbours' range. On a slope break that means a run of
cells all pinned to the same local extreme -- a flat bead. Strung along the
contour where the guard keeps firing, that is the necklace visible in detail2.

The prediction is specific and falsifiable: the diverted cells should be
concentrated on steep ground, and relaxing the bound should cut the error there
without hurting the flat ground, where the guard rarely fires and the two paths
agree anyway.

THE THREE ARMS
--------------
  current    reject when the estimate leaves [min, max] of the neighbours
  unguarded  reject only a non-finite estimate or a failed solve
  widened    reject when the estimate leaves the neighbour range by more than
             the neighbourhood's own spread, which still catches the 1e15
             blow-ups the guard was written for while allowing the honest
             extrapolation a gradient requires

Scored against the archive DEM, overall and by slope band, on the window the
artifact was spotted in.
"""
import json, os, sys
sys.path.insert(0, r"C:\Replicalm\src")
import numpy as np
from scipy.spatial import cKDTree
from replicalm import classify, grid as G, kriging as K, interpolate

GND = r"C:\Replicalm\render\work\l0s395_02000_00000_gnd.las"
REF = (r"D:\_Archive_EdgeFixed\Yuc_South\South_Glas\South_Glas"
       r"\South_GLAS_l0s395_DEM_0p5m_v1.tif")
LAS = (r"D:\GLiHT_LAS_orig\Yuc_South\South_Glas\South_Glas"
       r"\AMIGACarb_Yuc_South_GLAS_Apr2013_l0s395.las")
OUT = r"C:\Replicalm\render\guard"
os.makedirs(OUT, exist_ok=True)
X0, Y0, NW, NH = 2100, 0, 700, 700       # the detail2 window
RADIUS, MAXP, MINP = 14.4, 16, 3

full = G.grid_from_raster(REF)
sub = G.Grid(full.origin_x + X0 * full.cell, full.origin_y - Y0 * full.dy,
             full.cell, NW, NH, cell_y=full.cell_y)
wkt = interpolate.source_srs(LAS)

arr, _ = classify.read_points(GND)
g = arr[arr["Classification"] == 2]
x, y, z = g["X"].astype("f8"), g["Y"].astype("f8"), g["Z"].astype("f8")
mnx, mny, mxx, mxy = sub.bounds
pad = RADIUS + 2
k = ((x >= mnx - pad) & (x <= mxx + pad) & (y >= mny - pad) & (y <= mxy + pad))
x, y, z = x[k], y[k], z[k]
print("%d ground points around the window" % len(z))

v = K.fit_variogram(x, y, z, model="spherical", max_lag=10.0)
mdl = K.MODELS[v["model"]]
p = v["params"]
print("variogram: %s, range %.2f m, radius %.1f m\n" % (v["model"], p[2], RADIUS))

tree = cKDTree(np.c_[x, y])
gx, gy = sub.cell_centres()
flat = np.c_[gx.ravel(), gy.ravel()]
dist, idx = tree.query(flat, k=min(MAXP, len(z)), distance_upper_bound=RADIUS)


def run(mode):
    out = np.full(len(flat), -9999.0)
    fb = np.zeros(len(flat), bool)
    for c in range(len(flat)):
        d = dist[c]
        good = np.isfinite(d)
        n = int(good.sum())
        if n < MINP:
            continue
        ii = idx[c][good]
        pz = z[ii]
        px, py = x[ii], y[ii]
        dxy = np.hypot(px[:, None] - px[None, :], py[:, None] - py[None, :])
        A = np.empty((n + 1, n + 1))
        A[:n, :n] = mdl(dxy, *p)
        A[:n, n] = 1.0; A[n, :n] = 1.0; A[n, n] = 0.0
        b = np.empty(n + 1)
        b[:n] = mdl(d[good], *p); b[n] = 1.0
        try:
            w = np.linalg.solve(A, b)
            est = float(w[:n] @ pz)
            if not np.isfinite(est):
                bad = True
            elif mode == "current":
                bad = est < pz.min() or est > pz.max()
            elif mode == "unguarded":
                bad = False
            else:                      # widened
                spread = float(pz.max() - pz.min())
                bad = (est < pz.min() - spread) or (est > pz.max() + spread)
        except np.linalg.LinAlgError:
            bad = True
        if bad:
            ww = 1.0 / np.maximum(d[good], 1e-9) ** 2
            est = float((ww @ pz) / ww.sum())
            fb[c] = True
        out[c] = est
    return out.reshape(NH, NW), fb.reshape(NH, NW)


from osgeo import gdal
gdal.UseExceptions()
rds = gdal.Open(REF)
rb = rds.GetRasterBand(1)
ref = rb.ReadAsArray(X0, Y0, NW, NH).astype("f8")
nod = rb.GetNoDataValue()
rm = np.isfinite(ref) & (ref != 0.0)
if nod is not None:
    rm &= ref != nod
gyy, gxx = np.gradient(np.where(rm, ref, np.nan), full.dy, full.cell)
slope = np.degrees(np.arctan(np.hypot(gxx, gyy)))

BANDS = [(0, 5), (5, 10), (10, 20), (20, 30), (30, 90)]
print("%-11s %8s %9s %9s   %s" % ("arm", "fallback", "rmse", "median",
                                  "median |err| by slope band (deg)"))
print("%-11s %8s %9s %9s   %s" % ("", "", "", "",
      "  ".join("%5.0f-%-3.0f" % b for b in BANDS)))
print("-" * 92)
res = {}
for mode in ("current", "unguarded", "widened"):
    dem, fb = run(mode)
    K.write_geotiff(dem, sub, wkt, os.path.join(OUT, "guard_%s.tif" % mode))
    m = rm & (dem > -9000) & np.isfinite(slope)
    d = np.abs(dem[m] - ref[m])
    cells = []
    for lo, hi in BANDS:
        kk = m & (slope >= lo) & (slope < hi)
        cells.append(np.median(np.abs(dem[kk] - ref[kk])) if kk.sum() > 50
                     else float("nan"))
    res[mode] = {"fallback": float(fb[m].mean()),
                 "rmse": float(np.sqrt((d ** 2).mean())),
                 "median": float(np.median(d)),
                 "by_slope": [float(c) for c in cells],
                 "over_half_m": float((d > 0.5).mean())}
    print("%-11s %7.1f%% %9.3f %9.4f   %s"
          % (mode, 100 * fb[m].mean(), res[mode]["rmse"], res[mode]["median"],
             "  ".join("%9.3f" % c for c in cells)))

# where does the guard actually fire?
dem_c, fb_c = run("current")
m = rm & np.isfinite(slope)
print("\nwhere the current guard diverts to inverse distance:")
for lo, hi in BANDS:
    kk = m & (slope >= lo) & (slope < hi)
    if kk.sum() > 50:
        print("  %5.0f-%-3.0f deg: %6.2f%% of cells diverted (%d cells)"
              % (lo, hi, 100 * fb_c[kk].mean(), kk.sum()))

with open(os.path.join(OUT, "guard_test.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, indent=1)
print("\nrasters and guard_test.json in %s" % OUT)
