r"""Does the clutter rule hold on ground it was not fitted to?

The rule -- demote a ground return standing more than 0.20 m above the tenth
percentile of ground returns within 0.75 m -- was chosen on one 500 x 500 window
of South_GLAS_l0s395, against labels taken from the delivered cloud's own
TerraScan classification. Eight candidates on one patch of one tile is enough to
overfit comfortably, so nothing follows from that fit until it is checked
elsewhere.

This applies the same rule, unchanged, to windows it has never seen: another
part of l0s395, and the relief windows of three other tiles from the
recalibration. Scored against each reference surface on the measures that do not
need a G1 render -- roughness on flat ground, RMSE, and error above twenty
degrees, which is where a rule that cleans by eating architecture would show.

WHAT WOULD COUNT AS HOLDING UP
------------------------------
  roughness    moves toward 1.00x the reference on every window, not past it
  RMSE         falls
  >20 deg      falls, or at worst holds. A rise here means the rule is removing
               structure, and no improvement on flat ground would excuse it.

A rule that only works where it was fitted is a description of that window.
"""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import ROOT, SRC, DEM_ROOT
sys.path.insert(0, SRC)
import numpy as np
from osgeo import gdal
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
from replicalm import classify, grid as G, kriging as K
gdal.UseExceptions()

PATCH, HEIGHT, PCT = 0.75, 0.20, 10.0
TESTS = [
    # name, ground LAS, reference raster, window (None = whole clip)
    ("l0s395 mounds", os.path.join(ROOT, r"render_ncalm\work\l0s395_00000_00000_gnd.las"),
     os.path.join(DEM_ROOT, r"Yuc_South\South_Glas\South_Glas\South_GLAS_l0s395_DEM_0p5m_v1.tif"),
     (2100, 0, 500)),
    ("l0s395 relief", os.path.join(ROOT, r"tests\recal\l0s395_smrf_1pass.las"),
     os.path.join(ROOT, r"tests\recal\clips\l0s395_reference.tif"), None),
    ("l0s417 relief", os.path.join(ROOT, r"tests\recal\l0s417_smrf_1pass.las"),
     os.path.join(ROOT, r"tests\recal\clips\l0s417_reference.tif"), None),
    ("l2s444 relief", os.path.join(ROOT, r"tests\recal\l2s444_smrf_1pass.las"),
     os.path.join(ROOT, r"tests\recal\clips\l2s444_reference.tif"), None),
]


def rough(z, m):
    f = np.nan_to_num(np.where(m, z, np.nan), nan=np.nanmean(z[m]))
    return np.abs(f - ndi.uniform_filter(f, 5))


print("rule: demote ground returns more than %.2f m above the %gth percentile "
      "within %.2f m\n" % (HEIGHT, PCT, PATCH))
print("%-15s %9s %10s %9s %9s %9s %9s"
      % ("window", "removed", "roughness", "vs ref", "rmse", ">20deg", "verdict"))
print("-" * 82)
rows = []
for name, las, refp, win in TESTS:
    if not (os.path.exists(las) and os.path.exists(refp)):
        print("%-15s inputs missing" % name); continue
    rds = gdal.Open(refp)
    if win:
        X0, Y0, N = win
        ref = rds.GetRasterBand(1).ReadAsArray(X0, Y0, N, N).astype("f8")
        fullg = G.grid_from_raster(refp)
        g = G.Grid(fullg.origin_x + X0 * fullg.cell,
                   fullg.origin_y - Y0 * fullg.dy, fullg.cell, N, N,
                   cell_y=fullg.cell_y)
    else:
        ref = rds.GetRasterBand(1).ReadAsArray().astype("f8")
        g = G.grid_from_raster(refp)
    gt = rds.GetGeoTransform(); cy, cx = abs(gt[5]), abs(gt[1])
    valid = np.isfinite(ref) & (ref != 0)
    nod = rds.GetRasterBand(1).GetNoDataValue()
    if nod is not None:
        valid &= ref != nod
    gy, gx = np.gradient(np.where(valid, ref, np.nan), cy, cx)
    slope = np.degrees(np.arctan(np.hypot(gx, gy)))
    flat = valid & (slope < 5); steep = valid & (slope >= 20)
    ref_rough = rough(ref, valid)[flat].mean() if flat.sum() > 500 else np.nan

    arr, _ = classify.read_points(las)
    gnd = arr[arr["Classification"] == 2]
    mnx, mny, mxx, mxy = g.bounds
    k = ((gnd["X"] >= mnx - 25) & (gnd["X"] <= mxx + 25) &
         (gnd["Y"] >= mny - 25) & (gnd["Y"] <= mxy + 25))
    gnd = gnd[k]
    x = gnd["X"].astype("f8"); y = gnd["Y"].astype("f8"); z = gnd["Z"].astype("f8")
    if len(z) < 5000:
        print("%-15s only %d ground returns" % (name, len(z))); continue
    tree = cKDTree(np.c_[x, y])
    nb = tree.query_ball_point(np.c_[x, y], PATCH)
    floor = np.array([np.percentile(z[i], PCT) if len(i) >= 3 else z[j]
                      for j, i in enumerate(nb)])
    drop = (z - floor) > HEIGHT

    res = {}
    for tag, keep in (("before", np.ones(len(z), bool)), ("after", ~drop)):
        v = K.fit_variogram(x[keep], y[keep], z[keep], model="spherical",
                            max_lag=10.0)
        dem, _i = K.krige_grid(x[keep], y[keep], z[keep], g, radius=20.0,
                               max_points=16, variogram=v, verbose=False)
        m = valid & (dem > -9000); e = np.abs(dem - ref)
        res[tag] = {
            "roughness": float(rough(dem, m)[flat].mean()) if flat.sum() > 500 else np.nan,
            "rmse": float(np.sqrt((e[m] ** 2).mean())),
            "over20": float((e[steep & m] > 0.5).mean()) if (steep & m).sum() > 200 else None}
    b, a = res["before"], res["after"]
    better = (a["rmse"] < b["rmse"] and
              (a["over20"] is None or b["over20"] is None or a["over20"] <= b["over20"]))
    print("%-15s %8.2f%% %10s %9s %9s %9s %9s"
          % (name, 100 * drop.mean(),
             "%.4f" % a["roughness"] if np.isfinite(a["roughness"]) else "n/a",
             "%.2fx" % (a["roughness"] / ref_rough) if np.isfinite(ref_rough) else "n/a",
             "%.4f" % a["rmse"],
             "n/a" if a["over20"] is None else "%.2f%%" % (100 * a["over20"]),
             "holds" if better else "FAILS"))
    print("%-15s %8s %10s %9s %9s %9s"
          % ("  (before)", "-",
             "%.4f" % b["roughness"] if np.isfinite(b["roughness"]) else "n/a",
             "%.2fx" % (b["roughness"] / ref_rough) if np.isfinite(ref_rough) else "n/a",
             "%.4f" % b["rmse"],
             "n/a" if b["over20"] is None else "%.2f%%" % (100 * b["over20"])))
    rows.append({"window": name, "removed_pct": float(100 * drop.mean()),
                 "ref_roughness": float(ref_rough), "before": b, "after": a,
                 "holds": bool(better)})

with open(os.path.join(ROOT, "render_ncalm", "cleanup", "validation.json"),
          "w", encoding="utf-8") as fh:
    json.dump({"rule": {"patch_m": PATCH, "height_m": HEIGHT,
                        "percentile": PCT}, "windows": rows}, fh, indent=1,
              default=float)
print("\n%d of %d windows improved on both RMSE and steep-ground error"
      % (sum(r["holds"] for r in rows), len(rows)))
