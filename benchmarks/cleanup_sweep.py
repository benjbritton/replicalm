r"""Does demoting high ground returns clear the specks without costing terrain?

The nugget sweep established that a smoother suppresses the specks by blurring
them, at a strength (4% of sill) far beyond what measurement noise justifies
(about 0.6%). This tests the alternative: remove the offending returns instead.

Each arm demotes ground returns standing more than `tolerance` above the tenth
percentile of ground returns within 0.6 m, then kriges and renders through the
same recipe. The control is the untouched classification.

WHAT WOULD MEAN IT WORKS
------------------------
  specks fall        dark cells against the archive drop toward the archive's
                     own level
  roughness falls    toward 1.00x the archive rather than past it -- a filter
                     that drives roughness below 1.0 is removing terrain, not
                     clutter
  structures hold    RMSE and error above 20 degrees do not rise. If a
                     tolerance strips platform edges the steep-band figure
                     moves first.
  little is removed  a rule that demotes several per cent of ground returns is
                     not finding clutter, it is trimming the surface

The archive is rendered at the same window extent, because G1 normalizes each
raster by its own extremes and a full-tile stretch is not comparable with a
window one.
"""
import json, os, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import ROOT, SRC, DEM_ROOT
sys.path.insert(0, SRC)
import numpy as np
from osgeo import gdal
from scipy import ndimage as ndi
from replicalm import classify, cleanup, grid as G, kriging as K
gdal.UseExceptions()

X0, Y0, N = 1500, 5250, 500
RADIUS = 20.0
TOLERANCES = [None, 0.30, 0.20, 0.15, 0.10, 0.05]
REF_DEM = os.path.join(DEM_ROOT, r"Yuc_South\South_Glas\South_Glas"
                                 r"\South_GLAS_l0s395_DEM_0p5m_v1.tif")
GND = os.path.join(ROOT, r"render_ncalm\work\l0s395_00000_04000_gnd.las")
OUT = os.path.join(ROOT, "render_ncalm", "cleanup")
ARCH_G1 = os.path.join(ROOT, r"render_ncalm\nugget\_archive_window\rvt\tile"
                             r"\l0s395_arch_G1_0p5m_v1.tif")
G1ENV = os.path.dirname(os.environ.get(
    "REPLICALM_G1_PYTHON", r"C:\Users\benja\anaconda3\envs\g1\python.exe"))
G1SCRIPT = os.environ.get("REPLICALM_RVT_SCRIPT", r"C:\g1\tools\GLiHT_rvt.py")
os.makedirs(OUT, exist_ok=True)

_env = {k: v for k, v in os.environ.items() if k.upper() not in
        ("GDAL_DATA", "PROJ_LIB", "GDAL_DRIVER_PATH", "PYTHONPATH",
         "PYTHONHOME", "CONDA_PREFIX")}
_env["PATH"] = os.pathsep.join([
    G1ENV, os.path.join(G1ENV, "Library", "bin"),
    os.path.join(G1ENV, "Library", "usr", "bin"),
    os.path.join(G1ENV, "Scripts"), _env.get("PATH", "")])
_env["GDAL_DATA"] = os.path.join(G1ENV, "Library", "share", "gdal")
_env["PROJ_LIB"] = os.path.join(G1ENV, "Library", "share", "proj")

full = G.grid_from_raster(REF_DEM)
sub = G.Grid(full.origin_x + X0 * full.cell, full.origin_y - Y0 * full.dy,
             full.cell, N, N, cell_y=full.cell_y)
rds = gdal.Open(REF_DEM)
ref = rds.GetRasterBand(1).ReadAsArray(X0, Y0, N, N).astype("f8")
gt = rds.GetGeoTransform(); cy, cx = abs(gt[5]), abs(gt[1])
valid = np.isfinite(ref) & (ref != 0)
gy, gx = np.gradient(np.where(valid, ref, np.nan), cy, cx)
slope = np.degrees(np.arctan(np.hypot(gx, gy)))
flat = valid & (slope < 5)
steep = valid & (slope >= 20)
wkt = rds.GetProjection()


def lum(p):
    d = gdal.Open(p)
    return np.dstack([d.GetRasterBand(i + 1).ReadAsArray(0, 0, N, N)
                      for i in range(3)]).mean(axis=2).astype("f8")


def rough(z, m):
    f = np.nan_to_num(np.where(m, z, np.nan), nan=np.nanmean(z[m]))
    return np.abs(f - ndi.uniform_filter(f, 5))


arch = lum(ARCH_G1)
arch_med = np.median(arch[flat])
ref_rough = rough(ref, valid)[flat].mean()
print("archive at this window: roughness %.4f, flat median luminance %.0f\n"
      % (ref_rough, arch_med))

arr0, _ = classify.read_points(GND)
mnx, mny, mxx, mxy = sub.bounds
keep = ((arr0["X"] >= mnx - 25) & (arr0["X"] <= mxx + 25) &
        (arr0["Y"] >= mny - 25) & (arr0["Y"] <= mxy + 25))
arr0 = arr0[keep]
print("%d returns around the window, %d classified ground\n"
      % (len(arr0), int((arr0["Classification"] == 2).sum())))

print("%-9s %8s %9s %10s %9s %9s %9s"
      % ("tolerance", "demoted", "roughness", "vs arch", "dark", "rmse", ">20deg"))
print("-" * 70)
rows = []
for tol in TOLERANCES:
    tag = "raw" if tol is None else "tol%03d" % int(round(tol * 100))
    if tol is None:
        arr, rep = arr0, {"demoted": 0, "fraction": 0.0}
    else:
        arr, rep = cleanup.demote(arr0, radius=0.6, tolerance=tol,
                                  percentile=10.0, verbose=False)
    g = arr[arr["Classification"] == 2]
    x, y, z = g["X"].astype("f8"), g["Y"].astype("f8"), g["Z"].astype("f8")
    v = K.fit_variogram(x, y, z, model="spherical", max_lag=10.0)
    dem, _info = K.krige_grid(x, y, z, sub, radius=RADIUS, max_points=16,
                              variogram=v, verbose=False)
    stage = os.path.join(OUT, tag, "tile")
    os.makedirs(stage, exist_ok=True)
    tif = os.path.join(stage, "l0s395_%s_DEM_0p5m_v1.tif" % tag)
    K.write_geotiff(dem, sub, wkt, tif)
    m = valid & (dem > -9000)
    e = np.abs(dem - ref)
    r = rough(dem, m)[flat].mean()

    rvt = os.path.join(OUT, tag, "rvt")
    g1 = os.path.join(rvt, "tile", "l0s395_%s_G1_0p5m_v1.tif" % tag)
    if not os.path.exists(g1):
        subprocess.run([os.path.join(G1ENV, "python.exe"), G1SCRIPT,
                        "--res", "0p5m", "--no-metadata", "--workers", "1",
                        "--dem-root", os.path.join(OUT, tag),
                        "--rvt-root", rvt],
                       capture_output=True, text=True, env=_env)
    dark = None
    if os.path.exists(g1):
        ln = lum(g1)
        off = np.median(ln[flat]) - arch_med
        dark = int((flat & (((ln - off) - arch) < -35)).sum())
    rows.append({"tolerance": tol, "tag": tag, "demoted": rep["demoted"],
                 "demoted_fraction": rep["fraction"], "roughness": float(r),
                 "roughness_vs_archive": float(r / ref_rough),
                 "dark_cells": dark, "rmse": float(np.sqrt((e[m] ** 2).mean())),
                 "over20": float((e[steep & m] > 0.5).mean())
                 if (steep & m).sum() > 200 else None})
    print("%-9s %8d %9.4f %9.2fx %9s %9.4f %8s"
          % ("none" if tol is None else "%.2f m" % tol, rep["demoted"], r,
             r / ref_rough, "n/a" if dark is None else dark,
             rows[-1]["rmse"],
             "n/a" if rows[-1]["over20"] is None else "%.2f%%" % (100 * rows[-1]["over20"])))

with open(os.path.join(OUT, "cleanup.json"), "w", encoding="utf-8") as fh:
    json.dump({"archive_roughness": float(ref_rough), "arms": rows}, fh,
              indent=1, default=float)
print("\narchive roughness %.4f; 1.00x matches the original workflow, below it "
      "is smoother." % ref_rough)
