r"""Fit a clutter rule against NCALM's own labels, and score it on the surface.

THE LABELS
----------
The delivered G-LiHT clouds already carry TerraScan's classification. Matching
our SMRF ground returns to them by position shows that 3.88% of what we call
ground, NCALM calls unclassified -- and at the cells that render as specks the
figure is 12.2% against 0.8% on control ground.

Rebuilding the surface with exactly those points removed is the oracle:

    roughness   0.0367 -> 0.0292, which is 1.00x the archive
    specks       3,901 -> 253
    RMSE        0.0825 -> 0.0140

So this is entirely a classification residue, and perfect is known.

WHAT IS BEING FITTED
--------------------
A rule using only what is available without TerraScan. Separation against the
labels, as d-prime:

    height above the patch floor   1.82
    patch height spread            1.21
    intensity                      0.37
    scan angle, return structure   under 0.2, useless

Two features carry it. The rule is a pair of thresholds rather than a fitted
model, because a threshold pair can be read, argued with and re-tuned on a new
survey, and because the target is not classification accuracy.

WHY ACCURACY IS THE WRONG TARGET
--------------------------------
Precision against the labels is low at every threshold -- the rejected class is
under 4% of the population, so most points above any cut are genuine ground. But
discarding genuine ground is close to free here: at four returns per square
metre, losing 15% still leaves 3.4, and kriging averages what remains. What is
expensive is leaving clutter in.

So each candidate is scored on the surface it produces -- roughness against the
archive, speck count, RMSE, and error above 20 degrees to catch a rule that
cleans flat ground by eating architecture -- not on how well it reproduces the
labels.
"""
import json, os, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import ROOT, SRC, DEM_ROOT
sys.path.insert(0, SRC)
import numpy as np
from osgeo import gdal
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
from replicalm import classify, grid as G, kriging as K
gdal.UseExceptions()

X0, Y0, N = 1500, 5250, 500
PATCH = 0.75
OUT = os.path.join(ROOT, "render_ncalm", "cleanup")
NUG = os.path.join(ROOT, "render_ncalm", "nugget")
REF = os.path.join(DEM_ROOT, r"Yuc_South\South_Glas\South_Glas"
                             r"\South_GLAS_l0s395_DEM_0p5m_v1.tif")
ALL_LAS = os.path.join(OUT, "trench_all.las")
OUR_GND = os.path.join(ROOT, r"render_ncalm\work\l0s395_00000_04000_gnd.las")
G1ENV = os.path.dirname(os.environ.get(
    "REPLICALM_G1_PYTHON", r"C:\Users\benja\anaconda3\envs\g1\python.exe"))
G1SCRIPT = os.environ.get("REPLICALM_RVT_SCRIPT", r"C:\g1\tools\GLiHT_rvt.py")
CANDIDATES = [(None, None), (0.20, 0.00), (0.25, 0.40), (0.30, 0.00),
              (0.30, 0.45), (0.40, 0.45), (0.30, 0.60), (0.40, 0.60)]

_env = {k: v for k, v in os.environ.items() if k.upper() not in
        ("GDAL_DATA", "PROJ_LIB", "GDAL_DRIVER_PATH", "PYTHONPATH",
         "PYTHONHOME", "CONDA_PREFIX")}
_env["PATH"] = os.pathsep.join([
    G1ENV, os.path.join(G1ENV, "Library", "bin"),
    os.path.join(G1ENV, "Library", "usr", "bin"),
    os.path.join(G1ENV, "Scripts"), _env.get("PATH", "")])
_env["GDAL_DATA"] = os.path.join(G1ENV, "Library", "share", "gdal")
_env["PROJ_LIB"] = os.path.join(G1ENV, "Library", "share", "proj")

full = G.grid_from_raster(REF)
sub = G.Grid(full.origin_x + X0 * full.cell, full.origin_y - Y0 * full.dy,
             full.cell, N, N, cell_y=full.cell_y)
rds = gdal.Open(REF)
ref = rds.GetRasterBand(1).ReadAsArray(X0, Y0, N, N).astype("f8")
gt = rds.GetGeoTransform(); cy, cx = abs(gt[5]), abs(gt[1])
ox, oy = gt[0] + X0 * cx, gt[3] - Y0 * cy
valid = np.isfinite(ref) & (ref != 0)
gy, gx = np.gradient(np.where(valid, ref, np.nan), cy, cx)
slope = np.degrees(np.arctan(np.hypot(gx, gy)))
flat = valid & (slope < 5); steep = valid & (slope >= 20)
wkt = rds.GetProjection()


def lum(p):
    d = gdal.Open(p)
    return np.dstack([d.GetRasterBand(i + 1).ReadAsArray(0, 0, N, N)
                      for i in range(3)]).mean(axis=2).astype("f8")


def rough(z, m):
    f = np.nan_to_num(np.where(m, z, np.nan), nan=np.nanmean(z[m]))
    return np.abs(f - ndi.uniform_filter(f, 5))


arch = lum(os.path.join(NUG, "_archive_window", "rvt", "tile",
                        "l0s395_arch_G1_0p5m_v1.tif"))
arch_med = np.median(arch[flat])
ref_rough = rough(ref, valid)[flat].mean()

orig, _ = classify.read_points(ALL_LAS)
ours, _ = classify.read_points(OUR_GND)
og = ours[ours["Classification"] == 2]
k = ((og["X"] >= ox - 25) & (og["X"] <= ox + N * cx + 25) &
     (og["Y"] <= oy + 25) & (og["Y"] >= oy - N * cy - 25))
og = og[k]
t = cKDTree(np.c_[orig["X"], orig["Y"], orig["Z"]])
dist, idx = t.query(np.c_[og["X"], og["Y"], og["Z"]], k=1)
label = (dist < 0.01) & (orig["Classification"][idx] == 1)

x = og["X"].astype("f8"); y = og["Y"].astype("f8"); z = og["Z"].astype("f8")
tree = cKDTree(np.c_[x, y])
nb = tree.query_ball_point(np.c_[x, y], PATCH)
floor = np.array([np.percentile(z[i], 10) if len(i) >= 3 else z[j]
                  for j, i in enumerate(nb)])
spread = np.array([np.ptp(z[i]) if len(i) >= 3 else 0.0 for i in nb])
above = z - floor
print("%d ground returns, %d NCALM-rejected (%.2f%%)\n"
      % (len(z), label.sum(), 100 * label.mean()))

print("%-16s %8s %8s %9s %10s %9s %9s %8s"
      % ("rule", "removed", "recall", "roughness", "vs arch", "specks",
         "rmse", ">20deg"))
print("-" * 86)
rows = []
for h, s in CANDIDATES:
    if h is None:
        drop = np.zeros(len(z), bool); tag = "none"
    else:
        drop = (above > h) & (spread > s) if s else (above > h)
        tag = "h>%.2f" % h + (" s>%.2f" % s if s else "")
    keep = ~drop
    slug = "rule_" + tag.replace(">", "").replace(".", "").replace(" ", "_")
    xx, yy, zz = x[keep], y[keep], z[keep]
    v = K.fit_variogram(xx, yy, zz, model="spherical", max_lag=10.0)
    dem, _i = K.krige_grid(xx, yy, zz, sub, radius=20.0, max_points=16,
                           variogram=v, verbose=False)
    st = os.path.join(OUT, slug, "tile"); os.makedirs(st, exist_ok=True)
    tif = os.path.join(st, "l0s395_%s_DEM_0p5m_v1.tif" % slug)
    K.write_geotiff(dem, sub, wkt, tif)
    rvt = os.path.join(OUT, slug, "rvt")
    g1 = os.path.join(rvt, "tile", "l0s395_%s_G1_0p5m_v1.tif" % slug)
    if not os.path.exists(g1):
        subprocess.run([os.path.join(G1ENV, "python.exe"), G1SCRIPT,
                        "--res", "0p5m", "--no-metadata", "--workers", "1",
                        "--dem-root", os.path.join(OUT, slug),
                        "--rvt-root", rvt],
                       capture_output=True, text=True, env=_env)
    m = valid & (dem > -9000); e = np.abs(dem - ref)
    r = rough(dem, m)[flat].mean()
    specks = None
    if os.path.exists(g1):
        ln = lum(g1); off = np.median(ln[flat]) - arch_med
        specks = int((flat & (((ln - off) - arch) < -35)).sum())
    o20 = (float((e[steep & m] > 0.5).mean()) if (steep & m).sum() > 200
           else None)
    rows.append({"height": h, "spread": s, "removed": int(drop.sum()),
                 "removed_pct": float(100 * drop.mean()),
                 "recall": float(100 * drop[label].mean()),
                 "roughness": float(r), "vs_archive": float(r / ref_rough),
                 "specks": specks, "rmse": float(np.sqrt((e[m] ** 2).mean())),
                 "over20": o20})
    print("%-16s %7.2f%% %7.1f%% %9.4f %9.2fx %9s %9.4f %7s"
          % (tag, 100 * drop.mean(), 100 * drop[label].mean(), r,
             r / ref_rough, "n/a" if specks is None else specks,
             rows[-1]["rmse"],
             "n/a" if o20 is None else "%.2f%%" % (100 * o20)))

with open(os.path.join(OUT, "rule_fit.json"), "w", encoding="utf-8") as fh:
    json.dump({"archive_roughness": float(ref_rough), "arms": rows}, fh,
              indent=1, default=float)
print("\noracle, for reference: roughness 0.0292 (1.00x), specks 253, "
      "rmse 0.0140, removing 3.88%")
