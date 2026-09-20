r"""Does a nugget remove the freckles without removing the terrain?

THE OBSERVATION
---------------
On flat ground the output carries single-cell dark specks that the reference
does not. Measured on the window at column 1500, row 5250 of South_GLAS_l0s395:
2,127 clusters, median size one cell, our surface higher than the reference at
78.9% of them by a median of 2.6 cm, and local roughness 0.0768 m against the
reference's 0.0358 m.

They are not invented. The same local rise appears in the reference at 83.4% of
them; ours is 1.60 times its amplitude. These are real protuberances -- rocks,
karst -- rendered sharp where the reference renders them smooth.

WHY
---
`fit_variogram` searches three nugget values, 0, 5% and 10% of the sill, and
picks exactly zero for this tile at every fitting window; 34 of 60 fits across
ten tiles came back with a zero nugget. A variogram with no nugget makes
ordinary kriging an exact interpolator: it honors every observation, so
per-return ranging noise passes into the raster unaltered. At four points per
square meter with sixteen neighbors inside about 1.2 m, each cell is decided by
one or two returns.

The original workflow kriged over the documented 20 m neighborhood in Surfer,
pulling in hundreds of points per solve, which smooths by averaging.

WHAT A NUGGET DOES
------------------
It models measurement error. With one, the kriging weights no longer sum to a
solution that reproduces the data exactly, so the surface becomes a smoother
rather than an interpolator. The question this answers is whether that removes
the specks while leaving the structures legible, and at what fraction.

WHAT IS SCORED
--------------
  roughness       mean deviation from a 2.5 m local mean, on flat ground. The
                  reference sets the target; matching it is the goal, beating
                  it means over-smoothing.
  anomaly ratio   our local rise at the freckle cells divided by the
                  reference's. 1.0 means the protuberances are rendered at the
                  same amplitude the original workflow gave them.
  freckles        cells markedly darker in our G1 than in the reference's, on
                  flat ground, once each arm is rendered through the same
                  recipe.
  structures      RMSE and error above 20 degrees, so that a setting which
                  smooths the specks away by smoothing everything away is
                  visible as such.
"""
import json, os, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import ROOT, SRC, DEM_ROOT
sys.path.insert(0, SRC)
import numpy as np
from osgeo import gdal
from scipy import ndimage as ndi
from replicalm import grid as G, kriging as K, classify, interpolate
gdal.UseExceptions()

X0, Y0, N = 1500, 5250, 500
REF_DEM = os.path.join(DEM_ROOT, r"Yuc_South\South_Glas\South_Glas"
                                 r"\South_GLAS_l0s395_DEM_0p5m_v1.tif")
REF_G1 = os.path.join(DEM_ROOT, r"Yuc_South\South_Glas\South_Glas"
                                r"\South_GLAS_l0s395_G1_0p5m_v1.tif")
GND = os.path.join(ROOT, r"render_ncalm\work\l0s395_00000_04000_gnd.las")
LAS = os.path.join(ROOT, r"render_ncalm\work\l0s395_00000_04000_gnd.las")
OUT = os.path.join(ROOT, "render_ncalm", "nugget")
os.makedirs(OUT, exist_ok=True)
FRACTIONS = ([float(a) / 100.0 for a in sys.argv[1:]] if len(sys.argv) > 1
             else [0.0, 0.05, 0.10, 0.20])
RADIUS = 20.0
G1_PYTHON = os.environ.get("REPLICALM_G1_PYTHON",
                           r"C:\Users\benja\anaconda3\envs\g1\python.exe")
G1_SCRIPT = os.environ.get("REPLICALM_RVT_SCRIPT", r"C:\g1\tools\GLiHT_rvt.py")

full = G.grid_from_raster(REF_DEM)
sub = G.Grid(full.origin_x + X0 * full.cell, full.origin_y - Y0 * full.dy,
             full.cell, N, N, cell_y=full.cell_y)
rds = gdal.Open(REF_DEM)
ref = rds.GetRasterBand(1).ReadAsArray(X0, Y0, N, N).astype("f8")
cy, cx = abs(rds.GetGeoTransform()[5]), abs(rds.GetGeoTransform()[1])
valid = np.isfinite(ref) & (ref != 0)
gy, gx = np.gradient(np.where(valid, ref, np.nan), cy, cx)
slope = np.degrees(np.arctan(np.hypot(gx, gy)))
flat = valid & (slope < 5)
steep = valid & (slope >= 20)
print("window %d x %d: %d valid, %d flat (<5 deg), %d steep (>=20 deg)"
      % (N, N, valid.sum(), flat.sum(), steep.sum()))


def rough(z, m):
    f = np.where(m, z, np.nan)
    fill = np.nan_to_num(f, nan=np.nanmean(f))
    return np.abs(fill - ndi.uniform_filter(fill, 5))


def anomaly(z, m):
    f = np.where(m, z, np.nan)
    fill = np.nan_to_num(f, nan=np.nanmean(f))
    return fill - ndi.uniform_filter(fill, 11)


def lum(path, x=0, y=0, w=N, h=N):
    d = gdal.Open(path)
    return np.dstack([d.GetRasterBand(i + 1).ReadAsArray(x, y, w, h)
                      for i in range(3)]).mean(axis=2).astype("f8")


ref_rough = rough(ref, valid)[flat].mean()
ref_anom = anomaly(ref, valid)
ref_lum = lum(REF_G1, X0, Y0)
print("reference roughness on flat ground: %.4f m\n" % ref_rough)

arr, _ = classify.read_points(GND)
g = arr[arr["Classification"] == 2]
mnx, mny, mxx, mxy = sub.bounds
k = ((g["X"] >= mnx - 25) & (g["X"] <= mxx + 25) &
     (g["Y"] >= mny - 25) & (g["Y"] <= mxy + 25))
g = g[k]
x, y, z = g["X"].astype("f8"), g["Y"].astype("f8"), g["Z"].astype("f8")
wkt = rds.GetProjection()
print("%d ground returns around the window" % len(z))

base = K.fit_variogram(x, y, z, model="spherical", max_lag=10.0)
nug0, psill0, rng = base["params"]
total_sill = nug0 + psill0
print("fitted variogram: nugget %.5f, sill %.5f, range %.2f m\n"
      % (nug0, total_sill, rng))

rows = []
for f in FRACTIONS:
    v = dict(base)
    v["params"] = (f * total_sill, (1.0 - f) * total_sill, rng)
    dem, info = K.krige_grid(x, y, z, sub, radius=RADIUS, max_points=16,
                             variogram=v, verbose=False)
    tag = "nug%03d" % int(round(f * 100))
    stage = os.path.join(OUT, tag, "tile")
    os.makedirs(stage, exist_ok=True)
    tif = os.path.join(stage, "l0s395_%s_DEM_0p5m_v1.tif" % tag)
    K.write_geotiff(dem, sub, wkt, tif)
    m = valid & (dem > -9000)
    e = np.abs(dem - ref)
    ao = anomaly(dem, m)
    fr = flat & m
    rows.append({
        "nugget_fraction": f,
        "roughness": float(rough(dem, m)[fr].mean()),
        "roughness_vs_ref": float(rough(dem, m)[fr].mean() / ref_rough),
        "rmse": float(np.sqrt((e[m] ** 2).mean())),
        "over20": float((e[steep & m] > 0.5).mean()) if (steep & m).sum() > 200 else None,
        "dem": tif, "tag": tag, "stage_root": os.path.join(OUT, tag),
        "fallback": float(info["fallback_fraction"]),
    })
    print("nugget %4.0f%%  roughness %.4f (%.2fx ref)  rmse %.4f  over20 %s"
          % (100 * f, rows[-1]["roughness"], rows[-1]["roughness_vs_ref"],
             rows[-1]["rmse"],
             "n/a" if rows[-1]["over20"] is None else "%.2f%%" % (100 * rows[-1]["over20"])))

# The G1 recipe runs under a different interpreter, and it must not inherit
# this one's geospatial environment: ArcGIS Pro's GDAL_DATA and PROJ_LIB point
# osgeo at the wrong library and the import fails with "No module named _gdal".
G1_ENV_ROOT = os.path.dirname(G1_PYTHON)
_env = {k: v for k, v in os.environ.items()
        if k.upper() not in ("GDAL_DATA", "PROJ_LIB", "GDAL_DRIVER_PATH",
                             "PYTHONPATH", "PYTHONHOME", "CONDA_PREFIX")}
_env["PATH"] = os.pathsep.join([
    G1_ENV_ROOT,
    os.path.join(G1_ENV_ROOT, "Library", "bin"),
    os.path.join(G1_ENV_ROOT, "Library", "usr", "bin"),
    os.path.join(G1_ENV_ROOT, "Scripts"), _env.get("PATH", "")])
_env["GDAL_DATA"] = os.path.join(G1_ENV_ROOT, "Library", "share", "gdal")
_env["PROJ_LIB"] = os.path.join(G1_ENV_ROOT, "Library", "share", "proj")

print("\nrendering each arm through the G1 recipe...")
for r in rows:
    rvt_root = os.path.join(r["stage_root"], "rvt")
    cmd = [G1_PYTHON, G1_SCRIPT, "--res", "0p5m", "--no-metadata",
           "--workers", "1", "--dem-root", r["stage_root"],
           "--rvt-root", rvt_root]
    p = subprocess.run(cmd, capture_output=True, text=True, env=_env)
    g1 = os.path.join(rvt_root, "tile",
                      "l0s395_%s_G1_0p5m_v1.tif" % r["tag"])
    if p.returncode != 0 or not os.path.exists(g1):
        r["freckles"] = None
        print("  %s: G1 step failed (%s)" % (r["tag"], (p.stderr or "").strip()[-70:]))
        continue
    ln = lum(g1)
    fr = flat & (ln - ref_lum < -35)
    lab, n = ndi.label(fr)
    r["freckles"] = int(fr.sum())
    r["freckle_clusters"] = int(n)
    # hold the dataset: GDAL frees the band with it, and a temporary
    # gdal.Open(...).GetRasterBand(1) leaves the band pointing at freed memory
    _ds = gdal.Open(r["dem"])
    dem = _ds.GetRasterBand(1).ReadAsArray().astype("f8")
    ao = anomaly(dem, valid & (dem > -9000))
    sel = flat & (ref_lum - ln > 35)
    r["anomaly_ratio"] = float(
        np.median(np.abs(ao[sel])) / max(np.median(np.abs(ref_anom[sel])), 1e-9)
    ) if sel.sum() > 50 else None
    print("  %s: %d freckle cells in %d clusters" % (r["tag"], fr.sum(), n))

with open(os.path.join(OUT, "nugget.json"), "w", encoding="utf-8") as fh:
    json.dump({"reference_roughness": float(ref_rough), "arms": rows}, fh,
              indent=1, default=float)

print("\n%s\nSUMMARY\n%s" % ("=" * 74, "=" * 74))
print("%-9s %11s %11s %10s %10s %10s"
      % ("nugget", "roughness", "vs ref", "freckles", "rmse", "over20"))
print("-" * 68)
for r in rows:
    print("%8.0f%% %11.4f %10.2fx %10s %10.4f %9s"
          % (100 * r["nugget_fraction"], r["roughness"],
             r["roughness_vs_ref"],
             "n/a" if r["freckles"] is None else r["freckles"],
             r["rmse"],
             "n/a" if r["over20"] is None else "%.2f%%" % (100 * r["over20"])))
print("\nreference roughness %.4f m; vs ref of 1.00 matches the original "
      "workflow,\nbelow 1.00 is smoother than it." % ref_rough)
