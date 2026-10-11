"""Automatic support radius and void handling, tested on centro454.

Two things under test, neither of which has a user-visible parameter:

  the radius   chosen at the knee of the coverage-versus-radius curve, so it
               comes from this survey's own gap structure rather than from
               the 20 m that happened to fit the calibration tile
  the voids    filled in the DEM so the surface stays continuous, marked
               1,1,1 in the G1 afterwards, and recorded as polygons -- rather
               than written into the DEM as an elevation of 1, which would
               put a hundred-metre cliff around every hole

The G1 is computed on the closed footprint before the 21 m strip, so RVT's
horizon scans always have real data within reach of any pixel that survives,
and the per-layer stretch is not set by terrain that is about to be discarded.
"""
import json
import os
import sys

import numpy as np
import rasterio
from rasterio import features
from scipy import ndimage

sys.path.insert(0, r"C:\Replicalm\src")
from replicalm import classify, finalise
from replicalm import visualize as VIS
from replicalm.grid import Grid
from replicalm.preprocess import dedupe_xy

DEM = (r"C:\Replicalm\test100526\centro454_nofix50"
       r"\AMIGACarb_Yuc_Centro_GLAS_Apr2013_l0s454_DEM.tif")
GND = (r"C:\Replicalm\test100526\centro454_nofix50"
       r"\AMIGACarb_Yuc_Centro_GLAS_Apr2013_l0s454_ground_cleaned.laz")
ARCH = (r"D:\_Archive_EdgeFixed\Yuc_Centro\Centro_Glas\Centro_Glas"
        r"\Centro_GLAS_l0s454_DEM_0p5m_v1.tif")
OUT = r"C:\Replicalm\test100526\centro454_autoedge"
STRIP_M = 21.0
os.makedirs(OUT, exist_ok=True)

with rasterio.open(DEM) as r:
    Z = r.read(1).astype("f8")
    prof, at, cell, crs, nd = (r.profile, r.transform, abs(r.transform.a),
                               r.crs, r.nodata)
H, W = Z.shape
have = np.isfinite(Z) & (Z != 0.0)
if nd is not None:
    have &= Z != nd
g = Grid(origin_x=at.c, origin_y=at.f, cell=cell, width=W, height=H)

p = classify._pdal().Pipeline(json.dumps({"pipeline": [
    {"type": "readers.las", "filename": GND, "nosrs": True}]}))
p.execute()
a = p.arrays[0]
gx, gy, gz, _ = dedupe_xy(a["X"].astype("f8"), a["Y"].astype("f8"),
                          a["Z"].astype("f8"), keep="mean")
del a, p
print("%d ground returns" % len(gz))

# the footprint is the interpolator's reach, not a fitted threshold
from replicalm.config import PRESETS
R = float(PRESETS["clear"].search_radius_m)
print("  support radius %.1f m, the configured search radius" % R)
closed, voids, rep = finalise.data_footprint(gx, gy, g, R, axis="x")
final, voids_f, rep_f = finalise.data_footprint(gx, gy, g, R, axis="x",
                                                strip_m=STRIP_M)
print("  closed footprint %d cells in %d area(s), %d void cells"
      % (rep["cells"], rep["areas"], rep["void_cells"]))
print("  after a %.0f m strip: %d cells in %d area(s), %d void cells"
      % (STRIP_M, rep_f["cells"], rep_f["areas"], rep_f["void_cells"]))

with rasterio.open(ARCH) as r:
    A = r.read(1).astype("f8")
    andv = r.nodata
if andv is not None:
    A[A == andv] = np.nan
A[A > 1e30] = np.nan
A[A == 0] = np.nan
Am = np.isfinite(A)
print("  vs the archive footprint: %+.1f%%"
      % (100.0 * (rep_f["cells"] / Am.sum() - 1)))

# DEM continuous over the closed footprint, so the G1 has no cliffs
Zc, filled = finalise.fill_voids(Z, have, closed)
print("  voids filled in the DEM: %d cells" % int(filled.sum()))
dem_closed = np.where(closed, Zc, 0.0)
tmp = os.path.join(OUT, "_dem_closed_for_g1.tif")
pr = dict(prof)
pr.update(dtype="float32", nodata=0.0, compress="lzw", tiled=True,
          blockxsize=512, blockysize=512, BIGTIFF="IF_SAFER")
with rasterio.open(tmp, "w", **pr) as dst:
    dst.write(dem_closed.astype("float32"), 1)

print("")
print("G1 on the closed footprint, before the strip")
written = VIS.g1(tmp, OUT, keep_layers=True, stem="centro454_autoedge",
                 progress=lambda s, m, f=None: print("  [%s] %s" % (s, m),
                                                     flush=True))

# now apply the strip to the DEM and to every raster the G1 step produced
dem_final = np.where(final, Zc, 0.0)
dem_out = os.path.join(OUT, "dem_autoedge_%.0fm_strip%.0fm.tif"
                       % (R, STRIP_M))
with rasterio.open(dem_out, "w", **pr) as dst:
    dst.write(dem_final.astype("float32"), 1)
print("")
print("wrote %s" % dem_out)

stamped = 0
with rasterio.open(written["g1"]) as r:
    rgb = r.read()
    gprof = r.profile
out_rgb = rgb.copy()
for b in range(out_rgb.shape[0]):
    band = out_rgb[b]
    band[~final] = 0
    out_rgb[b] = band
# inside the survey but carrying no measurement: 1,1,1, never 0,0,0
mark = final & filled
if mark.any():
    for b in range(out_rgb.shape[0]):
        out_rgb[b][mark] = 1
    stamped = int(mark.sum())
# and no genuine active pixel may read as pure black
blk = final & (out_rgb.max(axis=0) == 0)
if blk.any():
    for b in range(out_rgb.shape[0]):
        out_rgb[b][blk] = 1
g1_out = os.path.join(OUT, "G1_autoedge_%.0fm_strip%.0fm.tif" % (R, STRIP_M))
with rasterio.open(g1_out, "w", **gprof) as dst:
    dst.write(out_rgb)
print("wrote %s" % g1_out)
print("  void pixels stamped 1,1,1: %d; active pixels lifted off black: %d"
      % (stamped, int(blk.sum())))

for nm, path in written.get("layers", {}).items():
    with rasterio.open(path) as r:
        lay = r.read(1)
        lp = r.profile
    lay = VIS._centre_into(lay, final.shape)
    lay = np.where(final, lay, lp.get("nodata") or 0)
    lp = dict(lp); lp.update(height=final.shape[0],
                             width=final.shape[1])
    with rasterio.open(path, "w", **lp) as dst:
        dst.write(lay, 1)
print("  %d RVT layers trimmed to the same mask" % len(written.get("layers", {})))


def vectors(mask, stem, extra=None):
    feats = []
    for i, (geom, v) in enumerate(features.shapes(
            mask.astype("uint8"), mask=mask, transform=at, connectivity=4), 1):
        if v != 1:
            continue
        ring = geom["coordinates"][0]
        xs = [q[0] for q in ring]
        ys = [q[1] for q in ring]
        ar = abs(sum(xs[j] * ys[j + 1] - xs[j + 1] * ys[j]
                     for j in range(len(xs) - 1))) / 2.0
        pr_ = {"poly_id": i, "area_m2": round(ar, 2),
               "area_ha": round(ar / 10000.0, 4)}
        if extra:
            pr_.update(extra)
        feats.append({"type": "Feature", "properties": pr_,
                      "geometry": {"type": "Polygon", "coordinates": [ring]}})
    fp = os.path.join(OUT, stem + ".geojson")
    with open(fp, "w") as fh:
        json.dump({"type": "FeatureCollection",
                   "crs": {"type": "name",
                           "properties": {"name": crs.to_string()}},
                   "features": feats}, fh)
    print("  %-26s %3d polygon(s)  %s" % (stem, len(feats), fp))
    return feats


print("")
print("VECTORS")
vectors(closed, "centro454_active_closed")
vectors(final, "centro454_active_strip%.0fm" % STRIP_M)
if filled.any():
    vectors(final & filled, "centro454_voids", {"kind": "no_measurement"})
else:
    print("  centro454_voids             none")

scored = final & have & Am
d2 = ndimage.distance_transform_edt(scored) * cell
print("")
print("ELEVATION AGAINST THE ARCHIVE")
print("%-10s %10s %9s %9s %9s" % ("setback", "cells", "RMSE", "MAE", "bias"))
for sb in (0.0, 25.0, 50.0):
    m = scored & (d2 > sb)
    e = (Z - A)[m]
    print("%-10s %10d %9.4f %9.4f %+9.4f"
          % ("%.0f m" % sb, int(m.sum()), float(np.sqrt(np.mean(e ** 2))),
             float(np.mean(np.abs(e))), float(np.mean(e))))


def shp(m):
    f = ndimage.binary_fill_holes(m)
    per = int((f ^ ndimage.binary_erosion(f)).sum())
    ar = int(f.sum())
    _l, k = ndimage.label(m)
    return ar, per, per / max(2.0 * np.sqrt(np.pi * ar), 1e-9), int(k)


print("")
print("%-22s %10s %10s %11s %7s" % ("mask", "area", "perimeter",
                                    "raggedness", "pieces"))
for nm, m in (("archive", Am), ("closed", closed),
              ("closed + %.0f m strip" % STRIP_M, final)):
    ar, per, rg, k = shp(m)
    print("%-22s %10d %10d %11.2f %7d" % (nm, ar, per, rg, k))

os.remove(tmp)
with open(os.path.join(OUT, "autoedge_report.json"), "w") as fh:
    json.dump({"support_radius_m": R, "strip_m": STRIP_M,
               "closed": rep, "final": rep_f,
               "voids_filled": int(filled.sum()),
               "void_pixels_stamped": stamped,
               "dem": dem_out, "g1": g1_out}, fh, indent=2)
print("")
print("wrote %s" % os.path.join(OUT, "autoedge_report.json"))
