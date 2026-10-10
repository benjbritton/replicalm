"""The polygon edge procedure, run on a mask derived from the returns.

The previous attempt traced the solver's 150 m block union: 65.8 percent of
boundary cells fell within three cells of a block line and 28.4 percent of the
"active" area had no ground return within 20 m. The archive has none beyond
19.73 m, so the reference's footprint is a distance-to-data footprint and that
is what the mask has to be built from.

The procedure itself is unchanged -- close each area into a polygon with no
interior rings, mask, 0 outside, 1 inside without a measurement, erode -- it
is only being given the right mask to work on.

Thresholds are swept because the right one is a property of the survey's point
spacing, not a constant. The archive's own distance distribution is the target
to match.
"""
import json
import os
import sys

import numpy as np
import rasterio
from rasterio import features
from scipy import ndimage

sys.path.insert(0, r"C:\Replicalm\src")
from replicalm import classify
from replicalm.preprocess import dedupe_xy

DEM = (r"C:\Replicalm\test100526\centro454_nofix50"
       r"\AMIGACarb_Yuc_Centro_GLAS_Apr2013_l0s454_DEM.tif")
GND = (r"C:\Replicalm\test100526\centro454_nofix50"
       r"\AMIGACarb_Yuc_Centro_GLAS_Apr2013_l0s454_ground_cleaned.laz")
ARCH = (r"D:\_Archive_EdgeFixed\Yuc_Centro\Centro_Glas\Centro_Glas"
        r"\Centro_GLAS_l0s454_DEM_0p5m_v1.tif")
OUT = r"C:\Replicalm\test100526\centro454_polyedge_data"
THRESHOLDS = (3.0, 5.0, 10.0, 20.0)
ERODE_M = 21.0
HOLE_VALUE = 1.0
os.makedirs(OUT, exist_ok=True)

with rasterio.open(DEM) as r:
    Z = r.read(1).astype("f8")
    prof, at, cell, crs, nd = (r.profile, r.transform, abs(r.transform.a),
                               r.crs, r.nodata)
H, W = Z.shape
has_value = np.isfinite(Z) & (Z != 0.0)
if nd is not None:
    has_value &= Z != nd

with rasterio.open(ARCH) as r:
    A = r.read(1).astype("f8")
    andv = r.nodata
if andv is not None:
    A[A == andv] = np.nan
A[A > 1e30] = np.nan
A[A == 0] = np.nan
Am = np.isfinite(A)

p = classify._pdal().Pipeline(json.dumps({"pipeline": [
    {"type": "readers.las", "filename": GND, "nosrs": True}]}))
p.execute()
a = p.arrays[0]
gx, gy, gz, _ = dedupe_xy(a["X"].astype("f8"), a["Y"].astype("f8"),
                          a["Z"].astype("f8"), keep="mean")
del a, p
col = np.floor((gx - at.c) / cell).astype(np.int64)
row = np.floor((at.f - gy) / cell).astype(np.int64)
ok = (col >= 0) & (col < W) & (row >= 0) & (row < H)
occ = np.zeros((H, W), bool)
occ[row[ok], col[ok]] = True
dist = ndimage.distance_transform_edt(~occ, sampling=cell)
print("%d ground returns in %d occupied cells; archive covers %d cells"
      % (len(gz), int(occ.sum()), int(Am.sum())))

er = int(round(ERODE_M / cell))
st = ndimage.generate_binary_structure(2, 1)


def shape(m):
    f = ndimage.binary_fill_holes(m)
    per = int((f ^ ndimage.binary_erosion(f)).sum())
    ar = int(f.sum())
    _l, k = ndimage.label(m)
    e = m ^ ndimage.binary_erosion(m)
    ys, xs = np.nonzero(e)
    pitch = int(round(150.0 / cell))
    blk = float(np.mean(((xs % pitch) < 3) | ((xs % pitch) > pitch - 4))) \
        if len(xs) else 0.0
    return {"area": ar, "perimeter": per,
            "raggedness": per / max(2.0 * np.sqrt(np.pi * ar), 1e-9),
            "pieces": int(k), "on_block_lines": blk}


ref = shape(Am)
print("")
print("%-16s %10s %10s %11s %7s %12s %10s"
      % ("mask", "area", "perimeter", "raggedness", "pieces", "on block line",
         "vs archive"))
print("%-16s %10d %10d %11.2f %7d %11.1f%% %10s"
      % ("archive", ref["area"], ref["perimeter"], ref["raggedness"],
         ref["pieces"], 100 * ref["on_block_lines"], "-"))
rows = []
for t in THRESHOLDS:
    near = dist <= t
    closed = np.zeros_like(near)
    lab, n = ndimage.label(near)
    holes = 0
    for k in range(1, n + 1):
        comp = lab == k
        f = ndimage.binary_fill_holes(comp)
        holes += int(f.sum() - comp.sum())
        closed |= f
    eroded = ndimage.binary_erosion(closed, structure=st, iterations=er,
                                    border_value=0)
    s = shape(eroded)
    _l, n_poly = ndimage.label(eroded)
    rows.append({"threshold_m": t, "near": int(near.sum()),
                 "closed": int(closed.sum()), "holes_closed": holes,
                 "eroded": int(eroded.sum()), "polygons": int(n_poly),
                 **{k: v for k, v in s.items()}})
    print("%-16s %10d %10d %11.2f %7d %11.1f%% %+9.1f%%"
          % ("data <= %.0f m" % t, s["area"], s["perimeter"],
             s["raggedness"], s["pieces"], 100 * s["on_block_lines"],
             100.0 * (s["area"] / ref["area"] - 1)))

best = min(rows, key=lambda r: abs(r["area"] / ref["area"] - 1))
T = best["threshold_m"]
print("")
print("closest to the archive's footprint: %.0f m (%+.1f%%), %d polygon(s), "
      "%d interior cells closed"
      % (T, 100.0 * (best["area"] / ref["area"] - 1), best["polygons"],
         best["holes_closed"]))

near = dist <= T
closed = np.zeros_like(near)
lab, n = ndimage.label(near)
for k in range(1, n + 1):
    closed |= ndimage.binary_fill_holes(lab == k)
final = ndimage.binary_erosion(closed, structure=st, iterations=er,
                               border_value=0)

out = np.where(final, Z, 0.0)
no_meas = final & ~has_value
out[no_meas] = HOLE_VALUE
out[final & has_value & (Z == 0.0)] = HOLE_VALUE
print("")
print("DEM: %d cells outside (0), %d inside without a measurement (%g), "
      "%d with a measurement"
      % (int((~final).sum()), int(no_meas.sum()), HOLE_VALUE,
         int((final & has_value & (Z != 0.0)).sum())))
pr = dict(prof)
pr.update(dtype="float32", nodata=0.0, compress="lzw", tiled=True,
          blockxsize=512, blockysize=512, BIGTIFF="IF_SAFER")
dem_out = os.path.join(OUT, "dem_polyedge_data%.0fm_erode%.0fm.tif"
                       % (T, ERODE_M))
with rasterio.open(dem_out, "w", **pr) as dst:
    dst.write(out.astype("float32"), 1)
print("  wrote %s" % dem_out)


def write_vectors(mask, stem):
    feats = []
    for i, (geom, v) in enumerate(
            features.shapes(mask.astype("uint8"), mask=mask, transform=at,
                            connectivity=4), 1):
        if v != 1:
            continue
        rings = geom["coordinates"]
        xs = [q[0] for q in rings[0]]
        ys = [q[1] for q in rings[0]]
        ar = abs(sum(xs[j] * ys[j + 1] - xs[j + 1] * ys[j]
                     for j in range(len(xs) - 1))) / 2.0
        feats.append({"type": "Feature",
                      "properties": {"poly_id": i, "area_m2": round(ar, 2),
                                     "area_ha": round(ar / 10000.0, 3),
                                     "rings_dropped": len(rings) - 1},
                      "geometry": {"type": "Polygon",
                                   "coordinates": [rings[0]]}})
    gj = os.path.join(OUT, stem + ".geojson")
    with open(gj, "w") as fh:
        json.dump({"type": "FeatureCollection",
                   "crs": {"type": "name",
                           "properties": {"name": crs.to_string()}},
                   "features": feats}, fh)
    feats.sort(key=lambda f: -f["properties"]["area_ha"])
    print("  %s: %d polygon(s); largest %.2f ha, smallest %.4f ha"
          % (stem, len(feats), feats[0]["properties"]["area_ha"],
             feats[-1]["properties"]["area_ha"]))
    print("    %s" % gj)
    return feats


print("")
print("VECTORS")
write_vectors(closed, "centro454_data%.0fm_closed" % T)
write_vectors(final, "centro454_data%.0fm_eroded%.0fm" % (T, ERODE_M))

scored = final & has_value & Am & (Z != 0.0)
d2 = ndimage.distance_transform_edt(scored) * cell
print("")
print("ELEVATION AGAINST THE ARCHIVE")
print("%-10s %10s %9s %9s %9s" % ("setback", "cells", "RMSE", "MAE", "bias"))
for sb in (0.0, 10.0, 25.0, 50.0):
    m = scored & (d2 > sb)
    if m.sum() < 1000:
        continue
    e = (Z - A)[m]
    print("%-10s %10d %9.4f %9.4f %+9.4f"
          % ("%.0f m" % sb, int(m.sum()), float(np.sqrt(np.mean(e ** 2))),
             float(np.mean(np.abs(e))), float(np.mean(e))))
dd = dist[final]
print("")
print("distance to the nearest return inside the final mask: p50 %.2f, "
       "p99 %.2f, max %.2f m" % (float(np.percentile(dd, 50)),
                                 float(np.percentile(dd, 99)), float(dd.max())))
with open(os.path.join(OUT, "polyedge_data_report.json"), "w") as fh:
    json.dump({"thresholds": rows, "chosen_m": T, "erode_m": ERODE_M,
               "dem": dem_out}, fh, indent=2)
print("wrote %s" % os.path.join(OUT, "polyedge_data_report.json"))
