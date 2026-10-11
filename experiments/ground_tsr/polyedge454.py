"""Polygon-based edge treatment, tested on the untrimmed centro454 DEM.

The procedure, as specified:

  1  interpolate to the full extent of the points, no trim   (already done)
  2  close each active-data area into a polygon with no interior rings
  3  use those polygons as the active / no-data mask
  4  no-data becomes 0 in the DEM, which makes it 0,0,0 in the G1
  5  any cell inside a polygon with no value becomes 1, so 1,1,1 in the G1
     marks "inside the survey, no measurement" and 0,0,0 means "outside"
  6  erode the active area by 21 m

Run against the existing no-trim DEM rather than by reprocessing, so the edge
treatment is tested in isolation. Nothing in the pipeline is changed yet.

Step 5 has a consequence worth watching: a cell set to 1 sits ~130 m below the
terrain around it, so an interior hole would appear in the G1 as a pit with a
cliff wall rather than as a flat marker. The count of such cells is reported;
if it is zero the question is moot on this tile, and if it is not, that is
worth seeing before the rule is adopted.
"""
import json
import os
import sys

import numpy as np
import rasterio
from rasterio import features
from scipy import ndimage

sys.path.insert(0, r"C:\Replicalm\src")

SRC = (r"C:\Replicalm\test100526\centro454_nofix50"
       r"\AMIGACarb_Yuc_Centro_GLAS_Apr2013_l0s454_DEM.tif")
ARCH = (r"D:\_Archive_EdgeFixed\Yuc_Centro\Centro_Glas\Centro_Glas"
        r"\Centro_GLAS_l0s454_DEM_0p5m_v1.tif")
OUT = r"C:\Replicalm\test100526\centro454_polyedge"
ERODE_M = 21.0
HOLE_VALUE = 1.0
os.makedirs(OUT, exist_ok=True)

with rasterio.open(SRC) as r:
    Z = r.read(1).astype("f8")
    prof, tr, cell, crs = r.profile, r.transform, abs(r.transform.a), r.crs
    nd = r.nodata
H, W = Z.shape
print("source %d x %d at %.6f m, nodata %s" % (H, W, cell, nd))

# step 1 input: what the interpolator actually produced
valid = np.isfinite(Z)
if nd is not None:
    valid &= Z != nd
valid &= Z != 0.0
print("active cells from the interpolator: %d (%.2f%% of the lattice)"
      % (int(valid.sum()), 100.0 * valid.mean()))

# ---- step 2: one closed polygon per area, no interior rings
lab, n_lab = ndimage.label(valid)
print("")
print("CONNECTED ACTIVE AREAS: %d" % n_lab)
closed = np.zeros_like(valid)
areas = []
for k in range(1, n_lab + 1):
    comp = lab == k
    filled = ndimage.binary_fill_holes(comp)
    closed |= filled
    areas.append({"id": int(k), "cells": int(comp.sum()),
                  "cells_closed": int(filled.sum()),
                  "holes_closed": int(filled.sum() - comp.sum()),
                  "hectares": float(filled.sum() * cell * cell / 10000.0)})
areas.sort(key=lambda a: -a["cells_closed"])
for a in areas[:8]:
    print("  area %-3d %10d cells  %9.2f ha  interior holes closed %7d"
          % (a["id"], a["cells_closed"], a["hectares"], a["holes_closed"]))
if len(areas) > 8:
    print("  ... and %d smaller" % (len(areas) - 8))
total_holes = sum(a["holes_closed"] for a in areas)
print("  interior cells enclosed by the fill, total: %d" % total_holes)

# ---- step 6 first, so both masks can be polygonised and compared
er_cells = int(round(ERODE_M / cell))
st = ndimage.generate_binary_structure(2, 1)
eroded = ndimage.binary_erosion(closed, structure=st, iterations=er_cells,
                                border_value=0)
print("")
print("EROSION %.1f m (%d cells): %d -> %d cells, %d removed (%.2f%%)"
      % (ERODE_M, er_cells, int(closed.sum()), int(eroded.sum()),
         int(closed.sum() - eroded.sum()),
         100.0 * (closed.sum() - eroded.sum()) / max(closed.sum(), 1)))
lab2, n2 = ndimage.label(eroded)
print("  areas surviving erosion: %d" % n2)


# ---- vectors
def write_vectors(mask, stem):
    shapes = [(geom, int(v)) for geom, v in
              features.shapes(mask.astype("uint8"), mask=mask,
                              transform=tr, connectivity=4) if v == 1]
    feats = []
    for i, (geom, _v) in enumerate(shapes, 1):
        rings = geom["coordinates"]
        # exterior ring only: the specification allows no interior rings
        geom2 = {"type": "Polygon", "coordinates": [rings[0]]}
        area_m2 = 0.0
        xs = [p[0] for p in rings[0]]
        ys = [p[1] for p in rings[0]]
        for j in range(len(xs) - 1):
            area_m2 += xs[j] * ys[j + 1] - xs[j + 1] * ys[j]
        area_m2 = abs(area_m2) / 2.0
        feats.append({"type": "Feature",
                      "properties": {"poly_id": i,
                                     "area_m2": round(area_m2, 2),
                                     "area_ha": round(area_m2 / 10000.0, 3),
                                     "rings_dropped": len(rings) - 1},
                      "geometry": geom2})
    fc = {"type": "FeatureCollection",
          "crs": {"type": "name",
                  "properties": {"name": crs.to_string() if crs else None}},
          "features": feats}
    gj = os.path.join(OUT, stem + ".geojson")
    with open(gj, "w") as fh:
        json.dump(fc, fh)
    written = [gj]
    try:
        import fiona
        from fiona.crs import from_wkt
        schema = {"geometry": "Polygon",
                  "properties": {"poly_id": "int", "area_m2": "float",
                                 "area_ha": "float", "rings_dropped": "int"}}
        gpkg = os.path.join(OUT, stem + ".gpkg")
        with fiona.open(gpkg, "w", driver="GPKG", schema=schema,
                        crs=from_wkt(crs.to_wkt()) if crs else None,
                        layer=stem) as dst:
            dst.writerecords(feats)
        written.append(gpkg)
    except Exception as exc:
        print("  (GeoPackage not written: %s)" % str(exc)[:60])
    print("  %s: %d polygon(s), largest %.2f ha"
          % (stem, len(feats),
             max([f["properties"]["area_ha"] for f in feats], default=0.0)))
    for p in written:
        print("    %s" % p)
    return feats


print("")
print("VECTORS")
f_closed = write_vectors(closed, "centro454_active_closed")
f_eroded = write_vectors(eroded, "centro454_active_eroded21m")

# ---- steps 4 and 5: write the DEM under the convention
out = np.where(eroded, Z, 0.0)
inside_novalue = eroded & ~valid
n_marked = int(inside_novalue.sum())
out[inside_novalue] = HOLE_VALUE
zero_valued = eroded & valid & (Z == 0.0)
out[zero_valued] = HOLE_VALUE
print("")
print("DEM UNDER THE CONVENTION")
print("  cells set to 0 (outside the polygons):      %d"
      % int((~eroded).sum()))
print("  cells set to %g (inside, no measurement):    %d"
      % (HOLE_VALUE, n_marked))
print("  active cells whose value was exactly 0:     %d"
      % int(zero_valued.sum()))
print("  active cells carrying a measurement:        %d"
      % int((eroded & valid & (Z != 0.0)).sum()))
p = dict(prof)
p.update(dtype="float32", nodata=0.0, compress="lzw", tiled=True,
         blockxsize=512, blockysize=512, BIGTIFF="IF_SAFER")
dem_out = os.path.join(OUT, "dem_polyedge_21m.tif")
with rasterio.open(dem_out, "w", **p) as dst:
    dst.write(out.astype("float32"), 1)
print("  wrote %s" % dem_out)

# ---- measurement against the archive
with rasterio.open(ARCH) as r:
    A = r.read(1).astype("f8")
    andv = r.nodata
if andv is not None:
    A[A == andv] = np.nan
A[A > 1e30] = np.nan
A[A == 0] = np.nan
Am = np.isfinite(A)


def shape(m):
    fill = ndimage.binary_fill_holes(m)
    per = int((fill ^ ndimage.binary_erosion(fill)).sum())
    ar = int(fill.sum())
    _l, k = ndimage.label(m)
    return per, ar, per / max(2.0 * np.sqrt(np.pi * ar), 1e-9), int(k)


print("")
print("%-22s %10s %10s %11s %8s" % ("mask", "area", "perimeter",
                                    "raggedness", "pieces"))
for nm, m in (("archive", Am), ("untrimmed", valid), ("closed", closed),
              ("closed + 21 m erode", eroded)):
    per, ar, rg, k = shape(m)
    print("%-22s %10d %10d %11.2f %8d" % (nm, ar, per, rg, k))

final = eroded & valid & (Z != 0.0)
common = Am & final
dist = ndimage.distance_transform_edt(common) * cell
print("")
print("ELEVATION AGAINST THE ARCHIVE on the eroded active area")
print("%-12s %10s %9s %9s %9s" % ("setback", "cells", "RMSE", "MAE", "bias"))
for sb in (0.0, 10.0, 25.0, 50.0):
    m = common & (dist > sb)
    if m.sum() < 1000:
        continue
    e = (Z - A)[m]
    print("%-12s %10d %9.4f %9.4f %+9.4f"
          % ("%.0f m" % sb, int(m.sum()), float(np.sqrt(np.mean(e ** 2))),
             float(np.mean(np.abs(e))), float(np.mean(e))))
print("")
print("coverage: archive %d, ours after erosion %d (%+.1f%%)"
      % (int(Am.sum()), int(final.sum()),
         100.0 * (final.sum() / max(Am.sum(), 1) - 1)))

with open(os.path.join(OUT, "polyedge_report.json"), "w") as fh:
    json.dump({"erode_m": ERODE_M, "hole_value": HOLE_VALUE,
               "areas": areas, "areas_after_erosion": int(n2),
               "polygons_closed": len(f_closed),
               "polygons_eroded": len(f_eroded),
               "interior_cells_closed": int(total_holes),
               "cells_marked_no_measurement": n_marked,
               "dem": dem_out}, fh, indent=2)
print("wrote %s" % os.path.join(OUT, "polyedge_report.json"))
