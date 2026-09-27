r"""The three-strip cell inventory as a GIS layer, so sites can be picked by eye.

Choosing test surfaces by reading coordinates out of a terminal table and typing
them into ArcGIS is slow and gets them wrong: three sites nominated that way
turned out to be closed canopy, and two of them had only two strips. The
inventory already holds everything needed to see that at a glance -- strip
coverage, slope, canopy height, understory fraction -- so this writes it out to
be loaded and looked at.

Four formats, because the obvious two are the least useful ones. A CSV carries
no coordinate system at all -- XY Table To Point has to be told -- and ArcGIS
Pro's Catalog pane does not list a `.geojson` extension, only `.json`. So the
ones to open are the shapefile and the GeoPackage, which carry EPSG:32616
themselves and need nothing said about them. The CSV and JSON are written too,
for anything that prefers them.

`suggests` is a hint, not a classification. It flags which rung a cell could
serve, on the same thresholds used to pick the existing ones, and a cell can
suit none of them.
"""
import csv
import gzip
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

RESULTS = os.environ.get(
    "REPLICALM_NR_RESULTS",
    os.path.join(os.path.dirname(os.path.abspath(__file__)),
                 "results", "nr_block"))
CRS = "EPSG:32616"

FIELDS = ["x", "y", "strips_n", "strips", "slope_pct", "band_pct",
          "canopy_p90_m", "over2m_pct", "ground_per_m2", "resid_sd_m",
          "tile", "suggests"]


def suggests(c):
    """Which rung this cell could serve, on the thresholds already in use."""
    out = []
    band = 100 * c["frac_band_015_035"]
    over2 = 100 * c["frac_over_2m"]
    if c["canopy_p90_m"] < 0.5 and over2 < 2 and band < 2:
        out.append("bare")
    if over2 < 5 and c["canopy_p90_m"] < 3 and band >= 15:
        out.append("low_understory")
    if over2 < 5 and c["canopy_p90_m"] < 3 and 5 <= band < 15:
        out.append("sparse_understory")
    if c["slope_pct"] > 12:
        out.append("riser")
    if over2 > 35 and c["canopy_p90_m"] > 8:
        out.append("forest")
    return "|".join(out)


# Shapefile field names are capped at ten characters and silently truncated,
# which turns canopy_p90_m and ground_per_m2 into collisions waiting to happen.
# Mapped explicitly instead.
SHP_NAMES = {"strips_n": "strips_n", "strips": "strips",
             "slope_pct": "slope_pct", "band_pct": "band_pct",
             "canopy_p90_m": "canopy_m", "over2m_pct": "over2m_pct",
             "ground_per_m2": "gnd_m2", "resid_sd_m": "resid_sd",
             "tile": "tile", "suggests": "suggests"}


def _write_vector(rows, shp_path, gpkg_path):
    """A shapefile and a GeoPackage, both carrying the CRS themselves."""
    try:
        from osgeo import ogr, osr
    except ImportError:
        print("  (no GDAL here, so no shapefile or GeoPackage)")
        return []
    ogr.UseExceptions()
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(32616)

    out = []
    for path, driver, layer_name in ((shp_path, "ESRI Shapefile", "nr_cells"),
                                     (gpkg_path, "GPKG", "nr_cells_3strip")):
        drv = ogr.GetDriverByName(driver)
        if drv is None:
            print("  (no %s driver)" % driver)
            continue
        if os.path.exists(path):
            drv.DeleteDataSource(path)
        ds = drv.CreateDataSource(path)
        lyr = ds.CreateLayer(layer_name, srs, ogr.wkbPoint)
        is_shp = driver == "ESRI Shapefile"
        for f in FIELDS:
            if f in ("x", "y"):
                continue
            name = SHP_NAMES[f] if is_shp else f
            if f in ("strips_n",):
                fd = ogr.FieldDefn(name, ogr.OFTInteger)
            elif f in ("strips", "tile", "suggests"):
                fd = ogr.FieldDefn(name, ogr.OFTString)
                fd.SetWidth(60)
            else:
                fd = ogr.FieldDefn(name, ogr.OFTReal)
            lyr.CreateField(fd)
        defn = lyr.GetLayerDefn()
        for r in rows:
            feat = ogr.Feature(defn)
            for f in FIELDS:
                if f in ("x", "y"):
                    continue
                v = r[f]
                if v is None:
                    continue
                feat.SetField(SHP_NAMES[f] if is_shp else f, v)
            pt = ogr.Geometry(ogr.wkbPoint)
            pt.AddPoint(float(r["x"]), float(r["y"]))
            feat.SetGeometry(pt)
            lyr.CreateFeature(feat)
            feat = None
        ds = None
        out.append(path)
    return out


def main(argv=None):
    src = os.path.join(RESULTS, "cells_3strip.json.gz")
    with gzip.open(src, "rt", encoding="utf-8") as fh:
        cells = json.load(fh)
    print("%d three-strip cells from %s" % (len(cells), src))

    rows = []
    for c in cells:
        rows.append({
            "x": round(c["x"], 2), "y": round(c["y"], 2),
            "strips_n": len(c["strips"]),
            "strips": "-".join(str(s) for s in c["strips"]),
            "slope_pct": round(c["slope_pct"], 2),
            "band_pct": round(100 * c["frac_band_015_035"], 2),
            "canopy_p90_m": round(c["canopy_p90_m"], 2),
            "over2m_pct": round(100 * c["frac_over_2m"], 1),
            "ground_per_m2": c["ground_per_m2"],
            "resid_sd_m": c.get("resid_sd_m"),
            "tile": c["tile"], "suggests": suggests(c)})

    csv_path = os.path.join(RESULTS, "nr_cells_3strip.csv")
    with open(csv_path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)

    gj_path = os.path.join(RESULTS, "nr_cells_3strip.geojson")
    with open(gj_path, "w", encoding="utf-8") as fh:
        json.dump({
            "type": "FeatureCollection",
            "crs": {"type": "name", "properties": {"name": CRS}},
            "features": [{"type": "Feature",
                          "geometry": {"type": "Point",
                                       "coordinates": [r["x"], r["y"]]},
                          "properties": {k: r[k] for k in FIELDS
                                         if k not in ("x", "y")}}
                         for r in rows]}, fh)

    shp_path = os.path.join(RESULTS, "nr_cells_3strip.shp")
    gpkg_path = os.path.join(RESULTS, "nr_cells_3strip.gpkg")
    written = _write_vector(rows, shp_path, gpkg_path)

    counts = {}
    for r in rows:
        for s in (r["suggests"].split("|") if r["suggests"] else ["none"]):
            counts[s] = counts.get(s, 0) + 1
    print("")
    for k in sorted(counts, key=lambda k: -counts[k]):
        print("  %-20s %6d" % (k, counts[k]))
    print("")
    for w in written:
        print("wrote %s" % w)
    print("wrote %s" % csv_path)
    print("wrote %s   (rename to .json if Catalog will not list it)" % gj_path)
    print("CRS is %s; cells are 20 m and the coordinates are their centres."
          % CRS)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
