r"""The ladder rungs as a point layer, so they can be opened rather than typed.

Entering an easting and a northing by hand in ArcGIS Pro is slow and easy to get
wrong, and three sites nominated that way turned out to be closed canopy. The
rungs are few and fixed, so they are better carried as a file: load it, zoom to
layer, and they are on the map with their labels.

The set itself is the record of which surfaces the Optimization work is measured
on, so it lives in the repository rather than in a chat log. Coordinates are
UTM 16N, EPSG:32616.
"""
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

RESULTS = os.environ.get(
    "REPLICALM_NR_RESULTS",
    os.path.join(os.path.dirname(os.path.abspath(__file__)),
                 "results", "nr_block"))
EPSG = 32616

# name, easting, northing, role, note
RUNGS = [
    ("bare", 317124.75, 1970163.18, "control",
     "dirt road; every return ground to 3 m, nothing above 0.15 m, slope 1.0%"),
    ("clutter_dense", 317130.60, 1970149.50, "clutter",
     "21.4% of returns in 0.15-0.35 m, canopy p90 0.32 m, 15 m from the bare rung"),
    ("clutter_sparse", 316698.41, 1970492.90, "clutter",
     "3.7% in 0.20-0.35 m; the sparse case observation 6 actually describes"),
    ("riser_open", 316810.00, 1970329.60, "real relief",
     "16.3% open slope, 11.3 ground returns per m2, no canopy. The clean "
     "false-positive test: nothing is standing here, so anything the threshold "
     "removes is a mistake"),
    ("mound_flank", 316622.34, 1970507.42, "architecture",
     "southeast flank of a mound, identified by Ben. Crest 10 m west at 79.2 m "
     "falling 15.7 m to a flat 30 m east; slope 86% within 3 m, 74% within "
     "10 m. Under 85% canopy at 5.0 ground returns per m2, so removals here "
     "are ambiguous between vegetation and terrain -- it shows behaviour on "
     "real architecture, it does not measure the false-positive rate"),
    ("forest", 316450.80, 1970610.00, "canopy",
     "Ka'Kabish platform complex, 18.7 m canopy, 77% of returns above 2 m"),
]

FIELDS = ["name", "role", "easting", "northing", "note"]


def main(argv=None):
    rows = [{"name": n, "role": r, "easting": x, "northing": y, "note": note}
            for n, x, y, r, note in RUNGS]

    base = os.path.join(RESULTS, "nr_ladder_rungs")
    with open(base + ".csv", "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)
    gj = {"type": "FeatureCollection",
          "crs": {"type": "name",
                  "properties": {"name": "EPSG:%d" % EPSG}},
          "features": [{"type": "Feature",
                        "geometry": {"type": "Point",
                                     "coordinates": [r["easting"], r["northing"]]},
                        "properties": {k: r[k] for k in ("name", "role", "note")}}
                       for r in rows]}
    for ext in (".geojson", ".json"):
        with open(base + ext, "w", encoding="utf-8") as fh:
            json.dump(gj, fh, indent=1)

    written = [base + e for e in (".csv", ".geojson", ".json")]
    try:
        from osgeo import ogr, osr
    except ImportError:
        print("no GDAL here, so no shapefile")
    else:
        ogr.UseExceptions()
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(EPSG)
        for path, driver, layer in ((base + ".shp", "ESRI Shapefile", "rungs"),
                                    (base + ".gpkg", "GPKG", "nr_ladder_rungs")):
            drv = ogr.GetDriverByName(driver)
            if drv is None:
                continue
            if os.path.exists(path):
                drv.DeleteDataSource(path)
            ds = drv.CreateDataSource(path)
            lyr = ds.CreateLayer(layer, srs, ogr.wkbPoint)
            for f in ("name", "role", "note"):
                fd = ogr.FieldDefn(f, ogr.OFTString)
                fd.SetWidth(120 if f == "note" else 24)
                lyr.CreateField(fd)
            defn = lyr.GetLayerDefn()
            for r in rows:
                feat = ogr.Feature(defn)
                for f in ("name", "role", "note"):
                    feat.SetField(f, r[f])
                pt = ogr.Geometry(ogr.wkbPoint)
                pt.AddPoint(r["easting"], r["northing"])
                feat.SetGeometry(pt)
                lyr.CreateFeature(feat)
                feat = None
            ds = None
            written.append(path)

    print("%d rungs, EPSG:%d" % (len(rows), EPSG))
    for r in rows:
        print("  %-15s %11.2f %12.2f  %s" % (r["name"], r["easting"],
                                             r["northing"], r["role"]))
    print("")
    for w in written:
        print("wrote %s" % w)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
