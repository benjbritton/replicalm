"""What is actually in Centro_GLAS_l0s441_GND_v1.laz, and what made the DEM?

Reports only what the files state: LAS header fields, the full Classification
histogram, and the LAS 1.4 classification flag bits, which matter here because
point data record format 6 carries Synthetic / KeyPoint / Withheld / Overlap
as separate bits rather than as class numbers. A TerraScan model keypoint can
therefore appear as the KeyPoint flag on a class 2 point rather than as class
8, and the two cases are distinguishable only by reading the bits.

Also reports the per-class planimetric spacing, because a keypoint set is a
thinned subset and shows a far larger nearest-neighbour distance than the full
ground class.
"""
import json
import sys

import numpy as np
import rasterio

sys.path.insert(0, r"C:\Replicalm\src")
from replicalm import classify
from scipy.spatial import cKDTree

GND = (r"D:\_LAZ\Yuc_Centro\Centro_Glas\Centro_Glas"
       r"\Centro_GLAS_l0s441_GND_v1.laz")
RAW = (r"D:\GLiHT_LAS_orig\Yuc_Centro\Centro_Glas\Centro_Glas"
       r"\AMIGACarb_Yuc_Centro_GLAS_Apr2013_l0s441.las")
DEMS = [
    (r"D:\_Archive_EdgeFixed\Yuc_Centro\Centro_Glas\Centro_Glas"
     r"\Centro_GLAS_l0s441_DEM_0p5m_v1.tif"),
]
CLSNAME = {0: "never classified", 1: "unclassified", 2: "ground",
           3: "low vegetation", 4: "medium vegetation", 5: "high vegetation",
           6: "building", 7: "low point (noise)", 8: "reserved / model key",
           9: "water", 10: "rail", 11: "road surface", 12: "overlap (1.1-1.3)",
           13: "wire guard", 18: "high noise"}


def header(path):
    p = classify._pdal().Pipeline(json.dumps({"pipeline": [
        {"type": "readers.las", "filename": path, "nosrs": True,
         "count": 1}]}))
    p.execute()
    md = p.metadata
    if isinstance(md, str):
        md = json.loads(md)
    md = md.get("metadata", md)
    las = md.get("readers.las", md)
    if isinstance(las, list):
        las = las[0]
    keys = ("major_version", "minor_version", "dataformat_id", "count",
            "software_id", "system_id", "creation_year", "creation_doy",
            "global_encoding", "scale_x", "offset_x", "vlr_count",
            "compressed", "project_id")
    print("  " + path.rsplit("\\", 1)[-1])
    for k in keys:
        if k in las:
            print("    %-18s %s" % (k, las[k]))
    extra = [k for k in las if k.startswith("vlr_")]
    for k in sorted(extra):
        v = las[k]
        if isinstance(v, dict):
            print("    %-18s userid=%s recordid=%s desc=%s"
                  % (k, v.get("user_id"), v.get("record_id"),
                     v.get("description")))
    return las


def load(path):
    p = classify._pdal().Pipeline(json.dumps({"pipeline": [
        {"type": "readers.las", "filename": path, "nosrs": True}]}))
    p.execute()
    return p.arrays[0]


print("LAS HEADERS AS WRITTEN")
h_gnd = header(GND)
print("")
h_raw = header(RAW)

print("")
print("CLASSIFICATION HISTOGRAM, ground delivery")
a = load(GND)
cl = a["Classification"].astype(np.int32)
tot = len(cl)
print("  %d points total" % tot)
for v in np.unique(cl):
    n = int((cl == v).sum())
    print("    class %-3d %-24s %10d  %6.2f%%"
          % (v, CLSNAME.get(int(v), "?"), n, 100.0 * n / tot))

print("")
print("LAS 1.4 CLASSIFICATION FLAG BITS, ground delivery")
names = [f for f in a.dtype.names]
for f in ("Synthetic", "KeyPoint", "Withheld", "Overlap", "ClassFlags"):
    if f in names:
        v = a[f].astype(np.int64)
        if f == "ClassFlags":
            for bit, lbl in ((0, "Synthetic"), (1, "KeyPoint"),
                             (2, "Withheld"), (3, "Overlap")):
                n = int(((v >> bit) & 1).sum())
                print("    ClassFlags bit %d %-12s %10d  %6.2f%%"
                      % (bit, lbl, n, 100.0 * n / tot))
        else:
            n = int((v != 0).sum())
            print("    %-12s set on %10d points  %6.2f%%"
                  % (f, n, 100.0 * n / tot))
    else:
        print("    %-12s dimension absent" % f)
print("  dimensions present: %s" % ", ".join(names))

print("")
print("PER-CLASS PLANIMETRIC SPACING, ground delivery")
for v in np.unique(cl):
    m = cl == v
    if m.sum() < 1000:
        continue
    xy = np.column_stack([a["X"][m], a["Y"][m]])
    k = min(200000, len(xy))
    sub = xy[np.random.default_rng(0).choice(len(xy), k, replace=False)]
    d, _ = cKDTree(xy).query(sub, k=2)
    print("    class %-3d n=%-9d median nearest neighbour %.3f m"
          % (v, int(m.sum()), float(np.median(d[:, 1]))))

print("")
print("CLASSIFICATION HISTOGRAM, original flightline LAS")
b = load(RAW)
cr = b["Classification"].astype(np.int32)
print("  %d points total" % len(cr))
for v in np.unique(cr):
    n = int((cr == v).sum())
    print("    class %-3d %-24s %10d  %6.2f%%"
          % (v, CLSNAME.get(int(v), "?"), n, 100.0 * n / len(cr)))

print("")
print("DEM TAGS AS WRITTEN")
for d in DEMS:
    with rasterio.open(d) as r:
        print("  " + d.rsplit("\\", 1)[-1])
        print("    driver %s, %dx%d, cell %.4f, dtype %s, nodata %s"
              % (r.driver, r.height, r.width, abs(r.transform.a),
                 r.dtypes[0], r.nodata))
        tags = r.tags()
        for k in sorted(tags):
            print("    %-22s %s" % (k, str(tags[k])[:110]))
        for ns in r.tag_namespaces():
            if ns:
                t = r.tags(ns=ns)
                if t:
                    print("    [%s] %s" % (ns, json.dumps(t)[:300]))
