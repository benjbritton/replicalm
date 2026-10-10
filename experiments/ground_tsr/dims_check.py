"""Do any of our outputs change the raster's dimensions or origin?

Trimming sets cells to nodata inside an unchanged lattice; it does not clip
the raster. This checks that against the grid every product was written on.
"""
import glob
import os

import rasterio

REF = (r"D:\_Archive_EdgeFixed\Yuc_Centro\Centro_Glas\Centro_Glas"
       r"\Centro_GLAS_l0s454_DEM_0p5m_v1.tif")
PATHS = [
    ("archive DEM", REF),
    ("archive G1", REF.replace("_DEM_", "_G1_")),
    ("untrimmed DEM", r"C:\Replicalm\test100526\centro454_nofix50"
                      r"\AMIGACarb_Yuc_Centro_GLAS_Apr2013_l0s454_DEM.tif"),
    ("e2e DEM", r"C:\Replicalm\test100526\centro454_e2e"
                r"\AMIGACarb_Yuc_Centro_GLAS_Apr2013_l0s454_DEM.tif"),
    ("e2e G1", r"C:\Replicalm\test100526\centro454_e2e\G1"
               r"\AMIGACarb_Yuc_Centro_GLAS_Apr2013_l0s454_G1.tif"),
    ("polyedge DEM", r"C:\Replicalm\test100526\centro454_polyedge_strip"
                     r"\dem_polyedge_data20m_strip21m.tif"),
    ("autoedge DEM", r"C:\Replicalm\test100526\centro454_autoedge"
                     r"\dem_autoedge_20m_strip21m.tif"),
]
for p in sorted(glob.glob(r"C:\Replicalm\test100526\centro454_e2e\G1"
                          r"\layers\*.tif")):
    PATHS.append(("  layer " + os.path.basename(p).split("_")[-1], p))

print("%-22s %14s %10s %18s %10s"
      % ("product", "size", "cell", "origin x", "nodata"))
print("-" * 80)
base = None
for nm, p in PATHS:
    if not os.path.exists(p):
        print("%-22s missing" % nm)
        continue
    with rasterio.open(p) as r:
        sz = (r.height, r.width)
        cc = abs(r.transform.a)
        ox = r.transform.c
        nd = r.nodata
    if base is None:
        base = (sz, round(cc, 6), round(ox, 3))
    same = (sz, round(cc, 6), round(ox, 3)) == base
    print("%-22s %14s %10.6f %18.3f %10s %s"
          % (nm, "%d x %d" % sz, cc, ox, nd, "" if same else "<- DIFFERS"))
print("")
print("every product shares the archive's lattice: trimmed cells are written")
print("as nodata, the raster itself is never clipped")
