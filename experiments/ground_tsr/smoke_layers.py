"""Do all four G1 layers now share the DEM's lattice?

Also checks the composite is unchanged: it already centred each layer before
blending, so the fix must move the sidecars into line without moving it.
"""
import os
import shutil
import sys
import tempfile

import numpy as np
import rasterio

sys.path.insert(0, r"C:\Replicalm\src")
from replicalm import visualize as VIS

SRC = (r"C:\Replicalm\test100526\centro454_e2e"
       r"\AMIGACarb_Yuc_Centro_GLAS_Apr2013_l0s454_DEM.tif")
OLD_G1 = (r"C:\Replicalm\test100526\centro454_e2e\G1"
          r"\AMIGACarb_Yuc_Centro_GLAS_Apr2013_l0s454_G1.tif")

with rasterio.open(SRC) as r:
    z = r.read(1)
    prof = r.profile
# a window is enough and keeps this to seconds
win = (slice(3000, 3800), slice(400, 1200))
sub = z[win]
p2 = dict(prof)
p2.update(height=sub.shape[0], width=sub.shape[1],
          transform=rasterio.windows.transform(
              rasterio.windows.Window(400, 3000, 800, 800), prof["transform"]))
d = tempfile.mkdtemp()
dem_p = os.path.join(d, "t_DEM.tif")
with rasterio.open(dem_p, "w", **p2) as dst:
    dst.write(sub, 1)

written = VIS.g1(dem_p, d, keep_layers=True)
with rasterio.open(dem_p) as r:
    base = (r.height, r.width, round(r.transform.c, 4), round(r.transform.f, 4))
print("DEM lattice: %d x %d at origin %.4f, %.4f" % base)
print("")
print("%-10s %12s %16s %16s %8s" % ("layer", "size", "origin x", "origin y",
                                    "matches"))
allok = True
for nm, path in sorted(written["layers"].items()):
    with rasterio.open(path) as r:
        got = (r.height, r.width, round(r.transform.c, 4),
               round(r.transform.f, 4))
    ok = got == base
    allok &= ok
    print("%-10s %12s %16.4f %16.4f %8s"
          % (nm, "%d x %d" % (got[0], got[1]), got[2], got[3],
             "yes" if ok else "NO"))
with rasterio.open(written["g1"]) as r:
    print("")
    print("composite: %d x %d at origin %.4f  (%s)"
          % (r.height, r.width, r.transform.c,
             "matches" if (r.height, r.width) == base[:2] else "DIFFERS"))
shutil.rmtree(d, ignore_errors=True)
print("")
print("all four layers on the DEM lattice: %s" % allok)
assert allok
