"""Does edge_fix_enabled=False actually keep every interpolated cell?"""
import sys

import numpy as np

sys.path.insert(0, r"C:\Replicalm\src")
from replicalm import finalise
from replicalm.config import PRESETS
from replicalm.grid import Grid

for n in ("clear", "baseline"):
    c = PRESETS[n]
    print("%-9s edge_fix_enabled=%-6s interpolator=%s"
          % (n, c.edge_fix_enabled, c.interpolator))

g = Grid(origin_x=0.0, origin_y=200.0, cell=1.0, width=200, height=200)
dem = np.full((200, 200), -9999.0)
dem[50:150, 50:150] = 100.0            # a 40 x 40 patch of "interpolated" cells
wkt = ('PROJCS["t",GEOGCS["g",DATUM["d",SPHEROID["s",6378137,298.257223563]],'
       'PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]],'
       'PROJECTION["Transverse_Mercator"],UNIT["metre",1]]')
import tempfile
import os

print("")
print("%-26s %10s %10s" % ("configuration", "kept", "of input"))
for label, kw in (("edge fixing on, strip 21 m",
                   dict(radius_m=4.0, strip_m=21.0, axis="x")),
                  ("edge fixing on, erosion 4 m",
                   dict(radius_m=4.0, strip_m=None)),
                  ("edge fixing OFF",
                   dict(radius_m=0.0, strip_m=None))):
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "t.tif")
        rep = finalise.finalise(dem.copy(), g, wkt, p, nodata_in=-9999.0,
                                verbose=False, **kw)
        print("%-26s %10d %10d" % (label, rep["cells_after_trim"], 10000))
print("")
print("10000 is the patch the interpolator produced; OFF must keep all of it")
