"""How long does the G1 stage actually add to a tile?

Timed on the v4 DEM, 7368 by 1699 at 0.5 m -- 12.5 million cells, a full
G-LiHT transect. Each of the four layers is timed separately, then the blend
and the writes, so the figure can be scaled to other tile sizes rather than
quoted as one opaque number.
"""
import os
import shutil
import sys
import tempfile
import time

sys.path.insert(0, r"C:\Replicalm\src")
import numpy as np
import rasterio

from replicalm import visualize as VIS

DEM = (r"C:\Replicalm\test100526\centro454_e2e_v4"
       r"\AMIGACarb_Yuc_Centro_GLAS_Apr2013_l0s454_DEM.tif")

with rasterio.open(DEM) as r:
    H, W = r.height, r.width
    cell = abs(r.transform.a)
print("tile %d x %d at %.3f m = %.1f million cells"
      % (H, W, cell, H * W / 1e6))

dem, c, profile, crs = VIS._read_dem(DEM)
print("")
print("%-10s %9s" % ("layer", "seconds"))
total = 0.0
raws = {}
for name, _t, _m, _b, _c, _g in VIS.G1_LAYERS:
    t0 = time.time()
    raw = VIS._centre_into(VIS.compute_layer(name, dem, c), dem.shape)
    dt = time.time() - t0
    raws[name] = raw
    total += dt
    print("%-10s %9.1f" % (name, dt))
print("%-10s %9.1f" % ("layers", total))

t0 = time.time()
canvas = np.zeros(dem.shape, dtype=np.float32)
for name, transparency, mode, brightness, contrast, gamma in VIS.G1_LAYERS:
    layer = VIS.normalise(raws[name])
    layer = VIS.adjust(layer, brightness, contrast, gamma)
    canvas = VIS.blend(canvas, layer, mode, 1.0 - transparency / 100.0)
blend_t = time.time() - t0
print("%-10s %9.1f" % ("blend", blend_t))

d = tempfile.mkdtemp()
t0 = time.time()
for name in raws:
    VIS._write_float(raws[name], os.path.join(d, "%s.tif" % name), profile,
                     crs)
wl = time.time() - t0
t0 = time.time()
VIS._write_rgb(canvas, os.path.join(d, "g1.tif"), profile, crs)
wc = time.time() - t0
size = sum(os.path.getsize(os.path.join(d, f)) for f in os.listdir(d))
shutil.rmtree(d, ignore_errors=True)
print("%-10s %9.1f" % ("write x4", wl))
print("%-10s %9.1f" % ("write G1", wc))
grand = total + blend_t + wl + wc
print("")
print("TOTAL %.1f s = %.1f min for %.1f million cells"
      % (grand, grand / 60.0, H * W / 1e6))
print("  %.2f s per million cells" % (grand / (H * W / 1e6)))
print("  output added: %.0f MB (composite plus four layers)"
      % (size / 1e6))
print("")
print("note: an installer build was compressing on the same machine during")
print("this measurement, so treat these as an upper bound")
