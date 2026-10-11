"""The five fixes, on small cases, before another 23-minute tile run."""
import json
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, r"C:\Replicalm\src")
from replicalm import finalise
from replicalm.grid import Grid

print("1. FRAME-CLIPPED RUNS ARE NOT STRIPPED")
m = np.zeros((40, 100), bool)
m[5:35, 0:60] = True                 # touches column 0, ends at 59 inside
out = finalise.transect_edge_mask(m, 1.0, 10.0, "x")
row = out[20]
c = np.flatnonzero(row)
print("   row spans columns %d..%d before, %d..%d after"
      % (0, 59, int(c[0]), int(c[-1])))
print("   west edge kept at the frame: %s ; east edge stripped 10 cells: %s"
      % (c[0] == 0, c[-1] == 49))
assert c[0] == 0 and c[-1] == 49

m2 = np.zeros((40, 100), bool)
m2[5:35, 10:60] = True               # both edges interior
c2 = np.flatnonzero(finalise.transect_edge_mask(m2, 1.0, 10.0, "x")[20])
print("   interior run 10..59 becomes %d..%d (both stripped): %s"
      % (int(c2[0]), int(c2[-1]), c2[0] == 20 and c2[-1] == 49))
assert c2[0] == 20 and c2[-1] == 49

print("")
print("2. FOOTPRINT VECTORS")
g = Grid(origin_x=500000.0, origin_y=2000000.0, cell=1.0, width=100,
         height=80)
mask = np.zeros((80, 100), bool)
mask[10:40, 10:40] = True
mask[50:70, 60:90] = True            # a second, detached area
voids = np.zeros_like(mask)
voids[20:24, 20:24] = True
with tempfile.TemporaryDirectory() as d:
    p = os.path.join(d, "fp.geojson")
    rep = finalise.write_footprint_vectors(mask, g, None, p, voids=voids)
    print("   %d features: %d active, %d void"
          % (rep["features"], rep["active"], rep["voids"]))
    with open(p) as fh:
        fc = json.load(fh)
    for f in fc["features"]:
        pr = f["properties"]
        n_rings = len(f["geometry"]["coordinates"])
        print("   %-16s %8.3f ha, %d ring(s)"
              % (pr["kind"], pr["area_ha"], n_rings))
        assert n_rings == 1, "interior rings must be dropped"
    xs = [q[0] for q in fc["features"][0]["geometry"]["coordinates"][0]]
    print("   georeferenced: x from %.1f (origin %.1f)"
          % (min(xs), g.origin_x))
    assert min(xs) >= g.origin_x - 1
assert rep["active"] == 2 and rep["voids"] == 1

print("")
print("3. VOID FILL KEEPS THE SURFACE CONTINUOUS")
dem = np.full((50, 50), 130.0)
dem[20:25, 20:25] = -9999.0
have = dem != -9999.0
want = np.ones_like(have)
out, filled = finalise.fill_voids(dem, have, want)
print("   %d cells filled; range inside the void %.3f to %.3f"
      % (int(filled.sum()), float(out[20:25, 20:25].min()),
         float(out[20:25, 20:25].max())))
print("   no cliff (all within 1 mm of the surround): %s"
      % bool(np.abs(out[20:25, 20:25] - 130.0).max() < 0.001))
assert np.abs(out[20:25, 20:25] - 130.0).max() < 0.001

print("")
print("4. CONFIG RESOLVES THE RADIUS FROM THE SEARCH RADIUS")
from replicalm.config import PRESETS
for n in ("clear", "baseline"):
    c = PRESETS[n]
    r = (getattr(c, "footprint_radius_m", None) or c.search_radius_m)
    print("   %-9s edge_fix_enabled=%-5s footprint %.1f m, strip %.1f m"
          % (n, c.edge_fix_enabled, r, c.edge_strip_m))
print("")
print("all checks passed")
