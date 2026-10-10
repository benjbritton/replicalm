"""The N/S strip, on a synthetic block, before another tile run."""
import sys

import numpy as np

sys.path.insert(0, r"C:\Replicalm\src")
from replicalm import finalise
from replicalm.config import PRESETS
from replicalm.grid import Grid

c = PRESETS["clear"]
print("config: strip across-track %.1f m, ends %.1f m, footprint %.1f m"
      % (c.edge_strip_m, c.edge_strip_ends_m,
         c.footprint_radius_m or c.search_radius_m))

g = Grid(origin_x=0.0, origin_y=300.0, cell=1.0, width=200, height=300)
rng = np.random.default_rng(0)
x = rng.uniform(40.0, 160.0, 80000)
y = rng.uniform(40.0, 260.0, 80000)

print("")
print("returns occupy x 40-160, y 40-260 (a 120 by 220 m block)")
print("%-28s %12s %12s %8s" % ("configuration", "columns", "rows", "cells"))
for nm, strip, ends in (("no strip", None, None),
                        ("across-track only", 21.0, None),
                        ("across-track and ends", 21.0, 21.0)):
    m, v, rep = finalise.data_footprint(x, y, g, 5.0, axis="x",
                                        strip_m=strip, ends_m=ends)
    rr = np.flatnonzero(m.any(axis=1))
    cc = np.flatnonzero(m.any(axis=0))
    print("%-28s %12s %12s %8d"
          % (nm, "%d-%d" % (cc[0], cc[-1]), "%d-%d" % (rr[0], rr[-1]),
             int(m.sum())))

m, v, rep = finalise.data_footprint(x, y, g, 5.0, axis="x", strip_m=21.0,
                                    ends_m=21.0)
rr = np.flatnonzero(m.any(axis=1))
cc = np.flatnonzero(m.any(axis=0))
print("")
print("expected columns ~56-143 and rows ~56-243:")
print("  5 m footprint pad outward, then 21 m stripped back on all four sides")
ok = (54 <= cc[0] <= 58 and 141 <= cc[-1] <= 145
      and 54 <= rr[0] <= 58 and 241 <= rr[-1] <= 245)
print("  got columns %d-%d, rows %d-%d : %s"
      % (cc[0], cc[-1], rr[0], rr[-1], "ok" if ok else "WRONG"))
print("  report: %s" % {k: rep[k] for k in ("strip_m", "ends_m", "areas",
                                            "cells", "void_cells")})
assert ok
print("")
print("frame safety: a block touching the lattice edge keeps that side")
g2 = Grid(origin_x=0.0, origin_y=300.0, cell=1.0, width=200, height=300)
x2 = rng.uniform(0.0, 160.0, 80000)
y2 = rng.uniform(40.0, 260.0, 80000)
m2, _v, _r = finalise.data_footprint(x2, y2, g2, 5.0, axis="x",
                                     strip_m=21.0, ends_m=21.0)
cc2 = np.flatnonzero(m2.any(axis=0))
print("  returns reach x=0; mask starts at column %d (0 = frame kept)"
      % int(cc2[0]))
assert cc2[0] == 0
print("")
print("all checks passed")
