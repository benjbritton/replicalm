"""Frame protection must depend on evidence, not on touching the border.

Two cases that look identical in the mask and must behave differently:

  grid sized to the cloud   the border IS the survey edge, strip it
  grid smaller than cloud   the border is a clip, leave it alone
"""
import sys

import numpy as np

sys.path.insert(0, r"C:\Replicalm\src")
from replicalm import finalise
from replicalm.grid import Grid

rng = np.random.default_rng(0)
print("%-34s %12s %12s %8s" % ("case", "columns", "rows", "cells"))

# 1. grid sized to the cloud: data touches every border, all sides are real
g = Grid(origin_x=0.0, origin_y=200.0, cell=1.0, width=200, height=200)
x = rng.uniform(0.0, 200.0, 120000)
y = rng.uniform(0.0, 200.0, 120000)
m, _v, rep = finalise.data_footprint(x, y, g, 5.0, axis="x", strip_m=21.0,
                                     ends_m=21.0)
cc = np.flatnonzero(m.any(axis=0))
rr = np.flatnonzero(m.any(axis=1))
print("%-34s %12s %12s %8d"
      % ("grid sized to cloud", "%d-%d" % (cc[0], cc[-1]),
         "%d-%d" % (rr[0], rr[-1]), int(m.sum())))
print("   past border: %s" % rep["cloud_past_border"])
assert cc[0] >= 20 and rr[0] >= 20, "a real edge must be stripped"

# 2. the same cloud on a grid that clips it: borders are not edges
g2 = Grid(origin_x=40.0, origin_y=160.0, cell=1.0, width=120, height=120)
m2, _v2, rep2 = finalise.data_footprint(x, y, g2, 5.0, axis="x",
                                        strip_m=21.0, ends_m=21.0)
cc2 = np.flatnonzero(m2.any(axis=0))
rr2 = np.flatnonzero(m2.any(axis=1))
print("%-34s %12s %12s %8d"
      % ("grid clips the cloud", "%d-%d" % (cc2[0], cc2[-1]),
         "%d-%d" % (rr2[0], rr2[-1]), int(m2.sum())))
print("   past border: %s" % rep2["cloud_past_border"])
assert cc2[0] == 0 and rr2[0] == 0, "a clipped border must be kept"
assert m2.sum() > 0.95 * m2.size, "nothing should be stripped here"

# 3. clipped on the west only: that side kept, the other three stripped
g3 = Grid(origin_x=40.0, origin_y=260.0, cell=1.0, width=240, height=320)
x3 = rng.uniform(0.0, 200.0, 120000)
y3 = rng.uniform(40.0, 200.0, 120000)
m3, _v3, rep3 = finalise.data_footprint(x3, y3, g3, 5.0, axis="x",
                                        strip_m=21.0, ends_m=21.0)
cc3 = np.flatnonzero(m3.any(axis=0))
rr3 = np.flatnonzero(m3.any(axis=1))
print("%-34s %12s %12s %8d"
      % ("clipped west only", "%d-%d" % (cc3[0], cc3[-1]),
         "%d-%d" % (rr3[0], rr3[-1]), int(m3.sum())))
print("   past border: %s" % rep3["cloud_past_border"])
assert cc3[0] == 0, "the clipped west border must be kept"
assert cc3[-1] < m3.shape[1] - 1, "the east edge is real and must be stripped"
print("")
print("all checks passed")
