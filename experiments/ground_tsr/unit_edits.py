"""Check both edits on small cases before spending an hour of compute."""
import sys

import numpy as np

sys.path.insert(0, r"C:\Replicalm\src")
from replicalm.preprocess import dedupe_xy
from replicalm import finalise
from replicalm.grid import Grid

print("DEDUPE")
x = np.array([0.0, 0.0, 1.0, 1.0, 1.0, 2.0])
y = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
z = np.array([10.0, 10.4, 5.0, 5.2, 5.4, 7.0])
for keep, want in (("lowest", [10.0, 5.0, 7.0]), ("mean", [10.2, 5.2, 7.0])):
    gx, gy, gz, rep = dedupe_xy(x, y, z, keep=keep)
    ok = np.allclose(np.sort(gz), np.sort(want))
    print("  %-7s -> %s  %s  (dropped %d, groups %d)"
          % (keep, np.round(gz, 4).tolist(), "ok" if ok else "WRONG",
             rep["dropped"], rep["collision_groups"]))
    assert ok

# the case that matters: coincident pairs straddling the surface
rng = np.random.default_rng(0)
n = 20000
bx = rng.uniform(0, 100, n)
by = rng.uniform(0, 100, n)
base = 50 + 0.01 * bx
lo = np.column_stack([bx, by, base - 0.016])
hi = np.column_stack([bx, by, base + 0.012])
P = np.vstack([lo, hi])
for keep in ("lowest", "mean"):
    gx, gy, gz, rep = dedupe_xy(P[:, 0], P[:, 1], P[:, 2], keep=keep)
    resid = gz - (50 + 0.01 * gx)
    print("  straddling pairs, keep=%-7s mean residual %+.5f m"
          % (keep, resid.mean()))

print("")
print("FOOTPRINT")
g = Grid(origin_x=0.0, origin_y=100.0, cell=1.0, width=100, height=100)
# a diagonal swath with a thinning margin and one far outlier
rng = np.random.default_rng(1)
pts = []
for r in range(5, 95):
    c0 = 20 + int(0.3 * r)
    w = 30 - int(abs(r - 50) * 0.25)
    cs = rng.choice(np.arange(c0, c0 + max(w, 3)), size=max(w // 2, 2),
                    replace=False)
    for c in cs:
        pts.append((c + 0.5, 100 - (r + 0.5)))
pts.append((95.5, 100 - 50.5))              # a detached outlier
P = np.array(pts)
disc = finalise.support_mask(P[:, 0], P[:, 1], g, 4.0)
span = finalise.return_footprint(P[:, 0], P[:, 1], g, axis="x")


def shape(m, name):
    from scipy import ndimage as ndi
    fill = ndi.binary_fill_holes(m)
    per = int((fill ^ ndi.binary_erosion(fill)).sum())
    area = int(fill.sum())
    rag = per / max(2.0 * np.sqrt(np.pi * area), 1e-9)
    # how many separate pieces?
    _lab, npieces = ndi.label(m)
    # per-row contiguity: a row-coherent mask has exactly one run per row
    runs = 0
    for r in range(m.shape[0]):
        d = np.diff(np.concatenate([[0], m[r].astype(np.int8), [0]]))
        runs += int((d == 1).sum())
    rows = int(m.any(axis=1).sum())
    print("  %-10s area %6d, perimeter %5d, raggedness %5.2f, pieces %3d, "
          "runs per occupied row %.2f"
          % (name, area, per, rag, npieces, runs / max(rows, 1)))


shape(disc, "disc 4 m")
shape(span, "row span")
print("")
print("a row-coherent footprint has exactly 1.00 runs per row by "
      "construction, so islands and scallops cannot be represented")
