r"""Does the fitted correlation range stabilise, or does it track the window?

A variogram fitted over too short a lag sees only the nugget and reports a range
that is an artefact of the window rather than a property of the terrain. The
test is simple: fit the same points at increasing lags. A real range settles; a
truncated one keeps climbing with the window that produced it.
"""
import json, os, sys
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import ROOT, SRC, DEM_ROOT, LAS_ROOT, TESTS
sys.path.insert(0, SRC)
import numpy as np
from replicalm import classify, kriging as K

OUT = TESTS + r"\\three_tiles"
LAGS = (5.0, 10.0, 20.0, 40.0, 80.0)

print("%-8s %8s" % ("tile", "lag m"), end="")
for L in LAGS:
    print("%10.0f" % L, end="")
print("     verdict")
print("-" * 78)

rows = []
for tile in ("l8s431", "l0s444", "l0s395"):
    p = os.path.join(OUT, "%s_ground.las" % tile)
    if not os.path.exists(p):
        print("%-8s missing %s" % (tile, p)); continue
    arr, _ = classify.read_points(p)
    gnd = arr[arr["Classification"] == 2]
    x, y, z = gnd["X"], gnd["Y"], gnd["Z"]
    ranges, sills = [], []
    for L in LAGS:
        try:
            v = K.fit_variogram(x, y, z, model="spherical", max_lag=L,
                                n_lags=14)
            ranges.append(v["params"][2]); sills.append(v["params"][1])
        except Exception:
            ranges.append(float("nan")); sills.append(float("nan"))
    # a range that settles changes little between the last two windows
    last2 = [r for r in ranges[-2:] if np.isfinite(r)]
    settled = (len(last2) == 2 and abs(last2[1] - last2[0]) / max(last2) < 0.15)
    tracking = all(ranges[i] < ranges[i + 1] for i in range(len(ranges) - 1)
                   if np.isfinite(ranges[i]) and np.isfinite(ranges[i + 1]))
    verdict = ("settles at %.1f m" % last2[-1] if settled
               else "tracks the window" if tracking else "unstable")
    print("%-8s %8s" % (tile, "range"), end="")
    for r in ranges:
        print("%10.2f" % r, end="")
    print("   %s" % verdict)
    print("%-8s %8s" % ("", "sill"), end="")
    for sv in sills:
        print("%10.4f" % sv, end="")
    print()
    rows.append({"tile": tile, "lags": list(LAGS), "ranges": ranges,
                 "sills": sills, "settled": bool(settled)})

with open(os.path.join(OUT, "lag_stability.json"), "w", encoding="utf-8") as fh:
    json.dump(rows, fh, indent=1, default=float)
print("\nif the range climbs with every window, the 5 m fit was truncated and")
print("the radius derived from it was too small for the terrain.")
