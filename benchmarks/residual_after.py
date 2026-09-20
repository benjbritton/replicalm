r"""Did the tighter filter dissolve the hybrid state, or is it terrain?

The residual diagnosis was run on DEMs built with threshold_m = 1.0, which the
filter sweep showed to be six to twelve times too loose. This repeats it on the
best filter for each tile, at two lag resolutions: 60 m for comparison with the
earlier numbers, and 15 m to see below four metres, where l0s395's structure was
previously invisible.
"""
import json, os, sys
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import ROOT, SRC, DEM_ROOT, LAS_ROOT, TESTS
sys.path.insert(0, SRC)
from replicalm import residual

with open(TESTS + r"\\clips\index.json", encoding="utf-8") as fh:
    clips = {c["tile"]: c for c in json.load(fh)}

S = TESTS + r"\\class_sweep"
BEST = {"l0s395": "csfrigid2", "l8s431": "smrfs05t05", "l0s444": "smrfs05t025"}
BEFORE = {"l0s395": 0.74, "l8s431": 0.69, "l0s444": 0.92}

print("%-8s %-12s %7s %9s %9s %9s %8s %s"
      % ("tile", "filter", "max lag", "nug/sill", "was", "bin m",
         "range m", "resolved"))
print("-" * 86)
rows = []
for tile, tag in BEST.items():
    p = os.path.join(S, "%s_%s.tif" % (tile, tag))
    if not os.path.exists(p):
        cands = [f for f in os.listdir(S)
                 if f.startswith(tile) and f.endswith(".tif")]
        print("%-8s missing %s; have %s" % (tile, os.path.basename(p),
                                            ", ".join(cands[:4])))
        continue
    ref = clips[tile]["clipped_reference"]
    for lag in (60.0, 15.0):
        try:
            r = residual.compare(p, ref, max_lag_m=lag)
            rows.append({"tile": tile, "filter": tag, "max_lag_m": lag, **r})
            print("%-8s %-12s %7.0f %9.2f %9.2f %9.2f %8s %s"
                  % (tile, tag, lag, r["nugget_to_sill"], BEFORE[tile],
                     r["lag_resolution_m"],
                     "%.1f" % r["range_m"] if r["range_m"] else "none",
                     r["range_resolved"]))
        except Exception as e:
            print("%-8s %-12s %7.0f FAILED %s" % (tile, tag, lag, str(e)[:40]))
    print()

with open(os.path.join(S, "residual_after.json"), "w", encoding="utf-8") as fh:
    json.dump(rows, fh, indent=1, default=float)

print("verdicts at the finer resolution:")
for r in rows:
    if r["max_lag_m"] == 15.0:
        print("  %-8s %s" % (r["tile"], r["verdict"]))
