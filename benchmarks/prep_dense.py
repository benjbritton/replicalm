r"""Clip the two hyper-dense tiles that exercise the ceiling rule.

WHY THESE TWO, AND WHY NOT AT RANDOM
------------------------------------
The five pilot tiles were drawn at random and all five came back well
conditioned. That tests the native-first default, which is the branch that runs
almost always, but it leaves the ceiling rule unexercised: a remedy that never
fired has not been shown to work.

Density here is all returns over the header bounding box, which is the only
density a header can give. It is not the ground density that produced l0s444's
89% fallback -- that tile reads 6.42 on this metric and 8.54 in ground points
per square metre after classification, because the bounding box includes the
empty margin either side of the flight strip and because only a third of its
returns reached the floor. The metric is therefore useful for ordering tiles and
useless as an absolute threshold, so the screen is anchored on l0s444's own
value rather than on a round number.

  l0s435   11.17, the densest in the pool, 1.74x l0s444, and steep at 6.6 deg
  l1s400    9.24, 1.44x, a different region and the least bounding-box padding
            in the top band

If the ground density of these lands near or above l0s444's 8.54 and the
fallback spikes, the ceiling rule gets its test. If it lands there and the
fallback does not spike, the rule's trigger is wrong and that is worth more than
a confirmation.
"""
import json, os, sys, time
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import ROOT, SRC, DEM_ROOT, LAS_ROOT, TESTS
sys.path.insert(0, SRC)
from replicalm import clip

WANT = ["l0s435", "l1s400"]
OUT = TESTS + r"\\dense_clips"

with open(TESTS + r"\\pool_traits.json", encoding="utf-8") as fh:
    pool = json.load(fh)
by = {r["tile"]: r for r in pool}

sel = []
for t in WANT:
    r = by[t]
    sel.append({"region": r["region"], "tile": t, "las": r["las"],
                "reference_dem": r["reference_dem"],
                "chosen": "hyper-dense screen: %.2f pts/m2 all returns over the "
                          "header bounding box, against l0s444's 6.42"
                          % r["density"]})
    print("%-8s %-14s %6.2f pts/m2 header  slope %5.2f deg  relief %5.1f m"
          % (t, r["region"], r["density"], r["mean_slope_deg"], r["relief_m"]))
print()

recs = []
for t in sel:
    t0 = time.time()
    try:
        rec = clip.prepare(t, OUT, win_m=400.0)
        rec["seconds"] = round(time.time() - t0, 1)
        recs.append(rec)
        print("  -> %d points, %.0f s\n" % (rec["points"], time.time() - t0))
    except Exception as e:
        print("  FAILED %s: %s\n" % (t["tile"], e))

with open(os.path.join(OUT, "index.json"), "w", encoding="utf-8") as fh:
    json.dump(recs, fh, indent=1)
for r in recs:
    area = r["window"]["window_m"] ** 2
    print("%-8s %d points in the clip, %.1f pts/m2 all returns"
          % (r["tile"], r["points"], r["points"] / area))
print("\nprepared %d of %d" % (len(recs), len(sel)))
