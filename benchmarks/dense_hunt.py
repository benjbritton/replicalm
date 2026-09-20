r"""Find a second tile that actually reaches l0s444's ground density.

WHY THE FIRST SCREEN MISSED
---------------------------
The header metric ordered tiles by all returns over the header bounding box, and
neither of the two it picked was dense where it counted:

  l0s435   header 11.17, clip 5.0 all returns, 2.96 ground
  l1s400   header  9.24, clip 11.8 all returns, 2.60 ground

Two separate things broke the prediction. The bounding box is not the flight
strip, so header density and clip density disagree in both directions -- l0s435
fell by half, l1s400 rose. And the share of returns that reach the floor varies
enormously: 88% on l2s443, 32% on l0s444, 22% on l1s400. All returns per square
metre therefore says very little about ground points per square metre, which is
the quantity that degenerates the kriging system.

WHAT THIS DOES INSTEAD
----------------------
Measures the two steps that the header cannot. Clips a 400 m window from each of
the top candidates and reports its real all-return density, then classifies the
best few and reports real ground density. Clipping is seconds per tile and
classification under a minute, so measuring is cheaper than modelling the
relationship.

The target is l0s444's 8.54 ground points per square metre, the density at which
89% of cells fell back to inverse distance. Anything near it exercises the
ceiling rule. If nothing in the pool reaches it, that is the finding: l0s444 is
an outlier and the rule guards a case that occurs once in 453 tiles.
"""
import json, os, sys, time
sys.path.insert(0, r"C:\Replicalm\src")
from dataclasses import replace
from replicalm import clip, classify
from replicalm.config import PRESETS, GroundPass

N_CLIP = 18          # candidates to clip and measure
N_CLASSIFY = 4       # of those, how many to classify
TARGET = 8.54        # l0s444's ground density, the condition being sought
ALREADY = {"l8s431", "l0s444", "l0s395", "l4s478", "l0s417", "l2s444",
           "l0s419", "l2s443", "l0s435", "l1s400"}
OUT = r"C:\Replicalm\tests\dense_hunt"
os.makedirs(OUT, exist_ok=True)

with open(r"C:\Replicalm\tests\pool_traits.json", encoding="utf-8") as fh:
    pool = json.load(fh)
cands = [r for r in pool if r["tile"] not in ALREADY]
cands.sort(key=lambda r: -r["density"])
cands = cands[:N_CLIP]

print("clipping %d candidates, ordered by header density\n" % len(cands))
print("%-8s %-14s %8s %8s %10s" % ("tile", "region", "header", "clip", "points"))
print("-" * 54)
recs = []
for r in cands:
    t = {"region": r["region"], "tile": r["tile"], "las": r["las"],
         "reference_dem": r["reference_dem"], "chosen": "hyper-dense hunt"}
    try:
        rec = clip.prepare(t, OUT, win_m=400.0, verbose=False)
    except TypeError:
        rec = clip.prepare(t, OUT, win_m=400.0)
    except Exception as e:
        print("%-8s %-14s FAILED %s" % (r["tile"], r["region"], str(e)[:30]))
        continue
    area = rec["window"]["window_m"] ** 2
    rec["clip_density"] = rec["points"] / area
    rec["header_density"] = r["density"]
    rec["multi_return_frac"] = r["multi_return_frac"]
    recs.append(rec)
    print("%-8s %-14s %8.2f %8.2f %10d"
          % (r["tile"], r["region"], r["density"], rec["clip_density"],
             rec["points"]))

recs.sort(key=lambda r: -r["clip_density"])
print("\nclassifying the top %d by clip density\n" % N_CLASSIFY)
print("%-8s %10s %10s %10s %8s" % ("tile", "clip", "ground", "ground frac",
                                   "vs 8.54"))
print("-" * 52)
cfg = replace(PRESETS["ncalm"],
              passes=[GroundPass(algorithm="smrf", slope=0.05, threshold_m=0.5)],
              remove_low_noise=False, remove_outliers=False)
found = []
for rec in recs[:N_CLASSIFY]:
    area = rec["window"]["window_m"] ** 2
    gp = os.path.join(OUT, "%s_ground.las" % rec["tile"])
    try:
        if not os.path.exists(gp):
            classify.classify_tile(rec["clipped_las"], gp, cfg, verbose=False)
        arr, _ = classify.read_points(gp)
        n = int((arr["Classification"] == 2).sum())
        gd = n / area
        rec["ground_density"] = gd
        found.append(rec)
        print("%-8s %10.2f %10.2f %10.0f%% %7.2fx"
              % (rec["tile"], rec["clip_density"], gd,
                 100 * gd / rec["clip_density"], gd / TARGET))
    except Exception as e:
        print("%-8s FAILED %s" % (rec["tile"], str(e)[:40]))

with open(os.path.join(OUT, "index.json"), "w", encoding="utf-8") as fh:
    json.dump(found, fh, indent=1, default=float)

best = max(found, key=lambda r: r["ground_density"]) if found else None
print()
if best and best["ground_density"] >= 6.0:
    print("%s reaches %.2f ground pts/m2: dense enough to exercise the ceiling."
          % (best["tile"], best["ground_density"]))
    print("run: pilot_sweep.py %s %s"
          % (OUT, r"C:\Replicalm\tests\dense2"))
elif best:
    print("best found is %s at %.2f ground pts/m2, against l0s444's %.2f."
          % (best["tile"], best["ground_density"], TARGET))
    print("The pool may simply not hold a second tile of that density.")
else:
    print("no candidate classified successfully")
