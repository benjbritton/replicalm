r"""Sweep the kriging: which radius and neighbour count preserve terrain detail?

Classification fixed at the best combination the first sweep found -- slope 0.10,
threshold 1.0 m -- since those are settled and re-running them costs minutes and
changes nothing.
"""
import json, sys, time
from dataclasses import replace
sys.path.insert(0, r"C:\Replicalm\src")
from replicalm import calibrate
from replicalm.config import PRESETS

with open(r"C:\Replicalm\tests\clips\index.json", encoding="utf-8") as fh:
    clips = {c["tile"]: c for c in json.load(fh)}
c = clips["l0s395"]

base = PRESETS["ncalm"]
cfg = replace(base, passes=[replace(p, slope=0.10, threshold_m=1.0)
                            for p in base.passes])
print("classification fixed: slope %.2f, threshold %.1f m\n"
      % (cfg.passes[0].slope, cfg.passes[0].threshold_m))

t0 = time.time()
res = calibrate.sweep_interpolation(
    las_path=c["clipped_las"],
    reference_dem=c["clipped_reference"],
    out_dir=r"C:\Replicalm\tests\sweep_interp",
    cfg=cfg,
    radii=(5.0, 10.0, 20.0),
    max_points=(8, 16, 32, 64))
print("\ntotal %.0f s" % (time.time() - t0))
