r"""First run of the classifier on a real tile, to prove the plumbing.

Uses the LAS sitting in C:\g1 rather than synthetic data, because the thing
most likely to be wrong is an assumption about the file, not the algorithm.
"""
import os, sys, time
sys.path.insert(0, r"C:\Replicalm\src")

from replicalm.config import PRESETS
from replicalm import classify

SRC = r"D:\GLiHT_LAS_orig\Yuc_Campeche\Campeche\AMIGACarb_Chiap_Campeche_NFI_Apr2013_l2s505.las"
OUT = r"C:\Replicalm\tests\out"
os.makedirs(OUT, exist_ok=True)

print("source: %s  (%.1f MB)" % (os.path.basename(SRC),
                                 os.path.getsize(SRC) / 1e6))
arr, meta = classify.read_points(SRC)
print("points: %d" % len(arr))
print("fields: %s" % ", ".join(arr.dtype.names[:12]))
import numpy as np
print("x %.1f-%.1f   y %.1f-%.1f   z %.2f-%.2f"
      % (arr["X"].min(), arr["X"].max(), arr["Y"].min(), arr["Y"].max(),
         arr["Z"].min(), arr["Z"].max()))
cls, cnt = np.unique(arr["Classification"], return_counts=True)
print("incoming classification: %s"
      % ", ".join("%d:%d" % (c, n) for c, n in zip(cls, cnt)))
print("extent: %.0f x %.0f m" % (arr["X"].max() - arr["X"].min(),
                                 arr["Y"].max() - arr["Y"].min()))

cfg = PRESETS["ncalm"]
print("\nrunning the NCALM translation, two SMRF passes:")
for p in cfg.passes:
    print("   %s window %.0f  slope %.3f  threshold %.1f   [%s]"
          % (p.algorithm, p.window_m, p.slope, p.threshold_m, p.stands_in_for))

t0 = time.time()
res = classify.classify_tile(SRC, os.path.join(OUT, "l2s505_ground.las"), cfg)
print("\nelapsed %.1f s" % (time.time() - t0))
for k, v in res.items():
    print("  %-12s %s" % (k, v))

