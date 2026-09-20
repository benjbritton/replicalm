r"""Diagnose the residuals already on disk: noise, or their smoothing?"""
import glob, json, os, sys, time
sys.path.insert(0, r"C:\Replicalm\src")
import numpy as np
from replicalm import residual

with open(r"C:\Replicalm\tests\clips\index.json", encoding="utf-8") as fh:
    clips = {c["tile"]: c for c in json.load(fh)}

M = r"C:\Replicalm\tests\methods"
print("%-8s %-9s %8s %8s %8s %9s %8s  %s"
      % ("tile", "method", "resid sd", "nugget", "sill", "nug/sill",
         "range m", "Moran I"))
print("-" * 92)
rows = []
for tile in ("l0s395", "l8s431", "l0s444"):
    ref = clips[tile]["clipped_reference"]
    for p in sorted(glob.glob(os.path.join(M, "%s_*.tif" % tile))):
        name = os.path.basename(p)[len(tile) + 1:-4]
        try:
            r = residual.compare(p, ref, max_lag_m=60.0)
            rows.append({"tile": tile, "method": name, **r})
            print("%-8s %-9s %8.3f %8.5f %8.5f %9.2f %8s  %6.3f"
                  % (tile, name, r["residual_sd_m"], r["nugget"], r["sill"],
                     r["nugget_to_sill"],
                     "%.1f" % r["range_m"] if r["range_m"] else "none",
                     r["morans_i"]))
        except Exception as e:
            print("%-8s %-9s FAILED %s" % (tile, name, str(e)[:44]))
    print()

with open(os.path.join(M, "residual_diagnosis.json"), "w",
          encoding="utf-8") as fh:
    json.dump(rows, fh, indent=1, default=float)

print("verdicts:")
for r in rows:
    print("  %-8s %-9s %s" % (r["tile"], r["method"], r["verdict"]))
