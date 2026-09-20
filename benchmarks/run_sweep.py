r"""First calibration sweep: does any SMRF setting reproduce the reference DEM?

Six combinations on South_GLAS_l0s395, scored on elevation RMSE combined with
loss of terrain complexity, against the DEM the TerraScan and Surfer pipeline
produced for the same ground.
"""
import json, sys, time
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import ROOT, SRC, DEM_ROOT, LAS_ROOT, TESTS
sys.path.insert(0, SRC)
from replicalm import calibrate
from replicalm.config import PRESETS

with open(TESTS + r"\\clips\index.json", encoding="utf-8") as fh:
    clips = {c["tile"]: c for c in json.load(fh)}

c = clips["l0s395"]
print("tile %s, %d points in a %.0f m window"
      % (c["tile"], c["points"], c["window"]["window_m"]))
print("reference: %s\n" % c["clipped_reference"])

GRID = {"slope": [0.10, 0.158, 0.213],
        "threshold_m": [1.0, 3.0]}

t0 = time.time()
res = calibrate.sweep(
    las_path=c["clipped_las"],
    reference_dem=c["clipped_reference"],
    out_dir=TESTS + r"\\sweep_l0s395",
    base_cfg=PRESETS["ncalm"],
    param_grid=GRID,
    cell=0.5, radius=20.0, max_points=32)
print("\ntotal %.0f s" % (time.time() - t0))

if res["ranked"]:
    print("\n%-26s %8s %9s %9s %9s %9s"
          % ("combination", "score", "rmse m", "bias m", "cplx loss", "pts/m2"))
    for r in res["ranked"]:
        print("%-26s %8.4f %9.3f %9.3f %9.3f %9.2f"
              % (r["tag"], r["score"], r["rmse_m"], r["bias_m"],
                 r["complexity_loss"], r["density_pts_m2"]))
