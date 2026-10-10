"""Centro GLAS l0s454 through the actual pipeline, Clear profile, 0.50 m.

This replaces my scratch harness with `replicalm.pipeline.process`, which is
what the installed product runs. Two faults in the earlier renders are removed
by doing that rather than by adjusting anything:

  the cleanup       my harness went assign -> outlier -> smrf -> class 2 and
                    stopped, so the low-vegetation stage never ran. That stage
                    exists to remove the fine speckle on flat ground, which is
                    what the earlier panel was covered in.
  the stretch       I had overridden visualize.normalise with a range shared
                    across three surfaces so panels would be comparable. That
                    is wrong for image quality: adjust()'s brightness, contrast
                    and gamma expect a per-layer min-max input. pipeline.process
                    calls visualize.g1, which normalises per layer per tile, as
                    the recipe intends.

`on_grid` adopts the archive's own 0.50 m lattice, so the output is cell for
cell comparable with the delivered product rather than merely the same
nominal resolution.

The comparison then renders three panels from the three G1 images as they are,
with no restretching by me at all, and measures the specific failures predicted
of the Clear profile: how much ground the cleanup removes, whether detail is
lost against the archive, and how the trimmed edge compares in shape.
"""
import json
import os
import sys
import time

import numpy as np
import rasterio

sys.path.insert(0, r"C:\Replicalm\src")
from replicalm import pipeline
from replicalm.config import PRESETS

STEM = "Centro_GLAS_l0s454"
RAW = (r"D:\GLiHT_LAS_orig\Yuc_Centro\Centro_Glas\Centro_Glas"
       r"\AMIGACarb_Yuc_Centro_GLAS_Apr2013_l0s454.las")
ARCH_DEM50 = (r"D:\_Archive_EdgeFixed\Yuc_Centro\Centro_Glas\Centro_Glas"
              r"\Centro_GLAS_l0s454_DEM_0p5m_v1.tif")
OUT = r"C:\Replicalm\test100526\centro454_clear50"
os.makedirs(OUT, exist_ok=True)

cfg = PRESETS["clear"]
p0 = cfg.passes[0]
print("profile Clear")
print("  smrf    window %.0f, slope %.3f, threshold %.2f, scalar %.1f, "
      "cell %.2f" % (p0.window_m, p0.slope, p0.threshold_m, p0.scalar,
                     p0.cell_m))
print("  cleanup height %.3f m over a %.2f m patch, plane-fitted"
      % (getattr(cfg, "clean_height_m", float("nan")),
         getattr(cfg, "clean_patch_m", float("nan"))))
print("  interpolator %s, tension %s, neighbours %s"
      % (cfg.interpolator, cfg.spline_tension, cfg.spline_neighbours))
with rasterio.open(ARCH_DEM50) as r:
    print("  grid adopted from the archive: %d x %d at %.6f m"
          % (r.height, r.width, abs(r.transform.a)))


def say(stage, message, fraction=None):
    print("  [%-10s] %s" % (stage, message), flush=True)


t0 = time.time()
out = pipeline.process(RAW, OUT, cfg=cfg, make_g1=True, progress=say,
                       keep_ground=True, on_grid=ARCH_DEM50,
                       ground_stage="after")
print("")
print("finished in %.1f min" % ((time.time() - t0) / 60.0))
for k in sorted(out):
    v = out[k]
    if isinstance(v, (str, int, float)) or v is None:
        print("  %-28s %s" % (k, v))
with open(os.path.join(OUT, "process_report.json"), "w") as fh:
    json.dump({k: (v if isinstance(v, (str, int, float, list, dict))
                   or v is None else str(v)) for k, v in out.items()},
              fh, indent=2, default=str)
print("")
print("wrote %s" % os.path.join(OUT, "process_report.json"))
