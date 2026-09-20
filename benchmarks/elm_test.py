r"""Does ELM collapse the nugget on l0s444?

Same tile, same classification parameters, same kriging. The only difference is
whether low noise is removed before the ground filter runs. If the residual is
noise we are introducing, the nugget-to-sill ratio should fall from 0.92.
"""
import json, os, sys, time
from dataclasses import replace
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import ROOT, SRC, DEM_ROOT, LAS_ROOT, TESTS
sys.path.insert(0, SRC)
import numpy as np
from replicalm import calibrate, classify, grid as G, kriging as K, residual, interpolate
from replicalm.config import PRESETS, GroundPass

OUT = TESTS + r"\\elm"
os.makedirs(OUT, exist_ok=True)
with open(TESTS + r"\\clips\index.json", encoding="utf-8") as fh:
    clips = {c["tile"]: c for c in json.load(fh)}

for tile in ("l0s444", "l0s395", "l8s431"):
    c = clips[tile]
    g = G.grid_from_raster(c["clipped_reference"])
    wkt = interpolate.source_srs(c["clipped_las"])
    area = (g.bounds[2] - g.bounds[0]) * (g.bounds[3] - g.bounds[1])
    radius = {"l0s395": 5.0, "l8s431": 20.0, "l0s444": 5.0}[tile]
    print("=== %s, radius %.0f m ===" % (tile, radius))
    print("%-14s %9s %8s %9s %9s %9s %8s %8s"
          % ("noise stages", "ground", "pts/m2", "rmse m", "cplx cand",
             "nug/sill", "range m", "fallback"))
    for label, on in (("off", False), ("ELM+outlier", True)):
        base = replace(PRESETS["ncalm"],
                       passes=[GroundPass(algorithm="smrf", slope=0.10,
                                          threshold_m=1.0)],
                       remove_low_noise=on, remove_outliers=on)
        tag = "%s_%s" % (tile, "on" if on else "off")
        las_p = os.path.join(OUT, "%s.las" % tag)
        dem_p = os.path.join(OUT, "%s.tif" % tag)
        try:
            classify.classify_tile(c["clipped_las"], las_p, base, verbose=False)
            arr, _ = classify.read_points(las_p)
            gnd = arr[arr["Classification"] == 2]
            v = K.fit_variogram(gnd["X"], gnd["Y"], gnd["Z"],
                                model="spherical", max_lag=radius)
            dem, info = K.krige_grid(gnd["X"], gnd["Y"], gnd["Z"], g,
                                     radius=radius, max_points=16,
                                     variogram=v, verbose=False)
            K.write_geotiff(dem, g, wkt, dem_p)
            sc = calibrate.compare_to_reference(dem_p, c["clipped_reference"])
            rs = residual.compare(dem_p, c["clipped_reference"], max_lag_m=60.0)
            print("%-14s %9d %8.2f %9.3f %9.4f %9.2f %8s %7.1f%%"
                  % (label, len(gnd), len(gnd)/area, sc["rmse_m"],
                     sc["complexity_candidate"], rs["nugget_to_sill"],
                     "%.1f" % rs["range_m"] if rs["range_m"] else "none",
                     100 * info.get("fallback_fraction", 0)))
        except Exception as e:
            print("%-14s FAILED %s" % (label, str(e)[:56]))
    print("   reference complexity %.4f\n"
          % calibrate.compare_to_reference(dem_p, c["clipped_reference"])["complexity_reference"])
