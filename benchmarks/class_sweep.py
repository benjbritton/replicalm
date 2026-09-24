r"""Is the excess terrain complexity coming from the ground filter?

Fifteen interpolation settings all produced 1.6 to 2.8 times the reference's
terrain complexity. If that is returns we are calling ground and the original
did not -- low vegetation, understorey, noise -- then tightening the filter
should bring the texture down, and no interpolator setting ever would.

Scored on complexity rather than RMSE, because RMSE is already close on every
tile and is not the quantity in question.
"""
import json, os, sys, time
from dataclasses import replace
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import ROOT, SRC, DEM_ROOT, LAS_ROOT, TESTS
sys.path.insert(0, SRC)
import numpy as np
from replicalm import calibrate, classify, grid as G, kriging as K, interpolate
from replicalm.config import PRESETS, GroundPass

OUT = TESTS + r"\\class_sweep"
os.makedirs(OUT, exist_ok=True)
with open(TESTS + r"\\clips\index.json", encoding="utf-8") as fh:
    clips = {c["tile"]: c for c in json.load(fh)}

# progressively stricter ground extraction
VARIANTS = [
    ("smrf s.10 t1.0", dict(algorithm="smrf", slope=0.10, threshold_m=1.0)),
    ("smrf s.10 t0.5", dict(algorithm="smrf", slope=0.10, threshold_m=0.5)),
    ("smrf s.05 t0.5", dict(algorithm="smrf", slope=0.05, threshold_m=0.5)),
    ("smrf s.05 t0.25", dict(algorithm="smrf", slope=0.05, threshold_m=0.25)),
    ("smrf s.02 t0.25", dict(algorithm="smrf", slope=0.02, threshold_m=0.25)),
    # PMF (Zhang et al. 2003) on the same slope/threshold ladder as SMRF, so the
    # two morphological families are compared at matched settings rather than at
    # each one's own defaults.
    ("pmf s.10 t1.0",  dict(algorithm="pmf", slope=0.10, threshold_m=1.0)),
    ("pmf s.10 t0.5",  dict(algorithm="pmf", slope=0.10, threshold_m=0.5)),
    ("pmf s.05 t0.5",  dict(algorithm="pmf", slope=0.05, threshold_m=0.5)),
    ("pmf s.05 t0.25", dict(algorithm="pmf", slope=0.05, threshold_m=0.25)),
    ("pmf s.02 t0.25", dict(algorithm="pmf", slope=0.02, threshold_m=0.25)),
    ("csf rigid2",     dict(algorithm="csf", csf_rigidness=2,
                            csf_threshold_m=0.5)),
    ("csf rigid3",     dict(algorithm="csf", csf_rigidness=3,
                            csf_threshold_m=0.3)),
]

RADIUS = {"l0s395": 5.0, "l8s431": 20.0, "l0s444": 5.0}   # each tile's own best

print("%-8s %-16s %9s %8s %9s %10s %8s"
      % ("tile", "filter", "ground", "pts/m2", "rmse m", "cplx cand", "ratio"))
print("-" * 78)
rows = []
for tile in ("l0s395", "l8s431", "l0s444"):
    c = clips[tile]
    g = G.grid_from_raster(c["clipped_reference"])
    wkt = interpolate.source_srs(c["clipped_las"])
    area = (g.bounds[2] - g.bounds[0]) * (g.bounds[3] - g.bounds[1])
    base = PRESETS["ncalm"]
    ref_c = None
    for name, kw in VARIANTS:
        tag = "%s_%s" % (tile, name.replace(" ", "").replace(".", ""))
        las_p = os.path.join(OUT, "%s.las" % tag)
        dem_p = os.path.join(OUT, "%s.tif" % tag)
        cfg = replace(base, passes=[GroundPass(**kw)])
        try:
            classify.classify_tile(c["clipped_las"], las_p, cfg, verbose=False)
            arr, _ = classify.read_points(las_p)
            gnd = arr[arr["Classification"] == 2]
            if len(gnd) < 5000:
                raise RuntimeError("only %d ground points" % len(gnd))
            v = K.fit_variogram(gnd["X"], gnd["Y"], gnd["Z"], model="spherical",
                                max_lag=RADIUS[tile])
            dem, _ = K.krige_grid(gnd["X"], gnd["Y"], gnd["Z"], g,
                                  radius=RADIUS[tile], max_points=16,
                                  variogram=v, verbose=False)
            K.write_geotiff(dem, g, wkt, dem_p)
            sc = calibrate.compare_to_reference(dem_p, c["clipped_reference"])
            ref_c = sc["complexity_reference"]
            ratio = sc["complexity_candidate"] / ref_c
            rows.append({"tile": tile, "filter": name,
                         "ground": int(len(gnd)), "density": len(gnd)/area,
                         "ratio": ratio, **sc})
            print("%-8s %-16s %9d %8.2f %9.3f %10.4f %7.2fx"
                  % (tile, name, len(gnd), len(gnd)/area, sc["rmse_m"],
                     sc["complexity_candidate"], ratio))
        except Exception as e:
            print("%-8s %-16s FAILED %s" % (tile, name, str(e)[:40]))
    if ref_c:
        print("%-8s %-16s reference complexity %.4f" % ("", "", ref_c))
    print()

with open(os.path.join(OUT, "class_sweep.json"), "w", encoding="utf-8") as fh:
    json.dump(rows, fh, indent=1, default=float)

print("closest texture match per tile:")
for tile in ("l0s395", "l8s431", "l0s444"):
    a = sorted([r for r in rows if r["tile"] == tile],
               key=lambda r: abs(r["ratio"] - 1.0))
    if a:
        print("  %-8s %-16s ratio %.2fx  rmse %.3f m  %.2f pts/m2"
              % (tile, a[0]["filter"], a[0]["ratio"], a[0]["rmse_m"],
                 a[0]["density"]))
