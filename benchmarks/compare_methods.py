r"""Pragmatic against correct: local ordinary kriging versus universal kriging.

The pragmatic route keeps ordinary kriging and stays local, which satisfies
stationarity by staying inside a small neighbourhood. The correct route models
the trend with a drift term and can afford a wider one. Both run on all three
tiles, against the same references, with the same classification.
"""
import json, os, sys, time
from dataclasses import replace
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import ROOT, SRC, DEM_ROOT, LAS_ROOT, TESTS
sys.path.insert(0, SRC)
import numpy as np
from replicalm import calibrate, classify, grid as G, kriging as K, universal as U, interpolate
from replicalm.config import PRESETS

OUT = TESTS + r"\\methods"
os.makedirs(OUT, exist_ok=True)
with open(TESTS + r"\\clips\index.json", encoding="utf-8") as fh:
    clips = json.load(fh)

base = PRESETS["ncalm"]
cfg = replace(base, passes=[replace(p, slope=0.10, threshold_m=1.0)
                            for p in base.passes])

# pragmatic: ordinary kriging, tight neighbourhood, variogram fitted to match
# correct: universal kriging with a linear drift, wider neighbourhood allowed
METHODS = [
    ("OK r5",    dict(kind="ok", radius=5.0,  fit=5.0,  k=16)),
    ("OK r20",   dict(kind="ok", radius=20.0, fit=20.0, k=16)),
    ("UK1 r5",   dict(kind="uk", radius=5.0,  fit=5.0,  k=16, drift=1)),
    ("UK1 r20",  dict(kind="uk", radius=20.0, fit=20.0, k=16, drift=1)),
    ("UK1 r40",  dict(kind="uk", radius=40.0, fit=40.0, k=32, drift=1)),
]

print("%-8s %-9s %8s %9s %10s %9s %8s"
      % ("tile", "method", "rmse m", "cplx ref", "cplx cand", "fallback", "filled%"))
print("-" * 74)
rows = []
for c in clips:
    tile = c["tile"]
    g = G.grid_from_raster(c["clipped_reference"])
    wkt = interpolate.source_srs(c["clipped_las"])
    gnd_path = os.path.join(TESTS + r"\\three_tiles", "%s_ground.las" % tile)
    if not os.path.exists(gnd_path):
        classify.classify_tile(c["clipped_las"], gnd_path, cfg, verbose=False)
    arr, _ = classify.read_points(gnd_path)
    gnd = arr[arr["Classification"] == 2]
    x, y, z = gnd["X"], gnd["Y"], gnd["Z"]

    for name, m in METHODS:
        tag = "%s_%s" % (tile, name.replace(" ", ""))
        p = os.path.join(OUT, "%s.tif" % tag)
        t0 = time.time()
        try:
            v = K.fit_variogram(x, y, z, model="spherical", max_lag=m["fit"])
            if m["kind"] == "ok":
                dem, info = K.krige_grid(x, y, z, g, radius=m["radius"],
                                         max_points=m["k"], variogram=v,
                                         verbose=False)
            else:
                dem, info = U.krige_universal(x, y, z, g, radius=m["radius"],
                                              max_points=m["k"], variogram=v,
                                              drift_order=m["drift"],
                                              verbose=False)
            K.write_geotiff(dem, g, wkt, p)
            sc = calibrate.compare_to_reference(p, c["clipped_reference"])
            fb = info.get("fallback_fraction", 0.0)
            rows.append({"tile": tile, "method": name, "seconds": round(time.time()-t0,1),
                         "fallback_fraction": fb, **sc})
            print("%-8s %-9s %8.3f %9.4f %10.4f %8.1f%% %7.1f%%"
                  % (tile, name, sc["rmse_m"], sc["complexity_reference"],
                     sc["complexity_candidate"], 100 * fb,
                     100 * sc["coverage_candidate"]))
        except Exception as e:
            print("%-8s %-9s FAILED %s" % (tile, name, str(e)[:44]))
    print()

with open(os.path.join(OUT, "methods.json"), "w", encoding="utf-8") as fh:
    json.dump(rows, fh, indent=1, default=float)

print("best per tile by rmse:")
for tile in {r["tile"] for r in rows}:
    a = sorted([r for r in rows if r["tile"] == tile], key=lambda r: r["rmse_m"])
    print("  %-8s %-9s rmse %.3f  complexity %.4f vs %.4f"
          % (tile, a[0]["method"], a[0]["rmse_m"],
             a[0]["complexity_candidate"], a[0]["complexity_reference"]))
