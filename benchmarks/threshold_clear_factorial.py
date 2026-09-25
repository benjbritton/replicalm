r"""Do the SMRF threshold and the Clear rule do the same job?

The class sweep showed that tightening the elevation threshold from 0.5 m to
0.25 m collapses the positive bias against the reference -- +0.0245 to +0.0022
on l0s444, +0.0294 to +0.0022 on l8s431 -- while dropping 9 to 13% of ground
returns. Clear removes 17.75% for what is described as the same reason: low
vegetation accepted as ground. One works at classification, the other after it.

If they are the same instrument, running both should either find nothing left
for the second to remove, or overshoot into negative bias by cutting terrain.

Held at the locked configuration throughout: slope 0.1584, working cell 0.50 m,
16 neighbours, minimum 3. An earlier run of this comparison built its passes by
hand and inherited GroundPass's field defaults instead, putting SMRF on a 1.0 m
working cell. Those numbers were wrong and are superseded by these.
"""
import json, os, sys
from dataclasses import replace
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, r"C:\Replicalm\benchmarks")
from paths import SRC, TESTS
sys.path.insert(0, SRC)
import numpy as np
from replicalm import calibrate, classify, cleanup, grid as G, kriging as K, interpolate
from replicalm.config import PRESETS, GroundPass

OUT = os.path.join(TESTS, "threshold_clear")
os.makedirs(OUT, exist_ok=True)
with open(os.path.join(TESTS, "clips", "index.json"), encoding="utf-8") as fh:
    clips = {c["tile"]: c for c in json.load(fh)}
RADIUS = {"l0s395": 5.0, "l8s431": 20.0, "l0s444": 5.0}
# slope comes from the locked pass; nothing here varies it

print("%-8s %-14s %9s %7s %9s %9s %9s %8s"
      % ("tile", "arm", "ground", "clear%", "rmse m", "bias m", "mae m", "cplx"))
print("-" * 82)
rows = []
for tile in ("l0s395", "l8s431", "l0s444"):
    c = clips[tile]
    g = G.grid_from_raster(c["clipped_reference"])
    wkt = interpolate.source_srs(c["clipped_las"])
    area = (g.bounds[2] - g.bounds[0]) * (g.bounds[3] - g.bounds[1])
    base = PRESETS["ncalm"]
    for thr in (0.5, 0.25):
        # Derived from the locked pass, not built fresh: GroundPass's own field
        # defaults are cell 1.0 m and slope 0.158, where the locked pass uses
        # 0.50 m and 0.1584. Constructing one by hand and naming only two fields
        # silently ran this comparison at the wrong working cell the first time.
        cfg = replace(base, passes=[replace(base.passes[0], threshold_m=thr)])
        tag = "%s_t%s" % (tile, str(thr).replace(".", ""))
        las_p = os.path.join(OUT, tag + ".las")
        try:
            classify.classify_tile(c["clipped_las"], las_p, cfg, verbose=False)
            arr, _ = classify.read_points(las_p)
        except Exception as e:
            print("%-8s t=%-12s FAILED %s" % (tile, thr, str(e)[:40])); continue
        for clear_on in (False, True):
            name = "t%.2f %s" % (thr, "clear" if clear_on else "raw")
            a, report = arr, None
            removed = 0.0
            if clear_on:
                # clear() returns (array, report); the report carries the
                # demoted fraction directly, which is the quantity in question
                a, report = cleanup.clear(arr.copy(), patch=base.clean_patch_m,
                                          height=base.clean_height_m,
                                          percentile=base.clean_percentile,
                                          verbose=False)
                removed = 100.0 * report["fraction"]
            gnd = a[a["Classification"] == 2]
            if len(gnd) < 5000:
                print("%-8s %-14s only %d ground" % (tile, name, len(gnd))); continue
            try:
                v = K.fit_variogram(gnd["X"], gnd["Y"], gnd["Z"], model="spherical",
                                    max_lag=RADIUS[tile])
                dem, _ = K.krige_grid(gnd["X"], gnd["Y"], gnd["Z"], g,
                                      radius=RADIUS[tile], max_points=16,
                                      variogram=v, verbose=False)
                dem_p = os.path.join(OUT, "%s_%s.tif" % (tag, "clear" if clear_on else "raw"))
                K.write_geotiff(dem, g, wkt, dem_p)
                sc = calibrate.compare_to_reference(dem_p, c["clipped_reference"])
            except Exception as e:
                print("%-8s %-14s FAILED %s" % (tile, name, str(e)[:40])); continue
            ratio = sc["complexity_candidate"] / sc["complexity_reference"]
            rows.append({"tile": tile, "arm": name, "threshold_m": thr,
                         "clear": clear_on,
                         "slope": cfg.passes[0].slope,
                         "smrf_cell_m": cfg.passes[0].cell_m,
                         "ground": int(len(gnd)), "clear_removed_pct": round(removed, 2),
                         "clear_report": report,
                         "density": len(gnd) / area, "ratio": ratio, **sc})
            print("%-8s %-14s %9d %6.2f%% %9.4f %+9.4f %9.4f %7.2fx"
                  % (tile, name, len(gnd), removed, sc["rmse_m"], sc["bias_m"],
                     sc["mae_m"], ratio))
    print()
with open(os.path.join(OUT, "threshold_clear.json"), "w", encoding="utf-8") as fh:
    json.dump(rows, fh, indent=1, default=float)
print("wrote %s (%d arms)" % (os.path.join(OUT, "threshold_clear.json"), len(rows)))
