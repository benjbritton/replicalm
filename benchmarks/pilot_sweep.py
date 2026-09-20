r"""Do the rules hold on tiles that had no part in making them?

WHAT IS BEING TESTED
--------------------
Three rules came out of the calibration set, and all three were read off three
tiles, two of which were drawn from one Centro survey:

  native-first        do not thin. Full density won on RMSE for every tile that
                      was not degenerate.
  the 4.4 ceiling     thin only when the kriging system fails, which on the
                      calibration set meant only l0s444 at 8.54 pts/m2.
  a derived radius    take the search radius from the fitted variogram range
                      rather than from the documented 20 m in the source method.

The pilot draws five South_NFI tiles, a survey neither calibration tile came
from. A rule that only holds where it was made is not a rule.

WHAT WOULD FALSIFY EACH
-----------------------
  native-first        a tile where thinning lowers RMSE without the fallback
                      having spiked first.
  the ceiling         a tile with a high fallback that thinning does not fix,
                      or a well-conditioned tile that still wants thinning.
  the range split     l8s431 fitted 3.0 m and l0s395 14.4 m. Two tiles is not a
                      distribution. If the five land anywhere across that
                      interval it is a continuum, and the split was an artefact
                      of a sample of two.

CLASSIFICATION IS NOT HELD FIXED
--------------------------------
It cannot be: the calibration set produced a different winner on each of its
three tiles, so there is no single setting to carry forward. Each pilot tile
runs all three leading candidates and the best by RMSE is the one whose
interpolation figures are reported. Whether that choice is still tile-dependent
in a new survey is itself part of what the pilot measures. If one candidate wins
five times, the deployment gains a default it does not currently have.
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

# argv lets the hyper-dense supplementary run take the identical code path
# rather than a copy of it, so a difference in its results cannot be a
# difference in the script.
CLIPS = sys.argv[1] if len(sys.argv) > 1 else TESTS + r"\\pilot_clips"
OUT = sys.argv[2] if len(sys.argv) > 2 else TESTS + r"\\pilot"
os.makedirs(OUT, exist_ok=True)
with open(os.path.join(CLIPS, "index.json"), encoding="utf-8") as fh:
    clips = json.load(fh)

CANDIDATES = [
    ("smrf_t050", dict(algorithm="smrf", slope=0.05, threshold_m=0.5)),
    ("smrf_t025", dict(algorithm="smrf", slope=0.05, threshold_m=0.25)),
    ("csf_r2",    dict(algorithm="csf", csf_rigidness=2, csf_threshold_m=0.5)),
]
# l0s444 fell back on 89% of cells; every well-conditioned tile was under 0.2%.
# Anything above this is a spike by any reading, and nothing observed so far
# lands between.
FALLBACK_SPIKE = 0.05
CEILING = 4.43


def interpolate_at(x, y, z, g, wkt, ref, dem_p, var_p):
    """Fit, krige and score one point set on the tile's own reference grid."""
    v = K.fit_variogram(x, y, z, model="spherical", max_lag=10.0)
    radius = K.radius_from_variogram(v)
    dem, info = K.krige_grid(x, y, z, g, radius=radius, max_points=16,
                             variogram=v, verbose=False)
    K.write_geotiff(dem, g, wkt, dem_p)
    var = info["variance"]
    K.write_geotiff(var, g, wkt, var_p)
    sc = calibrate.compare_to_reference(dem_p, ref)
    vmask = np.isfinite(var) & (var > -9000)
    return {"range_m": float(v["params"][2]), "radius_m": float(radius),
            "fallback": float(info["fallback_fraction"]),
            "rmse_m": float(sc["rmse_m"]),
            "ratio": float(sc["complexity_candidate"] /
                           sc["complexity_reference"]),
            "complexity_reference": float(sc["complexity_reference"]),
            "variance_coverage": float(vmask.mean())}


results = []
for c in clips:
    tile = c["tile"]
    g = G.grid_from_raster(c["clipped_reference"])
    wkt = interpolate.source_srs(c["clipped_las"])
    area = (g.bounds[2] - g.bounds[0]) * (g.bounds[3] - g.bounds[1])
    print("\n%s  %s" % ("=" * 78, tile))

    arms = []
    for name, params in CANDIDATES:
        t0 = time.time()
        cfg = replace(PRESETS["ncalm"], passes=[GroundPass(**params)],
                      remove_low_noise=False, remove_outliers=False)
        gnd_p = os.path.join(OUT, "%s_%s_ground.las" % (tile, name))
        try:
            if not os.path.exists(gnd_p):
                classify.classify_tile(c["clipped_las"], gnd_p, cfg,
                                       verbose=False)
            arr, _ = classify.read_points(gnd_p)
            gnd = arr[arr["Classification"] == 2]
            X, Y, Z = gnd["X"], gnd["Y"], gnd["Z"]
            dens = len(Z) / area
            r = interpolate_at(X, Y, Z, g, wkt, c["clipped_reference"],
                               os.path.join(OUT, "%s_%s.tif" % (tile, name)),
                               os.path.join(OUT, "%s_%s_var.tif" % (tile, name)))
            r.update({"tile": tile, "classifier": name, "params": params,
                      "ground_points": int(len(Z)), "density": float(dens),
                      "thinned": False, "seconds": round(time.time() - t0, 1)})
            arms.append(r)
            print("  %-10s %8.2f pts/m2  range %6.2f m  radius %5.1f  "
                  "fallback %5.1f%%  rmse %.3f  ratio %.2fx"
                  % (name, dens, r["range_m"], r["radius_m"],
                     100 * r["fallback"], r["rmse_m"], r["ratio"]))
        except Exception as e:
            print("  %-10s FAILED %s" % (name, str(e)[:60]))

    if not arms:
        continue
    best = min(arms, key=lambda r: r["rmse_m"])
    print("  best by rmse: %s" % best["classifier"])

    # The ceiling is a remedy, not a default. It is applied only where the
    # native-density solve actually failed, which is the rule as written.
    if best["fallback"] > FALLBACK_SPIKE:
        print("  fallback %.1f%% exceeds %.0f%%: applying the %.2f pts/m2 cap"
              % (100 * best["fallback"], 100 * FALLBACK_SPIKE, CEILING))
        gnd_p = os.path.join(OUT, "%s_%s_ground.las"
                             % (tile, best["classifier"]))
        arr, _ = classify.read_points(gnd_p)
        gnd = arr[arr["Classification"] == 2]
        X, Y, Z = gnd["X"], gnd["Y"], gnd["Z"]
        sel = K.thin_declustered(X, Y, CEILING)
        try:
            r = interpolate_at(
                X[sel], Y[sel], Z[sel], g, wkt, c["clipped_reference"],
                os.path.join(OUT, "%s_capped.tif" % tile),
                os.path.join(OUT, "%s_capped_var.tif" % tile))
            r.update({"tile": tile, "classifier": best["classifier"],
                      "params": best["params"], "ground_points": int(sel.size),
                      "density": float(sel.size / area), "thinned": True})
            arms.append(r)
            print("  capped     %8.2f pts/m2  range %6.2f m  radius %5.1f  "
                  "fallback %5.1f%%  rmse %.3f  ratio %.2fx"
                  % (r["density"], r["range_m"], r["radius_m"],
                     100 * r["fallback"], r["rmse_m"], r["ratio"]))
        except Exception as e:
            print("  capped     FAILED %s" % str(e)[:60])
    else:
        print("  fallback %.2f%% is within tolerance: no cap, native kept"
              % (100 * best["fallback"]))
    results.extend(arms)

with open(os.path.join(OUT, "pilot.json"), "w", encoding="utf-8") as fh:
    json.dump(results, fh, indent=1, default=float)

print("\n\n%s\nSUMMARY: best arm per tile, native density\n%s"
      % ("=" * 78, "=" * 78))
print("%-8s %-10s %9s %8s %8s %9s %8s %9s"
      % ("tile", "classifier", "pts/m2", "range m", "radius", "fallback",
         "rmse m", "cplx"))
print("-" * 78)
for c in clips:
    rs = [r for r in results if r["tile"] == c["tile"] and not r["thinned"]]
    if not rs:
        continue
    b = min(rs, key=lambda r: r["rmse_m"])
    print("%-8s %-10s %9.2f %8.2f %8.1f %8.2f%% %8.3f %8.2fx"
          % (b["tile"], b["classifier"], b["density"], b["range_m"],
             b["radius_m"], 100 * b["fallback"], b["rmse_m"], b["ratio"]))
print("\nwrote %s" % os.path.join(OUT, "pilot.json"))
