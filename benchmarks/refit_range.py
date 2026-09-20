r"""Where does the correlation range actually sit, once the fit is allowed to find it?

THE OBSERVATION
---------------
Several tiles reported a fitted range of 14.39 m to four significant figures --
l0s395 across a thinning ladder that cut its point count by 56%, l4s478 under
three different classifiers. A fitted parameter does not hold four figures while
its data changes that much.

`fit_variogram` searches the range over np.linspace(lag[1], lag[-1] * 1.5, 24).
With max_lag=10 and 12 bins, lag[-1] is about 9.59, so the search cannot return
more than about 14.39. Those tiles were pinned to the ceiling of the search, and
`radius_from_variogram` passed that straight through to the search radius. The
radius on those tiles was set by a fitting parameter, not by the terrain.

WHAT THIS PASS DOES
-------------------
Refits each tile's ground returns at max_lag of 5, 10, 20, 40, 60 and 80 m and
records the fitted range beside the ceiling the search allowed. Three outcomes
are distinguishable:

  range plateaus            a real correlation length exists and the 10 m fit
                            was simply too short to see it. The rule stands and
                            only its max_lag needs raising.
  range tracks the ceiling  there is no range to find. The semivariance is still
                            climbing at 80 m, which is the signature of a
                            regional gradient rather than of correlation, and
                            the derived radius is a tunable constant wearing a
                            measurement's clothes.
  mixed by tile             the split between tiles is real, and which case a
                            tile falls into is itself the thing to detect.

The third is what universal.py already suspects from the calibration set, where
ranges climbed from 1.76 to 59.99 m as the window went from 5 to 80 m. This pass
puts numbers on it for eight tiles instead of three, and does it on the ground
files the sweeps already classified, so nothing is reclassified.

SUBSAMPLING
-----------
query_pairs at 80 m on 600,000 points is on the order of 10^10 pairs. Every fit
here therefore runs on the same random subsample of 20,000 points, so the
comparison across max_lag is like for like; the absolute semivariances are
unaffected by subsampling, only the pair counts are.
"""
import glob, json, os, sys
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import ROOT, SRC, DEM_ROOT, LAS_ROOT, TESTS
sys.path.insert(0, SRC)
import numpy as np
from replicalm import classify, kriging as K

SUBSAMPLE = 20000
SEED = 20260919
LAGS = [5.0, 10.0, 20.0, 40.0, 60.0, 80.0]
OUT = TESTS + r"\\refit_range.json"

# the ground files the sweeps already wrote, calibration and pilot alike
SOURCES = [
    ("l8s431", "calibration", TESTS + r"\\thin_val_l8s431\l8s431_ground.las"),
    ("l0s395", "calibration", TESTS + r"\\thin_val_l0s395\l0s395_ground.las"),
    ("l0s444", "calibration", TESTS + r"\\thin\ground.las"),
]
for p in sorted(glob.glob(TESTS + r"\\pilot\*_ground.las")):
    b = os.path.basename(p)
    tile = b.split("_")[0]
    cls = b[len(tile) + 1:-len("_ground.las")]
    SOURCES.append((tile, "pilot/" + cls, p))
for p in sorted(glob.glob(TESTS + r"\\dense\*_ground.las")):
    b = os.path.basename(p)
    tile = b.split("_")[0]
    cls = b[len(tile) + 1:-len("_ground.las")]
    SOURCES.append((tile, "dense/" + cls, p))

rng = np.random.default_rng(SEED)
rows = []
for tile, origin, path in SOURCES:
    if not os.path.exists(path):
        print("%-8s %-18s missing: %s" % (tile, origin, path))
        continue
    # one classifier per tile is enough for a question about terrain; the pilot
    # showed the fitted range moving by under 0.01 m between classifiers
    if origin.startswith("pilot/") and not origin.endswith("smrf_t050"):
        continue
    if origin.startswith("dense/") and not origin.endswith("smrf_t050"):
        continue
    arr, _ = classify.read_points(path)
    g = arr[arr["Classification"] == 2]
    if len(g) < 1000:
        print("%-8s %-18s only %d ground points" % (tile, origin, len(g)))
        continue
    idx = rng.choice(len(g), size=min(SUBSAMPLE, len(g)), replace=False)
    x, y, z = g["X"][idx], g["Y"][idx], g["Z"][idx]

    print("\n%-8s %-18s %d ground points, %d subsampled"
          % (tile, origin, len(g), len(idx)))
    print("   %8s %10s %10s %10s %9s %9s"
          % ("max_lag", "range m", "ceiling m", "pegged", "nugget", "sill"))
    for ml in LAGS:
        try:
            v = K.fit_variogram(x, y, z, model="spherical", max_lag=ml)
            nug, sill, rg = v["params"]
            ceiling = float(v["lags"][-1]) * 1.5
            pegged = rg >= ceiling * 0.995
            rows.append({"tile": tile, "origin": origin, "max_lag": ml,
                         "range_m": float(rg), "ceiling_m": ceiling,
                         "pegged": bool(pegged), "nugget": float(nug),
                         "sill": float(sill),
                         "semivariance": v["semivariance"],
                         "lags": v["lags"]})
            print("   %8.0f %10.2f %10.2f %10s %9.4f %9.4f"
                  % (ml, rg, ceiling, "yes" if pegged else "no", nug, sill))
        except Exception as e:
            print("   %8.0f  FAILED %s" % (ml, str(e)[:50]))

with open(OUT, "w", encoding="utf-8") as fh:
    json.dump(rows, fh, indent=1, default=float)

print("\n\n%s\nVERDICT PER TILE\n%s" % ("=" * 72, "=" * 72))
print("%-8s %-18s %28s  %s" % ("tile", "origin", "range at 5/10/20/40/60/80 m",
                               "reading"))
print("-" * 100)
for tile, origin, _ in SOURCES:
    rs = [r for r in rows if r["tile"] == tile and r["origin"] == origin]
    if not rs:
        continue
    rs.sort(key=lambda r: r["max_lag"])
    seq = "/".join("%.1f" % r["range_m"] for r in rs)
    # Not hitting the search ceiling does not make a range real. A range that
    # climbs with every window is the window's, whether or not it stopped short
    # of the bound -- so the test is whether the range stops growing, and
    # whether the sill stops growing with it. A stationary field gives both.
    ratios = [r["range_m"] / r["max_lag"] for r in rs]
    grow_r = rs[-1]["range_m"] / rs[0]["range_m"]
    grow_s = (rs[-1]["sill"] / rs[0]["sill"]) if rs[0]["sill"] > 0 else float("inf")
    if grow_r < 2.0 and grow_s < 2.0:
        verdict = "range and sill both stable: a correlation length exists"
    elif grow_s < 2.0:
        verdict = "sill stable, range grows %.0fx: fit unstable, field bounded" % grow_r
    else:
        verdict = ("range grows %.0fx and sill %.0fx across the windows: "
                   "unbounded, no correlation length to find"
                   % (grow_r, grow_s))
    print("%-8s %-18s %28s  ratio %.2f-%.2f  %s"
          % (tile, origin, seq, min(ratios), max(ratios), verdict))
print("\nwrote %s" % OUT)
