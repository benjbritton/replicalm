r"""Does a filter's surface rise with understory over ground selected as flat?

The earlier regression compared the filters with each other and found PMF 14 cm
below SMRF and CSF, with the gap widening as understory thickens. That says they
diverge; it does not say which is moving. This asks the question that does.

Over a cell chosen for flatness, terrain has no reason to be higher precisely
where the scrub is thicker. So if a filter's surface rises with understory and
another's does not, the first is accepting vegetation as ground -- an attribution
that needs no ground truth, only the flatness of the cell.

Datum. Each filter's surface is measured against the 1st percentile of all
returns in the same cell, which references no classification and no commercial
product. That datum is not unbiased: where understory is thick fewer pulses
reach soil, so the lowest return sits higher and the datum rises with the very
quantity on the x axis. But it rises identically for all three filters, so the
bias is common mode and the contrast between their slopes survives it. Only the
absolute level is affected, and no absolute claim is made here.
"""
import gzip, json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import ensure_blas, SRC
sys.path.insert(0, SRC)
ensure_blas()
import pdal
from replicalm.classify import las_reader

D = os.environ.get(
    "REPLICALM_NR_ROOT",
    "C:/NR_Lidar/Original_LAS_Files_15March2023/Original_LAS_Files_15March2023")
# Results live in the repository now; McLellan gave permission to publish
# measurements derived from the block on 2026-09-24. The environment
# variable still overrides, for working outside a checkout.
RESULTS = os.environ.get(
    "REPLICALM_NR_RESULTS",
    os.path.join(os.path.dirname(os.path.abspath(__file__)),
                 "results", "nr_block"))
HALF = 10.0

IN = os.path.join(RESULTS, sys.argv[1] if len(sys.argv) > 1
                  else "band_regression_tuned.json")
OUT = IN.replace("band_regression", "band_offset")
R = json.load(open(IN, encoding="utf-8"))
# whichever filters this run carried, so the tuned and tight PMF can be
# regressed side by side on the same cells
FILT = [f for f in ("smrf", "pmf", "pmf_tight", "csf")
        if any("mean_z" in r["filters"].get(f, {}) for r in R)]
ok = [r for r in R if all("mean_z" in r["filters"].get(f, {}) for f in FILT)]
print("reading %s; filters %s" % (IN, ", ".join(FILT)))
by_tile = {}
for r in ok:
    by_tile.setdefault(r["tile"], []).append(r)
print("%d cells across %d tiles" % (len(ok), len(by_tile)), flush=True)

rows = []
for n, (tile, cells) in enumerate(sorted(by_tile.items()), 1):
    pl = pdal.Pipeline(json.dumps({"pipeline": [las_reader(os.path.join(D, tile + ".las"))]}))
    pl.execute()
    a = pl.arrays[0]
    x, y, z = a["X"], a["Y"], a["Z"]
    for c in cells:
        m = (np.abs(x - c["x"]) <= HALF) & (np.abs(y - c["y"]) <= HALF)
        if m.sum() < 500:
            continue
        zc = z[m]
        rec = {"x": c["x"], "y": c["y"], "band": c["band"],
               "slope_pct": c["slope_pct"], "canopy_p90_m": c["canopy_p90_m"],
               "p1": float(np.percentile(zc, 1)),
               "p5": float(np.percentile(zc, 5)),
               "returns": int(m.sum())}
        for f in FILT:
            rec["off_" + f] = float(c["filters"][f]["median_z"]) - rec["p1"]
            rec["den_" + f] = float(c["filters"][f]["density"])
        rows.append(rec)
    print("  tile %2d/%d  %d cells" % (n, len(by_tile), len(rows)), flush=True)

json.dump(rows, open(OUT, "w", encoding="utf-8"), indent=1)

b = np.array([r["band"] for r in rows])
print("\nsurface height above the 1st percentile of all returns, vs understory")
print("  filter  mean offset   corr with band   slope per 10% band   n")
fits = {}
for f in FILT:
    o = np.array([r["off_" + f] for r in rows])
    s = np.polyfit(b, o, 1)
    fits[f] = (0.1 * s[0], float(np.corrcoef(b, o)[0, 1]))
    print("  %-9s %+8.3f m %13.2f %18.3f m %6d"
          % (f, o.mean(), fits[f][1], fits[f][0], len(o)))
print("\ncontrast in slope (common-mode datum bias cancels):")
for i, a_ in enumerate(FILT):
    for b_ in FILT[i + 1:]:
        print("  %-9s minus %-9s  %+.3f m per 10%% band"
              % (a_, b_, fits[a_][0] - fits[b_][0]))
print("\nwrote %s" % OUT)
