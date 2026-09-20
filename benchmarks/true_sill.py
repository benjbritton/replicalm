r"""Rescale the nugget against a sill measured where the curve actually flattens.

The fine-lag run gives a trustworthy nugget, because its first bins sit near the
origin where the extrapolation belongs. Its sill is the mean of bins from 10 to
15 m, which is only a sill if the curve has plateaued by then -- and this terrain
has no sill within 80 m on the raw elevations.

So: take the nugget from the fine run, measure the far-field behaviour from the
coarse run, and report the ratio against a sill that has been checked rather than
assumed. If semivariance is still climbing at 30 m the ratio is not defined at
all, and the honest statement is the nugget as a fraction of variance at a named
lag.
"""
import json, os, sys
sys.path.insert(0, r"C:\Replicalm\src")
import numpy as np
from osgeo import gdal
from replicalm import residual
gdal.UseExceptions()

with open(r"C:\Replicalm\tests\clips\index.json", encoding="utf-8") as fh:
    clips = {c["tile"]: c for c in json.load(fh)}
S = r"C:\Replicalm\tests\class_sweep"
BEST = {"l0s395": "csfrigid2", "l8s431": "smrfs05t05", "l0s444": "smrfs05t025"}


def read(p):
    d = gdal.Open(str(p)); b = d.GetRasterBand(1)
    a = b.ReadAsArray().astype("f8"); nod = b.GetNoDataValue()
    m = np.isfinite(a) & (a != 0.0)
    if nod is not None:
        m &= a != nod
    return a, m, abs(d.GetGeoTransform()[1])


print("%-8s %8s %8s %9s %9s %9s %9s  %s"
      % ("tile", "nugget", "g(15)", "g(30)", "g(60)", "slope", "ratio*",
         "plateau by 30 m?"))
print("-" * 92)
for tile, tag in BEST.items():
    p = os.path.join(S, "%s_%s.tif" % (tile, tag))
    ca, cm, cell = read(p)
    ra, rm, _ = read(clips[tile]["clipped_reference"])
    both = cm & rm
    resid = np.where(both, ca - ra, 0.0)

    # fine: an honest nugget
    lf, gf, _ = residual.gridded_variogram(resid, both, cell=cell,
                                           max_lag_cells=int(15 / cell),
                                           n_bins=24)
    nug = residual.characterise(lf, gf, None)["nugget"]

    # coarse: the far field
    lc, gc, _ = residual.gridded_variogram(resid, both, cell=cell,
                                           max_lag_cells=int(90 / cell),
                                           n_bins=30)

    def at(target):
        i = int(np.argmin(np.abs(lc - target)))
        return float(gc[i])

    g15, g30, g60 = at(15.0), at(30.0), at(60.0)
    seg = (lc >= 15.0) & (lc <= 30.0)
    slope = float(np.polyfit(lc[seg], gc[seg], 1)[0]) if seg.sum() > 2 else float("nan")
    # flat enough if the 15-30 m rise is under 10% of the level there
    flat = abs(slope * 15.0) < 0.10 * g30
    sill = g30 if flat else g60
    ratio = nug / sill if sill > 0 else float("nan")
    print("%-8s %8.5f %8.5f %9.5f %9.5f %9.2e %9.2f  %s"
          % (tile, nug, g15, g30, g60, slope, ratio,
             "yes" if flat else "NO, still climbing"))

print()
print("ratio* is the fine-lag nugget over a sill taken at 30 m where the curve")
print("has flattened, or at 60 m where it has not. Where it has not flattened")
print("even at 60 m the figure remains a lower bound on the true ratio.")
