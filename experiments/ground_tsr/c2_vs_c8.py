"""What separates class 2 from class 8 in the delivered ground population?

A classifier cannot be said to replicate a two-class partition until the
partition's basis is known. Class 8 is not a thinned model-keypoint set -- it
is denser than class 2 and carries no KeyPoint flag -- so the usual reading
does not apply and the distinction has to be measured.

Every candidate basis is checked against the dimensions the file already
carries: which flightline pass a point came from, where it sat in the return
train, how bright it was, how far off nadir, and when it was collected. The
spatial test asks whether the two classes interleave within a cell or occupy
separate ground, which separates a per-point decision from a per-swath one.
"""
import json
import sys

import numpy as np

sys.path.insert(0, r"C:\Replicalm\src")
from replicalm import classify

GND = (r"D:\_LAZ\Yuc_Centro\Centro_Glas\Centro_Glas"
       r"\Centro_GLAS_l0s441_GND_v1.laz")

p = classify._pdal().Pipeline(json.dumps({"pipeline": [
    {"type": "readers.las", "filename": GND, "nosrs": True}]}))
p.execute()
a = p.arrays[0]
cl = a["Classification"].astype(np.int32)
m2, m8 = cl == 2, cl == 8
print("class 2 = %d, class 8 = %d, total %d"
      % (int(m2.sum()), int(m8.sum()), len(cl)), flush=True)

print("")
print("CATEGORICAL DIMENSIONS, share of each class")
for f in ("PointSourceId", "ReturnNumber", "NumberOfReturns", "ScanChannel",
          "ScanDirectionFlag", "EdgeOfFlightLine", "UserData"):
    if f not in a.dtype.names:
        print("  %s absent" % f)
        continue
    v = a[f].astype(np.int64)
    vals = np.unique(v)
    print("  %s, %d distinct value(s)" % (f, len(vals)))
    if len(vals) > 12:
        print("    too many to tabulate; range %d..%d" % (vals.min(),
                                                          vals.max()))
        continue
    print("    %-8s %12s %12s" % ("value", "class 2", "class 8"))
    for u in vals:
        print("    %-8d %11.2f%% %11.2f%%"
              % (u, 100 * float((v[m2] == u).mean()),
                 100 * float((v[m8] == u).mean())))

print("")
print("CONTINUOUS DIMENSIONS, percentiles per class")
print("%-16s %-8s %9s %9s %9s %9s %9s"
      % ("dimension", "class", "p5", "p25", "median", "p75", "p95"))
for f in ("Z", "Intensity", "ScanAngleRank", "GpsTime"):
    if f not in a.dtype.names:
        continue
    v = a[f].astype("f8")
    for lbl, mm in (("2", m2), ("8", m8)):
        q = np.percentile(v[mm], [5, 25, 50, 75, 95])
        print("%-16s %-8s %9.3f %9.3f %9.3f %9.3f %9.3f"
              % (f, lbl, q[0], q[1], q[2], q[3], q[4]))

print("")
print("RETURN TRAIN POSITION, is a point the last of its pulse?")
if "ReturnNumber" in a.dtype.names:
    rn = a["ReturnNumber"].astype(np.int64)
    nr = a["NumberOfReturns"].astype(np.int64)
    for lbl, mm in (("class 2", m2), ("class 8", m8)):
        print("  %-8s only return %.2f%%, last of several %.2f%%, "
              "intermediate %.2f%%, first of several %.2f%%"
              % (lbl, 100 * float(((nr == 1)[mm]).mean()),
                 100 * float(((nr > 1) & (rn == nr))[mm].mean()),
                 100 * float(((nr > 1) & (rn > 1) & (rn < nr))[mm].mean()),
                 100 * float(((nr > 1) & (rn == 1))[mm].mean())))

print("")
print("SPATIAL ARRANGEMENT")
x = a["X"].astype("f8")
y = a["Y"].astype("f8")
for cell in (0.5, 5.0, 20.0):
    ci = np.floor((x - x.min()) / cell).astype(np.int64)
    cj = np.floor((y - y.min()) / cell).astype(np.int64)
    key = ci * 1000003 + cj
    o = np.argsort(key, kind="stable")
    ks = key[o]
    first = np.ones(len(ks), bool)
    first[1:] = ks[1:] != ks[:-1]
    gid = np.cumsum(first) - 1
    n = gid[-1] + 1
    c2 = np.bincount(gid, weights=m2[o].astype("f8"), minlength=n)
    c8 = np.bincount(gid, weights=m8[o].astype("f8"), minlength=n)
    both = (c2 > 0) & (c8 > 0)
    only2 = (c2 > 0) & (c8 == 0)
    only8 = (c2 == 0) & (c8 > 0)
    print("  %5.1f m cells: %d occupied, both classes %.2f%%, "
          "class 2 only %.2f%%, class 8 only %.2f%%"
          % (cell, n, 100 * both.mean(), 100 * only2.mean(),
             100 * only8.mean()))

print("")
print("GPS TIME STRUCTURE, do the classes come from different passes?")
if "GpsTime" in a.dtype.names:
    t = a["GpsTime"].astype("f8")
    tt = t - t.min()
    edges = np.linspace(0, tt.max(), 41)
    h2, _ = np.histogram(tt[m2], bins=edges)
    h8, _ = np.histogram(tt[m8], bins=edges)
    print("  40 equal time bins over %.1f s; share of class 8 in each:"
          % tt.max())
    frac = np.where(h2 + h8 > 0, h8 / np.maximum(h2 + h8, 1), np.nan)
    print("  " + " ".join("%.2f" % f if np.isfinite(f) else " -- "
                          for f in frac))
    print("  class 8 share: min %.2f, median %.2f, max %.2f over bins"
          % (np.nanmin(frac), np.nanmedian(frac), np.nanmax(frac)))
