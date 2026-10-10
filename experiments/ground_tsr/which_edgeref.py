"""Which edge treatment is the refined one, and is it the reference I used?

Dates from the G1 reprocessing work:

  2026-03-28  GLiHTmast.py
  2026-03-31  D:\\_DEMs                     production DEMs
  2026-04-16  D:\\_Archive_EdgeFixed        the reference I have been using
  2026-06-13 06:44  clean_las_edges.py v3.0 trim the point cloud, banded
  2026-06-13 07:08  _south_v2_surfer        tested through Surfer
  2026-06-13 07:14  mask_dem_edges.py       trim the raster, per row
  2026-06-13 08:40  D:\\_DEMs_EdgeMasked     run over all 457 tiles

So the raster per-row method was written half an hour after the point-cloud
method had been tested, and it is the one that went to production. And
_Archive_EdgeFixed predates both by two months, which means the surface I
have been scoring every edge against is the earlier treatment, not the
refined one.

This compares the two products on the same tile: whether the interiors are
the same surface, and how the two boundaries differ.
"""
import numpy as np
import rasterio
from scipy import ndimage

A_PATH = (r"D:\_Archive_EdgeFixed\Yuc_Centro\Centro_Glas\Centro_Glas"
          r"\Centro_GLAS_l0s454_DEM_0p33m_v1.tif")
B_PATH = (r"D:\_DEMs_EdgeMasked\Yuc_Centro\Centro_Glas\Centro_Glas"
          r"\Centro_GLAS_l0s454_DEM_0p33m_v1.tif")
C_PATH = (r"D:\_DEMs\Yuc_Centro\Centro_Glas\Centro_Glas"
          r"\Centro_GLAS_l0s454_DEM_0p33m_v1.tif")


def rd(p):
    with rasterio.open(p) as r:
        z = r.read(1).astype("f8")
        nd, cc, tr = r.nodata, abs(r.transform.a), r.transform
        h, w = r.height, r.width
    if nd is not None:
        z[z == nd] = np.nan
    z[z > 1e30] = np.nan
    z[z == 0] = np.nan
    return z, cc, (h, w)


import os
loaded = {}
for nm, p in (("archive Apr", A_PATH), ("edgemasked Jun", B_PATH),
              ("production", C_PATH)):
    if not os.path.exists(p):
        print("%-16s MISSING  %s" % (nm, p))
        continue
    z, cc, sh = rd(p)
    loaded[nm] = z
    print("%-16s %s at %.6f m, %d valid cells"
          % (nm, sh, cc, int(np.isfinite(z).sum())))

if "archive Apr" in loaded and "edgemasked Jun" in loaded:
    A, B = loaded["archive Apr"], loaded["edgemasked Jun"]
    ma, mb = np.isfinite(A), np.isfinite(B)
    both = ma & mb
    print("")
    print("COVERAGE")
    print("  April only      %9d cells" % int((ma & ~mb).sum()))
    print("  June only       %9d cells" % int((mb & ~ma).sum()))
    print("  both            %9d cells" % int(both.sum()))
    print("  June vs April   %+8.2f%%"
          % (100.0 * (mb.sum() / max(ma.sum(), 1) - 1)))
    d = (B - A)[both]
    d = d[np.isfinite(d)]
    print("")
    print("ELEVATION WHERE BOTH HAVE DATA")
    print("  identical cells %9d of %d (%.4f%%)"
          % (int((d == 0).sum()), len(d), 100.0 * (d == 0).mean()))
    print("  max |difference| %.6f m" % float(np.abs(d).max()))

    def shape(m):
        f = ndimage.binary_fill_holes(m)
        per = int((f ^ ndimage.binary_erosion(f)).sum())
        ar = int(f.sum())
        _l, k = ndimage.label(m)
        return ar, per, per / max(2.0 * np.sqrt(np.pi * ar), 1e-9), int(k)

    print("")
    print("%-16s %10s %10s %11s %7s" % ("mask", "area", "perimeter",
                                        "raggedness", "pieces"))
    for nm in ("archive Apr", "edgemasked Jun", "production"):
        if nm in loaded:
            ar, per, rg, k = shape(np.isfinite(loaded[nm]))
            print("%-16s %10d %10d %11.2f %7d" % (nm, ar, per, rg, k))

    # per-row span width, which is what mask_dem_edges.py controls
    print("")
    print("PER-ROW SPAN, first to last valid cell")
    for nm in ("archive Apr", "edgemasked Jun", "production"):
        if nm not in loaded:
            continue
        m = np.isfinite(loaded[nm])
        rows = np.flatnonzero(m.any(axis=1))
        widths, runs = [], []
        for r in rows[::25]:
            c = np.flatnonzero(m[r])
            widths.append(c[-1] - c[0] + 1)
            dd = np.diff(np.concatenate([[0], m[r].astype(np.int8), [0]]))
            runs.append(int((dd == 1).sum()))
        print("  %-16s median span %6.0f cells, median runs per row %.2f"
              % (nm, float(np.median(widths)), float(np.median(runs))))
