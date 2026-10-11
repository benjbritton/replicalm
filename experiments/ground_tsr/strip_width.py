"""How wide is the strip the reference actually removed, at each resolution?

mask_dem_edges.py expresses it as 64 pixels, which is 21.1 m at 0.33 m and
would be 32 m at 0.5 m. clean_las_edges.py states 21.0 m outright. If the
archive's 0.5 m products were stripped in pixels the two resolutions will
disagree in metres; if they were stripped in metres they will agree.
"""
import numpy as np
import rasterio

PAIRS = [
    ("0p33m", r"D:\_DEMs\Yuc_Centro\Centro_Glas\Centro_Glas"
              r"\Centro_GLAS_l0s454_DEM_0p33m_v1.tif",
     r"D:\_Archive_EdgeFixed\Yuc_Centro\Centro_Glas\Centro_Glas"
     r"\Centro_GLAS_l0s454_DEM_0p33m_v1.tif"),
    ("0p5m", r"D:\_DEMs\Yuc_Centro\Centro_Glas\Centro_Glas"
             r"\Centro_GLAS_l0s454_DEM_0p5m_v1.tif",
     r"D:\_Archive_EdgeFixed\Yuc_Centro\Centro_Glas\Centro_Glas"
     r"\Centro_GLAS_l0s454_DEM_0p5m_v1.tif"),
    ("1m", r"D:\_DEMs\Yuc_Centro\Centro_Glas\Centro_Glas"
           r"\Centro_GLAS_l0s454_DEM_1m_v1.tif",
     r"D:\_Archive_EdgeFixed\Yuc_Centro\Centro_Glas\Centro_Glas"
     r"\Centro_GLAS_l0s454_DEM_1m_v1.tif"),
]


def mask_of(p):
    with rasterio.open(p) as r:
        z = r.read(1).astype("f8")
        nd, cc = r.nodata, abs(r.transform.a)
    if nd is not None:
        z[z == nd] = np.nan
    z[z > 1e30] = np.nan
    z[z == 0] = np.nan
    return np.isfinite(z), cc


import os
print("%-7s %7s %12s %12s %10s %12s %11s"
      % ("res", "cell", "prod span", "arch span", "cells cut",
         "per side", "metres/side"))
print("-" * 80)
for nm, prod, arch in PAIRS:
    if not (os.path.exists(prod) and os.path.exists(arch)):
        print("%-7s missing" % nm)
        continue
    mp, cc = mask_of(prod)
    ma, _ = mask_of(arch)
    sp, sa = [], []
    rows = np.flatnonzero(mp.any(axis=1))
    for r in rows[::20]:
        c = np.flatnonzero(mp[r])
        sp.append(c[-1] - c[0] + 1)
        ca = np.flatnonzero(ma[r]) if ma[r].any() else None
        if ca is not None and len(ca):
            sa.append(ca[-1] - ca[0] + 1)
    msp, msa = float(np.median(sp)), float(np.median(sa))
    cut = msp - msa
    print("%-7s %7.4f %12.0f %12.0f %10.0f %12.1f %11.2f"
          % (nm, cc, msp, msa, cut, cut / 2.0, cc * cut / 2.0))
print("")
print("64 px would be: %.2f m at 0.33 m, %.2f m at 0.5 m, %.2f m at 1 m"
      % (64 * 0.329981, 64 * 0.500024, 64 * 1.0))
