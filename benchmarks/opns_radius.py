r"""Positive openness at a wider search radius, over a whole tile.

Reviving the openness layer by fixing its stretch makes it contribute again,
and what it contributes at RVT's default reach is a hard outline around every
piece of micro-relief. That reach is `svf_r_max`, counted in PIXELS, so a tile
delivered at 0.35 m cells gets a 3.5 m radius where one at 0.5 m gets 5 m --
the same setting means a different distance on every tile.

This recomputes the layer at a wider radius from the finished DEM and writes it
beside the others, so the composite can be built from it without recomputing
anything else. Nothing existing is modified.
"""
import argparse
import os
import time

import numpy as np


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("dem")
    ap.add_argument("out_tif")
    ap.add_argument("--r-max", type=int, default=30, help="radius in PIXELS")
    ap.add_argument("--n-dir", type=int, default=16)
    a = ap.parse_args(argv)

    import rasterio
    import rvt.vis

    with rasterio.open(a.dem) as r:
        z = r.read(1).astype(np.float32)
        nd, profile = r.nodata, r.profile.copy()
        cell = abs(r.transform.a)
    if nd is not None:
        z[z == nd] = np.nan
    z[z == 0] = np.nan
    print("%s: %d x %d at %.3f m; r_max %d px = %.2f m"
          % (os.path.basename(a.dem), z.shape[0], z.shape[1], cell,
             a.r_max, a.r_max * cell), flush=True)

    t0 = time.time()
    out = rvt.vis.sky_view_factor(dem=z, resolution=cell, compute_svf=False,
                                  compute_asvf=False, compute_opns=True,
                                  svf_n_dir=a.n_dir, svf_r_max=a.r_max,
                                  svf_noise=0, no_data=np.nan)
    opns = out["opns"] if isinstance(out, dict) else out
    opns = np.asarray(opns, np.float32)
    v = opns[np.isfinite(opns)]
    print("%.0f s; %.2f to %.2f deg, p1 %.2f p50 %.2f p99 %.2f"
          % (time.time() - t0, v.min(), v.max(), np.percentile(v, 1),
             np.percentile(v, 50), np.percentile(v, 99)), flush=True)

    profile.update(dtype="float32", count=1, compress="lzw", nodata=np.nan)
    with rasterio.open(a.out_tif, "w", **profile) as dst:
        dst.write(opns, 1)
    print("wrote %s" % a.out_tif)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
