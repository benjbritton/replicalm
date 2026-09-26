r"""Two 50 m elevation profiles through a structure, north-south and west-east.

Samples the finalised DEM, not the G1 composite: the composite is a rendering
and its values are brightness, not height. Sampling is bilinear at a quarter
of the DEM cell, so the profile is not stepped by the grid.

Stacked north-south over west-east, as one image. Both panels share an
elevation axis so the two sections can be read against each other.
"""
import argparse
import os

import numpy as np


def sample_line(z, tf, nodata, x0, y0, x1, y1, n):
    """Bilinear elevation along a straight line, NaN where there is no data."""
    cell_x, cell_y = tf.a, -tf.e
    xs = np.linspace(x0, x1, n)
    ys = np.linspace(y0, y1, n)
    c = (xs - tf.c) / cell_x - 0.5
    r = (tf.f - ys) / cell_y - 0.5
    c0, r0 = np.floor(c).astype(int), np.floor(r).astype(int)
    fc, fr = c - c0, r - r0
    h, w = z.shape
    out = np.full(n, np.nan)
    good = (r0 >= 0) & (r0 < h - 1) & (c0 >= 0) & (c0 < w - 1)
    if good.any():
        i, j = r0[good], c0[good]
        a, b = fr[good], fc[good]
        q = np.stack([z[i, j], z[i, j + 1], z[i + 1, j], z[i + 1, j + 1]])
        if nodata is not None:
            q[q == nodata] = np.nan
        q[q == 0] = np.nan
        out[good] = ((1 - a) * ((1 - b) * q[0] + b * q[1])
                     + a * ((1 - b) * q[2] + b * q[3]))
    d = np.hypot(xs - x0, ys - y0) - np.hypot(x1 - x0, y1 - y0) / 2.0
    return d, out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("dem")
    ap.add_argument("out_tif")
    ap.add_argument("--centre", nargs=2, type=float, required=True,
                    metavar=("EASTING", "NORTHING"))
    ap.add_argument("--length", type=float, default=50.0)
    ap.add_argument("--label", default="")
    a = ap.parse_args(argv)

    import rasterio
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cx, cy = a.centre
    half = a.length / 2.0
    with rasterio.open(a.dem) as r:
        tf, nd = r.transform, r.nodata
        cell = abs(tf.a)
        pad = int(half / cell) + 4
        col = int((cx - tf.c) / cell)
        row = int((tf.f - cy) / cell)
        win = rasterio.windows.Window(col - pad, row - pad, 2 * pad, 2 * pad)
        z = r.read(1, window=win, boundless=True,
                   fill_value=float("nan")).astype(np.float64)
        wtf = r.window_transform(win)

    n = int(a.length / (cell / 4)) + 1
    # north to south: start at the north end, so distance runs N -> S
    dn, zn = sample_line(z, wtf, nd, cx, cy + half, cx, cy - half, n)
    # west to east
    dw, zw = sample_line(z, wtf, nd, cx - half, cy, cx + half, cy, n)
    for tag, v in (("N-S", zn), ("W-E", zw)):
        f = v[np.isfinite(v)]
        print("%s  %d/%d samples, %.2f to %.2f m, relief %.2f m"
              % (tag, f.size, v.size, f.min(), f.max(), f.max() - f.min()))

    allv = np.concatenate([zn, zw])
    allv = allv[np.isfinite(allv)]
    lo, hi = allv.min(), allv.max()
    m = max(0.25, 0.08 * (hi - lo))
    lo, hi = lo - m, hi + m

    fig, axes = plt.subplots(2, 1, figsize=(9.0, 7.2), sharex=True, sharey=True,
                             dpi=300)
    title = a.label or os.path.basename(a.dem)
    fig.suptitle("%s\nprofiles through %.2f E, %.2f N  (UTM 16N, %.0f m)"
                 % (title, cx, cy, a.length), fontsize=10)
    for ax, (d, v, name, left, right) in zip(axes, [
            (dn, zn, "north - south", "N", "S"),
            (dw, zw, "west - east", "W", "E")]):
        ax.plot(d, v, "-", color="#1a1a1a", lw=1.1)
        ax.fill_between(d, lo, v, where=np.isfinite(v), color="#c8c8c8",
                        linewidth=0)
        ax.axvline(0.0, color="#b03030", lw=0.8, ls="--")
        ax.set_ylim(lo, hi)
        ax.set_xlim(d.min(), d.max())
        ax.grid(True, lw=0.3, color="#bbbbbb")
        ax.set_ylabel("elevation (m)", fontsize=9)
        ax.text(0.01, 0.92, name, transform=ax.transAxes, fontsize=9,
                fontweight="bold")
        ax.text(0.005, 0.03, left, transform=ax.transAxes, fontsize=9,
                color="#666666")
        ax.text(0.99, 0.03, right, transform=ax.transAxes, fontsize=9,
                color="#666666", ha="right")
        f = v[np.isfinite(v)]
        ax.text(0.99, 0.92, "relief %.2f m" % (f.max() - f.min()),
                transform=ax.transAxes, fontsize=9, ha="right")
    axes[1].set_xlabel("distance from centre (m)", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.95))

    # Written through GDAL, not through Pillow. Both matplotlib's TIFF backend
    # and PIL.Image.save with a .tif target abort the g1 interpreter outright --
    # exit 127, no traceback, stdout lost with the buffer. Rendering to an array
    # and handing it to rasterio avoids Pillow entirely.
    fig.canvas.draw()
    rgba = np.asarray(fig.canvas.buffer_rgba())
    rgb = rgba[:, :, :3].transpose(2, 0, 1).copy()
    with rasterio.open(a.out_tif, "w", driver="GTiff", height=rgb.shape[1],
                       width=rgb.shape[2], count=3, dtype="uint8",
                       compress="lzw") as dst:
        dst.write(rgb)
    plt.close(fig)
    print("wrote %s  (%d x %d)" % (a.out_tif, rgb.shape[2], rgb.shape[1]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
