r"""Measure every tile that has both a cloud and a reference DEM, so the pilot
tiles can be chosen on traits instead of on assumption.

Nothing here reads a point. The LAS header carries the bounding box, the point
count and the count by return number, which between them give density and a
canopy proxy; the reference DEM gives relief and slope. Both are cheap enough to
run across the whole pool, which is the point -- a trait-based draw is only
meaningful if the traits were measured on everything available to draw from.

TRAITS AND WHAT THEY STAND IN FOR
---------------------------------
  density           all returns per square metre of the header bounding box.
                    The trait that produced l0s444's 89% fallback.
  returns_per_pulse total points divided by first returns. One means a bare
                    surface; higher means the pulse found several layers on the
                    way down, which is the closest thing to a forest measure
                    available without classifying.
  last_return_frac  the share of returns that terminated the pulse. Low values
                    on a multi-layer tile mean few pulses reached the floor,
                    which is the "low ground returns" condition.
  mean_slope        from the reference DEM, downsampled. Selection only: a
                    coarse grid understates local slope, but the ordering
                    between tiles survives.
  relief            98th minus 2nd percentile elevation, in metres.
"""
import glob, json, os, re, struct, sys
import numpy as np

DEM_ROOT = r"D:\_Archive_EdgeFixed"
LAS_ROOT = r"D:\GLiHT_LAS_orig"
OUT = r"C:\Replicalm\tests\pool_traits.json"
COARSE = 300          # cells per side for the slope read


def tile_key(name):
    m = re.search(r"(l\d+s\d+)", name)
    return m.group(1) if m else None


def las_header(path):
    """Point count, bounds and return histogram, straight from the header."""
    with open(path, "rb") as fh:
        h = fh.read(400)
    if h[:4] != b"LASF":
        return None
    ver = (h[24], h[25])
    legacy_n = struct.unpack_from("<I", h, 107)[0]
    legacy_ret = struct.unpack_from("<5I", h, 111)
    mx, nx, my, ny, mz, nz = struct.unpack_from("<6d", h, 179)
    n, ret = legacy_n, list(legacy_ret)
    if ver >= (1, 4):
        n14 = struct.unpack_from("<Q", h, 247)[0]
        if n14:
            n = n14
            ret = list(struct.unpack_from("<15Q", h, 255))
    if not n:
        return None
    area = (mx - nx) * (my - ny)
    if area <= 0:
        return None
    first = ret[0] if ret and ret[0] else n
    # A pulse's last return is the one that stopped it. Summed over return
    # numbers this is not directly available, so it is approximated by the
    # deepest populated return bin per pulse: total minus first gives the
    # subsequent returns, and the ratio below is what separates canopy from bare.
    return {"points": int(n), "area_m2": float(area),
            "density": float(n / area),
            "returns_per_pulse": float(n / first) if first else float("nan"),
            "multi_return_frac": float(1.0 - (ret[0] / n)) if ret[0] else 0.0,
            "relief_las_m": float(mz - nz), "version": "%d.%d" % ver}


def dem_traits(path):
    from osgeo import gdal
    gdal.UseExceptions()
    d = gdal.Open(path)
    b = d.GetRasterBand(1)
    gt = d.GetGeoTransform()
    w = min(COARSE, d.RasterXSize); hgt = min(COARSE, d.RasterYSize)
    a = b.ReadAsArray(buf_xsize=w, buf_ysize=hgt).astype("f8")
    nod = b.GetNoDataValue()
    m = np.isfinite(a)
    if nod is not None:
        m &= a != nod
    m &= a != 0.0
    if m.sum() < 100:
        return None
    # the downsampled cell size, so slope is rise over the run actually used
    cx = abs(gt[1]) * d.RasterXSize / w
    cy = abs(gt[5]) * d.RasterYSize / hgt
    z = np.where(m, a, np.nan)
    gy, gx = np.gradient(z, cy, cx)
    slope = np.degrees(np.arctan(np.hypot(gx, gy)))
    v = a[m]
    return {"nodata_frac": float(1.0 - m.mean()),
            "relief_m": float(np.nanpercentile(v, 98) - np.nanpercentile(v, 2)),
            "mean_slope_deg": float(np.nanmean(slope)),
            "p95_slope_deg": float(np.nanpercentile(slope[np.isfinite(slope)], 95)),
            "cell_m": float(abs(gt[1])),
            "size": [d.RasterXSize, d.RasterYSize]}


dems, lass = {}, {}
for p in glob.glob(os.path.join(DEM_ROOT, "**", "*DEM_0p5m*.tif"), recursive=True):
    k = tile_key(os.path.basename(p))
    if k:
        dems.setdefault((p.split(os.sep)[2], k), p)
for p in glob.glob(os.path.join(LAS_ROOT, "**", "*.la[sz]"), recursive=True):
    k = tile_key(os.path.basename(p))
    if k:
        lass.setdefault((p.split(os.sep)[2], k), p)

pairs = sorted(set(dems) & set(lass))
print("tiles with both a cloud and a reference DEM: %d" % len(pairs))

rows = []
for i, (region, k) in enumerate(pairs):
    try:
        lh = las_header(lass[(region, k)])
        if lh is None:
            continue
        dt = dem_traits(dems[(region, k)])
        if dt is None:
            continue
        r = {"region": region, "tile": k, "las": lass[(region, k)],
             "reference_dem": dems[(region, k)]}
        r.update(lh); r.update(dt)
        rows.append(r)
    except Exception as e:
        print("  %s %s: %s" % (region, k, str(e)[:60]))
    if (i + 1) % 50 == 0:
        print("  %d/%d" % (i + 1, len(pairs)))

with open(OUT, "w", encoding="utf-8") as fh:
    json.dump(rows, fh, indent=1)
print("\nmeasured %d tiles -> %s" % (len(rows), OUT))

d = np.array([r["density"] for r in rows])
s = np.array([r["mean_slope_deg"] for r in rows])
m = np.array([r["multi_return_frac"] for r in rows])
for name, arr, unit in (("density", d, "pts/m2"), ("mean slope", s, "deg"),
                        ("multi-return frac", m, "")):
    print("%-18s min %7.2f  p25 %7.2f  median %7.2f  p75 %7.2f  max %7.2f %s"
          % (name, arr.min(), np.percentile(arr, 25), np.median(arr),
             np.percentile(arr, 75), arr.max(), unit))
