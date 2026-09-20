r"""Draw five South_NFI tiles for the pilot sweep, excluding those already used.

The calibration set was two Centro tiles drawn with seed 20260919 plus a named
South tile. Those three tuned the pipeline, so they cannot also test it. This
draw uses the same seed and the same pool with the used tiles removed, which
keeps the procedure reproducible without redrawing the same names.
"""
import glob, json, os, random, re

SEED = 20260919
USED = {"l8s431", "l0s444", "l0s395"}
DEM_ROOT = r"D:\_Archive_EdgeFixed"
LAS_ROOT = r"D:\GLiHT_LAS_orig"
OUT = r"C:\Replicalm\tests\pilot_tiles.json"


def tile_key(name):
    m = re.search(r"(l\d+s\d+)", name)
    return m.group(1) if m else None


def index(root, region, patterns):
    out = {}
    for pat in patterns:
        for p in glob.glob(os.path.join(root, "**", pat), recursive=True):
            if region not in p:
                continue
            k = tile_key(os.path.basename(p))
            if k and k not in out:
                out[k] = p
    return out


dem = index(DEM_ROOT, "South_NFI", ["*DEM_0p5m*.tif"])
las = index(LAS_ROOT, "South_NFI", ["*.las", "*.laz"])
pool = sorted((set(dem) & set(las)) - USED)
print("South_NFI tiles with both a cloud and a reference DEM: %d"
      % len(set(dem) & set(las)))
print("after excluding %s: %d" % (", ".join(sorted(USED)), len(pool)))

random.seed(SEED)
picked = random.sample(pool, 5)
print("drawn with seed %d: %s\n" % (SEED, ", ".join(picked)))

sel = [{"region": "Yuc_South_NFI", "tile": k, "las": las[k],
        "reference_dem": dem[k],
        "chosen": "random draw, seed %d, pool excludes %s"
                  % (SEED, "/".join(sorted(USED)))} for k in picked]
for s in sel:
    mb = os.path.getsize(s["las"]) / 1e6
    print("  %-8s %7.0f MB  %s" % (s["tile"], mb, os.path.basename(s["las"])))

with open(OUT, "w", encoding="utf-8") as fh:
    json.dump({"seed": SEED, "pool": len(pool), "excluded": sorted(USED),
               "tiles": sel}, fh, indent=1)
print("\nwrote %s" % OUT)
