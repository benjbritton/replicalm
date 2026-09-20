r"""Choose the calibration tiles, recording the draw so it can be repeated.

Two Centro tiles drawn at random from the 87 that have both an original cloud
and a reprocessed DEM, plus South_GLAS_l0s395, which is named rather than drawn
because its structures are known and it anchors the set to ground already
examined by eye.
"""
import glob, json, os, random, re

SEED = 20260919          # the date, so the draw is reproducible and unchosen
random.seed(SEED)

DEM_ROOT = r"D:\_Archive_EdgeFixed"
LAS_ROOT = r"D:\GLiHT_LAS_orig"


def tile_key(name):
    m = re.search(r"(l\d+s\d+)", name)
    return m.group(1) if m else None


def index(root, region, patterns):
    out = {}
    for pat in patterns:
        for p in glob.glob(os.path.join(root, "**", pat), recursive=True):
            if region and region not in p:
                continue
            k = tile_key(os.path.basename(p))
            if k and k not in out:
                out[k] = p
    return out


centro_dem = index(DEM_ROOT, "Yuc_Centro", ["*DEM_0p5m*.tif"])
centro_las = index(LAS_ROOT, "Centro", ["*.las", "*.laz"])
pairs = sorted(set(centro_dem) & set(centro_las))
print("Centro tiles with both a cloud and a reference DEM: %d" % len(pairs))

picked = random.sample(pairs, 2)
print("drawn with seed %d: %s" % (SEED, ", ".join(picked)))

sel = []
for k in picked:
    sel.append({"region": "Yuc_Centro", "tile": k,
                "las": centro_las[k], "reference_dem": centro_dem[k],
                "chosen": "random draw, seed %d" % SEED})

south_dem = index(DEM_ROOT, "Yuc_South", ["*l0s395*DEM_0p5m*.tif"])
south_las = index(LAS_ROOT, "South", ["*l0s395*.las", "*l0s395*.laz"])
if "l0s395" in south_dem and "l0s395" in south_las:
    sel.append({"region": "Yuc_South", "tile": "l0s395",
                "las": south_las["l0s395"],
                "reference_dem": south_dem["l0s395"],
                "chosen": "named: known structures, used for visual review"})
else:
    print("  South_GLAS_l0s395: dem=%s las=%s"
          % ("l0s395" in south_dem, "l0s395" in south_las))

print()
for s in sel:
    mb = os.path.getsize(s["las"]) / 1e6
    dmb = os.path.getsize(s["reference_dem"]) / 1e6
    print("%-10s %-8s cloud %7.1f MB   reference DEM %6.1f MB"
          % (s["region"], s["tile"], mb, dmb))
    print("           %s" % s["las"])
    print("           %s" % s["reference_dem"])

os.makedirs(r"C:\Replicalm\tests", exist_ok=True)
with open(r"C:\Replicalm\tests\calibration_tiles.json", "w", encoding="utf-8") as fh:
    json.dump({"seed": SEED, "candidates": len(pairs), "tiles": sel}, fh, indent=1)
print("\nrecorded in C:\Replicalm\tests\calibration_tiles.json")
