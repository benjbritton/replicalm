"""Pair a G-LiHT tile with its archive reference, unambiguously.

THE BUG THIS REPLACES
---------------------
`variant_g1.locate()` found a tile's reference by globbing `*<id>*0p5m*.tif`,
where `<id>` is the `l0s441` part of the name. Tile ids repeat across regions:
there is a Centro GLAS l0s441 and a Chiaps-1 GLAS l0s441, a South GLAS l0s395
and a North-2 GLAS l0s395. On the 458-tile corpus that glob is ambiguous for
223 tiles -- nearly half -- and returns whichever the filesystem lists first.
Every per-tile measurement in this work used explicit paths and so is
unaffected, but a survey-scale comparison driven by that function would have
scored half the corpus against the wrong reference and produced a plausible
distribution of nonsense.

Nothing needs renaming. Both names already carry the region:

    AMIGACarb_Yuc_Centro_GLAS_Apr2013_l0s441.las
    Centro_GLAS_l0s441_DEM_0p5m_v1.tif

so the information is present and the glob simply discarded it. What follows
builds a key from region, sensor and tile id and requires exactly one match,
raising on none and on more than one rather than choosing.

THE TWO NAMING SCHEMES
----------------------
Delivered clouds are `AMIGACarb_<region>_<sensor>_<month><year>_<id>`, where the
region carries underscores of its own and the month differs between campaigns.
Archive rasters are `<prefix>_<sensor>_<id>_DEM_<res>_v1.tif` with a shorter
region prefix, the seven in REGION_PREFIX below. The mapping between them is
flat and is written out rather than inferred, because inferring it is what the
glob was doing.
"""
import os
import re

# delivered-cloud region -> archive filename prefix
REGION_PREFIX = {
    "Yuc_South": "South",
    "Yuc_Centro": "Centro",
    "Yuc_Norte": "Norte",
    "Yuc_North_2": "North2",
    "Out_of_the_Yuc": "OutYuc",
    "Chiaps_1": "Chiaps1",
    "Chiap_Campeche": "ChiapsCampeche",
}

SENSORS = ("GLAS", "NFI")

_TILE = re.compile(
    r"^AMIGACarb_(?P<region>.+?)_(?P<sensor>GLAS|NFI)_"
    r"(?P<campaign>[A-Z][a-z]{2}\d{4})_(?P<id>l\d+s\d+)$")

_ARCHIVE = re.compile(
    r"^(?P<prefix>[A-Za-z0-9]+)_(?P<sensor>GLAS|NFI)_(?P<id>l\d+s\d+)_"
    r"(?P<product>[A-Za-z0-9]+)_(?P<res>[0-9a-zA-Z]+)_v\d+$")


class PairingError(RuntimeError):
    pass


def tile_key(name):
    """`AMIGACarb_Yuc_Centro_GLAS_Apr2013_l0s441` -> (Centro, GLAS, l0s441)."""
    stem = os.path.splitext(os.path.basename(str(name)))[0]
    # a cloud may arrive as .las.gz, so strip twice
    stem = os.path.splitext(stem)[0] if stem.endswith(".las") else stem
    m = _TILE.match(stem)
    if not m:
        raise PairingError("not a G-LiHT tile name: %r" % stem)
    region = m.group("region")
    if region not in REGION_PREFIX:
        raise PairingError(
            "unknown region %r in %r. Add it to REGION_PREFIX; guessing is "
            "what the glob this replaces was doing." % (region, stem))
    return REGION_PREFIX[region], m.group("sensor"), m.group("id")


def archive_key(path):
    """`Centro_GLAS_l0s441_DEM_0p5m_v1.tif` -> (Centro, GLAS, l0s441, DEM, 0p5m)."""
    stem = os.path.splitext(os.path.basename(str(path)))[0]
    m = _ARCHIVE.match(stem)
    if not m:
        return None
    return (m.group("prefix"), m.group("sensor"), m.group("id"),
            m.group("product"), m.group("res"))


def index_archive(root, product="DEM", res="0p5m"):
    """Every archive raster of one product and resolution, keyed by tile.

    Built once and reused: walking `_Archive_EdgeFixed` costs seconds and a
    survey-scale comparison would otherwise do it 458 times.
    """
    index = {}
    for dirpath, _dirs, files in os.walk(root):
        for f in files:
            if not f.lower().endswith(".tif"):
                continue
            k = archive_key(f)
            if not k or k[3] != product or k[4] != res:
                continue
            index.setdefault(k[:3], []).append(os.path.join(dirpath, f))
    return index


def reference_for(tile, index=None, root=None, product="DEM", res="0p5m"):
    """The one archive raster for this tile. Raises rather than choosing.

    Pass `index` from `index_archive` when resolving many tiles; pass `root` to
    build one for a single lookup.
    """
    if index is None:
        if root is None:
            raise PairingError("give either an index or an archive root")
        index = index_archive(root, product=product, res=res)
    key = tile_key(tile)
    hits = index.get(key, [])
    if not hits:
        raise PairingError(
            "no %s %s reference for %s (looked for %s_%s_%s)"
            % (res, product, tile, key[0], key[1], key[2]))
    if len(hits) > 1:
        raise PairingError(
            "%d references match %s, which should be impossible since region, "
            "sensor and tile id are all in the key:\n  %s"
            % (len(hits), tile, "\n  ".join(hits)))
    return hits[0]


def pair_all(tiles, root, product="DEM", res="0p5m"):
    """Resolve many tiles. Returns (pairs, failures) rather than raising.

    A survey-scale run should report which tiles it cannot score and process the
    rest, not stop on the first gap.
    """
    index = index_archive(root, product=product, res=res)
    pairs, failures = {}, {}
    for t in tiles:
        name = t if isinstance(t, str) else str(t)
        try:
            pairs[name] = reference_for(name, index=index)
        except PairingError as exc:
            failures[name] = str(exc)
    return pairs, failures


if __name__ == "__main__":
    import json
    import sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from paths import DEM_ROOT, ROOT

    manifest = os.path.join(ROOT, "tests", "full_clear", "batch_manifest.json")
    with open(manifest, encoding="utf-8") as fh:
        rec = json.load(fh)
    tiles = sorted(v["tile"] for v in rec.values()
                   if isinstance(v, dict) and v.get("tile"))
    pairs, failures = pair_all(tiles, DEM_ROOT)
    print("%d tiles: %d paired, %d unpaired"
          % (len(tiles), len(pairs), len(failures)))
    if failures:
        print("\nunpaired:")
        for t, why in sorted(failures.items()):
            print("  %-46s %s" % (t, why))
    # the pairing must also be injective: two tiles sharing a reference would
    # mean the key is not specific enough after all
    seen = {}
    for t, a in pairs.items():
        seen.setdefault(a, []).append(t)
    clashes = {a: ts for a, ts in seen.items() if len(ts) > 1}
    print("\narchives claimed by more than one tile: %d" % len(clashes))
    for a, ts in sorted(clashes.items()):
        print("  %s\n    %s" % (a, "\n    ".join(ts)))
