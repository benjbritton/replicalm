# Benchmark results

Measurement output from the harness, kept so that the numbers quoted in
`docs/` and `README.md` can be traced to the run that produced them without
re-processing anything.

## Paths are tokenised

Absolute drive letters have been replaced with tokens:

| token | resolves to |
|---|---|
| `${DEM_ROOT}` | reference DEMs, default `<REPLICALM_DATA>/_Archive_EdgeFixed` |
| `${LAS_ROOT}` | point clouds, default `<REPLICALM_DATA>/GLiHT_LAS_orig` |
| `${ROOT}` | the repository |
| `${WORK}` | the working directory for intermediates |

`benchmarks/paths.py` resolves them:

```python
from paths import expand
expand(record["las"])
```

This is an archive, not an input. The harness writes its own working files with
real paths as it runs; nothing here needs to be edited to re-run anything.

## What is not here

The point clouds and reference surfaces themselves — NASA G-LiHT products from
the AMIGACarb and Yucatan campaigns, not redistributed. `pool_traits.json`
records the measured traits of all 453 tiles that have both a cloud and a
reference surface, which is what the tile selection drew from, so the selection
is inspectable without the data.
