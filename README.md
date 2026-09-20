# Replicalm

**Bare-earth lidar processing that reproduces a commercial workflow with open tools.**

Replicalm turns a raw airborne lidar point cloud into a bare-earth digital
elevation model, following the NCALM workflow described in the Estrada-Belli et
al. 2025 supplementary material — using only PDAL, GDAL, NumPy and SciPy.

The published method depends on TerraScan, ArcGIS Pro, Golden Surfer and paid
LAStools modules. Without those licences the results cannot be reproduced and
the method cannot be applied to new surveys. This closes that gap.

![Archive output, a flat-calibrated configuration, and the locked baseline](docs/figures/fig1_flank_artifact.png)

*Left: the reference surface from the original commercial workflow. Centre: a
configuration tuned on flat sample windows. Right: the locked baseline. Same
point cloud, same 0.5 m grid, same visualisation recipe.*

---

## Results

Tested on `South_GLAS_l0s395` against the reference surface from the original
workflow, over the 12.29 million cells every configuration delivers:

| configuration | RMSE | 20–30° | 30–90° | terracing |
|---|---:|---:|---:|---:|
| tuned on flat windows | 0.3370 m | 17.38% | 40.38% | 9.83% |
| SMRF, single pass | 0.0678 m | 1.10% | 4.85% | 0.85% |
| **locked baseline** | **0.0601 m** | **1.01%** | **2.25%** | **0.26%** |
| reference | — | — | — | 0.00% |

Percentages are cells more than 0.5 m from the reference, by slope band.
Terracing is the share of steep cells whose gradient has collapsed below a
quarter of the reference's — a step where the reference has a continuous face.

Median difference from the reference is **1.5 mm**. Coverage is 99.56% and does
not vary between configurations, so none of the improvement comes from dropping
difficult cells.

![Error concentrates on steep ground](docs/figures/fig3_error_by_slope.png)

## The finding that shaped this project

The pipeline scored well and still produced a visible artifact on mound flanks.
The parameters were not at fault — the sample windows were. The window selector
chose the 400 m square with the highest reference coverage, which by
construction finds flat, open terrain. Every calibration window contained
between 0.00% and 0.20% of cells steeper than twenty degrees.

On this tile, 97.7% of all error above half a metre sits on slopes steeper than
ten degrees, which are 15.5% of the ground. Below ten degrees every
configuration agrees to within 0.08% of cells; above thirty they range from
2.25% to 40.38%. The calibration was being scored almost entirely on terrain
where the answer does not matter — and archaeological features are not on flat
ground.

`clip.best_window` now gates on coverage and then maximises relief, with a
ground-return density floor so that a window with relief but no returns beneath
it fails loudly instead of quietly.

Full account: [docs/posts/2026-09-20-calibrating-on-the-wrong-ground.md](docs/posts/2026-09-20-calibrating-on-the-wrong-ground.md)

## Quick start

```bash
# PDAL and GDAL come with ArcGIS Pro's Python, or install via conda:
conda install -c conda-forge pdal python-pdal gdal numpy scipy
```

```python
import config
from replicalm import classify, grid, kriging, finalise, interpolate

cfg = config.load()                 # the locked baseline
config.verify_baseline(cfg)         # stops the run if anything has drifted

g = grid.grid_for_las("tile.las", cell=1.0)
classify.classify_tile("tile.las", "ground.las", cfg)

arr, _ = classify.read_points("ground.las")
gnd = arr[arr["Classification"] == 2]
x, y, z = gnd["X"], gnd["Y"], gnd["Z"]

density = len(z) / ((g.bounds[2]-g.bounds[0]) * (g.bounds[3]-g.bounds[1]))
radius = kriging.radius_for_density(density, cfg.max_points)
v = kriging.fit_variogram(x, y, z)
dem, info = kriging.krige_grid(x, y, z, g, radius=radius, variogram=v)

finalise.finalise(dem, g, interpolate.source_srs("tile.las"),
                  "dem.tif", radius_m=radius)
```

`python config.py` prints the baseline and the source method it translates.

## The locked baseline

| setting | value | why |
|---|---|---|
| ground filter | `filters.smrf`, one pass | second pass gave nothing across seven tiles and doubled runtime |
| slope | 0.1584 | tan 9°, the source's pass-1 iteration angle |
| threshold | 0.5 m | the source's 3.0 m admits almost everything: 17.68% vs 10.34% error above 30° |
| working grid | 0.5 m | tied to the output cell; at 1.0 m steep faces terrace |
| ELM filter | off | removed 0.00% of points on every window tested |
| outlier filter | on | helps on one tile, hurts on another; sweep it per survey |
| search radius | density-scaled, ≤ 20 m | the ceiling is the source's figure |
| neighbours | 16 | binds before the radius does on well-covered ground |

## Relation to the published method

**Carried over unchanged:** 1 km tiles with 10 m buffers, the −0.5 m / 600 m
height-above-ground cuts, the ±0.2 m class 8 near-ground band, ground-only
export, LAS 1.2 output, a 20 m kriging search radius as the maximum, a 1 m
output grid.

**Translated, not ported:** the ground classification. TerraScan's progressive
TIN densification (Axelsson 2000) is not implemented by any free tool, and its
four parameters have no equivalent in SMRF. Angles were converted to slope
tolerances by taking their tangent, and the result was checked against reference
surfaces rather than assumed correct.

**Genuinely different:** one classification pass instead of two; a 0.5 m
elevation tolerance instead of 3 m; a declared output grid on whole multiples of
the cell size, so neighbouring tiles mosaic without resampling; and reduced
noise filtering. Each departure carries the measurement behind it.

Detail: [docs/processing_report.md](docs/processing_report.md)

## Layout

```
README.md               this file
config.py               baseline enforcement and configuration defaults
src/replicalm/          the processing modules
docs/
  processing_report.md  plain-language description and parameter audit
  open_observations.md  measured but unexplained; five open items
  posts/                technical write-ups
  figures/
benchmarks/             the calibration and validation harness
  results/              measurement output as JSON
archive/                superseded renders and experiments (not tracked)
```

## What is not settled

Five items are recorded in [docs/open_observations.md](docs/open_observations.md)
rather than smoothed over. The largest: one comparison on a steep window
produced a 2× difference in RMSE that no isolated variable has yet accounted
for. It was attributed twice and both attributions were wrong.

Results obtained before 20 September 2026 were measured on sample windows now
known to contain almost no steep ground. They should be treated as untested.

## Data

Point clouds and reference surfaces are G-LiHT products from the NASA
AMIGACarb / Yucatán campaigns and are not redistributed here. The benchmark
harness records which tiles it used and how each sample window was selected, in
`benchmarks/results/`.

## Licence

MIT — see [LICENSE](LICENSE).

---

Benjamin Jay Britton, 2026
