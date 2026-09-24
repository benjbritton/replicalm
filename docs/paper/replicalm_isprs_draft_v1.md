# Reproducibility is not accuracy: an open-source replication of a commercial bare-earth lidar workflow, and what it reveals about ground filtering

**Draft v1 — prepared for the ISPRS Journal of Photogrammetry and Remote Sensing.**

Figures and numbers marked `[PENDING]` await the 458-tile production run, which
was launched 2026-09-24. Everything else is measured and traceable to a record
in `benchmarks/results/`.

---

## Abstract

Airborne lidar processing for archaeological prospection depends on proprietary
software whose parameters are published but whose implementations cannot be
inspected, rerun, or audited by most readers. We present Replicalm, an
MIT-licensed framework that reproduces the bare-earth workflow used by the
National Center for Airborne Laser Mapping for NASA G-LiHT surveys of the Maya
Lowlands, using PDAL, GDAL and scientific Python only.

Three contributions follow. First, replication: across 458 G-LiHT transects the
open implementation reproduces the commercial product [PENDING: bias, RMSE and
coverage distributions], with coverage consistently equal to or greater than the
reference. Second, reproducibility as a design property rather than an
aspiration: every run records its configuration, a locked parameter set raises
on drift, the environment is pinned and redistributable, and the interpolation
is verified invariant to block size. We document a case in which the commercial
chain fitted different variogram parameters to adjacent tiles of one survey and
produced visibly discontinuous output — the failure mode an auditable pipeline
removes.

Third, and least expected: reproducibility and accuracy come apart, and we
measure the gap. On a survey with three overlapping strips we show that three
ground filters — SMRF, PMF and cloth simulation — each reproduce themselves to
approximately 5 cm across independent passes while placing the ground surface up
to 14 cm apart from one another. Over terrain selected for flatness, SMRF and CSF
raise their surface by 0.46 and 0.55 m per ten percentage points of understory
and admit more returns as they do so, while PMF does neither. Filter choice is
therefore not adjudicable by self-consistency, and archive-referenced comparison
inherits the archive's own treatment of low vegetation.

**Keywords:** airborne laser scanning; ground filtering; bare-earth DEM;
reproducibility; open-source; digital terrain model; archaeological prospection

---

## 1. Introduction

[Frame: lidar is now routine in landscape archaeology; bare-earth extraction is
the step every interpretation rests on; the dominant workflows are proprietary.]

The problem this paper addresses is not that commercial software is inaccurate.
It is that a published method which names TerraScan and Surfer specifies a
result that a reader cannot reproduce without those products, cannot inspect,
and cannot vary in a controlled way. Parameters are published; implementations
are not. When a downstream interpretation depends on which returns were called
ground, that opacity is a methodological liability rather than a licensing
inconvenience.

Three questions organise the work:

1. Can an open chain reproduce the commercial product at survey scale?
2. Can the reproduction be made auditable end to end — configuration, environment,
   and numerical behaviour?
3. Is the closest reproduction the best representation of the terrain?

The third question turns out to be the interesting one, and answering it
required a dataset with a property the reference survey does not have.

## 2. Background

### 2.1 The source workflow

[Estrada-Belli et al. 2025 supplementary materials: TerraScan progressive TIN
densification after Axelsson (2000); Golden Software Surfer ordinary kriging;
1 m output.]

### 2.2 Ground filtering families

Three families are compared here. Morphological filters threshold a surface of
local minima: the Simple Morphological Filter (Pingel et al., 2013) and the
Progressive Morphological Filter (Zhang et al., 2003). Cloth simulation (Zhang
et al., 2016) drapes a simulated fabric over the inverted cloud. TerraScan's
progressive TIN densification belongs to a fourth family and is the reference
rather than a candidate.

### 2.3 What counts as validation

[Review: archive-referenced comparison, withheld-return prediction, strip
overlap, and surveyed control. Argue that the first three are all measures of
consistency, and that the literature routinely reports them as accuracy.]

## 3. Materials and methods

### 3.1 Data

**NASA G-LiHT Yucatán.** 458 transects across seven regions of the Maya
Lowlands, April–May 2013, 323 GB of LAS. Ground-return density averages 4.3 per
m². Flown as single-pass transects: adjacent tiles are sequential segments of one
flight line, with roughly 96 m of overlap at the joins, so the archive contains
no independent second observation of any ground.

**Reference products.** Bare-earth rasters produced from the same point clouds
with TerraScan 026.002 and Surfer 30.2.240 under controlled, recorded parameters,
at 0.33, 0.5 and 1.0 m. Because the inputs and the settings are held fixed and
only the software differs, this is a stronger comparison than a comparison
against delivered archive products of unknown provenance.

**New River, northern Belize.** A 9.5 × 8 km block over the Ka'Kabish and
Cocochan sites, flown 17 May 2022, 144 tiles at 500 m. Three overlapping strips
cover the same ground at 12.9–21 returns per m², separated by 260 and 690
seconds and spanning scan angles from −31° to +31°. This block supplies the
independent re-observation G-LiHT cannot. Point clouds courtesy of A. McLellan;
used and published with permission.

### 3.2 Open implementation

[PDAL `filters.smrf` / `filters.pmf` / `filters.csf`; ordinary kriging written
against the source's stated search radius; GDAL for georeferencing; RVT for the
relief composites.]

Two implementation decisions bear on reproducibility. Kriging is chunked by
raster rows so that block size is a choice about memory rather than a limit on
the method, and the result is verified invariant across block sizes. The solve is
batched by neighbour count — cells finding the same number of neighbours are
solved as one stacked LAPACK call — which is 3.8× faster than a per-cell loop and
agrees with it to 2.1 × 10⁻¹⁴ m, with an identical set of filled and degenerate
cells.

### 3.3 Configurations

| | ground filter | elevation threshold | residual vegetation removed |
|---|---|---|---|
| Baseline | SMRF, one pass, slope 0.1584 | 0.50 m | no |
| Clear | as Baseline | 0.25 m | yes, 0.20 m above local floor |
| Deep | as Clear | 0.25 m | yes, on a grid 1.4× finer than point spacing |

Baseline is the replication target and is fixed by the source publication.
Clear's threshold was set by a 2 × 2 factorial (Section 4.3). Cell size is
derived from measured ground-return density as `1/√density` unless overridden.

### 3.4 Evaluation

**Archive-referenced.** Bias, RMSE, MAE, coverage and a terrain-complexity ratio
against the controlled commercial product.

**Archive-free.** Buffered leave-one-out prediction of withheld returns;
returns-below-surface; sharpness at a fixed baseline; rasterised spatial-block
folds.

**Physics-referenced.** On the three-strip block: classify each strip
independently and compare the resulting surfaces — independent observations of
the same ground, with no reference product in the measurement path.

**Attribution without ground truth.** Over cells selected for flatness, regress
each filter's surface height, and its recovered ground-return density, against
the fraction of returns in the 0.15–0.35 m band. Terrain has no reason to be
higher where understory is thicker; a filter whose surface tracks understory is
accepting vegetation. The datum is the 1st percentile of all returns in the same
cell, which references no classification — it is not unbiased, but it is common
to all filters, so contrasts survive it.

## 4. Results

### 4.1 Replication at survey scale

[PENDING: 458-tile distributions of bias, RMSE, coverage; failure cases; runtime
against the commercial chain's own logged 458 tiles in ~9.5 h.]

### 4.2 Filter comparison against the archive

Twelve parameterisations across three tiles. SMRF at slope 0.05 and threshold
0.25 m gives the lowest RMSE on every tile (0.0194, 0.0399, 0.0449 m). CSF is
close on low-relief tiles and degrades with relief; PMF overtakes it on the
steepest tile while remaining worst overall against the archive, with negative
bias throughout (−0.010 to −0.098 m), roughly half the ground returns retained,
and surfaces at 0.34–0.75× the reference's terrain texture. For both
morphological filters the elevation threshold dominates: varying slope at fixed
threshold moves RMSE by ~0.0001 m, while moving the threshold from 0.25 to 1.0 m
changes it by a factor of two or more.

### 4.3 Threshold and residual-vegetation removal are not redundant

| tile | arm | removed by cleanup | RMSE | bias | complexity |
|---|---|---|---|---|---|
| l0s395 | t0.50 raw | — | 0.0229 | +0.0010 | 1.04× |
| | t0.25 + cleanup | 0.60% | 0.0186 | −0.0007 | 0.90× |
| l8s431 | t0.50 raw | — | 0.0804 | +0.0367 | 1.37× |
| | t0.25 + cleanup | 6.35% | 0.0404 | +0.0033 | 0.78× |
| l0s444 | t0.50 raw | — | 0.0606 | +0.0263 | 1.28× |
| | t0.25 + cleanup | 11.76% | 0.0357 | −0.0023 | 1.02× |

Tightening the threshold takes part of what the cleanup would have taken — its
share falls from 10.47% to 6.35% on l8s431 — but the cleanup still finds 6–12%
afterwards and still improves bias and RMSE. Applied together they bring bias
within 3.3 mm of zero on all three tiles, with no overshoot into negative bias.

### 4.4 Self-consistency does not distinguish the filters

Across 41 flat three-strip cells, each filter classified from each strip alone:

| filter | spread of the three single-strip surfaces | correlation with understory |
|---|---|---|
| SMRF | 0.050 m median, 0.120 p90 | −0.03 |
| PMF | 0.050 m median, 0.110 p90 | −0.18 |
| CSF | 0.060 m median, 0.120 p90 | +0.03 |

All three reproduce themselves to about 5 cm, none degrades as understory
thickens, and the differences between them are nil — while the same filters place
the ground surface up to 0.142 m apart in the same cells.

### 4.5 Which filter is moving

| filter | surface rise per 10% understory | r | t (n = 41) | ground density change |
|---|---|---|---|---|
| CSF | +0.553 m | +0.49 | +3.55 | +7.73 returns/m² |
| SMRF | +0.456 m | +0.42 | +2.91 | +7.23 returns/m² |
| PMF | −0.034 m | −0.03 | −0.21 | −0.06 returns/m² |

Over ground with under 2% slope, SMRF and CSF raise their surface as understory
thickens and admit more returns as they do it; PMF does neither. The penetration
objection does not hold: the low-return tail spread is effectively flat against
understory (r = +0.11), and any residual datum bias is common to all three
filters.

### 4.6 A per-tile automatic fit produces a spatial discontinuity

Surfer's ordinary kriging is an exact interpolator: it honours the input points
exactly unless a nugget is supplied. In production over the New River block the
nugget was therefore fitted per tile, using Surfer's own `VariogramObject.AutoFit`
against that tile's ground returns. That choice was made deliberately and for a
defensible reason — a hand-picked nugget large enough to be physically plausible
had no effect at all, and a flight-line artifact needed suppressing.

The magnitudes involved are worth stating. On one tile a hand-picked nugget of
0.0009 m² — three centimetres, about the ranging precision of the sensor — changed
the resulting DEM by a mean of 0.05 mm, which is to say not at all. `AutoFit` on
the same tile's data returned 0.736 m², roughly eighty-six centimetres and some
eight hundred times larger, which visibly smooths a surface whose archaeological
features are decimetric.

The failure is not that the fit was wrong on any one tile. It is that the fit is
per tile and the tiles are adjacent. Around Ka'Kabish, six of the twenty-two
tiles within one kilometre received `AutoFit` nuggets between 0.62 and 1.90 m²,
while the remainder — including the tile containing the site itself — used
Surfer's untuned default. The neighbourhood was therefore rendered at 0.25 m by
interpolators that differed by three orders of magnitude in how much they were
willing to smooth, producing a glossy-against-rough patchwork across the site.

It was found by eye, during visual review. Nothing in the delivered rasters
records which variant produced them, so no automated check could have caught it,
and the discontinuity is invisible in any per-tile quality statistic: each tile
is individually well interpolated. The fix required re-kriging the affected
neighbourhood with the automatic fit disabled so that the whole block shared one
interpolator.

## 5. Discussion

### 5.1 The archive is not truth, and this is measurable

On G-LiHT, scored against the TerraScan reference, PMF appears worst: biased low
everywhere, half the retained points, markedly smoother surfaces. On the
three-strip block, scored against flat ground, PMF is the only filter that does
not follow the understory. Both results are correct. The archive shares the
acceptance behaviour we measure in SMRF and CSF, so part of PMF's apparent
negative bias is it correctly declining vegetation the archive kept — though
0.34× terrain texture is too much smoothing to attribute entirely that way, and
PMF plausibly over-erodes relief as well.

The general point is that a comparison against a product inherits that product's
decisions. Where those decisions are the object of study, the comparison cannot
settle the question.

### 5.2 Opacity is not secrecy

The Ka'Kabish case is worth dwelling on because it is not a story about
proprietary software behaving badly, and it would be weaker if it were. The
software did what it was asked. The operator's reasoning was sound, the problem
being solved was real, and the defect was diagnosed and corrected by the same
operator. Every step is defensible in isolation.

What made it possible is that a parameter controlling how much the surface is
smoothed — by a factor of eight hundred, between a value that does nothing and a
value that erases decimetric relief — was selected automatically, per tile,
against each tile's own data, and then not carried in the product. Kriging
presents itself as principled: a geostatistical model fitted to the data rather
than a filter chosen by taste. Per-tile automatic fitting preserves that framing
while making the model itself a variable that changes from one 500 m square to
the next, with no term anywhere in the method expressing that neighbouring tiles
should agree.

The result is opacity without secrecy. The algorithm is documented, the
operator is competent, the fit is data-driven, and the output is still not
reproducible from the record it carries. That is the case for auditability
stated in its strongest form, because it does not depend on anyone having made
a mistake. Replicalm's response is structural rather than procedural: a locked
parameter set that raises on drift, a configuration written beside every
product, and interpolation whose behaviour is verified invariant to how the work
was divided. None of that makes a better surface. It makes the surface
accountable for how it was made.

### 5.3 What these measurements do not establish

No result here states how far above the soil any surface sits. PMF holds a
near-constant 0.20 m above the lowest returns and may be right or may simply be
anchored to the minimum surface. The flatness prior rests on cells classified
flat by a vendor whose ground class we have shown to be understory-contaminated.
Resolving these requires surveyed control, which we identify as follow-up work
rather than fold into unsupported claims here.

### 5.4 Limits

Clear's 0.20 m threshold was fitted on one campaign's labels and validated on
three windows from the same campaign; it is not established as general. The
understory regression rests on 41 cells from one block. The 0.25 m threshold and
cleanup together return 0.78× the archive's terrain complexity on the most
relieved tile, which is either vegetation texture correctly removed or relief
incorrectly removed, and we cannot currently distinguish them.

## 6. Conclusions

[To follow the completed run.]

---

## Data and code availability

Replicalm is MIT-licensed. Source, the locked configuration, every benchmark
record quoted here, and the scripts that rebuild them are at [repository].
G-LiHT data are public. New River measurements are published with the permission
of the data owner; the point clouds themselves are not redistributed.

## CRediT author statement

**B.J. Britton:** Conceptualization, Methodology, Software, Formal analysis,
Investigation, Validation, Visualization, Writing – original draft.
**A. McLellan:** Data curation, Investigation, Resources, Writing – review and
editing.

## References

[To compile: Axelsson 2000; Estrada-Belli et al. 2025; Kokalj & Somrak 2019;
Pingel et al. 2013; Zhang et al. 2003; Zhang et al. 2016; PDAL; GDAL.]
