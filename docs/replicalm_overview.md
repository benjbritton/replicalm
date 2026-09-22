# Replicalm

### What it is, why it exists, and how it compares to TerraScan and Golden Surfer at reproducing the NCALM method

Benjamin Jay Britton, September 2026

---

## What it is

Replicalm is a software tool that turns raw airborne lidar point clouds into
bare-earth elevation models and archaeological relief visualizations. It
reproduces the processing workflow that NCALM used on the G-LiHT surveys of the
Maya lowlands, as documented in the supplementary materials of Estrada-Belli et
al. 2025.

It runs on free and open software throughout: PDAL for point-cloud handling and
ground classification, GDAL for coordinate systems and raster output, NumPy and
SciPy for the interpolation, and the Relief Visualization Toolbox for the
optional image step. It is MIT licensed.

## Why it exists

The published method depends on four commercial products: TerraScan for ground
classification, Golden Software Surfer for interpolation, ArcGIS Pro, and paid
LAStools modules. Together these represent several thousand dollars of licensing
before a single tile is processed.

That creates two problems for the field. Published results cannot be
independently reproduced by anyone without those licenses — which is most
researchers, most students, and most institutions outside well-funded programs.
And the method cannot be applied to new surveys by those same people, so the
data accumulating from lidar campaigns can only be processed by those who can
pay to process it.

Replicalm closes that gap. The method is the same method; the tools underneath
it are ones anyone can install.

## What it is useful for

The immediate use is processing G-LiHT and comparable airborne lidar over
archaeological landscapes: producing a ground surface with vegetation removed,
and from it the shaded relief images in which platforms, plazas, causeways,
terraces and looter trenches become legible.

Beyond that, three things make it useful as infrastructure rather than as a
one-off script:

**Every run records its own settings.** A configuration file is written beside
each output, so any raster can be traced back to the parameters that produced
it. The tool refuses to run silently under settings that differ from the
validated baseline — it says so.

**The output grid is declared rather than derived.** Neighboring tiles land on
the same lattice and mosaic without resampling, and rerunning reproduces the
previous result exactly.

**Parameters are measured, not inherited.** Each setting in the baseline carries
the measurement that chose it, recorded in the code beside the value.

## How it works

Given a point cloud, five stages:

1. **Classify ground.** Decide which returns are the soil and which are canopy,
   structures or brush. This is the hard part and the part that determines
   everything after it.
2. **Remove outliers.** Flag returns that cannot be ground.
3. **Measure height above ground** and apply the source method's cuts, including
   its near-ground class band.
4. **Interpolate** the scattered ground returns onto a regular grid by ordinary
   kriging.
5. **Finish the raster** — fill enclosed gaps, trim the unreliable fringe where
   cells were calculated from data on one side only, and write it so that no
   data is unambiguous.

An optional sixth stage produces the G1 relief composite from the finished
surface.

---

## How it compares to TerraScan and Golden Surfer

### What cannot be ported, only translated

TerraScan's "Classify Ground" implements progressive TIN densification after
Axelsson (2000). The source method drives it with four parameters: a 70 m
window, an 88° maximum terrain angle, iteration angles of 9° and 12° on two
successive passes, and a 3 m iteration distance.

No free tool implements that algorithm with that parameter surface. PDAL offers
SMRF (Pingel et al. 2013), a morphological filter; PMF (Zhang et al. 2003); and
CSF (Zhang et al. 2016), a cloth simulation. None of them has a terrain angle or
an iteration angle. ArcGIS Pro's equivalent offers three presets and exposes
none of the four.

So the classification stage is a translation, not a port, and any claim that it
has been replicated exactly would be false. Angles were converted to slope
tolerances by taking their tangent — tan 9° = 0.1584 — and the translation was
then checked against reference surfaces produced by the original workflow rather
than assumed correct.

Golden Surfer's kriging is likewise proprietary. The source specifies a 20 m
search radius and says nothing about the variogram model, its parameters, or the
maximum number of points per solve. Those had to be chosen and are recorded.

### What is carried over unchanged

The height-above-ground cuts at −0.5 m and 600 m, the ±0.2 m near-ground class 8
band, ground-only export, LAS 1.2 output, and a 20 m kriging search radius as
the maximum.

![The Pixoyal group under three treatments](figures/fig10_clear_vs_baseline.png)

*The Pixoyal group on South_GLAS_l0s395. Left: the original TerraScan and Surfer
workflow. Center: Replicalm's baseline — note the speckle on the plaza floors,
residual low vegetation the ground filter accepted. Right: Replicalm with
`Clear`. Same point cloud, same visualization recipe, all three normalized over
the same extent.*

### How closely the output agrees

Measured on `South_GLAS_l0s395`, a tile containing known structures, against the
reference surface produced by the original commercial workflow, over the 12.29
million cells both deliver:

| measure | value |
|---|---|
| median difference | 1.5 mm |
| RMSE, whole tile | 0.0601 m |
| cells more than 0.5 m apart, below 5° slope | 0.09% |
| between 5° and 10° | 0.17% |
| between 10° and 20° | 0.39% |
| between 20° and 30° | 1.01% |
| above 30° | 2.25% |
| coverage delivered | 99.56% |

Half the surface agrees to within a millimeter and a half. The disagreement is
concentrated almost entirely on steep ground: 97.7% of all error above half a
meter sits on slopes steeper than 10°, which are 15.5% of the tile.

### Where the two approaches genuinely differ

**Ground classification is where the real difference lies.** TerraScan rejects
low vegetation more completely than a translated SMRF pass does. Comparing
against TerraScan's own classification, which the delivered clouds carry, 3.88%
of the returns Replicalm calls ground were explicitly rejected by the original
workflow — and at the cells where that shows up visually as speckle on flat
ground, the figure is 12.2% against 0.8% on clean ground.

Replicalm's `Clear` rule recovers 92.5% of that without TerraScan, by demoting
returns that stand more than 0.20 m above the local ground floor. It removes
17.75% of ground returns to do what TerraScan does by removing 3.88%, so it
takes genuine returns with the clutter, and platform edges come out slightly
softer as a result.

**Interpolation differs in how the neighborhood is chosen.** The source fixes
the search radius at 20 m. Replicalm scales it to the measured point density,
with 20 m as the ceiling, because a fixed radius means different things at
different densities — at high density it packs the neighborhood so tightly that
the kriging system loses rank.

**The output grid differs by design.** The reference products have cells of
0.500042 × 0.499988 m on origins that are multiples of nothing, so neighboring
tiles do not align without resampling. Replicalm puts cell edges on whole
multiples of the cell size. It also derives the cell size from measured ground
density rather than fixing it, since a survey at one return per square meter
cannot support a 0.5 m raster and one at sixteen is wasted on it.

**Tiling differs for reasons that are not methodological.** The source works in
1 km tiles because a whole transect will not fit in memory on modest hardware.
Replicalm's interpolation is internally chunked, so it processes whole extents
by default, with tiling available where classification memory still demands it.

### Which is better

On elevation agreement over ordinary terrain, the two are equivalent — the
median difference is smaller than the sensor's own precision.

On steep ground, the original workflow is smoother and Replicalm's baseline
carries more texture, some of which is real micro-relief and some of which is
residual vegetation. With `Clear` applied, Replicalm predicts held-out ground
returns 13% more accurately than its own baseline and leaves roughly half as
many returns stranded beneath the modeled surface — measures that do not
reference the commercial output at all, and so can be compared without assuming
it is correct.

On reproducibility Replicalm is better, and not marginally: settings are
recorded per run, the grid is deterministic, and a configuration that has
drifted from the validated baseline stops the run rather than quietly producing
something different.

On cost the comparison is not close.

What cannot honestly be claimed is that Replicalm is more *accurate* than the
commercial workflow. Establishing that would require ground truth neither has,
and the measurements here compare a surface against another surface, or against
the point cloud both were derived from.

---

## Limits, and what is not settled

Nine items are recorded in `open_observations.md` rather than smoothed over. The
substantial ones:

- One comparison on a steep window produced a large difference that no isolated
  variable has yet accounted for. It was attributed twice and both attributions
  were wrong.
- The statistical outlier filter helps on one tile and hurts on another, and
  what distinguishes them is not known.
- `Clear` leaves a 14-point gap against what TerraScan's own labels achieve, and
  closing it needs a discriminator that separates clutter from small real
  features. At four returns per square meter a 0.4 m bush and a 0.4 m rock are
  described by about five points each, with the same geometry.

Results obtained before 20 September 2026 were measured on sample windows now
known to contain almost no steep ground, and should be treated as untested.

## Status

The processing pipeline is complete and validated. A Windows installer wrapping
a pinned environment is assembled but not yet compiled. `Clear` is validated and
available but deliberately not the default, since adopting it is a change to the
locked baseline and that is a decision to be made rather than inherited.

---

— Benjamin Jay Britton, 2026
