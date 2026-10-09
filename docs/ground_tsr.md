# `ground_tsr` — ground classification from raw returns

A reusable ground classifier that produces a bare-earth point population from
raw LAS/LAZ alone. It was built to reproduce the ground population an NCALM
TerraScan/Surfer chain delivers, and it is calibrated against one, but it reads
no delivered classification while it runs and does not require any to be
present.

Module: `replicalm.ground_tsr`
Selected by: `Settings.ground_engine = "tsr"`
Status: opt-in. `ground_engine` defaults to `"smrf"`, the locked baseline.

---

## 1. Inputs and outputs

**Input.** Three arrays of equal length — easting, northing, elevation — in a
projected coordinate system with metres as the unit. Nothing else is read: not
return number, not intensity, not scan angle, not GPS time, not any existing
classification. A minimum of 32 returns is required and a block of under a few
thousand will usually fail stage 2 for want of context.

**Output.** `classify_ground(x, y, z, params=None, adapt_density=True)` returns
`(mask, report)`.

- `mask` — boolean array over the input returns, true where ground.
- `report` — dict: return count, measured density, the effective working cell
  and curvature scales after density adaptation, the per-stage counts, and the
  accepted count at each refinement iteration.

The report exists so a run can be audited without being repeated. A block whose
ground fraction departs from its neighbours' is the signature of a
classification that has gone wrong, and that is visible in the report before
anything is interpolated.

**Through the pipeline.** `classify.classify_tile` branches on
`cfg.ground_engine`. Under `"tsr"` the PDAL pipeline runs in two parts with the
numpy classification between them, because the routine is not a PDAL stage.
Steps 6 and 7 — height above ground, the near-ground band, the class filter,
the LAS writer — are the same stages the baseline path uses, including the
retry when the ground surface cannot be triangulated. Switching engines does
not silently switch the output format.

`cfg.passes` is ignored under `"tsr"`. The routine carries its own stage
parameters, so one object describes one classifier rather than two places
having to agree. Unknown keys in `cfg.ground_tsr_overrides` raise rather than
being ignored.

---

## 2. Why the design has three filters rather than one

A morphological opening of a minimum surface is conservative. At the Replicalm
baseline it keeps almost no vegetation but discards real ground on closed
landforms — measured tile-wide against the delivered population, 1.87% of
planar ground against 7.79% on mounds and 8.34% on terrace shoulders — and the
resulting surface sits high over depression floors.

A multiscale curvature test is the opposite. It recovers most of that ground,
reaching 1.28% on planar and 4.71% on mounds, and admits enough low vegetation
that the surface is worse overall: 0.0634 m against the archive where the
opening gives 0.0621 m.

These are different mistakes. One discards true ground, the other keeps false
ground, and neither is reachable by tuning the other filter. So the routine
separates the two decisions and adds a third stage that adjudicates between
them. The adjudication is the part no single filter can perform, because it is
the only stage that sees both answers.

---

## 3. Stages

| # | stage | what it does |
|---|---|---|
| 1 | density | median returns/m² on a 10 m grid; sets the working cell |
| 2 | confident | morphological opening, deliberately conservative |
| 3 | candidates | multiscale curvature test, deliberately permissive |
| 4 | reference | robust surface from the confident set alone |
| 5 | adjudicate | accept candidates inside a locally scaled, asymmetric gate |
| 6 | refine | rebuild the surface from the accepted set, adjudicate again |

**Stage 1.** The median rather than the mean, so a strip of overlap between
flight lines does not set the working cell for a whole block.

**Stage 4.** A low quantile per cell, not a minimum: the minimum of a cell
follows the lowest noise spike that survived denoising, while the tenth
percentile tracks the soil and is stable at these densities. Holes are closed
by normalised convolution — the surface and its mask are smoothed separately
and divided — so the smoothing never reaches across a boundary. Filling with a
scene constant first would put a step at every hole edge and feed it into
everything computed downstream.

**Stage 5.** A candidate is accepted when its residual to the reference
surface falls inside a gate whose width is set by the local scatter of the
confident points, plus an allowance proportional to surface slope. Where the
confident points scatter — rough ground, sparse coverage, steep faces — the
gate opens; on smooth well-sampled ground it closes to its floor.

The gate is asymmetric, by roughly a factor of five at the calibrated settings.
What stands above the ground is vegetation, structure or noise. What sits below
is either real microtopography or a dropout the noise filter should already
have taken. These are not symmetric situations and a symmetric tolerance treats
them as though they were.

The gate adds to the confident set and never subtracts from it. See §5.

---

## 4. Parameters

Defaults live in `GroundTSRParams`. Every length is metres; every tolerance is
metres of elevation.

| parameter | default | role |
|---|---|---|
| `density_cell_m` | 10.0 | grid the density is measured on |
| `density_ref` | 20.0 | returns/m² the defaults were calibrated at |
| `cell_min_m` / `cell_max_m` | 0.5 / 2.0 | clamp on the adapted cell |
| `smrf_window_m` | 32.0 | opening window |
| `smrf_slope` | 0.158 | 9°, as the baseline |
| `smrf_threshold_m` | 0.25 | elevation tolerance |
| `smrf_scalar` | 0.0 | slope does not scale with window |
| `smrf_cell_m` | 1.0 | working cell, validated full-tile |
| `mcc_scales_m` | (0.5, 1.0, 1.5) | curvature scales |
| `mcc_tolerance_m` | 0.30 | curvature tolerance |
| `mcc_iterations` | 4 | passes per scale |
| `surface_cell_m` | 2.0 | reference surface cell |
| `surface_quantile` | 0.10 | quantile taken per cell |
| `surface_smooth_m` | 3.0 | smoothing of the reference surface |
| `gate_scale_m` | 8.0 | window the local scatter is read over |
| `gate_k_above` | 0.75 | multiples of local scatter, upward |
| `gate_k_below` | 4.0 | and downward |
| `gate_floor_above_m` | 0.05 | floor on the upward gate |
| `gate_floor_below_m` | 0.15 | floor on the downward gate |
| `gate_ceiling_m` | 1.00 | ceiling on both |
| `slope_allowance` | 0.30 | metres of gate per unit of slope |
| `refine_iterations` | 1 | stage 6 repeats |
| `prune_confident` | False | measured and rejected; see §5 |
| `min_ground_fraction` | 0.02 | safety net under pruning |
| `outlier_neighbours` | 8 | statistical outlier filter |
| `outlier_multiplier` | 2.5 | statistical outlier filter |

### Density adaptation

`scale_for_density(density)` returns a copy with lengths adapted to measured
return density. Point spacing goes as the inverse square root of density, so
the working cell scales that way and is then clamped.

Clamping rather than extrapolating is deliberate. Beyond the calibrated band
the right behaviour is the nearest configuration that was actually measured,
not a projection of a trend that was never tested there. Measured behaviour:

| density | working cell | reference cell | curvature scales |
|---|---|---|---|
| 40 /m² | 0.71 | 1.41 | 0.35, 0.71, 1.06 |
| 20 /m² | 1.00 | 2.00 | 0.50, 1.00, 1.50 |
| 10 /m² | 1.41 | 2.83 | 0.71, 1.41, 2.12 |
| 5 /m² | 2.00 | 4.00 | 1.00, 2.00, 3.00 |
| 1 /m² | 2.00 (clamped) | 4.00 (clamped) | 2.24, 4.00, 4.00 |

Pass `adapt_density=False` to hold the parameters exactly as given.

`density_ref` is raw returns per square metre, not ground returns per square
metre. The two differ by roughly a factor of three on this coverage, and
setting it from the wrong one silently moves the working cell off its validated
value — which happened during development and cost a calibration pass.

---

## 5. What was measured, and what was rejected

**Stage contribution**, eight development blocks of 240 m, scored as DEM RMSE
after biharmonic interpolation against the reference surface defined in §6:

| configuration | RMSE vs reference | vs archive | recall | precision |
|---|---|---|---|---|
| opening alone | 0.0470 | 0.0629 | 0.980 | 0.977 |
| curvature alone, ungated | 0.0518 | 0.0669 | 0.990 | 0.960 |
| adjudicated | **0.0434** | **0.0607** | 0.980 | 0.977 |

The gain is about 8% over the better of the two filters, and it holds across
every gate setting in the additive family — which is what distinguishes it from
a lucky parameter.

**Precision binds, not recall.** Widening the gate moved recall from 0.980 to
0.995 and made the surface worse every time. Admitted vegetation costs more
than recovered ground returns gain. This is why `gate_k_above` is small.

**Pruning was implemented, swept and rejected.** Letting the gate reject
confident points as well as admit candidates collapses recall to between 0.39
and 0.74 while precision goes to 1.000: the confident points scatter about the
quantile surface by more than a gate tuned for candidates allows. Every pruning
configuration scored worse than every additive one. The switch is retained,
defaulting off, so the result is reproducible rather than merely asserted.

**The gate floor was not the active parameter.** In the first calibration pass
all three values of `gate_floor_above_m` gave identical scores, because the
computed gate was always wider than its floor and the slope allowance was
holding it open. The second pass swept the slope allowance and `gate_k_above`
below 1.0, which is where the usable settings turned out to be.

---

## 6. Validation

Two reference surfaces, because they answer different questions.

**The reference surface** is the delivered Class 2 + Class 8 ground population
passed through the Replicalm biharmonic engine. Both sides of the comparison
then go through the same interpolator, so the score moves only with
classification, and an exact reproduction of the delivered population would
read zero. This is the optimisation target.

**The archive** is the delivered 0.5 m DTM. It is the scale on which the
0.052 m benchmark was quoted, so it is reported alongside, but it conflates
classification and interpolation differences and is not optimised against.

Class 2 and Class 8 are treated as one population throughout. The delivery
carries no attribute by which they could be separated — `PointSourceId`,
`ReturnNumber`, `NumberOfReturns`, `Intensity`, `ScanAngleRank` and `GpsTime`
are all zeroed in the file — and at 5 m scale 99.5% of occupied cells contain
both. The split has no reproducible basis; the union does.

Point-wise F1, precision and recall are computed as diagnostics and are never
optimised. A classifier can raise agreement while making the surface worse, and
during development it repeatedly did.

Blocks are drawn at random from the archive's valid area under a fixed seed,
eight for development and six held out, with the two forensic sites forced into
development because they are the terrain the routine has to handle. The
held-out blocks are scored once, with the configuration chosen on development,
so the reported figure is not the figure that was optimised.

---

## 7. Operational limits

Calibrated on G-LiHT Riegl VQ-480 coverage over tropical karst with closed
canopy, at a measured 19 to 24 raw returns per square metre on the calibration
blocks. The density scaling in §4 is what carries the routine to other
acquisitions, and within the clamped band it interpolates between measured
configurations rather than extrapolating.

- **Buildings are not modelled.** This is a bare-earth routine for terrain. A
  roof large enough to defeat the opening window will survive into the ground
  class.
- **Water is not detected.**
- **Block context is required.** Stages 2 and 4 need surroundings. Run with a
  buffer of at least the opening window and keep only each block's core; the
  full-tile harness uses 400 m blocks with a 40 m buffer.
- **Density below about 5 returns/m² is outside the calibrated band.** The
  parameters clamp rather than extrapolate, so the routine will run, but its
  behaviour there is not measured.
- **No parameter was fitted per block.** One configuration runs everywhere.
  Per-block fitting would raise every score reported here and would not
  survive onto an unseen dataset.
