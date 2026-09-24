# Open observations

Things measured but not explained. Each entry records what was seen, what was
ruled out and what remains, so that a later run does not re-derive the same dead
ends. An entry leaves this file when it is explained, not when it stops being
inconvenient.

---

## 1. The trench-window result is unexplained

**Status:** open, as of 2026-09-20.

**What was measured.** On a 500 x 500 cell window at column 1500, row 5250 of
`South_GLAS_l0s395` — the window with the highest residual error in the tile,
containing looter trenches and steep mound edges — two configurations differed
by far more than anything else tried:

| arm | RMSE | terracing | error above 30 deg |
|---|---|---|---|
| two SMRF passes, filters on, slope 0.1584 / 0.2126 | 0.083 m | 0.25% | 2.33% |
| one SMRF pass, filters off, slope 0.05 | 0.161 m | 2.44% | 8.67% |

Both ran at `cell_m` 0.5 and `threshold_m` 0.5. Terracing is the share of cells
above twenty degrees whose own gradient has collapsed below a quarter of the
reference's — a step where the reference has a continuous face.

**What has been ruled out.**

*Pass count.* Recalibration across seven tiles with relief, holding the filters
on and varying only the number of passes, found the two indistinguishable:
steep support 98.1% either way, error above thirty degrees 10.34% against
10.18%, terracing 2.87% against 2.90%. One pass is marginally better on RMSE.
The second pass roughly doubles classification time and returns nothing.

*Noise filters.* A four-arm test on three relief windows (33.2%, 13.1% and
11.3% of cells above twenty degrees), one pass everywhere and only the filters
varying, found ELM removes 0.00% of ground points at any threshold — arms that
ran it are identical to arms that did not — and the statistical outlier filter
removes 0.05% to 0.36% with a tile-dependent effect that is not remotely large
enough to account for a 2x difference in RMSE and a 10x difference in terracing.

**What remains.** The SMRF `slope` term, 0.1584 against 0.05, is the only
variable of the three that has not been isolated. It was dismissed earlier on
the strength of a flank test that found 0.05, 0.20 and 0.40 within noise of each
other — but that test ran at `cell_m` 1.0 and never measured terracing, so it
cannot speak to this.

**Why it matters.** The accepted `South_GLAS_l0s395` G1 was rendered at
`slope=0.1584` and the locked baseline uses it. If slope is the active
parameter, the baseline is correct by inheritance rather than by measurement,
and the value has never been swept on terrain where terracing is visible.

**The test, not yet run.** One pass, filters fixed, `cell_m` 0.5, sweeping slope
across 0.05 / 0.10 / 0.1584 / 0.2126 / 0.40 on the three relief windows, scored
on terracing and error above thirty degrees. Roughly twenty minutes.

**Method note.** This was attributed twice and wrongly — first to the two-pass
structure, then to the noise filters — both times by eliminating variables from
a comparison whose variables had not been fully enumerated. Three differed, not
two. Enumerate first.

---

## 2. The statistical outlier filter helps or hurts depending on the tile

**Status:** open, low priority.

Removing 0.05% to 0.36% of ground points, measured on three relief windows:

| window | RMSE without | RMSE with | above 30 deg without | with |
|---|---|---|---|---|
| l0s395 | 0.1745 | 0.0913 | 5.22% | 2.89% |
| l0s417 | 0.1685 | 0.2195 | 8.20% | 8.43% |
| l2s444 | 0.0481 | 0.0489 | 2.13% | 1.79% |

It is on in the baseline because the aggregate favors it, but the aggregate is
carried by one window. What distinguishes `l0s395` from `l0s417` is not known.
Worth sweeping per survey rather than trusting.

---

## 3. `thin_declustered` does not deliver the density it is asked for

**Status:** open.

Asking for 6.00 points per square meter returned 2.38; asking for 2.08 returned
1.32. The undershoot is roughly 1.5x to 2.5x and scales with how clustered the
returns are, because the routine keeps one point per cell of a grid sized from
the target and many cells are empty where returns follow scan lines.

The parameter is therefore not a density and should either be solved for or
renamed. Nothing currently depends on its accuracy, because thinning turned out
to buy no accuracy (see below), but any future use of it does.

---

## 4. Thinning restores the kriging solve without improving the result

**Status:** understood, recorded so it is not mistaken for a quality setting.

On `l0s444`'s relief window at 8.18 ground points per square meter, 76% of cells
fell back to inverse distance because the sixteen nearest neighbors sit within
about 0.8 m and every pairwise semivariance is nearly identical, so the matrix
loses rank. Thinning drops the fallback to 0.3%. It moves RMSE by 0.002 m.

Reducing the search radius does not help — 76.6% at 20 m against 76.4% at
3.16 m — because `max_points` binds first and the same sixteen points are
selected either way.

So thinning buys the kriging variance raster and the right to call the output
kriged rather than inverse-distance-weighted. It does not buy elevation
accuracy. Describe it that way.

---

## 5. Calibration windows were selected against the terrain that matters

**Status:** fixed in `clip.best_window`, recorded because it invalidates
earlier conclusions.

Windows were chosen for maximum reference coverage, which finds flat, open
ground. Measured after the fact, every calibration window contained between
0.00% and 0.20% of cells steeper than twenty degrees. On `l0s395` the error
below ten degrees is 0.01% to 0.08% of cells past half a meter under every
configuration tried, and above thirty degrees it ranges from 2.25% to 40.38%;
97.7% of all large error sits on slopes above ten degrees, which are 15.5% of
the tile.

Every parameter chosen on those windows — CSF over SMRF, one pass over two,
filters off, the tightened threshold — was chosen where the alternatives are
indistinguishable. Any conclusion in this project dated before 2026-09-20 and
derived from a coverage-selected window should be treated as untested rather
than as established. That includes the point-density and thinning results.

---

## 6. The flat-ground speckle is residual vegetation, and only partly removable

**Status:** diagnosed, partly fixed by `cleanup.clear()`, gap quantified.

**What it is.** Dark specks on otherwise flat ground, absent from the reference.
Median cluster size one cell. Our surface sits above the reference at 78.9% of
them, by a median of 2.6 cm, and local roughness there is 0.0768 m against the
reference's 0.0358 m.

**What they are made of.** Each speck carries a median of five classified ground
returns. They sit 5.8 cm above the surrounding 3 m median, 76.2% positive — a
one-sided bias, so not ranging noise, which would be symmetric. But they scatter
among themselves by 8.8 cm, more than the feature they describe. Where two
source chunks both contribute they disagree by 0.47 m against 0.11 m on control
ground. Mean scan angle is identical to control (11.5 against 11.2 degrees), so
incidence geometry is ruled out, and the returns are 11% darker.

Something intercepts pulses at variable heights regardless of viewing angle and
reflects less than soil or limestone. Low vegetation does exactly that. The
hypothesis is Ben's; the measurements above were taken to test it and did not
refute it.

**The ceiling.** The delivered clouds carry TerraScan's own classification. At
speck cells 12.2% of our ground returns are unclassified there, against 0.8% on
control ground. Removing exactly those — 3.88% of returns — gives roughness
1.00x the reference, specks down 94%, RMSE 5.9x better. So it is wholly a
classification residue.

**What is open.** `cleanup.clear()` recovers 92.5% of that by removing 17.75% of
ground returns where the oracle removes 3.88%. The collateral softens platform
edges. Closing the 14-point gap needs a discriminator that separates clutter
from small real features, and the ones tried do not.

---

## 7. Three clutter rules that do not work

**Status:** closed negative. Kept so they are not re-derived. Code in
`cleanup.py`.

| rule | synthetic brush | synthetic rock | real specks | real steep ground |
|---|---|---|---|---|
| height above patch floor | 72% removed | 100% removed | — | — |
| spread, slope-compensated | 54% | 92% | — | — |
| source disagreement | — | — | 3.4% hit | 10.2% hit |

The first two destroy a 0.4 m rock as readily as a 0.4 m bush: a feature
narrower than the patch has neighbours on the flat ground around it, so the
local floor stays low and the whole feature reads as high. The third damages
architecture three times harder than it cleans flat ground, because on a slope
two source chunks sample different parts of the same real surface.

At four ground returns per square metre a bush and a rock of the same footprint
are each described by about five points with the same geometry. The statistical
signal is real; it is not available per point.

**Note on method.** The synthetic test that condemned height-above-floor used a
0.4 m feature against a 0.6 m radius -- the pathological case. On real data at a
0.75 m patch it is the best discriminator available, d-prime 1.82, and is what
`clear()` uses. A synthetic test can be unrepresentative in either direction.

---

## 8. A nugget suppresses the specks by blurring them, and is the wrong fix

**Status:** closed. `fit_variogram(min_nugget_fraction=)` exists, defaults to
0.0, and is not used by the baseline.

Fitted nuggets come back at exactly zero on 34 of 60 fits across ten tiles,
making ordinary kriging an exact interpolator, so per-return noise passes into
the raster. Flooring the nugget smooths it away:

    0%  1.26x    2%  1.09x    4%  0.99x    6%  0.92x   20%  0.72x
    1%  1.16x    3%  1.04x    5%  0.96x   10%  0.83x

4% of sill matches the reference's flat-ground roughness. Measurement noise
alone justifies about 0.6%, so matching the reference means smoothing roughly
seven times harder than noise suppression requires. What produced the
reference's texture is not documented and is not established here; what is
established is that tuning to match it blurs real micro-relief. And the
classification fix reaches RMSE 0.0140 where the nugget reaches 0.0679.

**A measurement trap worth remembering.** Speck counts across nugget arms were
not comparable: G1 normalises each raster by its own extremes, and one outlier
sky-view-factor cell in the 5% arm restretched the whole image, making it look
44% better. Two arms were scored on a different tonal mapping than their
neighbours before this was caught.

---

## 9. The resolution knee sits at the mean point spacing

**Status:** measured on one window at one density. Implemented as
`grid.cell_for_density`.

Rasterised block cross-validation on the Pixoyal window at 4.74 ground returns
per square metre, mean spacing 0.46 m. Residual between the raster and returns
that did not build it:

| cell | median | marginal gain | cost vs 0.50 m |
|---|---:|---:|---:|
| 1.00 m | 0.0865 | — | 0.25x |
| 0.50 m | 0.0803 | -7.2% | 1.0x |
| 0.33 m | 0.0787 | -2.0% | 2.3x |
| 0.25 m | 0.0781 | -0.8% | 3.9x |

The curve turns at the mean point spacing. NCALM's 0.5 m suits this density
rather than being arbitrary, and the suspicion that it was costing detail was
wrong.

**What is open.** One density, one window. That the knee sits at the mean
spacing is physically sensible and is one measurement, not a law. Denser and
sparser tiles would test whether the knee moves as `1/sqrt(density)` predicts.

**Two metrics that could not answer this, and why.** `buffered_loo` predicts
point to point and never touches the raster, so it returned identical figures at
every cell size. Returns-below-surface improves as cells shrink for a reason
unrelated to quality: a coarse cell averages over more ground, so returns at the
low end fall beneath its single value. Both are sound for comparing
classifications and useless for comparing resolutions.

---

## 10. The tightened threshold and Clear together smooth below the reference on steep ground

**Status:** open, and the first question the Ka'Kabish surfaces are meant to
answer.

The SMRF elevation threshold and the `cleanup.clear()` rule remove overlapping
populations. A 2 x 2 factorial -- threshold 0.50 m against 0.25 m, crossed with
Clear on and off, at the locked slope of 0.1584 on l0s395, l8s431 and l0s444 --
settled how they interact:

| tile | arm | Clear removed | RMSE | bias | complexity |
|---|---|---|---|---|---|
| l0s395 | t0.50 raw | -- | 0.0229 | +0.0010 | 1.04x |
| | t0.50 + Clear | 0.80% | 0.0188 | -0.0006 | 0.90x |
| | t0.25 raw | -- | 0.0193 | +0.0002 | 0.95x |
| | t0.25 + Clear | 0.60% | 0.0186 | -0.0007 | 0.90x |
| l8s431 | t0.50 raw | -- | 0.0804 | +0.0367 | 1.37x |
| | t0.50 + Clear | 10.47% | 0.0499 | +0.0107 | 0.93x |
| | t0.25 raw | -- | 0.0476 | +0.0157 | 0.91x |
| | t0.25 + Clear | 6.35% | 0.0404 | +0.0033 | **0.78x** |
| l0s444 | t0.50 raw | -- | 0.0606 | +0.0263 | 1.28x |
| | t0.50 + Clear | 14.74% | 0.0405 | +0.0017 | 1.09x |
| | t0.25 raw | -- | 0.0431 | +0.0151 | 1.10x |
| | t0.25 + Clear | 11.76% | 0.0357 | -0.0023 | 1.02x |

They overlap without being redundant. Tightening the threshold takes some of
what Clear would have taken -- Clear's share falls from 10.47% to 6.35% on
l8s431 and from 14.74% to 11.76% on l0s444 -- yet Clear still finds 6 to 12%
afterwards and still improves both RMSE and bias. Together they put bias within
3.3 mm of zero on all three tiles, the best of any arm, with none of the
overshoot that compounding two removals might have produced. That is why
`clear` and `deep` now carry a 0.25 m threshold.

**What is open.** On l8s431, the tile with the most relief, the pair returns
**0.78x the reference's terrain complexity** -- 22% less texture than the
archive, against 0.93x for Clear alone. Two readings fit equally well:

- the archive is carrying vegetation texture it should not, which is the
  premise the whole cleanup rests on; or
- the pair has begun eroding real relief, and the smoothing is loss.

Nothing in an archive-referenced measurement can separate them, because the
archive is the thing in question. Note the direction of the risk: every
internal metric available here rewards a smoother surface, so a measurement
that says "better" is not evidence against the second reading.

**How it gets settled.** A regional block flown over Ka'Kabish and Cocochan in
northern Belize on 17 May 2022 carries three overlapping strips at 12.9 to 21
returns per m2, so the same ground was observed independently three times at
different scan angles. An engineered plaza floor under closed canopy holds the
substrate constant while the understory varies across it, which lets filter
error be scored against vegetation height rather than against another product.

The test surface is at 316450.8 E, 1970610.0 N (UTM 16N), 148 m from the
Ka'Kabish centre: 0.21% slope, 18.7 m canopy, 77% of returns above 2 m, 4.8% in
the 0.15-0.35 m band, three strips at a 57 degree spread. It pairs with a
lower-understory run on the same platform complex at x = 316509, y 1970589 to
1970631 -- three contiguous cells, 0.93 to 1.87% slope, 1.3 to 1.7% band, the
same three passes at nearly the same incidence.

**What the block has already settled.** Across 41 flat cells (slope under 2%,
three strips, band fraction 3.1 to 13.2%), each filter's surface was measured
against the 1st percentile of all returns in its own cell -- a datum that
references no classification and no commercial product:

| filter | surface rise per 10% band | r | t |
|---|---|---|---|
| CSF | +0.553 m | +0.49 | +3.55 |
| SMRF | +0.456 m | +0.42 | +2.91 |
| PMF | -0.034 m | -0.03 | -0.21 |

Over ground with no reason to be higher where the scrub is thicker, SMRF and CSF
rise with understory and PMF does not. Recovered ground density says the same
independently: SMRF +7.23 and CSF +7.73 returns per m2 per 10% band, PMF -0.06.
So SMRF and CSF are following the understory, and the ambiguity in the
filter-to-filter comparison -- which could equally have been PMF eroding terrain
-- is resolved in that direction.

That also reconciles the two datasets. On G-LiHT, scored against the TerraScan
archive, PMF looked worst: biased low on every tile at 0.34 to 0.75x the
reference's texture. But the archive shares the understory acceptance now
measured in SMRF and CSF, so part of PMF's negative bias is it correctly
declining vegetation the archive kept. Not all of it -- 0.34x is too much
smoothing to explain that way -- so PMF plausibly over-erodes relief *and*
rejects understory properly.

**What is still open.** How far above the soil any of them sits; whether PMF is
right or merely anchored to the minimum surface, since it holds a near-constant
0.20 m above the lowest returns; and whether the l8s431 complexity drop is
vegetation texture in the archive or lost relief. Those need surveyed control,
which is a follow-up rather than a dependency.

Measurements: `benchmarks/results/threshold_clear/threshold_clear.json` (12
arms), `benchmarks/results/nr_block/cells_3strip.json.gz` (12,506 cells with
three-strip coverage, trimmed from 63,253 scored; denominators in
`cells_summary.json`), `band_regression.json` (41 cells x 3 filters),
`band_offset.json` (the table above). Rebuild with `benchmarks/scan_block.py`,
`band_regression.py` and `band_offset.py`.

Point clouds by Alec McLellan; published here with his permission.
