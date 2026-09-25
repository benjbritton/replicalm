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

## 10. The cleanup does not erode relief; the figure that said it did was wrong

**Status:** closed on the measurement that raised it. The collateral cost of the
cleanup remains open and is quantified below.

**What raised it.** A 2 x 2 factorial -- SMRF elevation threshold 0.50 m against
0.25 m, crossed with the cleanup on and off -- reported that the tightened pair
returned 0.78x the archetype's terrain texture on l8s431, the tile with most
relief. Either the archetype was carrying vegetation texture, which is the
premise the cleanup rests on, or the pair had begun removing real relief.

**Why that figure was wrong, twice over.** The factorial built its passes as
`GroundPass(algorithm="smrf", slope=..., threshold_m=...)`, naming two fields and
inheriting that class's defaults for the rest -- a 1.0 m working cell where the
locked pass uses 0.50 m. And it ran through a solve that fell back to inverse
distance on 89% of l0s444's cells, because the delivered clouds carry duplicate
locations that make the kriging system singular (see observation 11).

**Rerun at the locked configuration, with de-duplication and the covariance
solve:**

| tile | arm | removed | RMSE | bias | texture |
|---|---|---|---|---|---|
| l0s395 | t0.50 raw | -- | 0.0249 | +0.0014 | 1.09x |
| | t0.50 + cleanup | 0.85% | 0.0190 | -0.0006 | 0.91x |
| | t0.25 + cleanup | 0.68% | 0.0186 | -0.0007 | 0.90x |
| l8s431 | t0.50 raw | -- | 0.1002 | +0.0447 | 1.65x |
| | t0.50 + cleanup | 11.62% | 0.0572 | +0.0131 | **1.03x** |
| | t0.25 + cleanup | 8.82% | 0.0486 | +0.0088 | 0.90x |
| l0s444 | t0.50 raw | -- | 0.0659 | +0.0298 | 1.34x |
| | t0.50 + cleanup | 15.14% | 0.0353 | +0.0046 | **1.03x** |
| | t0.25 + cleanup | 13.86% | 0.0351 | +0.0046 | 1.06x |

At the locked threshold the cleanup lands at 1.03x on both harder tiles -- within
three points of the archetype's own texture -- and the 0.78x that prompted the
question does not occur at any arm. The tightened threshold is what smooths past
the target, taking l8s431 to 0.90x, which is why 0.50 m stays.

**What is still open.** The cleanup removes 11.6% and 15.1% of ground returns on
the two harder tiles where the agreement ceiling -- the returns TerraScan itself
rejected -- is 3.88%. It reaches the right surface by removing four times as much
as it needs to, and the excess softens platform edges. Closing that gap needs a
discriminator that separates clutter from small real features; the ones tried are
in observation 7, and the Lamanai three-strip work in observation 12 is the
attempt to derive one from physical invariance rather than from the archetype's
labels.

Measurements: `benchmarks/results/threshold_clear/threshold_clear.json`, rebuilt
by `benchmarks/threshold_clear_factorial.py`.

---

## 11. Duplicate locations in the delivered clouds made the kriging system singular

**Status:** closed. Fixed in `kriging.py`; the fix is verified against the
previous implementation on three tiles.

**What it is.** Some delivered G-LiHT tiles carry exact repeated returns. On the
l0s444 clip, 13.05% of all points and 13.49% of ground points are duplicates on
XYZ -- same return number, same GPS time, same source, 99.9% of duplicate groups
identical in elevation. The l0s395 clip has none. The repetition is in the
archive, not in the processing: the ratio going into our classifier and coming
out of it is the same to two decimals.

**Why it mattered.** Two returns at one location give the kriging matrix two
identical rows, so it is exactly singular. No nugget, no ridge and no neighbour
count can invert it, and every cell whose neighbourhood held such a pair
abandoned the solve for the inverse-distance fallback. That put 89% of l0s444's
cells, and 90 to 99% of the first production tiles at 64 neighbours, onto a path
that was still being described as ordinary kriging.

**Two wrong diagnoses, recorded so they are not re-derived.** The first was that
a short correlation range caused it: five of six production tiles had ranges of
14.38 m and fell back anyway. The second was that a nugget would condition the
system: a floor from 0 to 10% of the sill changed l0s444's result by nothing at
all to four decimals, because the fitted nugget was already 11.1% of sill and
the floor never bound -- and because with the system built from semivariances a
nugget lands off-diagonal, where it makes conditioning worse.

**The fix, in three parts.** De-duplicate on XY before the neighbour tree is
built, keeping the lowest elevation where one location holds several. Build the
system from covariances, `C(h) = (nugget + sill) - gamma(h)`, which is the same
estimator but puts the total sill on the diagonal, exceeding the off-diagonals by
exactly the nugget -- the dominance that conditions it, and the reason a nugget
belongs there. Add a ridge of 1e-6 of the sill for anything that survives.

**Verification.** l0s395, no duplicates: same cells filled, no cell differing by
more than a centimetre. l8s431, four duplicates: 62 cells of 625,962 move, and
RMSE, bias and texture are unchanged. l0s444, 213,937 duplicates: fallback falls
from 89.0% to zero and 37% of cells move, being the cells that were previously
weighted averages.

**What this means for the archetype.** The G-LiHT gridding ran through Surfer
with no `DupMethod` and no variogram specified -- `GLiHT_dem.py` sets only the
algorithm, cell size, sectors and the 20 m / 64 / 1 search -- so Surfer applied
its defaults for both and completed 1,368 grids without failing. How it coped is
not determinable from the g1 sources, and should not be guessed at.

---

## 12. Which filter is deceived by understory, measured without ground truth

**Status:** open, and the Optimization work rests on it.

Across 41 flat cells in the New River block (slope under 2%, three overlapping
strips, band fraction 3.1 to 13.2%), each filter's surface was measured against
the 1st percentile of all returns in its own cell -- a datum referencing no
classification and no commercial product:

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
-- resolves in that direction.

That also reconciles two datasets that appeared to disagree. On G-LiHT, scored
against the archive, PMF looked worst: biased low on every tile at 0.34 to 0.75x
its texture. But the archive shares the understory acceptance now measured in
SMRF and CSF, so part of PMF's negative bias is it correctly declining vegetation
the archive kept. Not all of it -- 0.34x is too much smoothing to explain that
way -- so PMF plausibly over-erodes relief *and* rejects understory properly.

**The surfaces this is built on.** Test surface at 316450.8 E, 1970610.0 N (UTM
16N), 148 m from the Ka'Kabish centre: 0.21% slope, 18.7 m canopy, 77% of returns
above 2 m, 4.8% in the 0.15-0.35 m band, three strips at a 57 degree spread. Its
low-understory pair sits at x = 316509, y 1970589 to 1970631 on the same platform
complex, 0.93 to 1.87% slope, 1.3 to 1.7% band, the same three passes. The
Cocochan plaza at 320630.1 / 1966971.9 is on strips 1, 2, 3, with its own bare
control 120 m away at 320750.0 / 1966969.9.

**Two known problems with the band.** It is measured 0.15 to 0.35 m above a plane
fitted across a 20 m cell, which is a different datum from the cleanup's 0.20 m
above a 0.75 m patch floor and from the specks' 5.8 cm above a 3 m median. The
three are not commensurable, and the band's bounds were chosen to bracket the
cleanup's threshold rather than the phenomenon. The lower bound should be derived
from the bare control -- the distribution of return heights above a fitted plane
where nothing is standing -- rather than asserted.

**What is still open.** How far above the soil any filter sits; whether PMF is
right or merely anchored to the minimum surface, holding a near-constant 0.20 m
above the lowest returns. Those need surveyed control.

Point clouds by Alec McLellan; published with his permission.
