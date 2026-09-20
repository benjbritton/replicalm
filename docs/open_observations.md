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
