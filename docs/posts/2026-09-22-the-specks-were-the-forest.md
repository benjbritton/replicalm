# The specks were the forest

### An archaeologist notices freckles on a plaza floor, and the fix turns out to be upstream of everything I was looking at

Benjamin Jay Britton, 22 September 2026

---

Replicalm reproduces the NCALM bare-earth lidar workflow with open tools — PDAL,
GDAL, NumPy, SciPy — so that results depending on TerraScan, ArcGIS Pro and
Golden Surfer can be reproduced by people who do not have those licenses. An
[earlier post](2026-09-20-calibrating-on-the-wrong-ground.md) covered the part
where the test set turned out to be selecting against the terrain that mattered.

This one starts after that was fixed, with a surface that matched the reference
to 0.0601 m RMSE and a visualization that looked right.

Then a closer look at the Pixoyal group: small dark specks scattered across
otherwise flat ground, absent from the reference. Freckles. Not obviously wrong,
easy to overlook — and exactly the sort of thing that makes an archaeologist
hesitate over whether a faint mark is a low platform or an artifact of
processing.

![The Pixoyal group: original workflow, baseline, Clear](../figures/fig10_clear_vs_baseline.png)

*Left: the original commercial workflow. Center: where this post starts — note
the plaza floors. Right: where it ends.*

## What they were made of

Measuring before theorizing, on a 500 × 500 cell window:

```
median cluster size                        1 cell
our surface above the reference            78.9% of them, median 2.6 cm
local roughness there                      0.0768 m vs reference 0.0358 m
ground returns within 0.6 m of each        median 5
their height above the 3 m local median    +5.8 cm, 76.2% positive
their scatter among themselves             8.8 cm
```

Two things stand out. The bias is **one-sided** — 76.2% positive — and ranging
noise is symmetric, so this is not sensor noise. And the returns defining each
speck **scatter by more than the feature they define**: 8.8 cm of disagreement
describing a 5.8 cm bump. Signal-to-noise of about 0.66.

So: real returns, sitting genuinely above their surroundings, disagreeing with
each other about by how much.

## Two wrong answers first

**The interpolator.** Our variogram fits returned a nugget of exactly zero on 34
of 60 fits across ten tiles, which makes ordinary kriging an exact interpolator —
it honors every observation, so anything in the ground class goes straight into
the raster. Flooring the nugget smooths it away, and at 4% of sill the surface
matched the reference's flat-ground texture precisely.

It was the wrong fix for a reason worth stating: **measurement noise alone
justifies a nugget of about 0.6%.** Matching the reference's flat-ground texture
took seven times that. What produced that texture is not documented and we have
not established it — but whatever it was, tuning our surface to match it would
have meant blurring real micro-relief along with the clutter, rather than
removing the clutter.

**A measurement trap, while we are here.** Speck counts across the nugget arms
came out non-monotonic — 5% and 10% scoring far better than their neighbors.
That was not the surfaces. G1 normalizes each raster by its own extremes, and
one outlier sky-view-factor cell in the 5% arm restretched the entire image.
Two arms had been scored against a different tonal mapping than the rest before
anyone noticed. Per-raster normalization means two tiles rendered separately are
not tonally comparable, which has implications well beyond this sweep.

## What they actually are

The decisive evidence came from the flight geometry. Where two source chunks
both contribute returns to a speck, they **disagree by 0.47 m** against 0.11 m
on control ground. Mean scan angle is identical to control — 11.5° against 11.2°
— so it is not an incidence-angle artifact. And the returns are 11% darker than
control ground.

Something at those locations stops pulses at variable heights, independently of
viewing angle, and reflects less than soil or limestone. Limestone returns
consistently and brightly. **Low vegetation — scrub, brush, root mass — returns
at whatever height it happens to occupy for a given pulse.**

That is the reading the numbers pointed to, and it survived every
test aimed at it. It also explains the single-chunk pattern: one pass hits a
dense patch of brush at an angle, the overlapping pass finds a gap and reaches
the soil, and the two disagree by half a meter about where the ground is.

## The labels were in the data all along

The delivered G-LiHT clouds carry TerraScan's own classification. Matching our
SMRF ground returns to them by position gives a direct label for every point:

```
our ground returns, by what NCALM called them
                          in specks    control      all
NCALM ground (2)              50.4%      66.2%     65.2%
NCALM near-ground (8)         37.4%      33.0%     29.7%
NCALM unclassified (1)        12.2%       0.8%      5.1%
```

**A 15× enrichment.** At speck cells, 12.2% of what we call ground, the original
workflow explicitly threw out.

That makes the ceiling computable. Rebuild the surface with exactly those points
removed — 3.88% of returns — and:

```
                        before    after
flat-ground roughness    1.26x    1.00x the reference
specks                   3,901      253
RMSE                    0.0825   0.0140
```

Wholly a classification residue. Not interpolation, not the nugget, not the
search radius, not the slope handling. The freckles were vegetation we were
calling ground.

## Finding them without TerraScan

Three rules were tried. Two failed in a way worth recording.

Demoting returns that stand above their neighborhood's floor removes 72% of
synthetic brush — and 100% of a synthetic rock of the same footprint. A feature
narrower than the search radius has neighbors on the flat ground around it, so
the floor stays low and the whole feature reads as high. Adding a
slope-compensated spread test did no better: 54% of the brush, 92% of the rock,
and 0% of brush on a 20° slope.

Using the flight-chunk disagreement directly — the signal that produced the
diagnosis — demotes 3.1% of ground returns, touches 3.4% of speck cells and
**10.2% of steep cells.** It damages architecture three times harder than it
cleans flat ground, because on a slope two chunks legitimately sample different
parts of the same real surface.

What worked is the rule the synthetic test had condemned. Against TerraScan's
labels, height above the local floor separates at d′ = 1.82 — rejected points sit
0.42 m above their patch floor against 0.08 m for kept points. The synthetic test
had used a 0.4 m feature against a 0.6 m radius, the pathological case. On real
data at a 0.75 m patch it is the best discriminator available.

**Clear:** demote a ground return standing more than 0.20 m above the tenth
percentile of ground returns within 0.75 m.

```
                                       baseline     Clear     oracle
recall of what TerraScan rejected             —     92.5%      100%
ground returns removed                        —    17.75%     3.88%
flat-ground roughness vs reference        1.26x     0.99x     1.00x
specks                                    3,901       426       253
error above 20 degrees                    1.86%     0.29%         —
```

Validated unchanged on three windows it was not fitted to. It improves steep
ground by 6×, which was the outcome I most expected to go the other way.

## Measuring without the answer key

Every number above is referenced to NCALM's output. That is the right frame for
a replication and the wrong one for asking whether a surface is *accurate* — if
ours were better, all of those figures would get worse.

So: withhold ground returns, build the surface without them, and measure the
residual where those returns actually are. No reference DEM involved. And count
returns lying *below* the modeled surface, since a pulse can be stopped early by
vegetation but cannot arrive below the soil.

```
                                        baseline     Clear
predicts held-out ground returns (MAE)   0.0492 m   0.0429 m
returns left below the surface             8.11%      4.80%
```

13% better prediction of real measurements, and the returns left stranded
beneath the surface nearly halved. Clear is not merely more archive-like; it is
a better account of the ground.

![Pixoyal: the original workflow and Clear](../figures/fig9_pixoyal_ncalm_vs_replicalm.png)

*Left: the original workflow. Right: Replicalm with Clear.*

## What it costs, and what is still open

Clear removes 17.75% of ground returns where the oracle removes 3.88%. The
difference is collateral — genuine returns, some of them on platform edges,
which is why the architecture in the right-hand panel is very slightly softer
than the reference. At 4.2 returns per square meter, losing 17.75% still leaves
3.4, which the density work says is ample; but the 14-point gap is the honest
measure of not having TerraScan's labels.

Closing it needs a discriminator that separates clutter from small real features.
At four returns per square meter, a 0.4 m bush and a 0.4 m rock are each
described by about five points with the same geometry. The statistical signal is
real; it is not available per point. One avenue remains untried: every test here
looked only at returns already classified as ground, and the full return stack —
what else that pulse hit, and at what height — is discarded before any rule sees
it.

## The general point

The specks looked like an interpolation artifact. They were measured as one,
theorized about as one, and two fixes were built for them as one. They were
vegetation that the ground filter had accepted, and no amount of work downstream
of that could remove them — only blur them. The reference surface is smoother
than measurement noise alone would account for; how it got that way is not
something these measurements can settle.

The thing that broke it open was already in the data: the delivered clouds carry
the original classification, and nobody had thought to compare against it. Not a
new measurement — a new look at one already sitting in the file.

---

*Full record of what remains unexplained, including three items from this work:
[docs/open_observations.md](../open_observations.md)*

— Benjamin Jay Britton, 2026
