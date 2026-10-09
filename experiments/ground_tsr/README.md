# Reproducing the `ground_tsr` experimental baseline

These scripts regenerate every number quoted in `docs/ground_tsr.md` and in the
commit that preserves this state. They are kept as they were run, not tidied
into a library: each one hardcodes the absolute paths it used, so the record of
which data produced which figure survives.

`ground_tsr` is experimental. It is integrated and functional but does not meet
the 0.052 m benchmark, and on Centro GLAS l0s441 it converges to within
0.0148 m RMSE of the SMRF baseline. Nothing here argues for replacing SMRF.

## Data this depends on

| role | path |
|---|---|
| raw flightline, l0s441 | `D:\GLiHT_LAS_orig\Yuc_Centro\Centro_Glas\Centro_Glas\AMIGACarb_Yuc_Centro_GLAS_Apr2013_l0s441.las` |
| delivered ground, l0s441 | `D:\_LAZ\Yuc_Centro\Centro_Glas\Centro_Glas\Centro_GLAS_l0s441_GND_v1.laz` |
| delivered DTM, l0s441 | `D:\_Archive_EdgeFixed\Yuc_Centro\Centro_Glas\Centro_Glas\Centro_GLAS_l0s441_DEM_0p5m_v1.tif` |
| raw flightline, South 395 | `D:\GLiHT_LAS_orig\Yuc_South\South_Glas\South_Glas\AMIGACarb_Yuc_South_GLAS_Apr2013_l0s395.las` |
| delivered DTM, South 395 | `D:\_Archive_EdgeFixed\Yuc_South\South_Glas\South_Glas\South_GLAS_l0s395_DEM_0p5m_v1.tif` |

Outputs land under `C:\Replicalm\test100526\`, which is gitignored: it runs to
5.6 GB of derived rasters and point sets and is local test output, not
published data.

Every script runs in the `replicalm` conda environment. The bare interpreter in
that environment fails to import scipy (`0xC06D007F`), so invoke them as
`conda run -n replicalm --no-capture-output python <script>`.

## Run order

The target and the benchmark have to exist before anything can be scored.

1. **`gnd_classes.py`** — what the delivery actually contains. Establishes that
   `Centro_GLAS_l0s441_GND_v1.laz` holds two classes, 2 and 8, and that every
   attribute except X, Y, Z and Classification is zeroed.
2. **`gnd_which_classes.py`** — which classes the delivered surface was fitted
   through. Class 2 alone sits a median 16 mm below it and class 8 12 mm above;
   their union is centred to 3 mm. This is why the target is the union.
3. **`c2_vs_c8.py`** — whether the 2/8 split could be reproduced at all. It
   could not: no attribute distinguishes them and at 5 m scale 99.5 % of
   occupied cells hold both.
4. **`ts_ceiling.py`** — builds `dem_ts_c2c8.tif`, the delivered union through
   the Replicalm biharmonic engine. This is the reference surface everything is
   scored against, and its 0.0524 m against the archive is the benchmark.

Then the classifier itself.

5. **`smoke_tsr.py`** — config defaults and the density-scaling table. Confirms
   `ground_engine` still defaults to `smrf` and that 20 returns/m² reproduces
   the validated 1.0 m working cell.
6. **`tsr_calibrate.py`** — first calibration pass. Found three faults: the
   gate floor never bound, precision rather than recall was binding, and the
   gate could only add to the confident set.
7. **`tsr_calibrate2.py`** — second pass, which chose the defaults now in
   `GroundTSRParams`. Eight development and six held-out blocks drawn under a
   fixed seed. Pruning was swept here and rejected.
8. **`tsr_fulltile.py`** — full tile with density adaptation on: 0.0817 m.
9. **`tsr_fulltile2.py`** — full tile with the cell pinned, emitting the
   confident, candidate and adjudicated sets from one pass: 0.0787 m. This is
   the run that separated the two causes. The density rule cost the 4 %; the
   gate moved RMSE by 0.0001 m.

Then the inspection packages and the generalisation test.

10. **`g1_package.py`** — G1 panels for the forensic mound and depression, with
    the stage-4 surface rendered as the gate actually samples it. Shows stage 4
    wrong by 0.195 m RMSE on the mound against a gate 0.05–0.30 m wide.
11. **`g1_fulltile.py`** — full-tile G1 at native 0.5 m for the three surfaces
    and three difference rasters, on one shared stretch.
12. **`south395_panels.py`** — South GLAS l0s395, a tile no parameter was
    fitted on. Both classifiers land at 0.055 m against the delivered surface
    and 0.0057 m apart, which is what establishes the adjudication stage as
    inert in general rather than only on its calibration tile.

## Supporting measurements

- **`rejection_union.py`** — SMRF and MCC rejection of the delivered union by
  landform, tile-wide and on the two 600 m cores. Every rejection figure
  predating this was measured against class 2 alone, which is under half the
  delivered ground population.
- **`rescore_common.py`** — rescores the surfaces on the cells all of them
  filled. Necessary because the TIN surface left 31 % of block cells outside a
  convex hull, and a surface scored only where it succeeded is scored on an
  easier subset.
- **`mcc_fulltile.py`**, **`ts_tin_surface.py`** — produce `dem_mcc.tif` and
  `dem_tin_c2c8.tif`, which `rescore_common.py` reads. Kept because the
  documentation cites their numbers, not because either is a candidate.

## Known open issue

Stage 4 builds its reference surface as a tenth quantile over 2 m cells,
smoothed at 3 m. On the forensic mound that surface is wrong by 0.195 m RMSE
with a 1.51 m worst case, and on the depression 0.139 m with 1.29 m, against a
gate that is roughly 0.05–0.30 m wide. A gate cannot adjudicate against a
reference that is wrong by several times its own width, and the consistent
negative bias (−0.103 m and −0.082 m) is what a low quantile does on a convex
form. `g1_package.py` renders this directly.

The density rule in `GroundTSRParams.scale_for_density` is the second open
issue. It was measured to cost 4 % at full-tile scale and has not been
replaced; the full-tile runs pin the cell with `adapt_density=False` instead.
