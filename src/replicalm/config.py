"""The NCALM method as parameters, and the honest translation into what we have.

THE METHOD BEING REPLICATED
---------------------------
Estrada-Belli et al. 2025, supplementary materials, quoting the workflow given in
Estrada-Belli et al. 2023. Ten steps, reproduced verbatim in docs/ncalm_method.md.
The classification is TerraScan's "Classify Ground", run twice:

    pass 1   window 70.0 m   max terrain angle 88.0   iteration angle  9.0   iteration distance 3.0 m
    pass 2   window 70.0 m   max terrain angle 88.0   iteration angle 12.0   iteration distance 3.0 m

WHY THIS IS A TRANSLATION AND NOT A PORT
----------------------------------------
Those four numbers are the parameter surface of TerraScan's progressive TIN
densification, after Axelsson 2000. No tool available here implements that
algorithm with that parameter surface:

  ArcGIS Pro   Classify LAS Ground offers three presets -- Standard,
               Conservative, Aggressive -- and exposes none of the four.
  LAStools     lasground is the closest commercial equivalent and is not
               installed here; its parameters (step, spike, offset, bulge) are
               also not the same four.
  PDAL         filters.smrf implements Pingel et al. 2013, a simple morphological
               filter; filters.pmf implements Zhang et al. 2003; filters.csf
               implements Zhang et al. 2016 cloth simulation. All three are
               parameterised, none of them by terrain angle and iteration angle.

So the classification step cannot be replicated exactly, and any claim that it
has been is false. What CAN be done, and is what this package does, is to
reproduce the SHAPE of the method -- reset, two-stage ground extraction, height
above ground, class 8 near-ground band, ground-only export, kriged interpolation
at a 20 m search radius, mosaicked at 1 m -- with a classification stage whose
parameters are stated, tunable, and calibratable against reference output.

Everything downstream of classification replicates exactly. The tiling, the
buffers, the height thresholds, the class assignments, the kriging radius and
the output grid are all specified numerically in the source and are all
reproducible.

PARAMETER TRANSLATION
---------------------
The mapping below is a starting point, not an equivalence. It is derived from
what each parameter controls rather than from any published correspondence:

  TerraScan window 70 m
      the largest structure the filter will step over. SMRF's `window` is the
      maximum morphological window and does the same job, but it does not
      carry across directly: measured on synthetic flat-topped boxes, SMRF's
      progressive opening stops changing its labels above 32 m at a 0.5 m
      working cell, so 70 and 32 produce identical ground and 70 costs two and
      a half times the time. The configured pass therefore uses 32. Note the
      sign of the parameter -- `window` is the largest object the opening will
      STRIP, not the largest it will preserve, so reducing it below the
      saturation point is conservative rather than aggressive.
  TerraScan iteration angle 9 deg / 12 deg
      how steeply the TIN may be pushed toward a candidate point. SMRF has no
      angular term; its `slope` is a rise-over-run tolerance. tan(9 deg) = 0.158
      and tan(12 deg) = 0.213, which is where the default slopes come from.
  TerraScan iteration distance 3.0 m
      how far below the TIN facet a point may sit and still be taken as ground.
      SMRF's `threshold` is an elevation tolerance and is the nearest analogue.
  TerraScan max terrain angle 88 deg
      a permissiveness ceiling, effectively "allow almost any slope". SMRF has
      no equivalent; the ceiling is implicit in `slope`.

The two-pass structure is preserved: a coarse pass establishes ground, a second
pass refines it with a steeper tolerance, exactly as the source describes.

CALIBRATION IS THE POINT
------------------------
Because the translation is inexact, the defaults here should not be trusted on
new terrain. `replicalm calibrate` scores a parameter sweep against a reference
DEM produced by the original workflow, and reports which settings come closest.
Use it before processing a survey, not after.
"""
from dataclasses import dataclass, field, asdict, replace
import json
import math


# --- the source method, recorded so it cannot drift -------------------------

NCALM_TERRASCAN = {
    "pass_1": {"window_m": 70.0, "max_terrain_angle_deg": 88.0,
               "iteration_angle_deg": 9.0, "iteration_distance_m": 3.0},
    "pass_2": {"window_m": 70.0, "max_terrain_angle_deg": 88.0,
               "iteration_angle_deg": 12.0, "iteration_distance_m": 3.0},
    "tile_km": 1.0,
    "macro_tile_km": 20.0,
    "buffer_m": 10.0,
    "hag_low_cut_m": -0.5,
    "hag_high_cut_m": 600.0,
    "near_ground_band_m": (-0.2, 0.2),
    "near_ground_class": 8,
    "kriging_search_radius_m": 20.0,
    "dem_cell_m": 1.0,
    "output_las_version": "1.2",
    "source": ("Estrada-Belli et al. 2025 supplementary materials, quoting "
               "Estrada-Belli et al. 2023, J. Archaeological Science 157:105835"),
}


def _slope_from_angle(deg):
    """Angular iteration tolerance -> rise over run, for SMRF's slope term."""
    return round(math.tan(math.radians(deg)), 4)


@dataclass
class GroundPass:
    """One ground-classification pass.

    `algorithm` selects which implementation runs. The parameters that follow
    are the union of what those implementations accept; each reads only its own.
    """
    algorithm: str = "smrf"              # smrf | pmf | csf
    window_m: float = 70.0               # smrf, pmf: max morphological window
    slope: float = 0.158                 # smrf, pmf: rise over run
    threshold_m: float = 3.0             # smrf: elevation tolerance
    scalar: float = 1.25                 # smrf: slope scaling with window
    cell_m: float = 1.0                  # smrf: working grid
    csf_resolution_m: float = 1.0        # csf: cloth grid
    csf_rigidness: int = 2               # csf: 1 steep, 2 relief, 3 flat
    csf_threshold_m: float = 0.5         # csf: cloth-to-point distance
    reuse_ground: bool = False           # pass 2 refines pass 1

    # provenance: which TerraScan setting this pass stands in for
    stands_in_for: str = ""


@dataclass
class ReplicalmConfig:
    """A complete run. Serialise it beside the output so a result is explainable."""

    # Tiling is a memory device, not a method parameter. The source works in
    # 1 km tiles because a whole transect will not fit in RAM on modest
    # hardware; the output does not depend on it, and the G-LiHT reprocessing
    # did not use it. None means process the whole extent.
    #
    # krige_grid now chunks its neighbour query internally, so interpolation is
    # bounded regardless of grid size. Classification is the remaining reason to
    # tile: PDAL's SMRF and CSF hold the cloud and their own working grids in
    # memory, so a very large cloud on a small machine may still need it. Set
    # tile_km to 1.0 for the source's behaviour.
    tile_km: float = None
    buffer_m: float = NCALM_TERRASCAN["buffer_m"]
    chunk_cells: int = 1000000           # cells per neighbour query batch

    # Why this run departed from its method's defaults, in the operator's own
    # words. The defaults are what the source specifies and what every
    # measurement in docs/ was taken at, so a departure is a claim about the
    # survey -- that its correlation structure, vegetation or relief differs in
    # some way that the default does not suit. Recording the claim beside the
    # product is what makes the departure auditable rather than merely visible.
    # Empty on a default run.
    note: str = ""

    # noise removal, before classification
    #
    # The source method removes outliers at step 6, by height above the ground
    # model, which is after classification. TerraScan's Classify Ground also
    # rejects low points internally, as part of the algorithm. SMRF does not:
    # handed a return three metres below true ground it takes that as the local
    # minimum and builds terrain down to it, leaving a pit.
    #
    # Omitting these stages was visible in the calibration as excess terrain
    # complexity. Fifteen interpolation settings all produced 1.6 to 2.8 times
    # the reference's texture while matching its elevations to a few
    # centimetres -- the signature of spurious pits and spikes rather than of
    # detail. The densest tile, at 10 points per square metre, was the worst.
    # ELM is OFF because it is inert on this data, not because it is unwanted.
    # Measured on four windows -- one flat, three with 11% to 33% of cells above
    # twenty degrees -- it removed 0.00% of ground points at thresholds of 0.5,
    # 1.0 and 2.0 m, and every arm that ran it produced output identical to the
    # arm that did not. G-LiHT returns evidently carry no extended local minima
    # of the kind it looks for. Turn it on for data that does.
    remove_low_noise: bool = False       # filters.elm, extended local minimum
    elm_cell_m: float = 10.0
    elm_threshold_m: float = 1.0

    # The statistical outlier filter removes 0.05% to 0.36% of ground points and
    # its effect is tile-dependent rather than uniform: on l0s395 it halved RMSE
    # (0.175 -> 0.091) and cut error above thirty degrees from 5.22% to 2.89%;
    # on l0s417 it made both worse (0.169 -> 0.220, 8.20% -> 8.43%). It is left
    # on because the aggregate favours it, but that aggregate rests on one tile
    # and this is a setting worth sweeping per survey rather than trusting.
    remove_outliers: bool = True         # filters.outlier, statistical
    outlier_neighbours: int = 8
    outlier_multiplier: float = 2.5
    noise_class: int = 7                 # LAS 1.2 low point / noise

    # classification, steps 3 to 5
    # ONE PASS, not the source's two.
    #
    # This was measured twice and the first measurement was wrong. A window
    # comparison appeared to show two passes beating one by 0.083 m RMSE against
    # 0.161 m, with terracing 0.25% against 2.44% -- but the two-pass arm had the
    # noise filters on and the one-pass arm had them off, so what it measured was
    # the filters. Re-run across seven tiles with the filters held on, the two
    # arms are indistinguishable: steep support 98.1% either way, error above 30
    # degrees 10.34% against 10.18%, terracing 2.87% against 2.90%. The second
    # pass roughly doubles classification time and returns nothing.
    #
    # Two SMRF passes were never TerraScan's refinement anyway -- SMRF has no
    # refinement mode, so the second pass simply reclassifies the whole cloud
    # with a steeper tolerance, superseding the first. That is recorded in
    # classify.py as a divergence, and the measurement now says the divergence
    # costs nothing to make complete.
    #
    # `threshold_m` is 0.25, not the source's 3.0. `cell_m` is 1.0 as of
    # 2026-10-07, raised from 0.5, and the two measurements behind that
    # disagree with each other -- both are recorded because the conflict is
    # not resolved.
    #
    # AGAINST 1.0 m, measured earlier: left coarser than the output cell the
    # working grid quantizes steep faces into steps, with terracing of 5.06%
    # at 1 m, 2.44% at 0.5 m and 0.68% at 0.25 m.
    #
    # FOR 1.0 m, measured on the whole of l0s441, 55.2 million raw returns
    # classified twice with nothing else changed, scored over 11.9 million
    # cells against the TerraScan/Surfer reference:
    #
    #                 class-2      RMSE     MAE    p99   roughness    nugget
    #     cell 0.50  17893142    0.1027  0.0459   21.0       1.055  0.001497
    #     cell 1.00  17702632    0.0789  0.0420   20.5       1.015  0.001402
    #
    # RMSE falls 23.1%, MAE 8.3%, nugget 6.4%, and roughness moves from 1.055
    # to 1.015 against the reference's 1.000, for 1.1% fewer ground returns.
    # The same direction appeared first on eighteen 40 m windows, so this is a
    # replication rather than a single result.
    #
    # The two are not necessarily contradictory: terracing is a defect of
    # steep faces and the RMSE is a whole-tile average, so 1.0 m could be
    # better overall and worse on scarps. That has not been tested. Revert by
    # setting cell_m back to 0.5 here and in BASELINE.
    #
    # `scalar` is 0.00 and `window_m` is 32, both changed from the values this
    # file carried until 2026-10-04 (1.25 and 70). The three parameters were
    # measured together on two full tiles at opposite ends of the corpus --
    # l0s426 at 0.86 ground returns per square metre and 40 degrees of mean
    # slope, l0s395 at 4.24 and 8.5 -- against each tile's archive surface on
    # the archive's own lattice.
    #
    # `scalar` multiplies the slope-dependent part of SMRF's elevation
    # tolerance, so on a steep face the effective threshold is far above
    # whatever `threshold_m` names. On l0s426 the classification arrived 0.30 m
    # above the archive on flat ground and 0.73 m above it past 40 degrees
    # before any cleanup ran, and over half of what it called ground stood
    # more than half a metre above the lowest return in its own cell. Zeroing
    # the term took cells within 25 cm of the archive from 17.3% to 44.4%. It
    # was then run on l0s395, where the classification already sat at
    # +0.0002 m and a change could only cost: it improved every measure there
    # too, halving the bias past 40 degrees.
    #
    # `threshold_m` 0.25 was the single largest remaining gain: flat-ground
    # bias on l0s426 from +0.1062 to +0.0059 m and cells within 25 cm from
    # 58.5% to 64.6%, and on l0s395 a semivariance of 0.00667 against the
    # archive's own 0.00663 -- the closest texture match any configuration has
    # produced on that tile. It removes 22% of the ground returns 0.5 m kept.
    #
    # `window_m` 32 is a runtime change and not an accuracy one. On synthetic
    # flat-topped boxes the progressive opening's effect saturates at 32 m with
    # a 0.5 m working cell: from 32 to 128 m the labels are identical, and the
    # three real clips swept over 6 to 70 m returned identical ground counts.
    # It cut l0s395's classification from 4,692 s to 1,888 s. Below 32 m the
    # parameter is not inert -- it strips objects narrower than itself -- so
    # this is a ceiling to stop at, not a direction to continue in.
    #
    # `slope` stays at tan 9 degrees. Raising it to tan 25 recovered 21% more
    # ground returns on l0s426, exactly as a clip-level test predicted, and
    # cost 10.9 points of accuracy doing it: the recovered returns were canopy
    # that happened to be the only return in its cell, which is the error a
    # per-cell test structurally cannot see. Measured in isolation and
    # rejected.
    reset_classification: bool = True
    passes: list = field(default_factory=lambda: [
        GroundPass(slope=_slope_from_angle(9.0), threshold_m=0.25, cell_m=1.0,
                   scalar=0.0, window_m=32.0, reuse_ground=False,
                   stands_in_for="TerraScan pass 1: 70 m, 88 deg, 9 deg, 3 m; "
                                 "pass 2 dropped, see the note above"),
    ])

    # Which ground classifier runs.
    #
    # "smrf" is the locked baseline: the single morphological pass defined in
    # `passes` above, which is what every published Replicalm figure was built
    # with and what the installer ships.
    #
    # "tsr" is `replicalm.ground_tsr`, which composes a conservative opening,
    # a permissive curvature test and an adjudication stage between them. It
    # was built to reproduce an NCALM TerraScan/Surfer ground population from
    # raw returns and is calibrated against one; see that module's docstring
    # for the measured basis of every default and for its operational limits.
    # It is opt-in because the baseline is locked, not because it scores worse.
    #
    # `passes` is ignored under "tsr": the routine carries its own stage
    # parameters so that one object describes one classifier, rather than two
    # places having to agree.
    ground_engine: str = "smrf"          # smrf | tsr
    ground_tsr_overrides: dict = field(default_factory=dict)

    # height above ground, step 6
    hag_low_cut_m: float = NCALM_TERRASCAN["hag_low_cut_m"]
    hag_high_cut_m: float = NCALM_TERRASCAN["hag_high_cut_m"]
    near_ground_band_m: tuple = NCALM_TERRASCAN["near_ground_band_m"]
    near_ground_class: int = NCALM_TERRASCAN["near_ground_class"]

    # Residual low vegetation -- scrub, brush, root mass that stopped the pulse
    # above the soil -- survives the ground filter and appears in the surface as
    # a fine speckle on flat ground. Demoting returns that stand more than
    # `clean_height_m` above the tenth percentile within `clean_patch_m`
    # recovers 92.5% of what TerraScan rejected, takes flat-ground roughness
    # from 1.26x the reference to 0.99x, and cuts error above twenty degrees
    # from 1.86% to 0.29%. Validated on three windows it was not fitted to.
    #
    # It removes 17.75% of ground returns where TerraScan's own labels need only
    # 3.88%, so it takes genuine returns with the clutter and platform edges come
    # out slightly softer.
    #
    # THE THRESHOLD IS NOT UNIVERSAL. It was fitted against one tile's labels and
    # validated on three more from the same G-LiHT campaign. A different sensor,
    # flight height or vegetation regime may want a different figure. Set
    # clean_vegetation False for the unfiltered surface, which is the pure
    # translation of the source method and what the published comparison figures
    # are measured against.
    clean_vegetation: bool = True
    clean_patch_m: float = 0.75
    clean_height_m: float = 0.20
    # Measure the height against a plane fitted to the neighbourhood instead of
    # against the 10th percentile of its raw elevations. The percentile sits
    # downhill on a slope, so the published rule removes ground for being on an
    # incline: measured on l0s395 it moves the surface not at all below 9.5
    # degrees and by 12 cm beyond 40, with 99% of the movement downward. The
    # plane cancels the terrain and leaves what stands above it. Off by default
    # -- it is a departure from the source method, not a translation of it.
    clean_slope_aware: bool = False
    # Radius the plane is fitted over, when `clean_slope_aware` is set. Separate
    # from `clean_patch_m` because the two do different jobs: the threshold wants
    # a small neighbourhood, while the plane must sit on the ground and will be
    # tilted onto itself by a clump of scrub wider than the disc it is fitted
    # over. None means use `clean_patch_m`.
    clean_plane_patch_m: float = None
    clean_percentile: float = 10.0

    # export, step 7
    #
    # LAS 1.4 with point format 6, written as LAZ. The source specifies 1.2,
    # which stores its coordinate system as GeoTIFF keys rather than WKT -- fine
    # for UTM, lossy for anything less common -- and carries classification in
    # five bits against 1.4's eight. LAZ costs nothing here and is five times
    # smaller: 25.4 MB of LAS becomes 4.8 MB.
    #
    # `forward` is limited to scale and offset so that nothing else carries over
    # from the input header. Clouds delivered from TerraScan bring proprietary
    # Terrasolid records that inflate the file and mean nothing outside that
    # software; those are dropped, and the coordinate system is written
    # explicitly rather than forwarded.
    las_version: str = "1.4"
    las_point_format: int = 6
    las_compression: bool = True
    las_forward: str = "scale,offset"
    keep_classes: tuple = (2, 8)

    # interpolation, steps 8 and 9
    #
    # Every value below is meant to be adjustable by whoever runs this, and the
    # defaults are measurements rather than preferences. They were chosen on
    # South_GLAS_l0s395 and should be re-checked on unfamiliar terrain, which is
    # what `replicalm.calibrate` is for.
    # WHICH SURFACE ENGINE RUNS
    #
    # "kriging" estimates each cell as a weighted average of its neighbours
    # under second-order stationarity. Where ground returns arrive in clumps
    # separated by canopy voids the weights concentrate on the nearest clump and
    # the surface relaxes toward the local mean between them, which on sloped
    # ground reads as curvilinear grooves. No neighbourhood, nugget or cell size
    # removes that: the mechanism is in the estimator, not its parameters.
    #
    # "spline" is minimum curvature with tension, which has no local mean to
    # relax toward. Measured on Centro GLAS l0s441, interior cells, identical
    # accepted ground and grid, against a TerraScan/Surfer reference of 21.7
    # degrees at the 99th percentile of local slope: kriging 24.1, spline 20.9
    # at tension 0.03. See `replicalm.spline` for the tension response.
    #
    # Tension 0.03 is the operating point. The response flattens below it --
    # 0.02 measures 21.1 against 0.03's 20.9 -- and 0.03 carries the smallest
    # local excursions of any surface compared, the reference included, so it is
    # not over-smoothing to reach that figure.
    # "biharmonic" is the production engine as of 2026-10-10.
    #
    # The local spline that held this slot approximates: a smoothing term of
    # 1e-2 over 48 neighbours lets the surface sit off the observations, and
    # measured against the delivered archive on Centro GLAS l0s454 it retained
    # 0.505 of the reference's coherent detail below 1 m and 0.607 at 1.5 m.
    # The global minimum-curvature solve retains 0.914 and 0.972 on the same
    # tile, halves the error inside the boundary -- 0.0964 m to 0.0513 m at a
    # 25 m setback -- and runs in a third of the time.
    #
    # It is the Surfer half of the reference chain and it had been in the tree
    # since a304e8f without being reachable from here, so every release before
    # 1.4.0 shipped the local spline while the measurements quoted against the
    # archive came from this one.
    # Edge fixing, back ON as of 2026-10-10 with the method corrected.
    #
    # Two stages now, with separate causes. The footprint keeps cells within
    # the interpolator's reach of a ground return -- `search_radius_m`, since
    # a cell beyond it was extrapolated rather than interpolated -- closes
    # each area into a polygon with no interior rings, and fills any void
    # inside the survey with a continuous surface rather than a sentinel.
    # The across-track strip then removes the flight-pass stitching seam.
    #
    # It replaces a 4 m disc union that drew the boundary as overlapping
    # circles with detached islands, and a blocked-solver artefact that cut
    # sawtooth teeth at 150 m. Measured on Centro GLAS l0s454: footprint
    # +0.5% of the archive's area, raggedness 1.91 against 1.89, one polygon.
    edge_fix_enabled: bool = True

    # Overrides the footprint radius. None means use `search_radius_m`, which
    # is what the reference's own footprint corresponds to.
    footprint_radius_m: float = None

    # The along-track ends, trimmed as well as the sides. The sides carry the
    # flight-pass stitching seam; the ends carry the end-of-pass taper, where
    # the scan narrows and the solver has support on one side only. The
    # reference treats both at the same width -- the archive is inset exactly
    # 42 rows, 21.0 m, at north and south as well as east and west. Set to 0
    # to leave the ends alone.
    edge_strip_ends_m: float = 21.0

    interpolator: str = "biharmonic"     # kriging | spline | biharmonic
    biharmonic_core_m: float = 150.0     # solved domain, declared in
    biharmonic_buffer_m: float = 20.0    # metres so it is resolution
    biharmonic_pool_sigma_m: float = 0.25  # independent
    spline_tension: float = 0.03         # inverse metres; see replicalm.spline
    spline_neighbours: int = 48          # observations per local solve
    spline_smoothing: float = 1e-2       # Tikhonov, as a fraction of the kernel
    variogram_model: str = "spherical"   # spherical | exponential | linear
    variogram_fit_lag_m: float = 5.0     # lag out to which the model is fitted

    # The radius scales with point density, and is NOT derived from a variogram.
    #
    # Fitting variograms at windows of 5 to 80 m on ten tiles gave a range that
    # tracked the window every time, at 0.30 to 1.44 times its width, with a
    # sill that grew by 3x to 140x and never plateaued. There is no correlation
    # length in this field to derive a radius from, so `radius_from_variogram`
    # was returning roughly 1.4x `variogram_fit_lag_m` and calling it a
    # measurement. It is kept in kriging.py for terrain that does have a range,
    # and is not used here.
    #
    # A fixed radius appeared not to be the answer either, because it means
    # different things at different densities. At 3 points per square metre
    # `max_points` binds long before a 20 m radius does and 5 m and 20 m give
    # identical results; at 8 points per square metre the same 20 m radius
    # packed the neighbourhood tightly enough to lose matrix rank, putting 77 to
    # 88% of l0s444's cells onto the inverse-distance fallback.
    #
    # That measurement was taken at `max_points` 16, and the cause was the
    # point count, not the radius. Sixteen neighbours at 8 points per square
    # metre are all inside 0.79 m whatever the radius is set to -- nearly
    # coincident, so the covariance matrix loses rank. Sixty-four reach 1.58 m
    # at the same density, twice the span, and the conditioning improves
    # accordingly: moving to 64 inside a fixed 20 m took the fallback from
    # 1.18% of cells to 0.53% on l0s426 and from 0.23% to 0.10% on l0s395. The
    # radius was never the problem.
    #
    # So the neighbourhood is now the source's own: up to 64 points, minimum 1,
    # inside a fixed 20 m (Estrada-Belli et al. 2025 supplementary,
    # GLiHT_Methods_Materials section 5.2), which is also what the `baseline`
    # profile already used. Measured against each tile's archive surface on the
    # archive's lattice, with classification held constant:
    #
    #   l0s426, 1.08 returns per m2   semivariance at one cell 0.14876 -> 0.12997
    #   (sparse, steep)               against an archive of 0.10984, so 35% above
    #                                 it to 18%; coverage 86.6% -> 90.9%;
    #                                 elevation unchanged within 6 mm
    #   l0s395, 4.47 returns per m2   semivariance 0.00667 -> 0.00668 against an
    #   (dense, gentle)               archive of 0.00663, i.e. no smoothing where
    #                                 there was nothing to gain; three of four
    #                                 slope bands to exactly zero, the fourth
    #                                 +0.0128 -> +0.0082 m; coverage 99.1% ->
    #                                 99.6%
    #
    # Better on the sparse tile, neutral-to-better on the dense one, so this is
    # not density-conditional. It costs 7 to 9 times the interpolation -- 2,400
    # to 4,700 s per tile against 350 to 630 s -- which is a scheduling cost and
    # not a methodological one.
    #
    # NOT YET MEASURED: the densest tiles. l0s444 at 8.18 returns per m2 is
    # where the rank failure above was found, and 64 points reach only 1.58 m
    # there. The mechanism says conditioning should improve, and both measured
    # tiles agree, but 8 points per square metre is an extrapolation until a run
    # confirms it. `radius_for_density` is kept in place for that case and for
    # surveys whose density differs from this one.
    #
    # The source's documented 20 m is preserved in NCALM_TERRASCAN and is the
    # ceiling here, so nothing ever searches further than the method specifies.
    # Below 1.0 the grid is deliberately finer than the point spacing. That
    # buys rendering rather than measurement: slope and sky-view factor computed
    # without stair-stepping at cell boundaries, for a 2.0% gain in rasterised
    # fold residual and 2.3x the compute. The Deep profile sets it.
    cell_factor: float = 1.0

    # 64 POINTS INSIDE A FIXED 20 M, WHICH IS THE SOURCE'S OWN FIGURE
    #
    # A smaller neighbourhood puts spurious narrow ridges and grooves on slopes.
    # On Centro GLAS l0s441 the 99th percentile of local slope is 22.9 and 22.0
    # degrees in the TerraScan/Surfer rasters at 0.33 and 0.50 m, and 32.9
    # degrees in a Replicalm raster built with 16 points at a density-derived
    # radius of 3.73 m. Rebuilt from the same classified ground with 64 points
    # inside 20 m it is 23.9 degrees, and the inverse-distance fallback falls
    # from 0.65% of cells to 0.12%. In a plain hillshade the difference is a
    # network of thin dark incisions across every flank, which is what makes the
    # hillocks read as inverted, and they are gone at 20 m.
    #
    # The defect is NOT the vegetation cleanup: with the cleanup off the 99th
    # percentile is 33.2 degrees, marginally worse. It is not the G1 blend
    # either, because a plain hillshade shows it. It is the neighbourhood.
    #
    # HOW THIS WAS NEARLY MISSED, which is the part worth keeping
    #
    # On 2026-10-05 the default was set to 16 points at a density-derived radius
    # on the grounds that 64 inside 20 m "bought half a point of coverage and a
    # tenth of a percent of texture for nine times the interpolation". Those
    # measures -- median displacement by slope band, cell-scale semivariance,
    # share of cells within 0.25 m, coverage -- all average over a whole tile. A
    # fault occupying about 1% of cells and arranged in thin lines moves none of
    # them appreciably, and their silence was read as absence. It was visible in
    # the imagery throughout.
    #
    # Evaluate a candidate neighbourhood on a tail statistic as well as a
    # central one. The 99th percentile of local slope against the reference's is
    # the cheapest one that sees this.
    #
    # These four values move together and are locked per profile for that
    # reason: 64 points inside a radius derived for 16 would be packed tightly
    # enough to lose matrix rank, and 16 inside a fixed 20 m was measured doing
    # exactly that on dense data. 16 with a density-derived radius remains
    # selectable for speed, and costs the slopes.
    search_radius_m: float = NCALM_TERRASCAN["kriging_search_radius_m"]
    search_radius_mode: str = "fixed"    # density | variogram | fixed
    search_radius_factor: float = 4.0    # multiplier on the k-neighbour spacing
    search_radius_floor_m: float = 3.0
    search_radius_ceiling_m: float = NCALM_TERRASCAN["kriging_search_radius_m"]
    max_points: int = 64                 # the real control on neighbourhood
                                         # size; see the note above
    min_points: int = 1                  # below this a cell is left as nodata
    dem_cell_m: float = NCALM_TERRASCAN["dem_cell_m"]

    # The one-sided fringe is one search radius wide, because that is how far a
    # cell reaches for support that exists on one side only. Measured on
    # l0s395 at a 20 m radius: median error 6.6 to 9.9 m within 20 cells of the
    # boundary, 99-100% of cells past half a metre, falling to 0.002 m and 0.20%
    # beyond 40 cells. A fixed cell count cannot express that, so the trim is in
    # metres and scales with whatever radius ran.
    erode_m: float = None                # None follows the search radius
    erode_factor: float = 1.0
    erode_cells: int = 0                 # kept for callers that set it directly

    # How far a cell may sit from the nearest accepted ground return and still
    # be written. Morphological erosion of the coverage mask cannot express
    # this: it trims inward from whatever the interpolator chose to fill, so an
    # engine that extrapolates keeps a fringe no erosion depth reaches, and an
    # engine that does not loses good ground.
    #
    # 4.0 m is measured, not chosen. On Centro GLAS l0s441, distance to the
    # nearest return was thresholded across 1 to 7 m and scored against the
    # TerraScan/Surfer extent for the same tile:
    #
    #     radius   1.5     2.5     3.0     4.0     5.0     7.0
    #     IoU      0.9422  0.9737  0.9782  0.9782  0.9749  0.9655
    #     lost     5.08%   1.54%   0.81%   0.33%   0.15%   0.04%
    #     fringe   0.77%   1.13%   1.40%   1.86%   2.37%   3.42%
    #
    # 3.0 and 4.0 tie on agreement; 4.0 keeps 99.67% of the reference's cells
    # against 3.0's 99.19%, for 0.46% more fringe. Set to None to skip the
    # support mask and trim by erosion alone, which is the pre-1.1 behaviour.
    #
    # Interior voids are filled before the mask is applied. Ground returns are
    # sparse under canopy, so a bare distance threshold punches holes through
    # the middle of good terrain, and a hole inside the survey is not an edge.
    edge_fix_m: float = 4.0

    # The across-track strip, which is the edge treatment the reference
    # products themselves use. From the G-LiHT reprocessing work: where
    # adjacent flight passes overlap and are merged, the ground classification
    # at the overlap differs from the interior and leaves a seam about 21 m
    # wide down both sides of every transect. `clean_las_edges.py` removes
    # STRIP_M = 21.0 m per side; `mask_dem_edges.py` removes 64 px at 0.33 m,
    # which is the same 21 m.
    #
    # This replaces the morphological erosion rather than adding to it, and it
    # is the less aggressive of the two in the direction that matters: erosion
    # takes a ring, so it strips the along-track ends where there is no seam,
    # and it widens every interior void by its own depth. The strip touches
    # only the first and last run of each line.
    #
    # `edge_axis` follows the G-LiHT naming -- _GLAS_ flies north-south so its
    # seams are left and right, _NFI_ flies east-west so they are top and
    # bottom. None infers from the name, then from the aspect ratio.
    edge_strip_m: float = 21.0
    edge_axis: str = None                # "x" | "y" | None to infer

    # the documented figure, kept so the source is never lost even though the
    # default no longer uses it directly
    ncalm_search_radius_m: float = NCALM_TERRASCAN["kriging_search_radius_m"]

    # calibration window selection
    #
    # A window is only worth calibrating on if it has ground to calibrate
    # against. Selecting for relief alone put l8s431 on a window carrying 0.92
    # ground points per square metre, where every configuration failed on 93%
    # of cells above 30 degrees and the result swung the seven-tile aggregate
    # from 10.34% to 22.20%. Relief without returns measures the survey, not the
    # filter.
    #
    # 2.5 is below every usable window measured here (2.79 to 4.57) and well
    # above the failure at 0.92.
    calib_min_coverage: float = 0.70
    calib_min_ground_density: float = 2.5
    calib_min_steep_fraction: float = 0.02   # below this a window cannot
                                             # discriminate; l0s419 had 0.0002

    # bookkeeping
    min_ground_points: int = 500         # below this, rasterising is not honest
    notes: str = ""

    def to_json(self, path):
        d = asdict(self)
        d["_source_method"] = NCALM_TERRASCAN
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(d, fh, indent=2)
        return path

    @classmethod
    def from_json(cls, path):
        with open(path, encoding="utf-8") as fh:
            d = json.load(fh)
        d.pop("_source_method", None)
        d["passes"] = [GroundPass(**p) for p in d.get("passes", [])]
        return cls(**d)


# --- the locked baseline ----------------------------------------------------
#
# These values are measured, not chosen, and are fixed as of 2026-09-20. Every
# one of them is defended by a figure in the comments above and in
# docs/open_observations.md. They are recorded here separately from the preset
# so that a change to the preset can be detected rather than merely noticed:
# `verify_baseline()` raises if the two drift apart.
#
# Changing a value here is a deliberate act and should come with the measurement
# that justifies it, on terrain that can show the difference. Three of the four
# parameters this replaces were chosen on windows containing 0.00% to 0.20%
# steep ground, and all three were wrong.

BASELINE = {
    "locked": "2026-10-07",
    "algorithm": "smrf",
    "passes": 1,
    "slope": 0.1584,            # tan 9 deg, the source's pass-1 iteration angle
    "cell_m": 1.0,              # raised 2026-10-07; see the note above
    "scalar": 0.0,              # the slope-coupled tolerance, measured off
    "window_m": 32.0,           # saturated: PDAL's opening stops responding
                                # above about 8 m on this data, so 8, 32, 100
                                # and 200 give identical output
    "remove_low_noise": False,  # ELM: inert on this data, 0.00% removed
    "remove_outliers": True,    # tile-dependent, see the note above
    "clean_height_m": 0.20,
    "clean_patch_m": 0.75,
    "search_radius_ceiling_m": 20.0,
    "dem_cell_m": 1.0,

    # The interpolation neighbourhood varies by profile, so it is locked per
    # profile rather than globally. The four values move together: 64 points
    # inside a radius derived for 16 would be packed tightly enough to lose
    # matrix rank, and 16 inside a fixed 20 m was measured doing exactly that
    # on dense data. Locking them separately would allow a combination neither
    # was measured at, and the window reports any mismatch as drift -- which it
    # did, on a default run, when the shipped default and this lock disagreed.
    "profiles": {
        "ncalm":    {"threshold_m": 0.25, "clean_vegetation": True,
                     "max_points": 64, "min_points": 1,
                     "search_radius_mode": "fixed", "search_radius_m": 20.0},
        # Every profile uses the source's own neighbourhood: 20 m radius, up to
        # 64 points, minimum 1 (Estrada-Belli et al. 2025 supplementary;
        # GLiHT_Methods_Materials section 5.2). The replication target always
        # did, because a configuration whose purpose is fidelity should not
        # differ from its target in the size of the neighbourhood it averages.
        # The others were briefly set to 16 at a density-derived radius on
        # 2026-10-05 and put back the next day, when a smaller neighbourhood was
        # found to cut spurious narrow grooves into every slope -- see the note
        # at `search_radius_m`, and note what it cost to find.
        "baseline": {"threshold_m": 0.25, "clean_vegetation": False,
                     "max_points": 64, "min_points": 1,
                     "search_radius_mode": "fixed", "search_radius_m": 20.0},
        "clear":    {"threshold_m": 0.25, "clean_vegetation": True,
                     "max_points": 64, "min_points": 1,
                     "search_radius_mode": "fixed", "search_radius_m": 20.0},
        "deep":     {"threshold_m": 0.25, "clean_vegetation": True,
                     "max_points": 64, "min_points": 1,
                     "search_radius_mode": "fixed", "search_radius_m": 20.0},
    },
}


def profile_of(cfg):
    """Which preset a config is, or None if it matches none of them."""
    for name in BASELINE["profiles"]:
        if name in PRESETS and PRESETS[name] == cfg:
            return name
    return None


def verify_baseline(cfg=None, profile=None):
    """Raise if a config has drifted from the locked baseline.

    Call it before a production run. It compares only the parameters that were
    settled by measurement; everything else is free to vary.

    Two of those parameters -- the SMRF elevation threshold and whether the
    cleanup runs -- legitimately differ between profiles, so they are checked
    against that profile's locked pair rather than against one global value.
    A config matching no known profile is checked on the shared parameters
    only; the caller chose those two deliberately and is not drifting.
    """
    cfg = cfg if cfg is not None else PRESETS["clear"]
    profile = profile or profile_of(cfg)
    bad = []
    if len(cfg.passes) != BASELINE["passes"]:
        bad.append("passes: %d, baseline %d"
                   % (len(cfg.passes), BASELINE["passes"]))
    if cfg.passes:
        p = cfg.passes[0]
        for k in ("algorithm", "slope", "cell_m", "scalar", "window_m"):
            if getattr(p, k) != BASELINE[k]:
                bad.append("%s: %r, baseline %r"
                           % (k, getattr(p, k), BASELINE[k]))
    for k in ("remove_low_noise", "remove_outliers", "clean_height_m",
              "clean_patch_m", "search_radius_ceiling_m"):
        if getattr(cfg, k) != BASELINE[k]:
            bad.append("%s: %r, baseline %r" % (k, getattr(cfg, k), BASELINE[k]))
    if profile:
        want = BASELINE["profiles"][profile]
        if cfg.passes and cfg.passes[0].threshold_m != want["threshold_m"]:
            bad.append("threshold_m: %r, %s baseline %r"
                       % (cfg.passes[0].threshold_m, profile, want["threshold_m"]))
        if cfg.clean_vegetation != want["clean_vegetation"]:
            bad.append("clean_vegetation: %r, %s baseline %r"
                       % (cfg.clean_vegetation, profile, want["clean_vegetation"]))
        # The neighbourhood is four values that only make sense together, and
        # it differs by profile: the replication target takes the source's 64
        # inside a fixed 20 m, the others ship 16 with the radius from density.
        for k in ("max_points", "min_points", "search_radius_mode",
                  "search_radius_m"):
            if getattr(cfg, k) != want[k]:
                bad.append("%s: %r, %s baseline %r"
                           % (k, getattr(cfg, k), profile, want[k]))
    if bad:
        raise ValueError("config has drifted from the %s baseline:\n  %s"
                         % (BASELINE["locked"], "\n  ".join(bad)))
    return True


# Presets. NCALM is the translation above; the others exist because the
# translation is inexact and the right settings are terrain-dependent.

# The default pass, so a preset that changes one parameter inherits the rest
# rather than respecifying them -- GroundPass's own field defaults are not the
# same as the configured pass, and listing fields by hand has silently changed
# cell size and slope before.
_PASS = ReplicalmConfig().passes[0]

PRESETS = {
    # The translation at the published 0.5 m elevation threshold. Note this
    # entry still carries clean_vegetation=True, so it is not the published
    # method alone -- "baseline" is. The two differ only in that flag.
    "ncalm": ReplicalmConfig(),

    # The delivered method: the translation at a 0.25 m threshold, plus removal
    # of the residual low vegetation the ground filter still accepts. The
    # threshold came from a 2 x 2 factorial against the 0.5 m value, crossed
    # with Clear on and off, on l0s395, l8s431 and l0s444. The two instruments
    # overlap without being redundant -- tightening the threshold takes some of
    # what Clear would have taken (Clear's share fell 10.47% to 6.35% on
    # l8s431) and Clear still finds 6 to 12% afterwards. Together they put bias
    # within 3.3 mm of zero on all three tiles, the best of any arm, with no
    # sign of the overshoot that compounding two removals might have produced.
    #
    # The cleanup measures against a fitted plane rather than the 10th
    # percentile of a disc, because the percentile sits at the downhill edge of
    # that disc and so removes ground for being on a slope: on synthetic bare
    # ground it took 90% of it at 25 degrees and 97% at 40. On l0s426 the
    # percentile rule removed 65% of the classified ground and the plane rule
    # at 1.20 m removes 14%, while landing closer to the archive on every
    # measure -- cells within 25 cm 64.6% against 59.9%, and semivariance
    # 0.14312 against 0.15292 on an archive of 0.10984. The 1.20 m fit radius
    # was swept over 0.75 to 2.00 m on l0s395 and confirmed on l3s457 and
    # l0s426; it is separate from `clean_patch_m` because the threshold wants a
    # small neighbourhood and the plane has to BE the local surface, and a
    # clump of scrub wider than the disc it is fitted over tilts the plane onto
    # itself and escapes.
    #
    # Worth knowing before this is taken as settled: at 4.24 ground returns per
    # square metre the cleanup has become nearly redundant. On l0s395 the
    # classification alone now matches the best cleanup arm at 99.9% of cells
    # within 25 cm and sits closer to the archive's texture than it does
    # (0.00667 against 0.00633, archive 0.00663), so adding the plane rule
    # there pushes the surface slightly below the reference rather than onto
    # it. It earns its place at 0.86 returns per square metre and not at 4.24,
    # which argues for making it density-conditional rather than a fixed
    # property of the profile.
    "clear": ReplicalmConfig(clean_slope_aware=True,
                             clean_plane_patch_m=1.20),

    # The translation alone, with nothing removed beyond what the source method
    # removes. This is what the published comparison figures are measured
    # against, and the right choice on terrain where the cleanup threshold has
    # not been checked.
    "baseline": ReplicalmConfig(clean_vegetation=False,
                                max_points=64, min_points=1,
                                search_radius_mode="fixed",
                                search_radius_m=20.0),

    # Clear on a grid finer than the point spacing. The extra resolution is for
    # rendering, not for accuracy -- see cell_factor.
    "deep": ReplicalmConfig(cell_factor=0.7, clean_slope_aware=True,
                            clean_plane_patch_m=1.20),
    "ncalm_csf": ReplicalmConfig(passes=[
        GroundPass(algorithm="csf", csf_rigidness=2, csf_threshold_m=0.5,
                   stands_in_for="TerraScan pass 1, via cloth simulation"),
        GroundPass(algorithm="csf", csf_rigidness=3, csf_threshold_m=0.3,
                   reuse_ground=True,
                   stands_in_for="TerraScan pass 2, via cloth simulation"),
    ]),
    "conservative": ReplicalmConfig(passes=[
        GroundPass(slope=0.10, threshold_m=2.0,
                   stands_in_for="keeps less ground; use where relief is low"),
    ]),
}
