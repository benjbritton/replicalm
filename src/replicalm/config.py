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
      maximum morphological window and does the same job, so it carries across
      directly.
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
    # `threshold_m` is 0.5, not the source's 3.0, and `cell_m` is 0.5 rather
    # than 1.0. The source's figures scored 17.68% error above 30 degrees
    # against 10.34% for these. `cell_m` is the one parameter the source cannot
    # specify, since TerraScan has no SMRF working grid; left coarser than the
    # output cell it quantizes steep faces into steps, with terracing of 5.06%
    # at 1 m, 2.44% at 0.5 m and 0.68% at 0.25 m. Tie it to `dem_cell_m`.
    reset_classification: bool = True
    passes: list = field(default_factory=lambda: [
        GroundPass(slope=_slope_from_angle(9.0), threshold_m=0.5, cell_m=0.5,
                   reuse_ground=False,
                   stands_in_for="TerraScan pass 1: 70 m, 88 deg, 9 deg, 3 m; "
                                 "pass 2 dropped, see the note above"),
    ])

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
    interpolator: str = "kriging"        # kriging | cloudcompare
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
    # A fixed radius is not the answer either, because it means different things
    # at different densities. At 3 points per square metre `max_points` binds
    # long before a 20 m radius does and 5 m and 20 m give identical results; at
    # 8 points per square metre the same 20 m radius packed the neighbourhood
    # tightly enough to lose matrix rank, putting 77 to 88% of l0s444's cells
    # onto the inverse-distance fallback. `radius_for_density` sets it from the
    # distance at which `max_points` neighbours are actually found.
    #
    # The source's documented 20 m is preserved in NCALM_TERRASCAN and is the
    # ceiling here, so nothing ever searches further than the method specifies.
    # Below 1.0 the grid is deliberately finer than the point spacing. That
    # buys rendering rather than measurement: slope and sky-view factor computed
    # without stair-stepping at cell boundaries, for a 2.0% gain in rasterised
    # fold residual and 2.3x the compute. The Deep profile sets it.
    cell_factor: float = 1.0

    search_radius_m: float = None        # None scales it from ground density
    search_radius_mode: str = "density"  # density | variogram | fixed
    search_radius_factor: float = 4.0    # multiplier on the k-neighbour spacing
    search_radius_floor_m: float = 3.0
    search_radius_ceiling_m: float = NCALM_TERRASCAN["kriging_search_radius_m"]
    max_points: int = 16                 # the real control on neighbourhood
                                         # size; see the note above
    min_points: int = 3                  # below this a cell is left as nodata
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
    "locked": "2026-09-23",
    "algorithm": "smrf",
    "passes": 1,
    "slope": 0.1584,            # tan 9 deg, the source's pass-1 iteration angle
    "cell_m": 0.5,
    "remove_low_noise": False,  # ELM: inert on this data, 0.00% removed
    "remove_outliers": True,    # tile-dependent, see the note above
    "clean_height_m": 0.20,
    "clean_patch_m": 0.75,
    "search_radius_mode": "density",
    "search_radius_ceiling_m": 20.0,
    "max_points": 16,
    "dem_cell_m": 1.0,

    # Two parameters vary by profile rather than being fixed for every run, so
    # they are locked per profile instead of globally. Everything above is the
    # same whichever profile runs.
    "profiles": {
        "ncalm":    {"threshold_m": 0.50, "clean_vegetation": True,
                     "max_points": 16, "min_points": 3},
        # The replication target uses the source's own neighbourhood: 20 m
        # radius, up to 64 points, minimum 1 (Estrada-Belli et al. 2025
        # supplementary; GLiHT_Methods_Materials section 5.2). Replicalm's 16
        # and 3 were chosen on their own merits and are kept for the other
        # profiles, but a configuration whose purpose is fidelity should not
        # differ from its target in the size of the neighbourhood it averages.
        "baseline": {"threshold_m": 0.50, "clean_vegetation": False,
                     "max_points": 64, "min_points": 1},
        "clear":    {"threshold_m": 0.50, "clean_vegetation": True,
                     "max_points": 16, "min_points": 3},
        "deep":     {"threshold_m": 0.50, "clean_vegetation": True,
                     "max_points": 16, "min_points": 3},
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
        for k in ("algorithm", "slope", "cell_m"):
            if getattr(p, k) != BASELINE[k]:
                bad.append("%s: %r, baseline %r"
                           % (k, getattr(p, k), BASELINE[k]))
    for k in ("remove_low_noise", "remove_outliers", "clean_height_m",
              "clean_patch_m", "search_radius_mode",
              "search_radius_ceiling_m"):
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
        for k in ("max_points", "min_points"):
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
    "clear": ReplicalmConfig(),

    # The translation alone, with nothing removed beyond what the source method
    # removes. This is what the published comparison figures are measured
    # against, and the right choice on terrain where the cleanup threshold has
    # not been checked.
    "baseline": ReplicalmConfig(clean_vegetation=False,
                                max_points=64, min_points=1),

    # Clear on a grid finer than the point spacing. The extra resolution is for
    # rendering, not for accuracy -- see cell_factor.
    "deep": ReplicalmConfig(cell_factor=0.7),
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
