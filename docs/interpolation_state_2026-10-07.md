# Interpolation engine — state as of 2026-10-07

A resumption note for the surface-engine work, written before switching to the
global-engine task in `surface.txt`. It records what is committed, what was
measured, and which parameters were deliberately left alone, so none of it has
to be re-derived.

## Committed state

Branch `slope-aware-cleanup`, HEAD `9107a1f`
*"feat: integrate tensioned spline interpolator (tau=0.03) to resolve kriging
slope corrugation"*. Not pushed.

`src/replicalm/config.py`, all six presets:

```
interpolator       = "spline"
spline_tension     = 0.03
spline_neighbours  = 48
spline_smoothing   = 1e-2
```

Also in that commit: `preprocess.py` (shared upstream `dedupe_xy`, which removed
13.56% coincident XY returns on the l0s441 clip), `spline.py` (regularised
spline with tension, Mitasova & Mitas 1993), the engine switch in
`pipeline.py`, and `kriging.py` rewired to the shared dedupe.

Verified on commit: all six presets resolve to the spline engine, and the
default path reproduces an analytic surface (sloped plane plus sinusoid,
irregular sampling, seeded duplicates) at **0.0126 m rms**, 0.0565 m max,
mean +0.0001 m, no singular cells.

## The artifact being chased

Linear slope corrugation at 88.9696867 W 19.2060904 N, Centro GLAS l0s441,
against a TerraScan/Surfer archive reading 22.8 deg at the 99th percentile of
local slope on that window.

Ruled out, with measurements, over this session and the previous one:

| candidate | result |
|---|---|
| G1 blend / visualisation | artifact present in the bare hillshade |
| vegetation cleanup | 33.2 deg off against 32.9 on |
| output cell size | persists at 1.00 m |
| variogram nugget | 25.9 -> 24.0 deg, not a fix |
| flight line / swath | same points give a clean Surfer surface |
| **SMRF bin lattice geometry** | eight bin phases 0 to -70 deg; the ribbon bearing stacks coherently in the WORLD frame (2.69x mean) and collapses in the BIN frame (1.58x). p99 across phases 20.6-23.1 deg, spread 2.5 |
| **SMRF minimum-pick noise** | within-bin vertical spread p50 0.000 m, p90 0.070, p99 0.210; median 2 returns per bin. The minimum is not a noisy estimator here |
| ordinary kriging across voids | the mechanism; RST substantially improves it |

The archive shares the ribbon direction at nearly the same angular strength
(2.54x mean at 15 deg, against our 2.69x at ~1-15 deg), so the direction is not
something the pipeline invented.

## Why no further parameter change was made

Radially averaged PSD against the archive, by wavelength, is the measurement
that matters. On the tuning window the shipped configuration carries 2.14x the
archive's power at 1.0-1.7 m and 0.68x at 3-5 m.

Three sweeps followed, and the third invalidated the first two:

1. **smoothing x tension**, 20 pairs. Diagonal floor at ~0.217 against the
   shipped 0.235 — the two terms are one knob, overall stiffness. Raising
   smoothing *increased* the short-wavelength excess (2.14 -> 3.45), which is
   how solver locality was identified: each cell solves its own system over its
   own k nearest points, so smoothing loosens each cell's tie to its data and
   lets neighbouring solutions drift apart.
2. **k x smoothing**, 32 pairs. k 64 / smoothing 7e-2 / tau 0.03 scored 0.1514
   against 0.2351, a 35.6% reduction, interior to the grid, with k 64 winning at
   both tensions.
3. **Seven-window validation** of that pair. It improved 1 window of 7 and lost
   3.9 to 12.0% on the other six.

The reason is in the accepted-ground density. Over an identical 100 m
footprint: tuning window **1.79 pts/m2**, the six others **4.21, 4.41, 4.89,
5.24, 6.45, 7.27**. The tuning window was chosen long ago *because* it showed
the artifact, which makes it the canopy-void case. At 1.0-1.7 m it holds 2.14x
the archive's power; the six dense windows hold **0.05 to 0.16x** — the deficit
has the opposite sign, so more smoothing makes them worse.

## The panel, and what it establishes

Seven 60 m windows: the tuning window, plus median-relief and steepest windows
on each of l0s439, l0s441, l0s443. Relief 1.5 to 17.6 m, archive p99 13.3 to
26.7 deg. Window definitions and clip paths are in the scratchpad as
`validation_windows.json`, accepted ground cached as `ground_cache.npz`.

k x smoothing over the panel, k in {16,24,32,48}, smoothing in {1e-4,1e-3,1e-2},
tau 0.03, screened on three guards — surface stays inside the archive's
elevation range, max departure from the 5 m local median no worse than the
archive's, every window within 5 deg of archive p99:

| k / smoothing | range | excursion | max dp99 | sparse | dense |
|---|---|---|---|---|---|
| 16 / 1e-2 | 4/7 fail | 7/7 fail | 16.9 | 1.3561 | 0.1666 |
| 32 / 1e-2 | 3/7 fail | 1/7 fail | 4.2 | 0.4395 | 0.5209 |
| 48 / 1e-3 | 4/7 fail | 4/7 fail | 3.2 | 0.6848 | 0.5563 |
| 48 / 1e-2 (shipped) | 1/7 | **0/7** | 5.6 | 0.2351 | 0.8793 |

**Nothing in the grid passes all three guards on all seven windows.** The
shipped configuration has the best guard record and the worst dense-window
spectrum. The configuration that minimises mean spectral error (k 16 / 1e-2)
produces p99 39.6 deg against 22.8 on the sparse window and 2.87 m local
excursions.

Two conclusions were reached and neither is implemented:

**A constant neighbour count is the wrong parameterisation.** Preferred k by
window, at smoothing 1e-3, converted to the radius it implies:

| window | pts/m2 | prefers k | radius |
|---|---|---|---|
| l0s441 tuning | 1.79 | 48 | 2.92 m |
| l0s443 median | 4.41 | 24 | 1.32 m |
| l0s439 median | 4.21 | 16 | 1.10 m |
| l0s439 steepest | 4.89 | 16 | 1.02 m |
| l0s443 steepest | 5.24 | 16 | 0.99 m |
| l0s441 steepest | 6.45 | 16 | 0.89 m |
| l0s441 median | 7.27 | 16 | 0.84 m |

The dense windows cluster at 0.84-1.32 m of radius; the sparse one needs 2.92 m
because 1 m there holds about five points. k 16 is the floor of the grid and
five of six dense windows sit on it, so that preference may continue lower.

**A spectral objective cannot be used alone.** Flattening the spectrum toward
the archive requires a rougher surface, and this engine's extra roughness
arrives as localised overshoot rather than as the archive's distributed
texture. The objective cannot tell the two apart, which is why it nominated a
surface with 2.87 m excursions. Any adaptive scheme needs the guards inside the
objective, not applied afterwards.

## Open items unrelated to the above

- PDAL at 3.2.3; raise to 3.4.x and rebuild the installer. The 3.4.5
  triangulation bug and 3.2.3 gave **identical ground counts (18,119,998)**, so
  parameter measurements transfer across the two.
- `min_nugget_fraction` is measured and documented but never wired into
  `pipeline.py`.
- `replicalm_isprs_draft_v3.md` 4.11 and the discussion should frame this as a
  development finding; "published" has to come off anything attached to the
  458-tile corpus, `ReplicalmFigs/README.md` included.
- Repo-root `README.md` still describes 1.0.
- Survey-scale bias and RMSE distributions need an archive-lattice rebuild
  (1.5-2 days) or a stratified sample.

## Figures from this work

All in `C:\Replicalm\test100526\`, all untracked:

`smrf_grid_phase.png`, `smrf_grid_phase_prediction.png`,
`smrf_grid_phase_stacks.png`, `roughness_psd.png`,
`smoothing_tension_grid.png`, `neighbour_sweep.png`, `k_smoothing_grid.png`,
`validation_k64_s007.png`, `multi_window_sweep.png`,
`zoom_comparison_t002_vs_t003.png`, `zoom_three_panel_t003.png`.

## Environment note

`python.exe` invoked directly from the `replicalm` conda env fails any scipy
import with `0xC06D007F`, a delay-load failure. Run through
`conda run -n replicalm --no-capture-output python <script>`.
