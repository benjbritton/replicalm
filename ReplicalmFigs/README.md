# ReplicalmFigs

Figures for the Replicalm article, separated by the version of the program whose
output they show, because the two versions do not produce the same surface and a
figure that does not say which one made it cannot be used.

    Figs_v1.0/    the figures as they stood through the 1.0 line
    Figs_v1.1/    figures for 1.1, to be produced

## Three configurations, not two

"Released" and "1.0" are not interchangeable, and neither is the same as what
produced the 458-tile corpus. Anything labelled by version needs to name which
of these it means:

| | scalar | threshold | window | neighbourhood | residual removal |
|---|---|---|---|---|---|
| **1.0.0 as installed** (tagged 2026-09-24) | 1.25 | 0.25 m | 70 m | 16 / density | 10th-percentile disc |
| **the 458-tile corpus** (built 2026-09-25 onward) | 1.25 | **0.50 m** | 70 m | 16 / density | 10th-percentile disc |
| **1.1.0** | **0.00** | 0.25 m | **32 m** | 16 / density | **plane fitted over 1.20 m** |

The corpus differs from the installer because the elevation threshold was
changed the day after the 1.0.0 tag and before the survey ran. So the published
raster library and the published installer were not built by the same
configuration, and a figure drawn from the corpus is not a figure of 1.0.0.

## What Figs_v1.0 actually contains

The 21 files copied here are `docs/figures/` as it stood at the end of the 1.0
line. They were produced over several weeks of development, and **not all of
them correspond to the 1.0.0 default**: several predate the tag and were made
with intermediate configurations while a parameter was being settled. The folder
name records when they were current, not that a single configuration produced
them all.

Before any of these appears in the article, each one needs its configuration
confirmed against the `_config.json` beside the raster it was drawn from, or
regenerating. That is a smaller job than it sounds for the figures that
illustrate a mechanism — `fig3_error_by_slope`, `fig4_terracing_vs_cell` and
`fig5_classifier_slice` are about behaviour rather than about a product — and
unavoidable for the ones that show a surface.

## What Figs_v1.1 needs

Nothing has been produced yet. The measurements are in
`docs/paper/replicalm_isprs_draft_v3.md` sections 4.7 to 4.13 and the figures
they imply are:

- displacement from the archive by slope band, the four parameter arms on
  l0s426 (§4.7)
- synthetic bare ground removed against slope, percentile rule versus plane
  (§4.8)
- the one-sided admissibility test: where each surface sits relative to the
  cloud's own floor (§4.9) -- this is the figure that carries the discussion's
  main claim and does not exist in any form yet
- neighbourhood against texture, coverage and elevation (§4.11)
- 1.0 against 1.1 on both tiles (§4.13)

The surfaces for all of these already exist under
`tests/clear_fulltile/l0s426/` and `tests/clear_fulltile/l0s395/`, so they are
plotting jobs rather than processing jobs.
