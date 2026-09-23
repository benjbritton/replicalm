# Replicalm handoff note

`src/replicalm/config.py` is the authority on the baseline, and
`config.verify_baseline()` will say whether a configuration matches it. Where
this note and the code disagree, the code is right.

## The locked baseline, as of 2026-09-20

| setting | value | note |
|---|---|---|
| ground filter | `filters.smrf`, **one pass** | a second pass gave nothing measurable across seven tiles and doubled runtime |
| slope | 0.1584 | tan 9°, the source's pass-1 iteration angle |
| threshold | 0.5 m | the source's 3.0 m admits nearly everything |
| working grid (`cell_m`) | 0.5 m | tied to the output cell; at 1.0 m steep faces terrace |
| ELM filter | **off** | removed 0.00% of points on every window tested — inert on this data |
| statistical outlier filter | on | effect is tile-dependent; worth sweeping per survey |
| search radius | **density-scaled**, 20 m ceiling | the ceiling is the source's figure; a fixed radius means different things at different densities |
| neighbors (`max_points`) | 16 | binds before the radius does on well-covered ground |
| `tile_km` | `None` — whole extent | tiling is a memory device; kriging now chunks internally |
| output cell | derived, `1/√density` | 0.49 m on l0s395 |

## Recent additions

- `evaluate.py` — archive-independent metrics: buffered leave-one-out,
  returns-below-surface, sharpness at a fixed baseline. Needed because every
  other measure in the project compares against NCALM's output, and all of them
  would worsen if a surface were genuinely better than it.
- `cleanup.py` — the `Clear` rule, plus two rules that failed, kept with their
  numbers so they are not re-derived.
- `finalise.py` — hole filling, radius-scaled edge trim, single-band output
  under the black-means-nodata convention.
- `grid.cell_for_density`, chunked `krige_grid`,
  `fit_variogram(min_nugget_fraction=)`.
- Tokenized benchmark records under `benchmarks/results/deep/`.

## Active tracks

**Clear** — validated on three windows it was not fitted to. Now the
application's default method, selectable against Baseline and Deep. The locked
`BASELINE` dict carries `clean_vegetation: True` as of the 2026-09-23 relock,
so drift detection covers it; the `baseline` preset turns it off.

**Deep** — built and measured on one window of one tile. Density-derived cell
size (`cell ≈ 1/√density`, the resolution knee sitting at mean point spacing)
became the default for all three, so Deep is now defined by `cell_factor=0.7`:
deliberate oversampling to about 1.4x finer than that spacing, 0.34 m at 4.3
returns per m². On the rasterised fold test it improves the median residual 2.0%
over Clear, which is rendering quality rather than accuracy. One window
and one density; the rule wants measuring on a survey of different density
before it is trusted there.

## Next steps

1. Decide whether to promote `Clear` into the baseline.
2. Test `packaging/dist/Replicalm-1.0.0-setup.exe` on a machine that has none of
   the dependencies. The staged build is proven; the install sequence, the
   `conda-unpack` step it runs, and the shortcuts are not.
3. Validate Deep on tiles of different point density.
4. The trench-window result, `open_observations.md` item 1, remains unexplained;
   the SMRF `slope` sweep that would settle it has not been run.

## Not yet published

The repository has no remote and nothing has been pushed.
