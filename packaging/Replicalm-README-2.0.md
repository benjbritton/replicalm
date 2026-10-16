# Replicalm 2.0

Bare-earth processing for airborne lidar. Give it point clouds and it returns
ground surfaces with the vegetation and buildings removed, and optionally the
relief visualization used to read archaeological landscapes.

It reproduces the NCALM processing workflow documented in the Estrada-Belli et
al. 2025 supplementary materials, using only open software. The published
method needs TerraScan and Golden Surfer; this needs neither.

**What it is for.** Replicalm did not exist when the terrain products for our
study of the G-LiHT Yucatán surveys were made, and it was not used to make
them: the DEMs in which the Q2000 model detected 57,181 candidate Maya
structures were produced with TerraScan and Golden Surfer. Replicalm
replicates that method with open-source software. It is provided as open
source for two purposes: so that interested researchers can reproduce our
method, and so that anyone can use a free tool to ground-filter raw LAS/LAZ
point clouds, produce filtered DEMs, and make GeoTIFF visualizations of lidar
data.

The defaults follow the published method as closely as open tools allow, and
every one of them is listed and explained below so you can change any of them
deliberately for your own work.
A run records what it resolved each setting to, so a result can always be
traced back to the configuration that produced it.

---

## Why Version 2.0 and not 1.0

Version 2.0 produces slightly better data output: the ground filter's working
grid, the interpolation engine and the coverage edge have all changed, and the
same point cloud yields a visibly different surface. Also, Replicalm v2.0 has
an improved user interface. For integrity of results, products from 1.0 and
2.0 should not be mixed in one study. A download for the Windows-compatible
installer is available on this page.

---

## The defaults, and when to change them

These are the settings that produced the published output. Everything here can
be changed; the window exposes the ones most often wanted, the rest live in the
configuration.

### Ground classification

| setting | default | why, and when to change it |
|---|---|---|
| **Method** | `Clear` | The translated workflow plus removal of low vegetation the ground filter accepts as terrain. Switch to `Base` for the published method with nothing extra removed, or if your terrain is unlike the tropical forest floor the cleanup was developed on |
| cleanup height | 0.20 m | How far a return may stand above a plane fitted to its neighbourhood and still count as ground. Measured on one survey campaign — not a constant of nature |
| cleanup patch | 0.75 m | The neighbourhood that plane is fitted over |
| filter window | 32 m | The largest opening the ground filter applies. Above about 32 m its output stops changing, so larger values cost time and buy nothing |
| slope | 0.158 (9°) | The terrain slope the filter tolerates. Raise it on steeper ground |
| **Ground threshold** | 0.25 m | How far a return may stand above the provisional surface and still be ground. The setting worth reaching for when a tile comes out wrong |
| scalar | 0.0 | How far the threshold grows with local slope. Zero because a non-zero value was admitting canopy as ground on steep tiles |
| working cell | 1.00 m | The grid the filter builds its provisional surface on, separate from output resolution. Validated over a whole tile; 0.5 m was measurably worse overall, though slightly better on steep terraced faces |

### The surface

| setting | default | why, and when to change it |
|---|---|---|
| interpolator | `biharmonic` | A global minimum-curvature solve: the surface of least curvature that honours the data. Alternatives are `kriging`, which is what the published chain used, and `spline`, a faster local fit that loses about half the detail below a metre |
| solve block | 150 m core, 20 m overlap | What makes a whole tile tractable. Fixed in metres, so it does not change with resolution |
| pooling | 0.25 m | Constraints within this distance are pooled before the solve |
| coincident returns | mean | Two returns at one location are averaged. Taking the lowest biases every collision downward, which matters when a ground population holds returns either side of the surface |
| **Resolution** | `auto` | Output cell size, chosen from measured return density — near 0.5 m on these surveys. Set a number to match an existing product |

### The coverage edge

| setting | default | why, and when to change it |
|---|---|---|
| footprint radius | the search radius, 20 m | A cell is part of the survey if a ground return lies within the interpolator's reach of it. Beyond that nothing was measured and the surface would be invention |
| side strip | 21 m | The flight-pass stitching seam, where overlapping passes were merged and the ground classification differs across the join. **Measured on G-LiHT; set it to 0 if your survey has no such seam** |
| end strip | 21 m | The end-of-pass taper, where the scan narrows and the surface is supported on one side only. A different defect from the seam, which is why it is a separate setting |

### Output

| setting | default | why |
|---|---|---|
| no-data | 0 | No real ground value is ever exactly zero, so a gap can never be mistaken for terrain even if the file passes through software that discards the usual marker |
| raster extent | unchanged | Trimmed cells are written as no-data inside the input grid. Nothing is ever cropped, so products register against each other and against anything else on the same lattice |

---

## What changed since 1.0

**The ground filter was wrong in two parameters.** Its elevation tolerance
grew with local slope, and the factor controlling that had been carried across
from the TerraScan mapping without being tested; on a steep forested tile it
was admitting canopy as ground. Correcting it and tightening the threshold took
flat-ground error from 29.9 cm to **3.1 mm**, and cells landing within 25 cm of
the reference from 17% to 58%.

**The low-vegetation cleanup no longer removes ground for being on a slope.**
It used to compare each return against the lowest few percent nearby, which on
a hillside sit downhill of the point being judged. On synthetic bare ground it
deleted 90% of it at 25 degrees. It now measures against a plane fitted to the
neighbourhood, so the terrain cancels; on a steep tile it removes 14% where the
old rule removed 65%.

**The surface is built by a global minimum-curvature solve.** Earlier versions
used a local spline that approximates — a smoothing term let the surface sit
away from the observations, rounding off small relief and reversing some closed
forms. Measured against a delivered reference DTM, 25 m inside the boundary,
detail retained below 1 m went from **0.505 to 0.914** of the reference's, and
error from 0.0964 m to 0.0513 m.

**The coverage edge is derived from the data.** The active area is the cells
within the interpolator's reach of a ground return, closed into polygons with
no interior holes, with the seam stripped from the sides and the taper from the
ends. Earlier versions used a union of 4 m discs, which drew the boundary as
overlapping circles with detached islands.

**The coverage polygons are written out.** `<name>_footprint.geojson` beside
the DEM: one polygon per connected area of real data, interior voids as
separate features. It opens in ArcGIS and QGIS and works directly as a clip.

**Voids keep a continuous surface.** A hole in coverage is filled by
interpolation and its extent recorded, rather than written into the elevation
as a marker — which would sit a hundred metres below the terrain and put a
cliff into every derived layer.

**The G1 image works, and is now on by default.** In 1.0 the checkbox called
a script that was not in the installer, so on every machine it quietly did
nothing. The recipe is built in, the four visualizations it blends are written
alongside it on the DEM's own grid, and a run produces the composite unless
you ask it not to.

**The window stays usable and can be stopped.** In 1.0 processing ran inside
the window's own process, so Windows marked it "Not Responding" — accurately.
It now runs separately. **Cancel** halts a run and restores the window with
your selections intact.

**Files and folders.** Any number of files; folders are searched recursively.

**Four methods became two.** The two removed produced no meaningfully different
output.

**The bundled lidar library is newer** — PDAL 2.7.2 with the 3.4.5 bindings and
GDAL 3.9.2. The ground classification is unaffected: compared point for point
across the two libraries on two tiles, it is identical.

---

## Known issue

**The Interpolation control does not affect the default surface.** It sets the
neighbourhood for the kriging path, and kriging is not the default engine.
Nothing is computed wrongly, but the control is inert at the default setting
and the log line about interpolation describes parameters that did not run.
Do not use either to explain a difference between two runs.

---

## Installing

Run `Replicalm-2.0.0-setup.exe`. It needs 64-bit Windows 10 or 11 and about
4 GB of disk.

Nothing else has to be installed first. Python, PDAL, GDAL, PROJ and the
visualization libraries travel inside the installer in their own private copy,
which cannot interfere with any existing Python, conda, ArcGIS or QGIS
installation. The last step configures that copy for wherever you put it; it
takes a few seconds and must be allowed to finish.

## Using it

Start **Replicalm** from the Start Menu.

| field | what it is |
|---|---|
| **Input** | the point clouds. **Add file(s)…** takes any number of LAS or LAZ; **Add folder…** finds every cloud inside a folder, including subfolders |
| **Output folder** | where results are written; created if it does not exist |
| **Method** | **Clear** or **Base** — see the defaults table |
| **Resolution** | pixels per square metre; `auto` matches measured density |
| **Interpolation** | the kriging neighbourhood. Inert at the default engine — see **Known issue** |
| **Ground threshold (m)** | `auto` uses the method's measured 0.25 m |
| **Build the G1 image** | on by default — the relief visualization is the point of the tool for most work. Untick it on the user interface, or pass `--no-g1` on the command line, to write only the elevation model. That saves about a minute a tile and a little hard drive space, but having the G1 and the four visualizations behind it might well be worth it |
| **Save ground-filtered point cloud (LAZ)** | off by default. Saves the ground points the filter kept, so you can look at them yourself or try a different cleanup without starting over. The box beside it picks before or after the cleanup |
| **Processing notes (optional)** | recorded beside the output in your own words |

Press **Run**. A full survey tile takes roughly twenty minutes to an hour; a
smaller area takes minutes. Files are processed one after another — **for many
files use the command line**, which runs several at a time.

```
replicalm-batch.cmd D:\survey -o D:\output
replicalm-cli.cmd survey.las -o C:\output --pixels-per-m2 4
replicalm-cli.cmd survey.las -o C:\output --profile baseline
replicalm-cli.cmd survey.las -o C:\output --no-g1
```

The G1 is built by default, so `--g1` is only there for clarity; `--no-g1`
writes the elevation model alone.

`--help` lists the options for either. Anything the window can do, these can
do, and a run from either is reproducible from the settings file it writes.

## What you get

| file | what it is |
|---|---|
| `<name>_DEM.tif` | the elevation model: single-band GeoTIFF of ground height in metres, with its coordinate system |
| `<name>_footprint.geojson` | the coverage polygons, with interior voids as separate features |
| `<name>_config.json` | every setting used and what the run resolved it to |
| `<name>_ground_cleaned.laz` | the ground points, if you asked for them. Called `_ground.laz` if you saved them before the cleanup |
| `G1\<name>_G1.tif` | the G1 composite, unless you turned it off: three-band 8-bit, for viewing rather than measuring |
| `G1\layers\<name>_*.tif` | the four visualizations G1 blends, each on the DEM's grid |

### A note on the images

Replicalm renders one composite so you can look at what the ground filter
produced. It is not a relief visualization package. The visualizations come
from the **Relief Visualization Toolbox**, called as a library; for the full
range of relief products use RVT directly. It is excellent and it is free.

## Not in this version

**Output tiling.** Large clouds are processed whole, so memory scales with the
file. Split the cloud and mosaic the DEMs if a file is too large.

**Settings are not remembered between sessions.** The window starts with its
defaults every time.

**A choice of interpolator in the window.** The engine is set in the
configuration file.

## If something goes wrong

**"Global encoding WKT flag not set for point format 6 - 10"** — the file
declares LAS 1.4 with a modern point format but omits a flag its own header
requires. Replicalm reads it anyway, without the coordinate system, and says
so. Supply the projection in your GIS.

**Too much trimmed at the edges** — the 21 m strips assume a G-LiHT-style
stitching seam. If your survey has none, that is 21 m of good data removed from
each side and each end. Set both strip widths to 0.

**A ragged or scalloped boundary** — a defect of earlier versions that should
not occur here. The footprint GeoJSON shows exactly what the coverage was
judged to be, which is the quickest way to report it.

**"only N ground points"** — too little ground to build a surface. Usually the
cloud is very sparse, covers mostly water or canopy, or has already been
filtered to non-ground returns.

**The G1 image fails while the elevation model succeeds** — the visualization
step needs components the elevation path does not; on a managed machine an
application-control policy may block them. The elevation model is unaffected.

**A tile has no class 8 points, and the log said height above ground could not
be computed** — the ground surface could not be triangulated. That feeds only
the class 8 near-ground band and a check on absurd elevations; the ground
classification and the DEM are unaffected.

**Low mounds or shallow depressions look flattened** — ground returns on small
closed landforms are discarded at several times the rate of flat ground, and no
threshold recovers them. Every ground filter tested behaves this way, including
ones working on entirely different principles. It matters most if you are
mapping small earthworks.

## Licensing

Replicalm is MIT licensed — free to use, modify and redistribute, including
commercially.

It includes components under their own licenses, all permissive: PDAL, GDAL,
PROJ, NumPy, SciPy, Python and the Relief Visualization Toolbox. Their license
texts are in the `licenses` folder, and `THIRD-PARTY-NOTICES.md` lists them.

The Relief Visualization Toolbox is Apache-2.0, by ZRC SAZU and the University
of Ljubljana (UL FGG), included unmodified and called as a library. If you
publish work using the visualizations, cite it:

> Kokalj, Ž., Somrak, M. (2019). Why Not a Single Image? Combining
> Visualizations to Facilitate Fieldwork and On-Screen Mapping.
> *Remote Sensing* 11(7), 747.

## Further reading

The `docs` folder carries the full account: how the method was translated, how
closely the output matches the commercial workflow, what was measured and what
remains unsettled.

---

Benjamin Jay Britton, 2026
