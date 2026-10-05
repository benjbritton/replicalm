# Replicalm 1.1

Bare-earth processing for airborne lidar. Give it point clouds and it returns
ground surfaces with the vegetation and buildings removed, and optionally the
relief visualization used to read archaeological landscapes.

It reproduces the NCALM processing workflow documented in the Estrada-Belli et
al. 2025 supplementary materials, using only open software. The published method
needs TerraScan and Golden Surfer; this needs neither.

---

## What changed since 1.0

**The ground filter was wrong in two parameters, and both are fixed.** SMRF's
elevation tolerance grows with local slope, and the factor controlling that
growth had been carried across from the TerraScan parameter mapping without
being tested. On a steep forested tile it was admitting canopy as ground: the
classified surface arrived 30 cm above the reference on flat ground and 73 cm
above it past 40 degrees, before anything else ran. Measured across two tiles at
opposite ends of a 458-tile survey, correcting it and tightening the elevation
threshold took flat-ground error from 29.9 cm to **3.1 mm** and cells landing
within 25 cm of the reference from 17% to 58%.

**The maximum filter window is 32 m, not 70.** Above 32 m the filter's output
stops changing, so the larger value was costing two and a half times the
classification time for identical results.

**Interpolation uses the published neighbourhood** — 64 ground points inside a
fixed 20 m, which is what the source method specifies. It makes the surface
texture closer to the reference and fills more of the map, at the cost of
slower interpolation. The faster 1.0 setting is still selectable.

**The G1 image now works.** In 1.0 the checkbox called a script that was not in
the installer, so on every machine it quietly did nothing and left the elevation
model alone. The recipe is now built in.

**Files and folders.** The window processes any number of files, and folders are
searched for point clouds recursively.

**Four methods became two.** The two that are gone did not produce meaningfully
different output.

---

## Installing

Run `Replicalm-1.1.0-setup.exe` and follow the prompts. It needs 64-bit Windows
10 or 11 and about 4 GB of disk.

Nothing else has to be installed first. Python, PDAL, GDAL, PROJ and the
visualization libraries all travel inside the installer, in their own private
copy that cannot interfere with other software on the machine — including any
existing Python, conda, ArcGIS or QGIS installation.

The last step of installation configures that bundled copy for wherever you put
it. It takes a few seconds and must be allowed to finish.

## Using it

Start **Replicalm** from the Start Menu.

| field | what it is |
|---|---|
| **Input** | the point clouds to process. **Add file(s)…** takes any number of LAS or LAZ files; **Add folder…** takes a folder and finds every point cloud inside it, including in subfolders |
| **Output folder** | where results are written; it is created if it does not exist |
| **Method** | **Clear** or **Base**. Leave it on Clear unless you have a reason not to; both are described below |
| **Resolution** | how much detail the output map carries, in pixels per square metre. Leave it on `auto` unless you have a reason not to |
| **Interpolation** | how each output cell is estimated from nearby ground points. Leave it on `published` unless you need the run to finish sooner |
| **Ground threshold (m)** | how far a return may stand above the provisional surface and still count as ground. Leave it on `auto` unless a tile comes out wrong |
| **Also build the G1 image** | produces the relief visualization as well as the elevation model |
| **Keep ground points** | writes out the ground returns as a point cloud, so the classification can be inspected |

Press **Run**. The progress bar tracks the whole job across however many files
were added, and the message pane reports each stage as it finishes.

A full survey tile takes roughly half an hour to two hours depending on its size
and settings; a smaller area takes minutes. Adding ten tiles takes ten times as
long — they are processed one after another, not at once. **For a large number
of files, use the command line instead**, which runs several at a time and sizes
that to the machine's memory. See below.

### About the method

**Clear** is the default and the recommended choice. It runs the translated
NCALM workflow and then removes the low vegetation that any ground filter
accepts as terrain — the scrub, grass tussocks and shrubs that sit a few
centimetres above the soil.

It measures that height against a plane fitted to each return's neighbourhood.
The older rule compared each return against the lowest few percent of returns
nearby, which on a hillside is downhill of the point being judged, so it removed
ground for being on a slope: on synthetic bare ground with nothing standing on
it, it deleted 90% at 25 degrees and 97% at 40. The plane cancels the terrain
and leaves only what stands above it. On a steep forested tile it removes 14% of
the returns where the old rule removed 65%, and lands closer to the reference on
every measure.

**Base** runs the translated workflow alone, with nothing removed beyond what
the published method removes. It is the right choice in two cases: when you want
output directly comparable to the commercial workflow, and when you are working
on terrain unlike the tropical forest floor the cleanup was developed on. The
cleanup threshold is a measured value from one survey campaign, not a constant
of nature, and Base is the option that does not assume it travels.

How much the cleanup matters depends on how dense your returns are. At under one
ground return per square metre it is worth 6 percentage points of accuracy. At
four or more it removes about 1% of the returns and changes the surface by
millimetres — the classification is already doing the work. If your survey is
dense, Base and Clear will look much the same.

Whichever you choose is recorded in the settings file, so a result can always be
traced back to the method that produced it.

### About resolution

Resolution is given in **pixels per square metre**, and the cell size each
choice means is shown beside it. The relationship is not linear — four pixels
per square metre is a 0.50 m cell, not a 0.25 m one — which is why both numbers
are always on screen.

`auto` measures how densely the survey's ground returns actually fall and
chooses a cell to match, roughly the average spacing between returns. On the
G-LiHT surveys this comes out near 0.5 m, or 4 px/m². A finer grid than the data
supports adds file size without adding detail; a coarser one throws detail away.

Choose a number instead if you need to match an existing product, or if you are
mosaicking with rasters made at a particular resolution.

### About interpolation

Each output cell is estimated by kriging from the ground returns around it. Two
things control that: how many returns are averaged, and how far out it looks.
They move together, so they are offered as one choice rather than two numbers.

**published** — 64 returns inside a fixed 20 m. This is the neighbourhood the
source method specifies. On a sparse tile it brings the surface texture from 35%
rougher than the reference to 18%, and fills 90.9% of the reference's cells
rather than 86.6%. On a dense tile it changes the texture by a tenth of a
percent and still fills more of the map.

**fast** — 16 returns, with the distance worked out from how dense the returns
are. This is what 1.0 used. It interpolates **seven to nine times faster** and
gives up some of the texture improvement. On a small survey the difference in
wall-clock time is minutes; on hundreds of tiles it is days.

Neither setting changes where the surface sits. Across four very different
neighbourhoods the elevation moved by less than 6 mm — the accuracy comes from
the classification, and this choice is about texture, coverage and time.

### About the ground threshold

This is the one setting worth reaching for when a tile comes out wrong, because
it decides what the filter is willing to call ground. Too loose and low
vegetation is accepted as terrain; too tight and real relief is cut away.

`auto` uses the method's own value, **0.25 m**, which is measured. Against 0.50 m
it took flat-ground error on a steep test tile from 10.6 cm to 3.1 mm and
brought the surface texture closer to the reference than any other setting
tried. On a gentler tile it matched the reference's own texture to within 0.6%.

Raise it if a tile has had real terrain removed — terracing on steep faces,
flattened platform edges. Lower it if low vegetation is coming through as
ground. Change it because a tile came out wrong, not as a matter of routine.

### Recording why you changed something

The **Reason for changing** box records that claim in the settings file beside
the output, where anyone reading the result later can see both what was changed
and why. A result produced with altered settings and no stated reason cannot be
defended afterwards, including by you.

## What you get

In the output folder:

| file | what it is |
|---|---|
| `<name>_DEM.tif` | the elevation model: a single-band GeoTIFF of ground height in metres, carrying its coordinate system so it lands correctly in any GIS |
| `<name>_config.json` | every setting the run used, plus how many cells were filled across gaps and how far the edge was trimmed |
| `<name>_ground_cleaned.laz` | the ground returns, if you asked for them: LAS 1.4, compressed, carrying its coordinate system and nothing proprietary |
| `G1\<name>_G1.tif` | the G1 composite, if you asked for it: a three-band 8-bit image for viewing, not for measuring |
| `G1\layers\<name>_*.tif` | the four visualizations G1 is blended from — multi-directional hillshade, slope, positive openness and sky-view factor — each a single-band float GeoTIFF |

In the elevation model, **zero means no data**, and no real ground value is ever
exactly zero. So a gap can never be mistaken for terrain, and terrain can never
be mistaken for a gap, even if the file passes through software that discards
the usual no-data marker.

### A note on the images

Replicalm renders one composite so you can look at what the ground filter
produced. It is not a relief visualization package and is not trying to be one.
The visualizations come from the **Relief Visualization Toolbox**, which is
called as a library and does the actual work; if you want the full range of
relief products and control over their parameters, use RVT directly. It is
excellent and it is free.

## Processing a lot of files

The window processes files one after another. For a survey of many tiles, the
command line runs several at once and sizes how many to the machine's memory:

```
replicalm-batch.cmd D:\survey -o D:\output
replicalm-batch.cmd D:\survey -o D:\output --interpolation fast --pixels-per-m2 4
replicalm-batch.cmd D:\survey\*.las -o D:\output --memory-gb 32
```

It records what it has finished, so interrupting it and starting again continues
rather than repeating. A folder is searched recursively, exactly as the window
searches it.

For a single file, `replicalm-cli.cmd` does the same work as the window:

```
replicalm-cli.cmd survey.las -o C:\output
replicalm-cli.cmd survey.las -o C:\output --g1 --pixels-per-m2 4
replicalm-cli.cmd survey.las -o C:\output --profile baseline
```

`--help` lists the options for either. Anything the window can do, these can do,
and a run from either is reproducible from the settings file it writes.

## Not in this version

**Output tiling.** Very large point clouds are processed whole, which means the
memory needed scales with the file. A user-selectable output tile size, with a
default sized so that modest machines can process large surveys, is planned for
the next version. It needs more than a size setting to be correct: tiles have to
be interpolated with an overlap and then joined, or the seams show, and the
no-data edges of neighbouring tiles have to be reconciled rather than simply
butted together. Doing that carelessly would produce a map with visible joins,
which is worse than a slow one without them.

For now, if a file is too large for the available memory, split the point cloud
before processing and mosaic the elevation models afterwards in your GIS.

## If something goes wrong

**"Global encoding WKT flag not set for point format 6 - 10"** — the file
declares LAS 1.4 with a modern point format but omits a flag its own header
requires. Replicalm reads it anyway, without the coordinate system, and says so.
The elevation model will be correct but unprojected; supply the projection in
your GIS, or ask whoever produced the file to correct the header.

**The output has no coordinate system** — the source file carried none that
could be read. See above.

**The G1 image fails while the elevation model succeeds** — the visualization
step needs components the elevation path does not. On a managed machine, an
application-control policy may block them. The elevation model is unaffected and
the message pane says what went wrong.

**"only N ground points"** — the classifier found too little ground to build a
surface from. Usually the cloud is very sparse, covers mostly water or canopy,
or has already been filtered to non-ground returns.

**It is slower than 1.0** — the published interpolation neighbourhood is seven
to nine times more work than 1.0's. Classification is two and a half times
faster, so a whole tile comes out roughly 25% slower overall. Choose
`interpolation: fast` to get 1.0's speed; you will give up some of the texture
improvement and none of the accuracy.

## Licensing

Replicalm is MIT licensed — free to use, modify and redistribute, including
commercially.

It includes components under their own licenses, all permissive: PDAL, GDAL,
PROJ, NumPy, SciPy, Python and the Relief Visualization Toolbox. Their license
texts are installed in the `licenses` folder, and `THIRD-PARTY-NOTICES.md` lists
them.

The Relief Visualization Toolbox is Apache-2.0, by ZRC SAZU and the University
of Ljubljana (UL FGG). It is included unmodified and called as a library. If you
publish work using the visualizations, cite it:

> Kokalj, Ž., Somrak, M. (2019). Why Not a Single Image? Combining
> Visualizations to Facilitate Fieldwork and On-Screen Mapping.
> *Remote Sensing* 11(7), 747.

## Further reading

The `docs` folder in the installation carries the full account: how the method
was translated, how closely the output matches the commercial workflow, what was
measured and what remains unsettled.

---

Benjamin Jay Britton, 2026
