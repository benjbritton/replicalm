# Replicalm 1.3

Bare-earth processing for airborne lidar. Give it point clouds and it returns
ground surfaces with the vegetation and buildings removed, and optionally the
relief visualization used to read archaeological landscapes.

It reproduces the NCALM processing workflow documented in the Estrada-Belli et
al. 2025 supplementary materials, using only open software. The published method
needs TerraScan and Golden Surfer; this needs neither.

---

## What changed since 1.1

**The surface is now built by a tensioned spline rather than by kriging.**
Kriging put a fine corrugation on sloping ground — narrow ridges and grooves
running across the slope, visible in a plain hillshade and not a property of the
terrain. The default interpolator is now a tensioned spline at a tension of
0.03 over 48 neighbouring returns, which does not produce it. This changes the
texture of every output surface, so a 1.3 result and a 1.1 result over the same
tile will differ in fine detail even at identical settings.

**The ground filter's working grid went from 0.5 m to 1.0 m.** The filter builds
a provisional minimum-elevation surface on a grid of its own before deciding
what is ground, and that grid is separate from the output resolution. Measured
over a whole survey tile, raising it took error against the reference surface
from 10.3 cm to **7.9 cm**, a 23% reduction. One caveat worth recording: on
steep faces specifically, a terracing measurement went the other way — 2.44% at
0.5 m against 5.06% at 1.0 m — so the gain is a tile-wide average and not an
improvement everywhere. The 1.0 m value is the one validated at full-tile scale
and is what ships.

**Edge handling was rebuilt.** The old trim eroded inward from whatever the
interpolator had filled, using the kriging search radius — 20 m — even when
kriging was not the engine that ran. That removed real ground along swath edges
and could leave holes in the interior surrounded by valid data. The trim is now
two separate things: a support test that keeps cells within 4 m of an actual
ground return, with interior holes filled rather than opened, and a directional
21 m strip across the flight line where coverage genuinely thins. On the survey
tiles tested, no interior holes remain.

**The bundled lidar library is newer.** PDAL moves from 2.6.3 to **2.7.2**, its
Python bindings from 3.2.3 to **3.4.5**, and GDAL from 3.8.4 to **3.9.2**. The
ground classification is unaffected: SMRF's output was compared point for point
across the two libraries on two survey tiles, 105,521 and 49,071 ground returns,
and is identical. So results produced under 1.1 and under 1.3 are comparable on
classification, and differ only through the interpolation and edge changes
above.

**An experimental second ground classifier is included but switched off.**
`ground_tsr` works from raw point clouds alone and was built to reproduce a
commercial TerraScan ground population without using any commercial
classification as input. It is in the installation and it runs, but it does not
yet beat the filter that ships as the default: over a full tile it scored 7.87 cm
against the reference where the default scored 7.88 cm, and on a second tile
that no part of it was tuned on, 5.47 cm against 5.50 cm. Those differences are
not meaningful. It is present for continuity of the work, not as a recommended
setting, and reaching it requires editing a configuration file rather than
ticking anything.

**1.2.0, if you have it, is this same application on the previous lidar
library.** It carries the interpolation, working-grid and edge changes above but
PDAL 2.6.3. There is no reason to prefer it unless you need to match output
produced under the older library exactly.

---

## Known issue in this version

**The Interpolation control does not affect the default surface.** It sets the
neighbourhood for the kriging path — 64 returns within 20 m, or 16 with the
distance taken from density — and kriging is no longer the default engine. Under
the tensioned spline the surface is built from 48 neighbours regardless of what
this control says, and the line the log prints about interpolation describes
parameters that did not run.

Nothing is computed wrongly and no output is invalid: the spline uses its own
documented neighbourhood either way. But the control is inert at the default
setting and the log line is misleading, so do not use either to explain a
difference between two runs. If you need the kriging path and its
neighbourhood choice, set the interpolator back to kriging in the configuration
file; otherwise leave this control alone.

---

## Installing

Run `Replicalm-1.3.0-setup.exe` and follow the prompts. It needs 64-bit Windows
10 or 11 and about 4 GB of disk.

Nothing else has to be installed first. Python, PDAL, GDAL, PROJ and the
visualization libraries all travel inside the installer, in their own private
copy that cannot interfere with other software on the machine — including any
existing Python, conda, ArcGIS or QGIS installation.

The last step of installation configures that bundled copy for wherever you put
it. It takes a few seconds and must be allowed to finish.

Installing 1.3 alongside an earlier version is not supported in the sense that
they share an installation directory by default; point one of them elsewhere if
you want both.

## Using it

Start **Replicalm** from the Start Menu.

| field | what it is |
|---|---|
| **Input** | the point clouds to process. **Add file(s)…** takes any number of LAS or LAZ files; **Add folder…** takes a folder and finds every point cloud inside it, including in subfolders; **Remove** drops whatever is selected in the list. Folders are shown as `[folder]` and the count of clouds found is reported when you add one |
| **Output folder** | where results are written; it is created if it does not exist |
| **Method** | **Clear** or **Base**. Leave it on Clear unless you have a reason not to; both are described below |
| **Resolution** | how much detail the output map carries, in pixels per square metre. Leave it on `auto` unless you have a reason not to |
| **Interpolation** | the kriging neighbourhood. Inert at the default interpolator — see **Known issue** above |
| **Ground threshold (m)** | how far a return may stand above the provisional surface and still count as ground. Leave it on `auto` unless a tile comes out wrong |
| **Also build the G1 image** | produces the relief visualization as well as the elevation model |
| **Save ground-filtered point cloud (LAZ)** | writes out the ground returns, so the classification can be inspected or the cleanup re-derived without classifying again. Ticking it enables the selector beside it, which chooses whether to save the ground before or after the low-vegetation cleanup |
| **Processing notes (optional)** | anything you want recorded beside the output, in your own words |

Press **Run**. The progress bar tracks the whole job across however many files
were added, and the message pane reports each stage as it finishes.

**Cancel** stops the run. The window returns to how it was before Run, with
every selection still filled in, so you can change one thing and start again.
Whatever had already finished stays on disk; a later run over the same file
replaces it rather than reusing it, so a part-written intermediate cannot be
mistaken for a finished one.

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

This control is separate from the ground filter's own working grid, which is
1.0 m and is not exposed here. A 4 px/m² output is still classified on a 1.0 m
provisional surface.

### About interpolation

Each output cell is estimated from the ground returns around it by a **tensioned
spline**: tension 0.03, 48 neighbouring returns. That is the default and in 1.3
it is effectively the only path the window can reach, for the reason given under
**Known issue**.

The spline replaced kriging because kriging left a fine corrugation on slopes —
narrow ridges and grooves across the fall line, present in the surface and not
in the terrain. Adding a smoothing term to kriging was tried and rejected: it
removed the corrugation by removing real microrelief with it. The spline does
not produce the artefact in the first place.

Neither engine changes where the surface sits to any degree that matters. Across
four very different neighbourhoods the elevation moved by less than 6 mm — the
accuracy comes from the classification, and the choice of interpolator is about
texture, coverage and time.

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

One thing the threshold cannot fix: ground returns on small closed landforms —
low mounds and shallow depressions — are discarded at several times the rate of
flat ground, and no value of this setting recovers them. Measured against a
commercial ground classification over a whole tile, the filter misses 1.9% of
returns on planar ground, 7.8% on mounds and 8.3% on terrace shoulders. Every
other ground filter tested behaves the same way, including ones that work on
entirely different principles, so this is a property of ground filtering at
these densities rather than of this implementation. It matters most if you are
mapping small earthworks.

### Keeping a note of what you did

Every run writes a `_config.json` beside its output listing every setting it
used, so a result can always be traced back to the configuration that produced
it, whether or not anyone wrote anything down.

**Processing notes (optional)** adds your own account to that file. It is useful
when a survey needed something other than the defaults and the reason will not
be obvious to whoever reads the result later — including you, six months on.
Nothing requires it and nothing will ask you for it.

## What you get

In the output folder:

| file | what it is |
|---|---|
| `<name>_DEM.tif` | the elevation model: a single-band GeoTIFF of ground height in metres, carrying its coordinate system so it lands correctly in any GIS |
| `<name>_config.json` | every setting the run used, plus how many cells were filled across gaps, how far the edge was trimmed, and how many interior cells the support test restored |
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
replicalm-batch.cmd D:\survey -o D:\output --pixels-per-m2 4
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
default sized so that modest machines can process large surveys, is planned. It
needs more than a size setting to be correct: tiles have to be interpolated with
an overlap and then joined, or the seams show, and the no-data edges of
neighbouring tiles have to be reconciled rather than simply butted together.
Doing that carelessly would produce a map with visible joins, which is worse
than a slow one without them.

For now, if a file is too large for the available memory, split the point cloud
before processing and mosaic the elevation models afterwards in your GIS.

**A choice of interpolator in the window.** The engine is set in the
configuration file, not on screen. See **Known issue**.

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

**A tile has no class 8 points, and the log said height above ground could not
be computed** — the ground surface could not be triangulated on that tile.
Height above ground feeds only the class 8 near-ground band and a check on
absurd elevations; the ground classification and the elevation model are built
without it and are unaffected. In 1.1 this was traced to the bundled lidar
library on certain large clouds, seen on a tile of 55 million returns. 1.3 ships
a newer library, which may or may not resolve it — **that has not been
retested**, so treat the 1.1 behaviour as the expectation until it is.

**A tile's edges look over-trimmed, or there is a hole in the middle of good
data** — this was a defect through 1.1 and is what the edge rebuild addressed.
If you see it in 1.3 output, it is worth reporting with the tile, because the
support test is supposed to make the second case impossible.

**Low mounds or shallow depressions look flattened** — see the note under
**About the ground threshold**. This is a known limit of ground filtering at
these return densities and is not specific to Replicalm.

**It is slower than you expected** — the spline is not a cheap interpolator, and
the 1.1 timing figures were measured on the kriging path. Expect a full survey
tile to take between half an hour and two hours.

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
measured and what remains unsettled. `docs/ground_tsr.md` documents the
experimental classifier, including what it was measured against and where it
fails.

---

Benjamin Jay Britton, 2026
