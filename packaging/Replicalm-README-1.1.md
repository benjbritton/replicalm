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

**The low-vegetation cleanup no longer removes ground for being on a slope.**
It used to compare each return against the lowest few percent of returns nearby,
and on a hillside those sit downhill of the point being judged. On synthetic
bare ground carrying nothing at all it deleted 90% of it at 25 degrees and 97%
at 40. On a real steep tile it removed **60.8%** of the classified ground where
the present rule removes 14%, and the surface was still too high, because the
removal was compensating for the filter error above rather than correcting it.
The rule now measures against a plane fitted to each neighbourhood, so the
terrain cancels and only what stands above it is removed.

**It is faster than 1.0 as well as better.** Measured end to end on a tile of 14
million ground returns: 1.0 took 6,198 seconds, this takes 2,513 — **2.5 times
quicker**. The maximum filter window did it. Above 32 m the ground filter's
output stops changing, so 1.0's 70 m was buying nothing for two and a half times
the classification time.

**The source's larger interpolation neighbourhood is selectable, not the
default** — 64 ground points inside a fixed 20 m, against the standard 16. It
makes the surface texture closer to the reference and fills more of the map, but
on a tile near the usual return density it bought half a point of coverage for
nine times the interpolation time. Worth choosing on a sparse survey; not worth
paying for on every tile.

**The G1 image now works.** In 1.0 the checkbox called a script that was not in
the installer, so on every machine it quietly did nothing and left the elevation
model alone. The recipe is now built in, and the four visualizations it blends
are written alongside it.

**Files and folders.** The window takes any number of files, and folders are
searched for point clouds recursively, including subfolders.

**The window stays usable while it works, and you can stop it.** In 1.0
processing ran inside the window's own process, which stopped it repainting:
Windows added "(Not Responding)" to the title, and that was accurate rather than
a false alarm. Processing now runs separately, so the window answers throughout
— and because it is separate it can be stopped. **Cancel** halts the run and
returns everything to how it was before Run, with your selections still filled
in. There was no way to offer that in 1.0.

**Four methods became two.** The two that are gone did not produce meaningfully
different output: one was another method on a finer grid, which the resolution
control does directly, and the other differed only in using the superseded
cleanup.

**Two labels, and one dialog that is gone.** A box used to interrupt you if
anything differed from the defaults, to say that the defaults were measured.
The field it was guarding is now **Processing notes (optional)** and still
records whatever you type beside the output. **Keep ground points** is now
**Save ground-filtered point cloud (LAZ)**, which says what it saves.

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
| **Input** | the point clouds to process. **Add file(s)…** takes any number of LAS or LAZ files; **Add folder…** takes a folder and finds every point cloud inside it, including in subfolders; **Remove** drops whatever is selected in the list. Folders are shown as `[folder]` and the count of clouds found is reported when you add one |
| **Output folder** | where results are written; it is created if it does not exist |
| **Method** | **Clear** or **Base**. Leave it on Clear unless you have a reason not to; both are described below |
| **Resolution** | how much detail the output map carries, in pixels per square metre. Leave it on `auto` unless you have a reason not to |
| **Interpolation** | how each output cell is estimated from nearby ground points. Leave it on `standard` unless the survey is sparse |
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

### About interpolation

Each output cell is estimated by kriging from the ground returns around it. Two
things control that: how many returns are averaged, and how far out it looks.
They move together, so they are offered as one choice rather than two numbers.

**standard** — 16 returns, with the distance worked out from how dense the
returns are. This is the default and what 1.0 used.

**published** — 64 returns inside a fixed 20 m, the neighbourhood the source
method specifies. On a sparse tile, under one return per square metre, it brings
the surface texture from 35% rougher than the reference to 18% and fills 90.9%
of the reference's cells rather than 86.6%. On a tile at four returns per square
metre it changes the texture by a tenth of a percent and fills half a point
more — for **seven to nine times the interpolation time**. That is why it is not
the default: the gain is real but concentrated on sparse surveys, and the cost
applies to every tile. On hundreds of tiles the difference is days.

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

**A tile has no class 8 points, and the log said height above ground could not
be computed** — the ground surface could not be triangulated on that tile, which
is a defect in the bundled lidar library on certain large clouds. It has been
seen on a tile of 55 million returns. Height above ground feeds only the class 8
near-ground band and a check on absurd elevations; the ground classification and
the elevation model are built without it and are unaffected. The same tile
classified 18,119,998 ground returns either way, identical to the point. Nothing
needs doing unless you want the class 8 band, in which case processing that tile
with a newer build of the library produces it.

**It is slower than you expected** — if you chose `interpolation: published`,
that is seven to nine times the interpolation work of the standard setting. On
the standard setting a tile is about two and a half times **faster** than 1.0,
because the ground filter's window no longer does work that changes nothing.

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
