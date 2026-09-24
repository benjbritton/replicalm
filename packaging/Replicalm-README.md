# Replicalm 1.0

Bare-earth processing for airborne lidar. Give it a point cloud and it returns a
ground surface with the vegetation and buildings removed, and optionally the
relief visualizations used to read archaeological landscapes.

It reproduces the NCALM processing workflow documented in the Estrada-Belli et
al. 2025 supplementary materials, using only open software. The published method
needs TerraScan and Golden Surfer; this needs neither.

---

## Installing

Run `Replicalm-1.0.0-setup.exe` and follow the prompts. It needs 64-bit Windows
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
| **Point cloud** | the LAS or LAZ file to process |
| **Output folder** | where results are written; it is created if it does not exist |
| **Method** | which processing method to run. Leave it on **Clear** unless you have a reason not to; the three are described below |
| **Cell size (m)** | metres per pixel of the output map. Leave it on `auto` unless you have a reason not to |
| **Ground threshold (m)** | how far a return may stand above the provisional surface and still count as ground. Leave it on `auto` unless a tile comes out wrong |
| **Also build the G1 image** | produces the relief visualizations as well as the elevation model |
| **Keep classified points** | writes out the point cloud with ground returns labelled, so the classification can be inspected |

Press **Run**. The progress bar tracks the run and the message pane reports each
stage as it finishes. Processing a full survey tile takes tens of minutes to a
couple of hours depending on its size and whether the images are built; a
smaller area takes minutes.

### About the method

**Clear** is the default and the recommended choice. It runs the translated
NCALM workflow and then removes the low vegetation that any ground filter
accepts as terrain — the scrub, grass tussocks and shrubs that sit a few
centimetres above the soil. On the surveys it was developed against this
removes a small percentage of the returns and takes the speckled texture out of
the relief images, leaving platform surfaces and field boundaries legible.

**Baseline** runs the translated workflow alone, with nothing removed beyond
what the published method removes. It is the right choice in two cases: when
you want output directly comparable to the commercial workflow, and when you
are working on terrain unlike the tropical forest floor the cleanup was
developed on. The cleanup threshold is a measured value from one survey
campaign, not a constant of nature, and Baseline is the option that does not
assume it travels.

**Deep** is Clear on a grid about 1.4 times finer than the point spacing. The
extra resolution goes into rendering rather than measurement: slope and
sky-view factor come out without stair-stepping at cell boundaries, which reads
better at high zoom. It roughly doubles the run time and the file size, and it
does not make the surface more accurate. Use it for figures.

Whichever you choose is recorded in the settings file, so a result can always
be traced back to the method that produced it.

### About the ground threshold

This is the one setting worth reaching for when a tile comes out wrong, because
it decides what the classifier is willing to call ground. A return standing
higher than the threshold above the provisional surface is treated as something
else — vegetation, a building, a wire.

`auto` uses the method's own value: **0.25 m for Clear and Deep**, 0.50 m for
Baseline, which is the figure the published workflow specifies.

Raise it — 0.5 m, or higher on rugged ground — if the model looks scraped, with
ridge crests flattened or small rises missing. Lower it if flat ground comes out
speckled, or if low scrub is showing up as terrain. The two failures look quite
different once you have seen them: too loose leaves texture that should not be
there, too tight removes relief that should.

The 0.25 m default was measured. Against 0.50 m on three test tiles it cut the
elevation bias from +0.037 m to +0.016 m before any cleanup ran, and with the
cleanup applied it brought bias within 3 mm of zero on all three. On one steep
tile, though, the result also came out smoother than the reference product, so
on strongly relieved terrain it is worth comparing both.

### About cell size

`auto` measures how densely the survey's ground returns actually fall and
chooses a cell to match — roughly the average spacing between returns. On the
G-LiHT surveys this produces about 0.5 m. A finer grid than the data supports
adds file size without adding detail; a coarser one throws detail away.

Set a number instead if you need to match an existing product, or if you are
mosaicking with rasters made at a particular resolution.

## What you get

In the output folder:

| file | what it is |
|---|---|
| `<name>_DEM.tif` | the elevation model: a single-band GeoTIFF of ground height in metres, carrying its coordinate system so it lands correctly in any GIS |
| `<name>_config.json` | every setting the run used, plus how many cells were filled across gaps and how far the edge was trimmed |
| `<name>_ground.laz` | the classified point cloud, if you asked for it: LAS 1.4, compressed, carrying its coordinate system and nothing proprietary |
| `rvt\` | the G1 composite and the five visualizations it is blended from, if you asked for them |

In the elevation model, **zero means no data**, and no real ground value is ever
exactly zero. So a gap can never be mistaken for terrain, and terrain can never
be mistaken for a gap, even if the file passes through software that discards
the usual no-data marker.

## From a command line

For scripting, `replicalm-cli.cmd` in the installation folder does the same work:

```
replicalm-cli.cmd survey.las -o C:\output
replicalm-cli.cmd survey.las -o C:\output --g1 --cell 0.5
replicalm-cli.cmd survey.las -o C:\output --profile baseline
```

`--help` lists the options. Anything the window can do, this can do, and a run
from either is reproducible from the settings file it writes.

## If something goes wrong

**"Global encoding WKT flag not set for point format 6 - 10"** — the file
declares LAS 1.4 with a modern point format but omits a flag its own header
requires. Replicalm reads it anyway, without the coordinate system, and says so.
The elevation model will be correct but unprojected; supply the projection in
your GIS, or ask whoever produced the file to correct the header.

**The output has no coordinate system** — the source file carried none that
could be read. See above.

**The G1 images fail while the elevation model succeeds** — the visualization
step needs components the elevation path does not. On a managed machine, an
application-control policy may block them. The elevation model is unaffected.

**"only N ground points"** — the classifier found too little ground to build a
surface from. Usually the cloud is very sparse, covers mostly water or canopy,
or has already been filtered to non-ground returns.

## Licensing

Replicalm is MIT licensed — free to use, modify and redistribute, including
commercially.

It includes components under their own licenses, all permissive: PDAL, GDAL,
PROJ, NumPy, SciPy, Python and the Relief Visualization Toolbox. Their license
texts are installed in the `licenses` folder, and `THIRD-PARTY-NOTICES.md` lists
them.

If you publish work using the visualizations, cite the Relief Visualization
Toolbox:

> Kokalj, Ž., Somrak, M. (2019). Why Not a Single Image? Combining
> Visualizations to Facilitate Fieldwork and On-Screen Mapping.
> *Remote Sensing* 11(7), 747.

## Further reading

The `docs` folder in the installation carries the full account: how the method
was translated, how closely the output matches the commercial workflow, what was
measured and what remains unsettled.

---

Benjamin Jay Britton, 2026
