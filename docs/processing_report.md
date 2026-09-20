# Replicalm

*A plain-speech description of what this is, what it does, and how it relates to
the NCALM workflow described in the Estrada-Belli et al. 2025 supplementary
material.*

Benjamin Jay Britton, 2026

---

## What it is

Replicalm is a set of Python tools that turns a raw airborne lidar point cloud
into a bare-earth elevation model — a picture of the ground with the trees and
buildings taken out. It is built to follow the method that NCALM used to process
the G-LiHT surveys of the Maya lowlands, but using only software that is free
and open.

The reason for building it is simple. The published method depends on commercial
software — TerraScan, ArcGIS Pro, Golden Surfer, LAStools. If you do not have
those licences, you cannot reproduce the published results, and you cannot apply
the same processing to new data. Replicalm is an attempt to close that gap.

## What it does

You give it a point cloud and it gives you back a ground surface. In between it
does five things:

1. **Cuts the survey into manageable pieces.** Lidar tiles are large — a single
   G-LiHT tile can be two gigabytes. The work is done in one-kilometre blocks
   with a ten-metre overlap, so that nothing at a block edge is processed
   without its surroundings.

2. **Sorts the returns into ground and not-ground.** This is the hard part. Each
   laser pulse may bounce off a leaf, a branch, a rooftop, and finally the soil.
   The software has to decide which returns are the floor. Replicalm uses a
   filter called SMRF, which works by sliding a shape across the data and asking
   which points sit close enough to the lowest surface to be ground.

3. **Removes obvious errors.** A few returns land below the real ground —
   reflections, birds, noise in the sensor. These are flagged so they do not
   drag the surface down into a pit.

4. **Fills in the gaps between the ground points.** Ground returns are scattered
   irregularly; a map needs a value at every cell of a regular grid. Replicalm
   uses kriging, a statistical interpolation method, falling back to a simpler
   distance-weighted average where the mathematics becomes unstable.

5. **Tidies the edges.** Cells around the outside of the surveyed area are
   calculated from data on one side only, so they are unreliable and get
   trimmed. Small holes inside the covered area get filled. Areas with no data
   are set to black, and nothing inside the data area is allowed to be black, so
   that a hole can never be mistaken for terrain or the other way round.

The output is a GeoTIFF with two layers: the elevations, and a mask saying where
the data is real. From there the existing G1 visualisation recipe produces the
image archaeologists actually look at.

## The tools it uses

Everything is free and open:

- **PDAL** reads point clouds and runs the ground filter.
- **GDAL** reads and writes the map files and handles coordinate systems.
- **NumPy and SciPy** do the kriging and the arithmetic.
- **QGIS and CloudCompare** are useful for looking at results but are not
  required by the pipeline.

Nothing in the chain needs TerraScan, ArcGIS Pro, Surfer, or the paid LAStools
modules.

## How it compares to the published method

### What is the same

The overall shape of the workflow is the same, and several numbers are taken
directly from the source and not changed:

- one-kilometre processing tiles with ten-metre buffers
- height-above-ground limits of −0.5 m and 600 m
- the near-ground band of ±0.2 m, written to class 8
- exporting only the ground and near-ground classes
- LAS version 1.2 output
- a twenty-metre kriging search radius as the maximum
- a one-metre output grid as the standard product

### What is similar but not identical

**The ground filter.** This is the largest difference and it cannot be avoided.
The published method uses TerraScan's "Classify Ground", which builds a
triangulated surface and grows it upward — an approach with four settings:
window size, maximum terrain angle, iteration angle and iteration distance. No
free tool implements that algorithm with those settings. Replicalm uses SMRF
instead, which does the same job by different means and has different knobs.

The settings were translated by reasoning about what each one controls. The
source's 9° iteration angle becomes a slope tolerance of 0.1584, which is the
tangent of 9°. That is a defensible translation, not an equivalence, and the
result was checked against reference surfaces produced by the original workflow
rather than assumed to be right.

**Interpolation.** Both krige. The published method uses Golden Surfer;
Replicalm uses its own implementation. The search radius is the source's
twenty metres as an upper limit, but it is scaled down where the points are
densely packed, because a wide search on dense data makes the mathematics
unstable.

### What is different

**Two passes became one.** The source runs its ground classification twice, the
second pass refining the first. SMRF has no refining mode — each run starts
over — so a second pass simply redoes the work with a steeper tolerance.
Measured across seven areas, running it twice changed nothing and doubled the
time, so Replicalm runs one pass. This is a real departure from the source and
is recorded as such.

**The elevation tolerance is tighter.** The source allows a point to sit three
metres below the surface and still count as ground. Carried across directly,
that setting lets in almost anything: it produced four times the error of a
half-metre tolerance. Replicalm uses 0.5 m. The source's figure describes
TerraScan's behaviour, not SMRF's.

**There is one setting the source could not have.** SMRF works on its own
internal grid, and TerraScan has no equivalent. Left at one metre while
producing a half-metre map, it turns steep slopes into staircases. Tying it to
the output cell size removes most of that.

**The output grid is declared rather than derived.** The original products have
cells of 0.500042 by 0.499988 metres on origins that are multiples of nothing,
which means neighbouring tiles do not line up without resampling. Replicalm puts
cell edges on whole multiples of the cell size, so tiles join exactly and a
re-run reproduces the previous result. For direct comparison against an original
product, it can adopt that product's grid instead.

**Noise removal is reduced.** The source's workflow includes a low-point filter.
Tested here, one of the two filters available removed no points at all, on flat
and steep ground alike, so it is switched off by default with that measurement
recorded. The other removes a fraction of a percent and its effect varies by
tile.

## How well it works

Tested on South_GLAS_l0s395, a tile with known structures, against the
reference surface from the original workflow:

- Half the cells agree to within 1.5 millimetres.
- Across the whole tile the root-mean-square difference is 6 centimetres.
- On flat ground, fewer than one cell in ten thousand differs by more than half
  a metre.
- On slopes steeper than thirty degrees, about one cell in forty-five does.
- Platforms, plazas, range structures and looter trenches read clearly in the
  resulting visualisation.

The remaining disagreement is almost entirely on steep ground. Below ten
degrees of slope the two surfaces are effectively the same; above thirty degrees
they are not. Slopes above ten degrees make up about 15% of this tile but carry
about 98% of the large differences.

## What is not settled

Three things are recorded in `docs/open_observations.md` rather than glossed:

- One comparison on a steep window produced a much larger difference than any
  other setting has explained. The cause was attributed twice and both
  attributions turned out to be wrong. One candidate remains untested.
- The outlier filter helps on one tile and hurts on another, and it is not known
  what distinguishes them.
- Results obtained before 20 September 2026 were measured on sample areas that
  turned out to contain almost no steep ground, which is where the differences
  live. Those results should be treated as untested rather than as established.

## The settings, as locked on 20 September 2026

| setting | value |
|---|---|
| ground filter | SMRF, one pass |
| slope tolerance | 0.1584 (tangent of 9°) |
| elevation tolerance | 0.5 m |
| filter working grid | 0.5 m |
| low-point filter (ELM) | off — removed nothing when tested |
| outlier filter | on — effect varies by tile |
| kriging search radius | scaled to point density, never above 20 m |
| neighbours per cell | 16 |

These are recorded in code as well as in this table, and the software will
refuse a production run if the configuration has drifted from them.
