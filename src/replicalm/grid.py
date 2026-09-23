"""The grid, declared rather than derived.

WHY THIS EXISTS
---------------
CloudCompare derives its raster origin from the data bounding box and offers no
way to set it. Two runs on the same file, differing only in a global-shift flag
that should have been immaterial, produced origins 1 cm and 2.5 cm apart -- and
because returns cluster along scan lines rather than spreading evenly, that
reshuffled which points fell in which cell almost everywhere: 99.6% of cells
changed, by 0.071 m RMS and up to 1.50 m. At the magnitude of the sensor's own
vertical accuracy, and not reproducible between runs.

A grid whose origin is declared has none of that. The same rule applied to any
tile in a transect puts cell edges in the same places, so neighbouring tiles
mosaic without resampling and a rerun reproduces the previous result exactly.

THE RULE
--------
Origin snapped outward to a whole multiple of the cell size, in the projected
coordinate system. For a 1 m cell that is whole metres; for 0.5 m it is half
metres. Snapping outward rather than to nearest guarantees the grid covers the
data, so no return falls outside the raster it belongs to.
"""
import math
from dataclasses import dataclass


@dataclass(frozen=True)
class Grid:
    """A raster definition: where it starts, how big the cells are, how many.

    Frozen because a grid that changes after cells have been filled is the
    defect this module exists to prevent.
    """
    origin_x: float          # west edge of the first column
    origin_y: float          # north edge of the first row
    cell: float              # east-west cell size
    width: int
    height: int
    cell_y: float = None     # north-south size, when it differs from cell

    @property
    def dy(self):
        """North-south cell size. Square unless a reference says otherwise.

        Production grids are square. Reference DEMs from the original workflow
        are not: kriging to a data-derived extent in Surfer produced cells of
        0.500042 by 0.499988 m. The difference is 0.054 mm, which is nothing on
        one cell and four centimetres across eight hundred of them, so it is
        carried rather than rounded away.
        """
        return self.cell if self.cell_y is None else self.cell_y

    @property
    def bounds(self):
        """(minx, miny, maxx, maxy) of the raster's outer edges."""
        return (self.origin_x,
                self.origin_y - self.height * self.dy,
                self.origin_x + self.width * self.cell,
                self.origin_y)

    @property
    def geotransform(self):
        """GDAL's six-element transform, north-up."""
        return (self.origin_x, self.cell, 0.0, self.origin_y, 0.0, -self.dy)

    def cell_centres(self):
        """Coordinates of every cell centre, as two 2-D arrays."""
        import numpy as np
        x = self.origin_x + (np.arange(self.width) + 0.5) * self.cell
        y = self.origin_y - (np.arange(self.height) + 0.5) * self.dy
        return np.meshgrid(x, y)

    def describe(self):
        mnx, mny, mxx, mxy = self.bounds
        return ("%d x %d cells of %.3f m, origin %.3f %.3f, "
                "covering %.1f x %.1f m"
                % (self.width, self.height, self.cell, self.origin_x,
                   self.origin_y, mxx - mnx, mxy - mny))

    def aligns_with(self, other, tol=1e-6):
        """True if two grids share cell edges, so they mosaic without resampling."""
        if abs(self.cell - other.cell) > tol:
            return False
        dx = (self.origin_x - other.origin_x) / self.cell
        dy = (self.origin_y - other.origin_y) / self.dy
        return (abs(dx - round(dx)) < tol) and (abs(dy - round(dy)) < tol)


def covered_density(x, y, probe_cell=5.0):
    """Ground-return density over the area actually covered, not the extent.

    WHY THE BOUNDING BOX IS THE WRONG DENOMINATOR
    ---------------------------------------------
    A G-LiHT tile is a flight strip crossing its own bounding box on the
    diagonal, so most of the box holds no returns at all. Measured over the box,
    a full l0s395 tile reads 1.24 returns per square metre; measured where the
    returns are, about 4.2. A cell size derived from the first figure comes out
    at 0.90 m when 0.5 m is right -- nearly twice too coarse, throwing away
    detail the survey paid for.

    The covered area is counted by occupancy: bin the returns into coarse cells
    and total the cells that contain any. `probe_cell` is deliberately much
    larger than any plausible point spacing, so an occupied cell means "the
    survey covers this ground" rather than "a point landed exactly here". At
    5 m it holds 2.5 returns at the sparsest density worth processing and over a
    thousand at the densest, so it is not sensitive to the choice.

    Returns (density, covered_area_m2).
    """
    import numpy as np
    x = np.asarray(x, float); y = np.asarray(y, float)
    if len(x) == 0:
        return 0.0, 0.0
    ix = np.floor((x - x.min()) / probe_cell).astype(np.int64)
    iy = np.floor((y - y.min()) / probe_cell).astype(np.int64)
    occupied = len(np.unique(iy * (ix.max() + 1) + ix))
    area = occupied * probe_cell * probe_cell
    return (len(x) / area if area > 0 else 0.0), area


def cell_for_density(density, factor=1.0, floor=0.10, ceiling=2.0):
    """A cell size matched to how far apart the ground returns actually are.

    WHY NOT A FIXED NUMBER
    ----------------------
    The source method produces 0.5 m and 1 m products, reasonable for the
    surveys it was built on. They are not properties of the method. A cloud at
    one return per square meter cannot support a 0.5 m raster, and one at
    sixteen is wasted on it.

    WHERE THE RULE COMES FROM
    -------------------------
    Rasterised block cross-validation on the Pixoyal window of South_GLAS_l0s395
    at 4.74 ground returns per square meter, a mean spacing of 0.46 m. Residual
    between the raster and returns that did not build it:

        cell    median    marginal gain    cost vs 0.50 m
        1.00 m  0.0865          --              0.25x
        0.50 m  0.0803       -7.2%              1.0x
        0.33 m  0.0787       -2.0%              2.3x
        0.25 m  0.0781       -0.8%              3.9x

    The curve turns at 0.5 m, which is the mean point spacing. Finer grids keep
    improving, by fractions of a per cent for several times the compute and
    storage.

    So the cell tracks 1/sqrt(density). A `factor` below 1 oversamples
    deliberately, buying smoother slope and sky-view-factor rendering rather
    than accuracy -- about 2% on that window at 0.7.

    CAVEAT
    ------
    The knee was located at one density on one window. That it sits at the mean
    spacing is a physically sensible place for it to sit, and it is one
    measurement, not a law. On a survey of markedly different density, measure
    before trusting it.
    """
    d = float(density)
    if not (d > 0) or d != d or d == float("inf"):
        return float(ceiling)
    return float(min(max(factor / math.sqrt(d), floor), ceiling))


def snap_outward(minx, miny, maxx, maxy, cell):
    """Expand a bounding box to the next whole multiple of the cell size.

    Outward, not nearest: a grid that snapped inward would exclude returns that
    sit in the margin, and a point dropped for being one centimetre outside a
    raster is a bug that shows up as a hole much later.
    """
    ox = math.floor(minx / cell) * cell
    oy = math.ceil(maxy / cell) * cell
    ex = math.ceil(maxx / cell) * cell
    ey = math.floor(miny / cell) * cell
    width = int(round((ex - ox) / cell))
    height = int(round((oy - ey) / cell))
    return Grid(ox, oy, cell, width, height)


def grid_for_las(las_path, cell=1.0, pad_cells=0):
    """The canonical grid for one point cloud, from its header alone.

    Reading the header rather than the points means this is instant and gives
    the same answer whether or not the file has been filtered, so a tile's grid
    is fixed before any processing decision is made.
    """
    import json
    import pdal
    pl = pdal.Pipeline(json.dumps({"pipeline": [str(las_path)]}))
    pl.execute()
    md = pl.metadata if isinstance(pl.metadata, dict) else json.loads(pl.metadata)
    las = md["metadata"]["readers.las"]
    if isinstance(las, list):
        las = las[0]
    g = snap_outward(las["minx"], las["miny"], las["maxx"], las["maxy"], cell)
    if pad_cells:
        g = Grid(g.origin_x - pad_cells * cell, g.origin_y + pad_cells * cell,
                 cell, g.width + 2 * pad_cells, g.height + 2 * pad_cells)
    return g


def common_grid(grids):
    """One grid covering several, on the shared cell alignment.

    Used when a transect's tiles are mosaicked: every tile already sits on the
    same lattice, so the union is exact and nothing is resampled.
    """
    if not grids:
        raise ValueError("no grids to combine")
    cell = grids[0].cell
    for g in grids[1:]:
        if not g.aligns_with(grids[0]):
            raise ValueError("grids do not share a lattice; they cannot be "
                             "mosaicked without resampling")
    minx = min(g.bounds[0] for g in grids)
    miny = min(g.bounds[1] for g in grids)
    maxx = max(g.bounds[2] for g in grids)
    maxy = max(g.bounds[3] for g in grids)
    return Grid(minx, maxy, cell,
                int(round((maxx - minx) / cell)),
                int(round((maxy - miny) / cell)))


def grid_from_raster(path):
    """Adopt an existing raster's exact lattice, irregular or not.

    Production output uses `snap_outward`, which puts cell edges on whole
    multiples of the cell size so that neighbouring tiles share a lattice and
    mosaic without resampling. Reference DEMs from the original workflow do not
    have that property: kriging to a data-derived extent in Surfer produced
    cells of 0.500042 m on origins that are multiples of nothing, so no two of
    those tiles align with each other either.

    Comparing against one of them therefore has to be done on its grid. Building
    the candidate here rather than resampling afterwards keeps the comparison
    cell for cell, which is the only way a difference means what it says.
    """
    from osgeo import gdal
    gdal.UseExceptions()
    d = gdal.Open(str(path))
    gt = d.GetGeoTransform()
    if abs(gt[2]) > 1e-9 or abs(gt[4]) > 1e-9:
        raise ValueError("raster is rotated; this assumes north-up")
    return Grid(gt[0], gt[3], abs(gt[1]), d.RasterXSize, d.RasterYSize,
                cell_y=abs(gt[5]))
