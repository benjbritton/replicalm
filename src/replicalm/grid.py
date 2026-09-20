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
