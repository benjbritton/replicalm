"""Is the straight vertical cut the lattice boundary rather than a data edge?

Column 42 at 0.5 m is exactly 21 m from column 0 -- the strip width. If the
footprint reaches the raster's own left border in those rows, then the "first
valid cell" the strip measures from is the edge of the grid, not the edge of
the survey, and 21 m of good data is removed along a perfectly straight line
for a seam that is not there.
"""
import numpy as np
import rasterio

E2E = (r"C:\Replicalm\test100526\centro454_e2e"
       r"\AMIGACarb_Yuc_Centro_GLAS_Apr2013_l0s454_DEM.tif")
PRE = (r"C:\Replicalm\test100526\centro454_nofix50"
       r"\AMIGACarb_Yuc_Centro_GLAS_Apr2013_l0s454_DEM.tif")
ARCH = (r"D:\_Archive_EdgeFixed\Yuc_Centro\Centro_Glas\Centro_Glas"
        r"\Centro_GLAS_l0s454_DEM_0p5m_v1.tif")


def rd(p):
    with rasterio.open(p) as r:
        z = r.read(1).astype("f8")
        nd, cc = r.nodata, abs(r.transform.a)
    if nd is not None:
        z[z == nd] = np.nan
    z[z > 1e30] = np.nan
    z[z == 0] = np.nan
    return np.isfinite(z), cc


M, cell = rd(E2E)
P, _ = rd(PRE)
A, _ = rd(ARCH)
H, W = M.shape
strip_cells = int(round(21.0 / cell))
print("strip = %d cells = %.1f m; lattice is %d columns wide"
      % (strip_cells, strip_cells * cell, W))

touch_w = P[:, 0]
touch_e = P[:, W - 1]
print("")
print("ROWS WHOSE PRE-STRIP COVERAGE REACHES THE LATTICE BORDER")
print("  west border (column 0)      %5d rows of %d (%.1f%%)"
      % (int(touch_w.sum()), H, 100.0 * touch_w.mean()))
print("  east border (column %d)   %5d rows (%.1f%%)"
      % (W - 1, int(touch_e.sum()), 100.0 * touch_e.mean()))
print("  archive reaches west border %5d rows" % int(A[:, 0].sum()))
print("  archive reaches east border %5d rows" % int(A[:, W - 1].sum()))

wedge = np.flatnonzero(M.any(axis=1))
west = np.array([np.flatnonzero(M[r])[0] for r in wedge])
at_strip = west == strip_cells
print("")
print("  rows whose final west edge sits exactly at the strip width: %d"
      % int(at_strip.sum()))
rows_cut = wedge[at_strip]
if len(rows_cut):
    print("  of those, %d had coverage at column 0 before the strip (%.1f%%)"
          % (int(touch_w[rows_cut].sum()),
             100.0 * touch_w[rows_cut].mean()))
    print("  they run rows %d to %d" % (int(rows_cut.min()),
                                        int(rows_cut.max())))

lost = int((P[:, :strip_cells] & ~M[:, :strip_cells]).sum())
lost_real = int((P[:, :strip_cells] & ~M[:, :strip_cells]
                 & A[:, :strip_cells]).sum())
print("")
print("CELLS REMOVED INSIDE THE FIRST %d COLUMNS" % strip_cells)
print("  removed by the strip            %8d" % lost)
print("  of those, present in the archive %8d" % lost_real)
print("  i.e. real data the archive kept and we cut: %.1f%% of the band"
      % (100.0 * lost_real / max(lost, 1)))

# same question on the east side
loste = int((P[:, W - strip_cells:] & ~M[:, W - strip_cells:]).sum())
loste_r = int((P[:, W - strip_cells:] & ~M[:, W - strip_cells:]
               & A[:, W - strip_cells:]).sum())
print("  east band: removed %d, of which archive kept %d" % (loste, loste_r))

print("")
print("WHERE WE CUT DATA THE ARCHIVE KEPT, ANYWHERE")
cut_real = P & ~M & A
print("  %d cells (%.2f%% of the archive)"
      % (int(cut_real.sum()), 100.0 * cut_real.sum() / A.sum()))
rr = np.flatnonzero(cut_real.any(axis=1))
if len(rr):
    print("  spread over rows %d to %d" % (int(rr.min()), int(rr.max())))
    per_row = cut_real.sum(axis=1)
    worst = np.argsort(per_row)[::-1][:5]
    for w in worst:
        print("    row %5d: %4d cells" % (int(w), int(per_row[w])))
