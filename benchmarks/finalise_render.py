import json, os, sys
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import ROOT, SRC, DEM_ROOT, LAS_ROOT, TESTS
sys.path.insert(0, SRC)
import numpy as np
from osgeo import gdal
from replicalm import finalise, grid as G, interpolate
gdal.UseExceptions()

SRC = (ROOT + r"\\render\dem\Yuc_South\South_Glas\South_Glas"
       r"\South_GLAS_l0s395_DEM_0p5m_v1.tif")
DST = (ROOT + r"\\render\dem_final\Yuc_South\South_Glas\South_Glas"
       r"\South_GLAS_l0s395_DEM_0p5m_v1.tif")
LAS = (LAS_ROOT + r"\\Yuc_South\South_Glas\South_Glas"
       r"\AMIGACarb_Yuc_South_GLAS_Apr2013_l0s395.las")
RADIUS = float(sys.argv[1]) if len(sys.argv) > 1 else 14.4
os.makedirs(os.path.dirname(DST), exist_ok=True)

g = G.grid_from_raster(SRC)
ds = gdal.Open(SRC)
dem = ds.GetRasterBand(1).ReadAsArray().astype("f8")
wkt = ds.GetProjection() or interpolate.source_srs(LAS)
print("source %s" % g.describe())
rep = finalise.finalise(dem, g, wkt, DST, radius_m=RADIUS, nodata_in=-9999.0)
print(json.dumps(rep, indent=1))
