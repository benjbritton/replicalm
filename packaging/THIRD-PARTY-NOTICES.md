# Third-party notices

Replicalm is MIT licensed. The installer redistributes the components below,
each under its own license. Nothing here is modified; all are used as published.

| component | license | holder |
|---|---|---|
| PDAL | BSD 3-Clause | Howard Butler and PDAL contributors |
| GDAL / OGR | MIT / X11-style | Open Source Geospatial Foundation |
| PROJ | MIT | PROJ contributors |
| NumPy | BSD 3-Clause | NumPy Developers |
| SciPy | BSD 3-Clause | SciPy Developers |
| matplotlib | matplotlib (BSD-compatible) | Matplotlib Development Team |
| Python | PSF License 2.0 | Python Software Foundation |
| rvt-py 2.2.1 | Apache License 2.0 | ZRC SAZU and University of Ljubljana (UL FGG) |

## Apache 2.0 components

`rvt-py` is distributed under the Apache License, Version 2.0. A full copy of
that license is installed as `licenses/Apache-2.0.txt`. Its source is at
https://github.com/EarthObservation/RVT_py. It is redistributed unmodified.

The Apache 2.0 obligations met here: the license text accompanies the binary
distribution, copyright and attribution notices are retained, and no
modifications have been made that would require a change statement. The
distribution contains no NOTICE file, so none is propagated.

## Not redistributed

**ArcGIS Pro's Python environment.** It ships PDAL and GDAL, and it is used for
development on the author's machine, but it contains Esri-licensed components
and is not redistributable. The installer builds its own conda environment from
conda-forge for exactly this reason.

**G-LiHT data.** The point clouds and reference surfaces are NASA products from
the AMIGACarb and Yucatan campaigns and are not included.

**The G1 recipe.** `GLiHT_rvt.py` implements the composite published as Table 3
of Britton et al. 2025 and belongs to that work. Replicalm calls it if it is
present and reports its absence otherwise; it is not bundled.

## Citation

If the visualizations are used in published work, cite RVT:

> Kokalj, Ž., Somrak, M. (2019). Why Not a Single Image? Combining
> Visualizations to Facilitate Fieldwork and On-Screen Mapping.
> Remote Sensing 11(7), 747.
