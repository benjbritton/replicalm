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

**The G1 recipe.** The composite published as Table 3 of Britton et al. 2025.
Until 1.1 Replicalm called an external script, `GLiHT_rvt.py`, which belonged to
that work and was not bundled -- so on any machine but the author's the image
step did nothing. From 1.1 the recipe is implemented in
`replicalm/visualize.py`, which reproduces it exactly: on a full 2,890 x 15,950
tile the two agree on 100.0000% of pixels, with a maximum difference of zero
levels.

That implementation computes four visualizations through `rvt.vis` and blends
them. The blending arithmetic -- per-layer min-max normalisation, brightness,
contrast, gamma as power(x, 1/gamma), and the normal, difference, overlay and
multiply modes -- is Replicalm's own code following the published table. The
visualizations themselves are RVT's, called as a library and unmodified: no part
of `rvt-py` is copied, patched or vendored, and Replicalm does not reimplement
hillshade, slope, positive openness or sky-view factor.

Replicalm writes only the G1 composite and those four source layers. It does not
produce RVT's own product set and is not an alternative to it; for the full range
of relief visualizations and control over their parameters, use RVT directly.

## Citation

If the visualizations are used in published work, cite RVT:

> Kokalj, Ž., Somrak, M. (2019). Why Not a Single Image? Combining
> Visualizations to Facilitate Fieldwork and On-Screen Mapping.
> Remote Sensing 11(7), 747.
