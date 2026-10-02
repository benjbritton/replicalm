Bare-earth lidar processing that reproduces a commercial workflow with open
tools. Replicalm turns a raw airborne lidar point cloud into a bare-earth
digital elevation model and the G1 relief visualization, following the NCALM
workflow described in the Estrada-Belli et al. 2025 supplementary material.

The published method depends on two commercial products — TerraScan for ground
classification and Golden Surfer for interpolation and rasterization. Replicalm
replaces both with PDAL, GDAL, NumPy and SciPy, plus the Relief Visualization
Toolbox for the optional image step.

## What is in this release

`Replicalm-1.0.0-setup.exe` — a Windows installer carrying the application and a
packed conda environment, so it needs no existing Python, PDAL or GDAL on the
target machine.

    size    1,211,198,530 bytes (1.21 GB)
    sha256  9590b7cfc064b99b6210c557c2e2bf6f2418340ca71624fea938889d5db325bb
    built   2026-09-24

The source this installer was built from is the commit this tag points at.

## Using it

Run the installer, then either the desktop application or `replicalm-cli` for
batch work. Input is LAS or LAZ; output is a GeoTIFF digital elevation model and,
optionally, the six relief products.

The default method is **Clear**. **Baseline**, the strict translation of the
published workflow, is selectable alongside it. Cell size is derived from
measured ground-return density unless overridden.

## What this version is, and is not

It is a translation of a published method, with every departure from it
recorded. `docs/processing_report.md` states what is the same, what is similar
and what is different. `docs/open_observations.md` carries what is not settled,
including the negative results and the measurements that overturned earlier
conclusions.

It is not a validated instrument. The configuration was locked on terrain chosen
to discriminate between settings, and the open-problems record names the places
where that evidence is thin or where a result rests on a single window.

## Citation

`CITATION.cff` in the repository carries the machine-readable form. Cite the
archived release rather than the repository, so the version is unambiguous.
