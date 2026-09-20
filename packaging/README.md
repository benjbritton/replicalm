# Packaging

A Windows installer that wraps a pinned conda environment, rather than a frozen
one-file executable.

## Why not PyInstaller

PyInstaller can freeze this stack and the licensing is fine — its bootloader
exception permits bundled applications under any license, so MIT code stays MIT.
The problem is the native geospatial libraries.

GDAL needs `proj.db` and its `GDAL_DATA` directory at runtime. Frozen builds
routinely ship without them, and the failure arrives *late*: classification and
kriging complete, then the write produces a raster with no coordinate system.
PDAL compounds it — the Python bindings wrap a C++ library whose plugins load
dynamically, which dependency walkers miss.

A packed environment keeps the directory layout those libraries expect. It costs
size (roughly 1 GB installed) and buys a build that works on someone else's
machine.

## Requirements

| tool | why | install |
|---|---|---|
| conda | builds the environment | Miniconda or Anaconda |
| conda-pack | relocates it | `conda install -n base -c conda-forge conda-pack` |
| Inno Setup 6 | compiles the installer | https://jrsoftware.org/isdl.php |

## Build

```powershell
cd packaging
.\build.ps1                    # full build -> dist\Replicalm-1.0.0-setup.exe
.\build.ps1 -SkipEnv           # reuse an existing replicalm environment
.\build.ps1 -SkipInstaller     # stop after staging, for testing
```

Four steps: create the environment from `environment.yml`, `conda-pack` it into
`stage\env`, stage the package and launchers into `stage\app`, then call ISCC on
`replicalm.iss`.

## What gets installed

```
Replicalm\
  env\                  the packed conda environment
  src\replicalm\        the processing modules
  gui.py                the launcher
  Replicalm.cmd         start the GUI (pythonw, no console window)
  replicalm-cli.cmd     the same work from a command line
  LICENSE               MIT
  THIRD-PARTY-NOTICES.md
  licenses\             Apache-2.0.txt and the rest
```

Both `.cmd` files set `GDAL_DATA` and `PROJ_LIB` into the bundled environment
before starting Python, so nothing is required on the user's `PATH` and an
existing conda or OSGeo install on the machine cannot interfere.

`conda-unpack.exe` runs once at install time. It rewrites the absolute paths
baked into the environment at build time; without it every entry point still
points at the build machine. The installer treats its absence as fatal.

## The launcher

`gui.py` is tkinter — bundled with Python, nothing added to the installer, and
this needs four widgets and a log pane. It contains no processing logic; every
decision lives in `replicalm.pipeline`, so any result from the GUI can be
reproduced from the command line. Work runs on a thread because a full tile
takes minutes and a window that may or may not be working is worse than a slow
one that reports what it is doing.

## The optional image step

G1 needs `rvt-py` (Apache 2.0) and `GLiHT_rvt.py`, which implements the
composite published as Table 3 of Britton et al. 2025. The recipe is called, not
copied, so one definition of it exists. Set `REPLICALM_RVT_SCRIPT` to point at
the script; without it the DEM is produced normally and the image step reports
that it was skipped.

## Not redistributable

ArcGIS Pro's Python ships PDAL and GDAL and is what development uses here, but
it contains Esri-licensed components. The installer builds its own environment
from conda-forge for that reason.

## Known issue: Application Control and freshly downloaded binaries

On a workstation running Windows Defender Application Control or AppLocker, a
newly created conda environment can have its native extension modules blocked:

```
ImportError: DLL load failed while importing _cext:
An Application Control policy has blocked this file.
```

Observed here on `kiwisolver._cext`, which matplotlib needs and which the G1
step needs in turn — `rvt.blend` imports `rvt.default`, which imports
`rvt.blend_func`, which imports matplotlib. The DEM path is unaffected: numpy,
scipy, PDAL, GDAL, rasterio and tkinter all load.

It is a policy effect rather than a build fault. The same package imports
normally from an environment that has been installed for a while, and the
blocked file is a legitimate conda-forge artifact. Whatever the local rule is,
it distinguishes newly introduced unsigned binaries from established ones.

Consequences for building and shipping:

- The build machine may not be able to run the G1 step from the freshly built
  environment even though the environment is correct.
- The packed environment is unaffected as a *build* artifact; `conda-pack`
  copies files rather than importing them.
- A target machine under the same policy will need the installed application
  allowlisted. Signing the installer is the durable answer if this is going to
  be distributed inside a managed estate.

Test the G1 step on a machine without the policy, or from an allowlisted
location, before concluding that a build is broken.
