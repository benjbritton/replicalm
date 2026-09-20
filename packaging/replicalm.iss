; Replicalm installer
;
; Wraps a conda-packed environment rather than a frozen one-file executable.
; PyInstaller can freeze this stack, but GDAL and PDAL carry native data
; directories -- proj.db above all -- and a frozen build that loses them fails
; at write time, after the work is done, with an unprojected raster. A packed
; conda environment keeps the directory layout those libraries expect.
;
; Build with build.ps1, which creates the environment, packs it, and calls
; ISCC on this file. Inno Setup 6 must be installed:
;   https://jrsoftware.org/isdl.php

#define AppName      "Replicalm"
#define AppVersion   "1.0.0"
#define AppPublisher "Benjamin Jay Britton"
#define AppURL       "https://github.com/benjbritton/Replicalm"
#define StageDir     "stage"

[Setup]
AppId={{7A2C1E64-5B3D-4F90-9C2A-REPLICALM0001}
AppName={#AppName}
AppVersion={#AppVersion}
AppPublisher={#AppPublisher}
AppPublisherURL={#AppURL}
DefaultDirName={autopf}\{#AppName}
DefaultGroupName={#AppName}
OutputDir=dist
OutputBaseFilename=Replicalm-{#AppVersion}-setup
Compression=lzma2/max
SolidCompression=yes
WizardStyle=modern
; The packed environment is large; refuse to start rather than fail midway.
ExtraDiskSpaceRequired=1200000000
PrivilegesRequiredOverridesAllowed=dialog commandline
ArchitecturesInstallIn64BitMode=x64compatible
LicenseFile=..\LICENSE

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Tasks]
Name: "desktopicon"; Description: "Create a desktop shortcut"; \
    GroupDescription: "Additional shortcuts:"; Flags: unchecked

[Files]
; The conda-packed environment, already unpacked into stage\env by build.ps1.
Source: "{#StageDir}\env\*"; DestDir: "{app}\env"; \
    Flags: ignoreversion recursesubdirs createallsubdirs
; The package, the launcher and the entry points.
Source: "{#StageDir}\app\*"; DestDir: "{app}"; \
    Flags: ignoreversion recursesubdirs createallsubdirs
Source: "..\LICENSE"; DestDir: "{app}"; Flags: ignoreversion
Source: "THIRD-PARTY-NOTICES.md"; DestDir: "{app}"; Flags: ignoreversion
Source: "{#StageDir}\licenses\*"; DestDir: "{app}\licenses"; \
    Flags: ignoreversion recursesubdirs

[Icons]
Name: "{group}\{#AppName}"; Filename: "{app}\Replicalm.cmd"; \
    WorkingDir: "{app}"; IconFilename: "{app}\replicalm.ico"
Name: "{group}\Third-party notices"; Filename: "{app}\THIRD-PARTY-NOTICES.md"
Name: "{group}\Uninstall {#AppName}"; Filename: "{uninstallexe}"
Name: "{autodesktop}\{#AppName}"; Filename: "{app}\Replicalm.cmd"; \
    WorkingDir: "{app}"; IconFilename: "{app}\replicalm.ico"; \
    Tasks: desktopicon

[Run]
; conda-unpack rewrites the absolute paths baked into the environment at build
; time. Without it every script entry point still points at the build machine.
Filename: "{app}\env\Scripts\conda-unpack.exe"; \
    StatusMsg: "Configuring the environment for this machine..."; \
    Flags: runhidden waituntilterminated
Filename: "{app}\Replicalm.cmd"; Description: "Launch {#AppName}"; \
    Flags: postinstall nowait skipifsilent

[UninstallDelete]
Type: filesandordirs; Name: "{app}\env"
Type: filesandordirs; Name: "{app}\src"
