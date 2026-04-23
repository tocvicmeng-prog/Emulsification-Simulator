EmulSim 9.1.2 -- Polysaccharide Microsphere Emulsification Simulator
Windows 11 x64 Release -- First Edition
=====================================================================

QUICK START

  1. Install Python 3.11 or 3.12 from
     https://www.python.org/downloads/windows/
     (Python 3.13+ is not yet supported -- see ADR-001 for why.)
     Make sure "Add python.exe to PATH" is ticked.

  2. Double-click  install.bat
     (Takes 3-8 minutes; downloads ~1 GB of scientific Python deps.)

  3. Double-click  launch_ui.bat
     Your default browser opens at http://localhost:8501 with the
     EmulSim dashboard. Click "Manual (PDF)" in the upper-right
     corner for the full First Edition instruction manual.

WHAT IS IN THIS PACKAGE

  install.bat              Installer (creates a local .venv)
  launch_ui.bat            Launch the web UI
  launch_cli.bat           Open a CMD with `emulsim` on PATH
  uninstall.bat            Remove the local .venv
  INSTALL.md               Detailed install and troubleshooting guide
  LICENSE.txt              Software license
  RELEASE_NOTES.md         What's new in 9.1.2
  wheels/                  The emulsim Python wheel
  configs/                 Example TOML configurations
  docs/                    User manual (PDF + Markdown)

SYSTEM REQUIREMENTS

  * Windows 11 x64 (Windows 10 x64 also supported)
  * Python 3.11 or 3.12 (3.13+ currently unsupported)
  * 2 GB of free disk space for the .venv after install
  * 8 GB RAM recommended (4 GB minimum)
  * Internet connection during first install (downloads deps)

SUPPORT

  Read docs\User_Manual_First_Edition.pdf first -- it covers input
  ranges, the four polymer platforms, the full wet-lab protocols,
  troubleshooting, and the physics/chemistry behind every model.

  Report issues to:
  https://github.com/tocvicmeng-prog/Emulsification-Simulator/issues
