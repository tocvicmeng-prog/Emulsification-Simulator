# EmulSim Windows Installer Sources

This directory holds the Inno Setup build assets that produce
`EmulSim-<version>-Setup.exe` — the one-click Windows 11 x64
installer attached to each GitHub release.

## What is tracked here

| File | Role |
|---|---|
| `EmulSim.iss` | Inno Setup script: installer metadata, EULA page, file layout, Start-Menu / desktop shortcuts, post-install hook, uninstall steps, Python-presence check. |
| `LICENSE_AND_IP.txt` | End-user licence agreement shown on the installer's first page. States that intellectual property rights belong to Holocyte Pty Ltd, that the software is licensed under GPL-3.0, and that the canonical source is the GitHub repository. |
| `build_installer.bat` | Build helper. Rebuilds the wheel, stages runtime assets, and compiles the installer via `ISCC.exe`. |
| `README.md` | This file. |

## What is **not** tracked (gitignored)

- `stage/` — transient build directory, recreated each time
  `build_installer.bat` runs. Contains the runtime assets that get
  compressed into the installer payload: wheel, configs, docs,
  launcher batch files, LICENSE.
- Any build output (`release/*Setup.exe`, wheels in `dist/`, etc.).
  The installer ships as a GitHub Release asset, not a committed
  binary.

## Building from scratch

### Prerequisites

- **Python 3.11+** on `PATH` (for building the wheel).
- **Inno Setup 6** — install with
  `winget install -e --id JRSoftware.InnoSetup` or download from
  <https://jrsoftware.org/isdl.php>.

### Build

From the repo root:

```
installer\build_installer.bat
```

Output: `release\EmulSim-<version>-Setup.exe` (≈ 2.5 MB).

## What the installer does when run

1. **EULA page** — displays `LICENSE_AND_IP.txt` (Holocyte Pty Ltd
   IP, GPL-3.0, GitHub source URL). User must click
   *I accept the agreement* to continue.
2. **README page** — shows the quickstart `README.txt`.
3. **Install location** — default `%LOCALAPPDATA%\Programs\EmulSim`
   (per-user). Picked deliberately so the post-install `install.bat`
   step can create `.venv\` without requiring admin elevation; an
   admin-scoped Program Files install would make the post-install
   step fail with Access Denied on venv creation.
4. **Shortcuts options** — Start-Menu group (default), optional
   desktop shortcut.
5. **Python presence check** — if Python 3.11+ is not on `PATH`,
   the installer offers to open <https://www.python.org/downloads/windows/>
   in the default browser before proceeding.
6. **File extraction** — wheel, configs, docs, launchers, EULA,
   LICENSE, manual PDF.
7. **Post-install step** — runs the bundled `install.bat`, which:
   - creates a local `.venv\` in the install directory,
   - installs the wheel with the `[ui,optimization]` extras,
   - runs a smoke pipeline to verify.
   All bytes stay inside the install directory; no registry,
   system Python, or global state is touched beyond Inno Setup's
   standard uninstall record.
8. **Optional final action** — offer to open the First Edition
   PDF manual immediately.

## What the uninstaller does

- Runs `rmdir /s /q "{app}\.venv"` to purge the virtual env.
- Removes every file the installer placed.
- Removes Start-Menu and desktop shortcuts.
- Leaves your user data and any files outside the install
  directory untouched.

## Release process

1. Tag the feature set and bump `pyproject.toml` + `__init__.py`.
2. Run `installer\build_installer.bat` to produce the `.exe`.
3. `gh release create vX.Y.Z release/EmulSim-X.Y.Z-Setup.exe
    release/EmulSim-X.Y.Z-Windows-x64.zip --title "…"
    --notes-file RELEASE_NOTES.md`.
4. Announce.
