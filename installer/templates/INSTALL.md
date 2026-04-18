# EmulSim 9.1.1 — Windows 11 x64 Installation Guide

## 1. System Requirements

| Item | Supported |
|---|---|
| Operating System | Windows 11 x64 (primary); Windows 10 x64 (supported) |
| Python | 3.11 or 3.12 (pinned in pyproject.toml: `>=3.11,<3.13`; Python 3.10 and earlier are refused, Python 3.13+ is currently unsupported) |
| Disk | ~2 GB free for `.venv` after install |
| RAM | 8 GB recommended (4 GB minimum) |
| Internet | Required during first install (downloads scientific Python deps) |

## 2. Prerequisite: Install Python

1. Go to <https://www.python.org/downloads/windows/>.
2. Download the latest Python 3.11.x or 3.12.x **64-bit** installer.
   (Python 3.13+ is not yet supported; see ADR-001.)
3. Run the installer.
4. **Tick the box "Add python.exe to PATH"** before clicking Install.
5. After install completes, open a new Command Prompt and run
   `python --version` to confirm Python is on PATH.

If you skip the PATH step the EmulSim installer will not find Python
and will exit with an error directing you to re-install.

## 3. Install EmulSim

### 3.1 Default install (recommended)

Double-click `install.bat`. The installer will:

1. Verify a supported Python (3.11 or 3.12) is on PATH.
2. Create a local virtual environment at `.venv\` inside this folder.
3. Install the EmulSim wheel with the **UI** and **optimisation** extras
   (numpy, scipy, matplotlib, streamlit, plotly, PyTorch, BoTorch).
4. Run a fast smoke pipeline to confirm the install works.

Total time: 3–8 minutes on a typical broadband connection.

### 3.2 Install variants

Open a Command Prompt, `cd` to this folder, and run:

| Command | Effect |
|---|---|
| `install.bat` | Full install (UI + BO). Default. |
| `install.bat --no-opt` | UI only, no PyTorch / BoTorch. Faster, smaller. |
| `install.bat --core` | Core wheel only. No UI, no BO. Scripting only. |
| `install.bat --no-test` | Skip the smoke-pipeline verification. |

Flags can be combined: `install.bat --no-opt --no-test`.

### 3.3 What the installer creates

- `.venv\` — the virtual environment. Everything EmulSim needs lives
  here. Your system Python is left untouched.
- `configs\` — example TOML configurations (already shipped; the
  installer does not modify them).
- No registry entries, no Start-Menu shortcuts, no services. Fully
  self-contained.

## 4. First Run

### 4.1 Launch the web UI

Double-click `launch_ui.bat`. A command window opens and reports

```
[EmulSim 9.1.1] Launching web UI at http://localhost:8501
```

Your default browser should open automatically at that URL. Close
the command window to stop the server.

Click **Manual (PDF)** in the upper-right corner of the dashboard
to download the First Edition manual (or open the copy at
`docs\User_Manual_First_Edition.pdf`).

### 4.2 Use the CLI

Double-click `launch_cli.bat`. A Command Prompt opens with the
EmulSim virtual environment activated. Try:

```
emulsim run configs\default.toml
emulsim sweep --rpm-min 3000 --rpm-max 15000 --rpm-steps 5
emulsim uncertainty --n-samples 50
emulsim info
```

`emulsim --help` lists every subcommand.

### 4.3 Programmatic use

```
.venv\Scripts\python.exe
>>> import emulsim
>>> emulsim.__version__
'9.1.1'
>>> r = emulsim.run_pipeline()
>>> r.emulsification.d32
```

## 5. Removing EmulSim

Run `uninstall.bat`. This deletes the local `.venv\` directory only;
the rest of the release tree (wheel, configs, docs) is left in place
so you can reinstall later without re-downloading.

If you also want to remove the release tree, just delete the
`EmulSim-9.1.1-Windows-x64` folder. No other files are touched.

## 6. Troubleshooting

### "python not found on PATH"

You did not tick the "Add python.exe to PATH" checkbox. Re-run the
Python installer, choose "Modify", and tick the box — or uninstall
Python and start over.

### "Python 3.10 is below the required 3.11" / "Python 3.13 is above the supported range"

EmulSim 9.1.1 requires Python **3.11 or 3.12** (pinned by `requires-python =
">=3.11,<3.13"` in pyproject.toml; see `docs/decisions/ADR-001-python-version-policy.md`
for the why). Install a supported version alongside your existing one. If
multiple versions are on PATH the installer picks the first one — verify with
`python --version` before running install.bat.

### "pip install failed" during install.bat

The most common cause is a network interruption. Re-run
`install.bat`. pip will resume from cache for already-downloaded
packages.

If you are behind a corporate proxy, set `HTTP_PROXY` and
`HTTPS_PROXY` environment variables before running `install.bat`.

### Antivirus flags the .bat files

Batch installers that create `.venv` directories sometimes trigger
heuristic AV scanners. The four `.bat` files in this release are
plain text — open them in Notepad to inspect. Add an exclusion for
the release folder if needed.

### UI starts but shows a blank page

First load can take 5–20 seconds while Streamlit imports the
matplotlib / scipy graph. Wait, then refresh the browser.

### CUDA / GPU torch download

The `--no-opt` install avoids the PyTorch / BoTorch download. Use
this if you only need the forward simulator and not the inverse-
design BO feature. You can re-run `install.bat` later with the full
extras when you want BO.

## 7. Offline install (air-gapped system)

If the target machine has no Internet:

1. On an internet-connected machine with matching Python and
   Windows architecture, run
   ```
   pip download emulsim[ui,optimization] -d wheels_cache
   ```
   to download all deps into `wheels_cache\`.
2. Copy `wheels_cache\` onto the air-gapped machine alongside this
   release.
3. Run
   ```
   install.bat --no-test
   ```
   but first edit it to pass `--find-links wheels_cache --no-index`
   to the pip install step. (Contact support if you need a
   pre-built offline bundle.)

## 8. Verifying the install

The final step of `install.bat` runs

```
python -m emulsim run configs\fast_smoke.toml --quiet
```

which exercises the full L1 → L2 → L3 → L4 pipeline with a small
grid. A successful exit code (0) means the install is healthy. The
generated `output/` folder contains the smoke-run artefacts and can
be deleted if you don't need them.
