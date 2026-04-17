@echo off
REM ---------------------------------------------------------------------
REM EmulSim 9.0.0 -- Windows 11 x64 installer
REM
REM Creates a self-contained Python virtual environment in .venv\ next to
REM this script, installs the EmulSim wheel plus the UI and optimisation
REM extras, and runs a fast smoke pipeline to confirm the install works.
REM
REM PREREQUISITE: Python 3.11 or newer, from https://www.python.org/
REM               with "Add python.exe to PATH" ticked during install.
REM
REM USAGE
REM     install.bat                Full install (UI + optimisation extras)
REM     install.bat --core         Minimum core (no UI, no BO/torch)
REM     install.bat --no-opt       Core + UI (no BO/torch)
REM     install.bat --no-test      Skip the smoke-pipeline verification
REM ---------------------------------------------------------------------

setlocal EnableExtensions EnableDelayedExpansion

cd /d "%~dp0"

echo [EmulSim 9.0.0] Installer -- Windows 11 x64
echo.

REM ---- Locate Python on PATH -----------------------------------------
set "PY=python"
where %PY% >nul 2>&1
if errorlevel 1 (
    echo [install] ERROR: 'python' not found on PATH.
    echo          Install Python 3.11+ from
    echo          https://www.python.org/downloads/windows/
    echo          Tick "Add python.exe to PATH" during setup,
    echo          then re-run this installer.
    pause
    exit /b 1
)

%PY% --version
for /f "tokens=2" %%v in ('%PY% --version 2^>^&1') do set "PYVER=%%v"
for /f "tokens=1,2 delims=." %%a in ("%PYVER%") do (
    set "PYMAJ=%%a"
    set "PYMIN=%%b"
)
if !PYMAJ! lss 3 goto bad_python
if !PYMAJ! equ 3 if !PYMIN! lss 11 goto bad_python
goto pyok

:bad_python
echo [install] ERROR: Python !PYVER! is below the required 3.11.
echo          Install Python 3.11, 3.12, or 3.13 from python.org.
pause
exit /b 2

:pyok

REM ---- Parse args ----------------------------------------------------
set "EXTRAS=[ui,optimization]"
set "RUN_TEST=1"

:parse_args
if "%~1"=="" goto args_done
if /i "%~1"=="--core"     set "EXTRAS="
if /i "%~1"=="--no-opt"   set "EXTRAS=[ui]"
if /i "%~1"=="--no-test"  set "RUN_TEST=0"
shift
goto parse_args
:args_done

REM ---- Create / refresh the virtual environment ----------------------
if exist ".venv\Scripts\python.exe" (
    echo [install] Found existing .venv -- re-using.
) else (
    echo [install] Creating virtual environment at .venv\
    %PY% -m venv .venv
    if errorlevel 1 (
        echo.
        echo [install] ERROR: venv creation failed.
        echo.
        echo The install directory
        echo     %CD%
        echo is not writable by your user account. Most common cause:
        echo EmulSim was installed into C:\Program Files\EmulSim which
        echo requires administrator rights, but this post-install script
        echo runs as your regular user.
        echo.
        echo Fix:
        echo   1. Uninstall the current EmulSim (Control Panel, or the
        echo      "Uninstall EmulSim" shortcut in the Start Menu).
        echo   2. Download the latest installer from
        echo      https://github.com/tocvicmeng-prog/Emulsification-Simulator/releases/latest
        echo      Version 9.0.0 and newer install per-user by default
        echo      (into %%LOCALAPPDATA%%\Programs\EmulSim) so no admin
        echo      is required.
        echo.
        echo Alternative (advanced): right-click install.bat and pick
        echo   "Run as administrator" from the install directory.
        echo.
        pause
        exit /b 3
    )
)

set "VPY=.venv\Scripts\python.exe"

echo [install] Upgrading pip + wheel inside .venv
"%VPY%" -m pip install --quiet --upgrade pip wheel
if errorlevel 1 (
    echo [install] ERROR: pip/wheel upgrade failed. Check your network.
    pause
    exit /b 4
)

REM ---- Install the wheel --------------------------------------------
REM Auto-discover the wheel so version bumps don't require script edits.
REM Every shipped release has exactly one emulsim-*-py3-none-any.whl in
REM the wheels\ directory.
set "WHEEL="
for %%W in (wheels\emulsim-*-py3-none-any.whl) do set "WHEEL=%%W"
if "%WHEEL%"=="" (
    echo [install] ERROR: no wheel found in wheels\ directory.
    echo          Expected: wheels\emulsim-^<version^>-py3-none-any.whl
    pause
    exit /b 5
)
if not exist "%WHEEL%" (
    echo [install] ERROR: wheel not found at %WHEEL%
    pause
    exit /b 5
)
echo [install] Found wheel: %WHEEL%

if "%EXTRAS%"=="" (
    echo [install] Installing core wheel only (no UI, no optimisation)
    "%VPY%" -m pip install --force-reinstall "%WHEEL%"
) else (
    echo [install] Installing wheel with extras %EXTRAS%
    "%VPY%" -m pip install --force-reinstall "%WHEEL%%EXTRAS%"
)
if errorlevel 1 (
    echo [install] ERROR: pip install failed. See messages above.
    pause
    exit /b 6
)

REM ---- Verify -------------------------------------------------------
echo.
echo [install] Verifying import...
"%VPY%" -c "import emulsim; print('  emulsim', emulsim.__version__, 'imported OK')"
if errorlevel 1 (
    echo [install] ERROR: cannot import emulsim after install.
    pause
    exit /b 7
)

if "%RUN_TEST%"=="1" (
    if exist "configs\fast_smoke.toml" (
        echo [install] Running fast smoke pipeline...
        "%VPY%" -m emulsim run configs\fast_smoke.toml --quiet
        if errorlevel 1 (
            echo [install] WARNING: smoke run non-zero. Install is OK but
            echo           the pipeline did not complete cleanly. Review log.
        ) else (
            echo [install] Smoke pipeline: OK
        )
    )
)

echo.
echo =====================================================================
echo [install] DONE. EmulSim 9.0.0 is installed.
echo.
echo Program files were extracted from the wheel into:
echo     %CD%\.venv\Lib\site-packages\emulsim\
echo.
echo Next steps:
echo     Launch the web UI:         launch_ui.bat
echo     Command-line pipeline:     launch_cli.bat
echo     Read the manual:           docs\User_Manual_First_Edition.pdf
echo.
echo To remove:                    uninstall.bat
echo =====================================================================
echo.

REM Pause so the user sees success / error before the window closes,
REM UNLESS the installer is running us non-interactively (in which
REM case stdin is redirected and the pause would hang forever).
if "%NONINTERACTIVE%"=="" pause

endlocal
exit /b 0
