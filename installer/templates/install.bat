@echo off
REM ---------------------------------------------------------------------
REM EmulSim 9.2.0 -- Windows 11 x64 installer
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

echo [EmulSim 9.2.0] Installer -- Windows 11 x64
echo [install] Step 1: locating Python.
echo.

REM ---- Locate Python on PATH -----------------------------------------
REM Windows 11 trap: 'python' on PATH may be the Microsoft Store launcher
REM stub (%LOCALAPPDATA%\Microsoft\WindowsApps\python.exe). The stub exits
REM with errorlevel 0 and prints nothing, which used to silently break
REM the version parse below and crash the whole script at the next
REM 'if !PYMAJ! lss 3' (empty-variable syntax error). We now detect
REM an empty version string and give a clear remediation message.
set "PY=python"
where %PY% >nul 2>&1
if errorlevel 1 (
    echo [install] ERROR: 'python' not found on PATH.
    echo          Install Python 3.11 or 3.12 from
    echo          https://www.python.org/downloads/windows/
    echo          (3.13+ is not yet supported -- see ADR-001)
    echo          Tick "Add python.exe to PATH" during setup,
    echo          then re-run this installer.
    pause
    exit /b 1
)

REM Robust version probe via a Python one-liner. Handles RC/dev suffixes,
REM and an empty result tells us we likely hit the Store launcher stub.
set "PYVER="
for /f "delims=" %%v in ('%PY% -c "import sys; print(f'{sys.version_info[0]}.{sys.version_info[1]}.{sys.version_info[2]}')" 2^>nul') do set "PYVER=%%v"

if not defined PYVER (
    echo [install] ERROR: 'python' returned no version string.
    echo.
    echo          This almost always means 'python' on your PATH is the
    echo          Microsoft Store launcher stub at
    echo              %LOCALAPPDATA%\Microsoft\WindowsApps\python.exe
    echo          which opens the Store instead of running Python.
    echo.
    echo          Fix:
    echo          1. Install real Python from
    echo             https://www.python.org/downloads/windows/
    echo             Tick "Add python.exe to PATH" on the first install screen.
    echo          2. Windows Settings -^> Apps -^> Advanced app settings -^>
    echo             App execution aliases, turn OFF "python.exe" and
    echo             "python3.exe".
    echo          3. Open a NEW Command Prompt and run:    python --version
    echo             If it prints a version, re-run this installer.
    pause
    exit /b 1
)

echo [install] Detected Python version: !PYVER!
for /f "tokens=1,2 delims=." %%a in ("!PYVER!") do (
    set "PYMAJ=%%a"
    set "PYMIN=%%b"
)
if not defined PYMAJ goto bad_python
if not defined PYMIN goto bad_python
if !PYMAJ! lss 3 goto bad_python
if !PYMAJ! equ 3 if !PYMIN! lss 11 goto bad_python
REM v9.1.x: pyproject pins requires-python = ">=3.11,<3.13". Reject 3.13+
REM and any major-version >= 4 so the wheel install does not fail later
REM with a confusing pip error. See docs/decisions/ADR-001.
if !PYMAJ! gtr 3 goto bad_python_too_new
if !PYMAJ! equ 3 if !PYMIN! geq 13 goto bad_python_too_new
goto pyok

:bad_python
echo [install] ERROR: Python !PYVER! is below the required 3.11.
echo          Install Python 3.11 or 3.12 from python.org.
pause
exit /b 2

:bad_python_too_new
echo [install] ERROR: Python !PYVER! is above the supported 3.12.
echo          EmulSim requires Python ^>=3.11 and ^<3.13 (ADR-001).
echo          Install Python 3.11 or 3.12 from python.org and re-run.
pause
exit /b 3

:pyok
echo [install] Python check passed.
echo [install] Step 2: parsing install arguments.

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
REM IMPORTANT: this block uses a goto-based flow on purpose. An earlier
REM version nested an "if errorlevel 1 (...)" with several echo lines
REM that contained LITERAL "(" and ")" characters (Windows path examples,
REM "Run as administrator", etc.). cmd.exe pre-parses the entire
REM parenthesised if-block at parse time and counts parens — unescaped
REM parens inside echoed text threw off the counter and the whole file
REM failed to parse with ". was unexpected at this time." ALWAYS use
REM goto + labels for error-message blocks so the echoes never sit
REM inside an if (...) structure.
echo [install] Step 3: creating or re-using virtual environment.
if exist ".venv\Scripts\python.exe" goto venv_ready

echo [install] Creating virtual environment at .venv\
%PY% -m venv .venv
if errorlevel 1 goto err_venv
goto venv_ready

:venv_ready
set "VPY=.venv\Scripts\python.exe"
echo [install] Virtual environment ready: %VPY%
echo [install] Step 4: upgrading pip + wheel inside .venv.

"%VPY%" -m pip install --quiet --upgrade pip wheel
if errorlevel 1 goto err_pip_upgrade

REM ---- Install the wheel --------------------------------------------
echo [install] Step 5: locating wheel.
set "WHEEL="
for %%W in (wheels\emulsim-*-py3-none-any.whl) do set "WHEEL=%%W"
if not defined WHEEL goto err_wheel_missing
if not exist "%WHEEL%" goto err_wheel_missing
echo [install] Found wheel: %WHEEL%

echo [install] Step 6: installing EmulSim wheel with extras %EXTRAS%.
if not defined EXTRAS goto wheel_core
"%VPY%" -m pip install --force-reinstall "%WHEEL%%EXTRAS%"
goto wheel_done
:wheel_core
echo [install] Installing core wheel only, without extras.
"%VPY%" -m pip install --force-reinstall "%WHEEL%"
:wheel_done
if errorlevel 1 goto err_pip_install

REM ---- Verify -------------------------------------------------------
echo [install] Step 7: verifying import + smoke run.
"%VPY%" -c "import emulsim; print('  emulsim', emulsim.__version__, 'imported OK')"
if errorlevel 1 goto err_import

if not "%RUN_TEST%"=="1" goto skip_smoke
if not exist "configs\fast_smoke.toml" goto skip_smoke
echo [install] Running fast smoke pipeline...
"%VPY%" -m emulsim run configs\fast_smoke.toml --quiet
if errorlevel 1 goto warn_smoke
echo [install] Smoke pipeline: OK
goto smoke_done
:warn_smoke
echo [install] WARNING: smoke pipeline returned non-zero.
echo [install]          Install is usable but the pipeline did not complete cleanly.
:skip_smoke
:smoke_done

echo.
echo =====================================================================
echo [install] DONE. EmulSim 9.2.0 is installed.
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
REM UNLESS the installer is running us non-interactively.
if "%NONINTERACTIVE%"=="" pause
endlocal
exit /b 0


REM ======================================================================
REM Error handlers. Kept at file scope (not inside any if-block) so echoed
REM text can contain parens, dots, colons, URLs without tripping cmd.exe's
REM paren-counting parser.
REM ======================================================================

:err_venv
echo.
echo [install] ERROR: Python venv creation failed.
echo           Install directory: %CD%
echo           Common cause: non-writable directory, e.g. C:\Program Files.
echo           Fix: uninstall, then re-install. v9.2.0 per-user install goes
echo           into %%LOCALAPPDATA%%\Programs\EmulSim which needs no admin.
echo           Alternative: right-click install.bat, Run as administrator.
if "%NONINTERACTIVE%"=="" pause
endlocal
exit /b 3

:err_pip_upgrade
echo.
echo [install] ERROR: pip/wheel upgrade failed. Check your network connection
echo           and any corporate proxy / firewall rules for pypi.org.
if "%NONINTERACTIVE%"=="" pause
endlocal
exit /b 4

:err_wheel_missing
echo.
echo [install] ERROR: no EmulSim wheel found.
echo           Expected: wheels\emulsim-^<version^>-py3-none-any.whl
echo           The installer may have been truncated or extracted incorrectly.
echo           Re-download the Setup.exe from
echo           https://github.com/tocvicmeng-prog/Emulsification-Simulator/releases/latest
if "%NONINTERACTIVE%"=="" pause
endlocal
exit /b 5

:err_pip_install
echo.
echo [install] ERROR: pip install of the EmulSim wheel failed.
echo           See the messages above from pip for the exact cause.
echo           Common causes:
echo             - A dependency has no wheel for your Python version.
echo             - pypi.org is unreachable from this machine.
echo             - A required C runtime is missing on this Windows build.
if "%NONINTERACTIVE%"=="" pause
endlocal
exit /b 6

:err_import
echo.
echo [install] ERROR: cannot import emulsim after install.
echo           The wheel installed but Python cannot load it. Check the
echo           pip output above for warnings about missing shared libraries.
if "%NONINTERACTIVE%"=="" pause
endlocal
exit /b 7
