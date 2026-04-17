@echo off
REM ---------------------------------------------------------------------
REM EmulSim 9.0.0 -- Open a Command Prompt with emulsim on PATH (self-heal)
REM Drops into an interactive cmd.exe with the .venv activated, so you
REM can run emulsim subcommands directly:
REM     emulsim run configs\default.toml
REM     emulsim sweep --rpm-min 3000 --rpm-max 15000
REM     emulsim uncertainty --n-samples 50
REM     emulsim design --d32 50e-6 --d32-tol 10e-6
REM Type 'exit' to close.
REM
REM If the virtual environment is missing, this script offers to run
REM install.bat automatically so the first launch "just works".
REM ---------------------------------------------------------------------

setlocal EnableExtensions
cd /d "%~dp0"

REM Check both the venv activation script AND that emulsim is actually
REM installed inside the venv. A partial install (venv created, pip
REM install failed) would otherwise pass the first check and drop the
REM user into a shell without the emulsim command available.
if not exist ".venv\Scripts\activate.bat" goto setup
if not exist ".venv\Scripts\python.exe"  goto setup
".venv\Scripts\python.exe" -c "import emulsim" 1>nul 2>nul
if errorlevel 1 goto setup

:launch
echo [EmulSim 9.0.0] Command-line shell. Type 'exit' to close.
echo                 Example: emulsim run configs\default.toml
echo.
cmd /k ".venv\Scripts\activate.bat && cd /d %~dp0"
endlocal
exit /b %ERRORLEVEL%

:setup
echo [EmulSim 9.0.0] Virtual environment not found at:
echo                 %CD%\.venv
echo.
echo The installer's post-install step appears not to have completed.
echo (Most common cause: Python 3.11+ was not on PATH at install time.)
echo.

where python >nul 2>&1
if errorlevel 1 (
    echo [EmulSim] Python is not on PATH either.
    echo           1. Install Python 3.11 or newer from
    echo              https://www.python.org/downloads/windows/
    echo              (tick "Add python.exe to PATH" during setup).
    echo           2. Re-run launch_cli.bat and setup will continue.
    echo.
    pause
    endlocal
    exit /b 1
)

echo [EmulSim] Python was found:
for /f "delims=" %%v in ('python --version 2^>^&1') do echo           %%v
echo.
echo           Press any key to run setup now (creates .venv and
echo           installs the runtime -- approx. 3-8 minutes), or
echo           close this window to abort.
pause >nul

call "%~dp0install.bat" --no-test
if errorlevel 1 (
    echo.
    echo [EmulSim] Setup failed with error code %ERRORLEVEL%.
    echo           Scroll up to see the cause, or rerun
    echo           install.bat for a full re-install.
    pause
    endlocal
    exit /b %ERRORLEVEL%
)

echo.
echo [EmulSim] Setup finished successfully. Starting the CLI shell...
echo.
goto launch
