@echo off
REM ---------------------------------------------------------------------
REM EmulSim 9.2.2 -- Launch the Streamlit web UI (self-healing, diag-safe)
REM Opens http://localhost:8501 in your default browser.
REM Close the terminal window to stop the server.
REM
REM If the virtual environment is missing (e.g., the installer's post-
REM install step was skipped or failed quietly), this script offers to
REM run install.bat automatically so the first launch "just works".
REM
REM If the UI itself exits abnormally, the window is kept open with a
REM diagnostic block so the Python/Streamlit traceback is visible.
REM ---------------------------------------------------------------------

setlocal EnableExtensions
cd /d "%~dp0"

REM The .venv may exist as a partial directory if an earlier install.bat
REM run failed between venv creation and pip install. Check BOTH that
REM python.exe is present AND that the emulsim package is actually
REM importable; either miss drops into the self-heal setup branch.
if not exist ".venv\Scripts\python.exe" goto setup
".venv\Scripts\python.exe" -c "import emulsim" 1>nul 2>nul
if errorlevel 1 goto setup

:launch
echo [EmulSim 9.2.2] Starting the web UI. Your default browser should open
echo                 automatically at http://localhost:8501.
echo                 Close this window to stop the server.
echo.
".venv\Scripts\python.exe" -m emulsim ui
set "EXITCODE=%ERRORLEVEL%"
if not "%EXITCODE%"=="0" (
    echo.
    echo ================================================================
    echo [EmulSim] The UI exited unexpectedly with error code %EXITCODE%.
    echo.
    echo Scroll up to read the Python / Streamlit error output that
    echo appeared before this box. That is the real cause.
    echo.
    echo Common causes and fixes:
    echo   - Port 8501 already in use: close the other app and retry.
    echo   - Dependency import error: re-run install.bat to reinstall.
    echo   - Python version incompatibility: ensure Python 3.11 or 3.12.
    echo.
    echo For manual diagnosis, open a Command Prompt and run:
    echo   cd /d "%%LOCALAPPDATA%%\Programs\EmulSim"
    echo   .venv\Scripts\python.exe -m emulsim ui
    echo ================================================================
    pause
)
endlocal
exit /b %EXITCODE%

:setup
echo [EmulSim 9.2.2] Virtual environment not found at:
echo                 %CD%\.venv
echo.
echo This usually means the installer's post-install step was skipped
echo or could not complete (typically because a supported Python (3.11 or 3.12)
echo was not on PATH when the installer ran).
echo.

where python >nul 2>&1
if errorlevel 1 (
    echo [EmulSim] Python is not on PATH either.
    echo           1. Install Python 3.11 or newer from
    echo              https://www.python.org/downloads/windows/
    echo              (tick "Add python.exe to PATH" during setup).
    echo           2. Re-run launch_ui.bat and setup will continue.
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
echo [EmulSim] Setup finished successfully. Starting the UI...
echo.
goto launch
