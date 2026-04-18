@echo off
REM ---------------------------------------------------------------------
REM EmulSim 9.1.1 -- Uninstall
REM Removes the local .venv directory. Leaves installed Python, configs,
REM manual, and the wheel archive alone so you can re-install later.
REM ---------------------------------------------------------------------

setlocal EnableExtensions
cd /d "%~dp0"

echo [EmulSim 9.1.1] Uninstaller
echo.

if not exist ".venv" (
    echo [uninstall] No .venv directory found -- nothing to remove.
    pause
    exit /b 0
)

echo This will delete the .venv directory in
echo     %CD%
echo Your Python install and configs are NOT touched.
echo.
set /p CONFIRM="Proceed? (y/N): "
if /i not "%CONFIRM%"=="y" (
    echo [uninstall] Cancelled.
    exit /b 0
)

echo [uninstall] Removing .venv ...
rmdir /s /q ".venv"
if exist ".venv" (
    echo [uninstall] ERROR: could not fully remove .venv -- close any open
    echo             processes using it and try again.
    pause
    exit /b 2
)

echo [uninstall] Done. To re-install, run install.bat.
endlocal
exit /b 0
