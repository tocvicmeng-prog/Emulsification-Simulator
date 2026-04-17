@echo off
REM ---------------------------------------------------------------------
REM Build the EmulSim Windows installer (.exe) with Inno Setup.
REM
REM Prerequisites:
REM   * Python 3.11+ on PATH (for building the wheel)
REM   * Inno Setup 6 (winget install -e --id JRSoftware.InnoSetup)
REM
REM Usage (from repo root):
REM   installer\build_installer.bat
REM
REM Produces:
REM   release\EmulSim-<version>-Setup.exe
REM ---------------------------------------------------------------------

setlocal EnableExtensions EnableDelayedExpansion
cd /d "%~dp0\.."

echo [build-installer] 1/4  Building wheel + sdist
python -m pip install --quiet --upgrade build wheel || exit /b 1
if exist build rmdir /s /q build
if exist dist rmdir /s /q dist
if exist src\emulsim.egg-info rmdir /s /q src\emulsim.egg-info
python -m build --wheel --sdist || exit /b 2

echo [build-installer] 2/4  Staging runtime assets
if exist installer\stage rmdir /s /q installer\stage
mkdir installer\stage\wheels
mkdir installer\stage\configs
mkdir installer\stage\docs

copy /y dist\emulsim-*-py3-none-any.whl installer\stage\wheels\       > nul
copy /y configs\default.toml            installer\stage\configs\       > nul
copy /y configs\fast_smoke.toml         installer\stage\configs\       > nul
copy /y configs\stirred_vessel.toml     installer\stage\configs\       > nul

REM Manual: rebuild the PDF if Markdown is newer.
if exist docs\user_manual\build_pdf.py (
    python docs\user_manual\build_pdf.py > nul
)
copy /y docs\user_manual\polysaccharide_microsphere_simulator_first_edition.pdf ^
        installer\stage\docs\User_Manual_First_Edition.pdf > nul
copy /y docs\user_manual\polysaccharide_microsphere_simulator_first_edition.md ^
        installer\stage\docs\User_Manual_First_Edition.md > nul

REM Bundled launcher + installer scripts from the release tree.
for %%F in (install.bat launch_ui.bat launch_cli.bat uninstall.bat
            README.txt INSTALL.md RELEASE_NOTES.md) do (
    copy /y release\EmulSim-8.3.4-Windows-x64\%%F installer\stage\ > nul
)
copy /y LICENSE installer\stage\LICENSE.txt > nul

echo [build-installer] 3/4  Locating ISCC.exe
set "ISCC=%LOCALAPPDATA%\Programs\Inno Setup 6\ISCC.exe"
if not exist "%ISCC%" set "ISCC=C:\Program Files (x86)\Inno Setup 6\ISCC.exe"
if not exist "%ISCC%" set "ISCC=C:\Program Files\Inno Setup 6\ISCC.exe"
if not exist "%ISCC%" (
    echo [build-installer] ERROR: ISCC.exe not found.
    echo                  Install Inno Setup:
    echo                  winget install -e --id JRSoftware.InnoSetup
    exit /b 3
)

echo [build-installer] 4/4  Compiling installer
"%ISCC%" installer\EmulSim.iss || exit /b 4

echo.
echo [build-installer] DONE.
echo Output: release\EmulSim-8.3.4-Setup.exe
endlocal
exit /b 0
