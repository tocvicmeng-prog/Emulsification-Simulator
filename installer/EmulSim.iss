; EmulSim 8.3.2 -- Windows 11 x64 Inno Setup installer
; Build with: ISCC.exe installer\EmulSim.iss
; Produces:  release\EmulSim-8.3.2-Setup.exe

#define MyAppName       "EmulSim"
#define MyAppVersion    "8.3.2"
#define MyAppPublisher  "Holocyte Pty Ltd"
#define MyAppURL        "https://github.com/tocvicmeng-prog/Emulsification-Simulator"
#define MyAppExeLauncher "launch_ui.bat"

[Setup]
AppId={{9C7E0DFB-2A95-4F31-AE4B-61F7A9F1A2C0}
AppName={#MyAppName}
AppVersion={#MyAppVersion}
AppVerName={#MyAppName} {#MyAppVersion}
AppPublisher={#MyAppPublisher}
AppPublisherURL={#MyAppURL}
AppSupportURL={#MyAppURL}/issues
AppUpdatesURL={#MyAppURL}/releases
DefaultDirName={autopf}\{#MyAppName}
DefaultGroupName={#MyAppName}
DisableProgramGroupPage=auto
LicenseFile=LICENSE_AND_IP.txt
InfoBeforeFile=stage\README.txt
OutputDir=..\release
OutputBaseFilename=EmulSim-{#MyAppVersion}-Setup
SetupIconFile=
Compression=lzma2/max
SolidCompression=yes
WizardStyle=modern
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
PrivilegesRequired=lowest
PrivilegesRequiredOverridesAllowed=dialog
UninstallDisplayIcon={app}\docs\User_Manual_First_Edition.pdf
UninstallDisplayName={#MyAppName} {#MyAppVersion}
MinVersion=10.0.17763

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Tasks]
Name: "desktopicon";   Description: "Create a &desktop shortcut"; GroupDescription: "Additional shortcuts:"
Name: "startmenuicon"; Description: "Create a &Start Menu group";  GroupDescription: "Additional shortcuts:"; Flags: checkedonce

[Files]
Source: "stage\install.bat";            DestDir: "{app}"; Flags: ignoreversion
Source: "stage\launch_ui.bat";          DestDir: "{app}"; Flags: ignoreversion
Source: "stage\launch_cli.bat";         DestDir: "{app}"; Flags: ignoreversion
Source: "stage\uninstall.bat";          DestDir: "{app}"; Flags: ignoreversion
Source: "stage\README.txt";             DestDir: "{app}"; Flags: ignoreversion
Source: "stage\INSTALL.md";             DestDir: "{app}"; Flags: ignoreversion
Source: "stage\RELEASE_NOTES.md";       DestDir: "{app}"; Flags: ignoreversion
Source: "stage\LICENSE.txt";            DestDir: "{app}"; Flags: ignoreversion
Source: "LICENSE_AND_IP.txt";           DestDir: "{app}"; Flags: ignoreversion
Source: "stage\wheels\*";               DestDir: "{app}\wheels";  Flags: ignoreversion recursesubdirs createallsubdirs
Source: "stage\configs\*";              DestDir: "{app}\configs"; Flags: ignoreversion recursesubdirs createallsubdirs
Source: "stage\docs\*";                 DestDir: "{app}\docs";    Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
Name: "{group}\EmulSim (Web UI)"; Filename: "{app}\launch_ui.bat"; WorkingDir: "{app}"; Comment: "Launch the EmulSim web dashboard"; Tasks: startmenuicon
Name: "{group}\EmulSim (Command line)"; Filename: "{app}\launch_cli.bat"; WorkingDir: "{app}"; Comment: "Open a Command Prompt with emulsim on PATH"; Tasks: startmenuicon
Name: "{group}\User Manual (PDF)"; Filename: "{app}\docs\User_Manual_First_Edition.pdf"; WorkingDir: "{app}\docs"; Comment: "First Edition instruction manual"; Tasks: startmenuicon
Name: "{group}\Release Notes"; Filename: "{app}\RELEASE_NOTES.md"; WorkingDir: "{app}"; Tasks: startmenuicon
Name: "{group}\{cm:UninstallProgram,{#MyAppName}}"; Filename: "{uninstallexe}"; Tasks: startmenuicon
Name: "{autodesktop}\{#MyAppName} (Web UI)"; Filename: "{app}\launch_ui.bat"; WorkingDir: "{app}"; Comment: "Launch the EmulSim web dashboard"; Tasks: desktopicon

[Run]
; Post-install: set up the Python venv + pip-install the wheel. The
; existing install.bat does all of this and prints progress; running
; it inside the installer gives users one-click end-to-end setup.
Filename: "{app}\install.bat"; Description: "Finish setup (create Python environment and install runtime)"; Flags: postinstall shellexec skipifsilent
; Offer to open the manual at the end
Filename: "{app}\docs\User_Manual_First_Edition.pdf"; Description: "Open the User Manual (PDF)"; Flags: postinstall shellexec skipifsilent unchecked nowait

[UninstallRun]
; Let our uninstall.bat drop the .venv BEFORE Inno Setup wipes files.
Filename: "{cmd}"; Parameters: "/c rmdir /s /q ""{app}\.venv"""; Flags: runhidden; RunOnceId: "PurgeVenv"

[Code]
function IsPythonOnPath(): Boolean;
var
  ResultCode: Integer;
begin
  Result := Exec(ExpandConstant('{cmd}'), '/c where python',
                 '', SW_HIDE, ewWaitUntilTerminated, ResultCode)
            and (ResultCode = 0);
end;

procedure InitializeWizard();
begin
  { Leave default; the Python-availability check happens at install time. }
end;

function NextButtonClick(CurPageID: Integer): Boolean;
var
  ErrCode: Integer;
begin
  Result := True;
  if CurPageID = wpReady then
  begin
    if not IsPythonOnPath() then
    begin
      if MsgBox(
        'Python 3.11 or newer was not found on the system PATH.' + #13#10#13#10 +
        'EmulSim needs Python to run. You can:' + #13#10 +
        '  - Click YES to open the Python download page now.' + #13#10 +
        '    Install Python (tick "Add python.exe to PATH"),' + #13#10 +
        '    then re-run this installer.' + #13#10 +
        '  - Click NO to continue installing anyway. The program' + #13#10 +
        '    will refuse to launch until Python is available.',
        mbConfirmation, MB_YESNO) = IDYES then
      begin
        ShellExec('', 'https://www.python.org/downloads/windows/',
                  '', '', SW_SHOW, ewNoWait, ErrCode);
        Result := False;
      end;
    end;
  end;
end;
