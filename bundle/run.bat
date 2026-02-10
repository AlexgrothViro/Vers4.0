@echo off
setlocal EnableExtensions EnableDelayedExpansion

REM Pasta onde está este .bat (bundle/)
set "BUNDLE_DIR=%~dp0"
if "%BUNDLE_DIR:~-1%"=="\" set "BUNDLE_DIR=%BUNDLE_DIR:~0,-1%"

where wsl.exe >nul 2>nul
if errorlevel 1 (
  echo [ERRO] WSL nao encontrado no PATH. Instale com: wsl --install
  exit /b 1
)

for /f "tokens=2 delims=:" %%i in ('dism.exe /online /Get-FeatureInfo /FeatureName:Microsoft-Windows-Subsystem-Linux ^| findstr /i "State"') do set "WSL_FEATURE_STATE=%%i"
if defined WSL_FEATURE_STATE (
  set "WSL_FEATURE_STATE=!WSL_FEATURE_STATE: =!"
  if /i not "!WSL_FEATURE_STATE!"=="Enabled" (
    echo [ERRO] O recurso "Windows Subsystem for Linux" nao esta habilitado.
    echo [DICA] Abra o PowerShell como Administrador e execute: wsl --install
    echo [DICA] Reinicie o computador e tente novamente.
    exit /b 1
  )
)

set "WSL_DISTRO=%WSL_DISTRO%"
if not defined WSL_DISTRO set "WSL_DISTRO=Ubuntu"

wsl.exe -l -q | findstr /i /x "%WSL_DISTRO%" >nul
if errorlevel 1 (
  echo [ERRO] Distro WSL "%WSL_DISTRO%" nao encontrada.
  echo [DICA] Distros disponiveis:
  wsl.exe -l -q
  exit /b 1
)

for /f "delims=" %%i in ('wsl.exe -d "%WSL_DISTRO%" wslpath -a "%BUNDLE_DIR%"') do set "WSL_BUNDLE_DIR=%%i"
if not defined WSL_BUNDLE_DIR (
  echo [ERRO] Falha ao converter caminho Windows para WSL: "%BUNDLE_DIR%"
  exit /b 1
)

set "WSL_CMD=cd \"%WSL_BUNDLE_DIR%/..\" && bash bundle/run.sh %*"
wsl.exe -d "%WSL_DISTRO%" bash -lc "%WSL_CMD%"
