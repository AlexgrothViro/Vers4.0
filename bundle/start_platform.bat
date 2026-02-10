@echo off
setlocal EnableExtensions EnableDelayedExpansion

set "BUNDLE_DIR=%~dp0"
if "%BUNDLE_DIR:~-1%"=="\" set "BUNDLE_DIR=%BUNDLE_DIR:~0,-1%"

where wsl.exe >nul 2>nul
if errorlevel 1 (
  echo [ERRO] wsl.exe nao encontrado. Instale com: wsl --install
  exit /b 1
)

set "WSL_DISTRO=%WSL_DISTRO%"
if not defined WSL_DISTRO set "WSL_DISTRO=Ubuntu"

wsl.exe -l -q | findstr /i /x "%WSL_DISTRO%" >nul
if errorlevel 1 (
  echo [ERRO] Distro "%WSL_DISTRO%" nao encontrada.
  echo Distros disponiveis:
  wsl.exe -l -q
  exit /b 1
)

for /f "delims=" %%i in ('wsl.exe -d "%WSL_DISTRO%" wslpath -a "%BUNDLE_DIR%"') do set "WSL_BUNDLE_DIR=%%i"
if not defined WSL_BUNDLE_DIR (
  echo [ERRO] Falha ao converter path Windows para WSL.
  exit /b 1
)

set "WSL_CMD=cd \"%WSL_BUNDLE_DIR%/..\" && nohup bash bundle/run.sh ux > logs/wsl_ux_boot.log 2>&1 &"
wsl.exe -d "%WSL_DISTRO%" bash -lc "%WSL_CMD%"
if errorlevel 1 (
  echo [ERRO] Falha ao iniciar a plataforma no WSL.
  exit /b 1
)

start "" "http://localhost:8000"
echo [OK] Plataforma iniciada. Acesse: http://localhost:8000
