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

wsl.exe --status >nul 2>nul
if errorlevel 1 (
  echo [ERRO] O recurso WSL parece desabilitado neste Windows.
  echo         Habilite "Subsistema do Windows para Linux" e reinicie a maquina.
  echo         Guia rapido (PowerShell como Admin): wsl --install
  exit /b 1
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
