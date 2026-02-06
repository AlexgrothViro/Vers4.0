@echo off
setlocal
REM Pasta onde está este .bat (bundle/)
set "BUNDLE_DIR=%~dp0"
if "%BUNDLE_DIR:~-1%"=="\" set "BUNDLE_DIR=%BUNDLE_DIR:~0,-1%"
REM Converte caminho Windows -> WSL
for /f "delims=" %%i in ('wsl wslpath -a "%BUNDLE_DIR%"') do set "WSL_BUNDLE_DIR=%%i"

REM Sobe um nível para a raiz do projeto e roda o run.sh
wsl bash -lc "cd \"%WSL_BUNDLE_DIR%/..\" && bash bundle/run.sh %*"
