@echo off
setlocal EnableDelayedExpansion EnableExtensions

rem -------------------------
rem bundle/start_platform.bat - robust launcher WSL/Windows
rem - determines BUNDLE_DIR via %~dp0
rem - finds free port between 8000 and 8010, warns if 8000 busy
rem - exports PLATFORM_PORT into WSL and runs bundle/run.sh ux
rem - opens browser at selected port
rem -------------------------

rem diretório onde este .bat está (remove trailing backslash)
set "BUNDLE_DIR=%~dp0"
if "%BUNDLE_DIR:~-1%"=="\" set "BUNDLE_DIR=%BUNDLE_DIR:~0,-1%"

rem checar wsl.exe
where wsl.exe >nul 2>nul
if errorlevel 1 (
  echo [ERRO] wsl.exe nao encontrado. Instale com: wsl --install
  exit /b 1
)

rem distro WSL (pode ser sobrescrita por variavel de ambiente)
set "WSL_DISTRO=%WSL_DISTRO%"
if not defined WSL_DISTRO set "WSL_DISTRO=Ubuntu"

wsl.exe -l -q | findstr /i /x "%WSL_DISTRO%" >nul
if errorlevel 1 (
  echo [ERRO] Distro "%WSL_DISTRO%" nao encontrada.
  echo Distros disponiveis:
  wsl.exe -l -q
  exit /b 1
)

rem converter caminho Windows -> WSL (mantendo quoting)
for /f "delims=" %%i in ('wsl.exe -d "%WSL_DISTRO%" wslpath -a "%BUNDLE_DIR%"') do set "WSL_BUNDLE_DIR=%%i"
if not defined WSL_BUNDLE_DIR (
  echo [ERRO] Falha ao converter path Windows para WSL.
  exit /b 1
)

rem procurar porta livre entre 8000 e 8010
set "DEFAULT_PORT=8000"
set "PORT=%DEFAULT_PORT%"
for /L %%P in (8000,1,8010) do (
  netstat -ano -p tcp | findstr /C:":%%P" >nul 2>&1
  if errorlevel 1 (
    set "PORT=%%P"
    goto :PORT_CHOSEN
  )
)
:PORT_CHOSEN

if "%PORT%"=="%DEFAULT_PORT%" (
  echo Usando porta padrao %PORT%.
) else (
  echo Aviso: porta %DEFAULT_PORT% esta em uso. Usando fallback porta %PORT%.
)

rem Exportar a porta para o ambiente WSL via variável PLATFORM_PORT e iniciar o launcher no WSL.
set "WSL_CMD=cd \"%WSL_BUNDLE_DIR%/..\" && PLATFORM_PORT=%PORT% nohup bash bundle/run.sh ux > logs/wsl_ux_boot.log 2>&1 &"
wsl.exe -d "%WSL_DISTRO%" bash -lc "%WSL_CMD%"
if errorlevel 1 (
  echo [ERRO] Falha ao iniciar a plataforma no WSL.
  exit /b 1
)

rem abrir navegador na porta escolhida
start "" "http://localhost:%PORT%"
echo [OK] Plataforma iniciada. Acesse: http://localhost:%PORT%
exit /b 0
