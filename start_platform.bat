@echo off
setlocal EnableDelayedExpansion

rem -------------------------
rem start_platform.bat - robust launcher para Windows/WSL
rem - determina diretório do script em %%~dp0
rem - tenta porta 8000, se ocupada procura 8001..8010 e avisa
rem - exporta PLATFORM_PORT e delega para run_platform.* se existir
rem - fallback: python -m http.server (apenas para facilitar testes/UX)
rem -------------------------

rem diretório do script (com trailing backslash)
set "SCRIPT_DIR=%~dp0"

rem portas
set "DEFAULT_PORT=8000"
set "PORT=%DEFAULT_PORT%"

rem encontrar porta livre entre 8000 e 8010
for /L %%P in (8000,1,8010) do (
  rem verificar se existe uma linha contendo :PORT (TCP)
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

rem exportar para filhos (variavel de ambiente)
set "PLATFORM_PORT=%PORT%"

rem tentar delegar para um launcher local (mantem logica do projeto)
if exist "%SCRIPT_DIR%run_platform.bat" (
  pushd "%SCRIPT_DIR%"
  echo Chamando run_platform.bat com porta %PORT%...
  call "%SCRIPT_DIR%run_platform.bat" "%PORT%"
  popd
  exit /b %ERRORLEVEL%
)

if exist "%SCRIPT_DIR%run_platform.ps1" (
  echo Chamando run_platform.ps1 com porta %PORT%...
  powershell -ExecutionPolicy Bypass -File "%SCRIPT_DIR%run_platform.ps1" "%PORT%"
  exit /b %ERRORLEVEL%
)

rem Fallback: tentar python -m http.server para melhorar UX (NAO altera logica cientifica)
where python >nul 2>&1
if %ERRORLEVEL%==0 (
  pushd "%SCRIPT_DIR%"
  echo Iniciando servidor estatico de fallback: python -m http.server %PORT% (servindo "%SCRIPT_DIR%")...
  start "" "http://localhost:%PORT%/"
  python -m http.server %PORT%
  popd
  exit /b %ERRORLEVEL%
)

echo Nao encontrou launcher local nem Python. Inicie a plataforma manualmente na porta %PORT%.
exit /b 1
