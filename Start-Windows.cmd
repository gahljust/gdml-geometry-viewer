@echo off
where py >nul 2>nul
if errorlevel 1 (
  python -B "%~dp0launch.py" %*
) else (
  py -3 -B "%~dp0launch.py" %*
)
if errorlevel 1 pause
