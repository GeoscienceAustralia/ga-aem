@echo OFF

REM Add Python3.x installation directory to the PATH variable
set path=%PYTHON3_ROOT%;%PATH%

REM Add the ga-aem python code to the PYTHONPATH variable
set PYTHONPATH=%CD%\..;%PYTHONPATH%
echo %PYTHONPATH%

REM Run the example
python skytem_example.py

pause

