echo off

REM execute (Call) "ga-aem_vars.bat" batch script to add executables and dependencies to your search path
CALL %GA-AEM_ROOT%\scripts\ga-aem_vars.bat

ctlinedata2georefimage.exe ctlinedata2georefimage-log-stretch.con
ctlinedata2georefimage.exe ctlinedata2georefimage-linear-stretch.con

pause
