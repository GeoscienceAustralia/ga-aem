echo off

REM execute (Call) "ga-aem_vars.bat" batch script to add executables and dependencies to your search path
CALL %GA-AEM_ROOT%\scripts\ga-aem_vars.bat

REM using column numbers to define fields 
ctlinedata2curtainimage.exe ctlinedata2curtainimage-log-stretch-columns.con

REM using .hdr headerfile and field names to define fields 
REM ctlinedata2curtainimage.exe ctlinedata2curtainimage-log-stretch-hdr.con

REM using ASEGGDF2 .dfn headerfile and field names to define fields 
REM ctlinedata2curtainimage.exe ctlinedata2curtainimage-log-stretch-dfn.con

pause
