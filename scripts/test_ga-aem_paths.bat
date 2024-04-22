@echo off

REM Ideally you would set GA-AEM_ROOT as an environment variable in your user environment (e.g. Start | Edit environment variable for your account). After setting, be sure to open a fresh command window or Windows Explorer window.
set GA-AEM_ROOT=%LocalAppData%\GA-AEM

REM Call batch script to setup the GA-AEM programs and dependency paths
CALL %GA-AEM_ROOT%\scripts\ga-aem_vars.bat

ECHO ==========================================================================
ECHO GA-AEM_ROOT=%GA-AEM_ROOT%
ECHO GA-AEM_PATH=%GA-AEM_PATH%
ECHO ==========================================================================
ECHO PATH=%PATH%
ECHO ==========================================================================
galeisbstdem.exe
ECHO ==========================================================================
galeiallatonce.exe
ECHO ==========================================================================
ctlinedata2curtainimage.exe
ECHO ==========================================================================
ctlinedata2slicegrids.exe
ECHO ==========================================================================

pause


