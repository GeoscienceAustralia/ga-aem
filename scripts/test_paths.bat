@echo off

REM Reset the path to the very minimum
set path=C:\Windows\System32

REM Ideally you would set GA-AEM_ROOT as a variable in you environment (e.g. Start | Edit environment variable for your account)
set GA-AEM_ROOT=C:\Users\rossc\AppData\Local\GA-AEM

REM Call batch script to setup the GA-AEM programs and dependency paths
CALL %GA-AEM_ROOT%\scripts\ga-aem_vars.bat

galeisbstdem.exe
ctlinedata2sgrid.exe
ctlinedata2slicegrids.exe


pause


