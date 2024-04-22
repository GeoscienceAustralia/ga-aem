@echo off

REM Ideally you would set GA-AEM_ROOT as a variable in you user environment (e.g. Start | Edit environment variable for your account)
REM set GA-AEM_ROOT=C:\Users\[YourUserName]\AppData\Local\GA-AEM
REM set GA-AEM_ROOT=%LocalAppData%\GA-AEM

REM Call "ga-aem_vars.bat" batch script to add executables and dependencies to your search path
CALL %GA-AEM_ROOT%\scripts\ga-aem_vars.bat

REM Standalone
REM galeisbstdem.exe galeisbstdem.con

REM Use 4 MPI processes
REM mpiexec -np 4 galeisbstdem.exe galeisbstdem.con

REM Use 4 OpenMP threads
galeisbstdem.exe galeisbstdem.con 4

REM If you don't have MPI installed use the non-MPI version
REM Use 4 OpenMP threads 
REM galeisbstdem-nompi.exe galeisbstdem.con 4

pause


