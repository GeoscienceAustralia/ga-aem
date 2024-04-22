@echo off

REM Ideally you would set GA-AEM_ROOT as a variable in you user environment (e.g. Start | Edit environment variable for your account)
REM set GA-AEM_ROOT=C:\Users\[YourUserName]\AppData\Local\GA-AEM
REM set GA-AEM_ROOT=%LocalAppData%\GA-AEM

REM Call "ga-aem_vars.bat" batch script to add executables and dependencies to your search path
CALL %GA-AEM_ROOT%\scripts\ga-aem_vars.bat

REM Run standalone on a single CPU
REM galeisbstdem.exe galeisbstdem.con

REM Run using 4 OpenMP threads
REM galeisbstdem.exe galeisbstdem.con 4	

REM Use 4 MPI processes
mpiexec -np 4 galeisbstdem.exe galeisbstdem.con	

REM If you do not have MPI installed substitute galeisbstdem-nompi.exe in place of galeisbstdem.exe
REM galeisbstdem-nompi.exe galeisbstdem.con

pause
