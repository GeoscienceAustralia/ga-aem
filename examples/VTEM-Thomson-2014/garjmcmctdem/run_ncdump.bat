@echo off

REM This script uses NetCDF's ncdump program to dump out a text file header for each of the probability map files which are stored in stored in NetCDF format.

REM Ideally you would set GA-AEM_ROOT as a variable in you user environment (e.g. Start | Edit environment variable for your account)
REM set GA-AEM_ROOT=C:\Users\[YourUserName]\AppData\Local\GA-AEM
REM set GA-AEM_ROOT=%LocalAppData%\GA-AEM

REM Call "ga-aem_vars.bat" batch script to add executables and dependencies to your search path
CALL %GA-AEM_ROOT%\scripts\ga-aem_vars.bat

cd output\pmaps
FOR /f %%F IN ('dir /s/b *.nc') DO (
  echo %%F
  ncdump.exe -h %%F > %%F.txt  
)

pause