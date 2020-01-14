@echo off
set path=C:\Program Files\netCDF 4.6.1\bin;%pasth%

cd output\pmaps
FOR /f %%F IN ('dir /s/b *.nc') DO (
  echo %%F
  ncdump.exe -h %%F > %%F.txt  
)

pause