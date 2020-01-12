@echo off
rem set path=Z:\code\third_party\netCDF-4.3.3.1\bin;%pasth%
set path=C:\Program Files\netCDF 4.6.1\bin;%pasth%


rem cd ncfiles
FOR /f %%F IN ('dir /s/b *.nc') DO (
  echo %%F
  ncdump.exe -h %%F > %%F.txt  
)

pause