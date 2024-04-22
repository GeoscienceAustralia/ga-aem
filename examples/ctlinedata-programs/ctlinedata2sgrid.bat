echo off

REM Ideally you would set GA-AEM_ROOT as a variable in you user environment (e.g. Start | Edit environment variable for your account)
REM set GA-AEM_ROOT=C:\Users\[YourUserName]\AppData\Local\GA-AEM
REM set GA-AEM_ROOT=%LocalAppData%\GA-AEM

REM Call "ga-aem_vars.bat" batch script to add executables and dependencies to your search path
CALL %GA-AEM_ROOT%\scripts\ga-aem_vars.bat

REM using column numbers to define fields 
ctlinedata2sgrid.exe ctlinedata2sgrid-columns.con

REM using .hdr headerfile and field names to define fields 
REM ctlinedata2sgrid.exe ctlinedata2sgrid-hdr.con

REM using ASEGGDF2 .dfn headerfile and field names to define fields 
REM ctlinedata2sgrid.exe ctlinedata2sgrid-dfn.con


PAUSE





