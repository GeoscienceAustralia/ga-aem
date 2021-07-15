@ECHO OFF

SET target=C:\Users\rossc\OneDrive\CSIRO\ga-aem
MKDIR %target%
XCOPY /I/Y/D x64\Release\*.exe %target%\bin

pause
