@ECHO OFF

REM INSTALL_DIR is the directory for installing the built executables, libraries and examples package after successfully building in the Visual Studio IDE
SET INSTALL_DIR=%LocalAppData%\GA-AEM
REM SET INSTALL_DIR=%LocalAppData%\GA-AEM-vs2022-ide
REM SET INSTALL_DIR=C:\myprograms\GA-AEM
REM SET INSTALL_DIR=..\install-windows-vs2022-ide

MKDIR %INSTALL_DIR%
XCOPY /I/Y/D bin\x64\Release\*.exe %INSTALL_DIR%\bin
XCOPY /I/Y/D bin\x64\Release\*.dll %INSTALL_DIR%\bin
XCOPY /I/Y/D lib\x64\Release\*.lib %INSTALL_DIR%\lib

XCOPY /E/I/Y/D ..\docs %INSTALL_DIR%\docs
XCOPY /E/I/Y/D ..\scripts %INSTALL_DIR%\scripts
XCOPY /E/I/Y/D ..\examples %INSTALL_DIR%\examples
XCOPY /E/I/Y/D ..\matlab %INSTALL_DIR%\matlab
XCOPY /E/I/Y/D ..\python %INSTALL_DIR%\python
XCOPY /-I/Y/D ..\src\gatdaem1d.h %INSTALL_DIR%\include\gatdaem1d.h

pause
