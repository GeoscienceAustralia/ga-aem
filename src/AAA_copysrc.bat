@echo off
del *.cpp
del *.h
del *.c

set src=z:\code\C++\src

xcopy %src%\general_constants.h
xcopy %src%\general_types.h
xcopy %src%\file_formats.h

xcopy %src%\blocklanguage.*
xcopy %src%\fielddefinition.*
xcopy %src%\file_utils.*
xcopy %src%\gaforwardmodeltdem.*
xcopy %src%\galeisbstdem.*
xcopy %src%\garjmcmctdem.*
xcopy %src%\general_utils.*
xcopy %src%\geometry3d.*
xcopy %src%\le.*
xcopy %src%\matrix_ops.*
xcopy %src%\random.*
xcopy %src%\rjmcmc1d.*
xcopy %src%\rjmcmc1dtdeminverter.*
xcopy %src%\tdem1dmodel.*
xcopy %src%\tdemsystem.*
xcopy %src%\vector_utils.*
pause
