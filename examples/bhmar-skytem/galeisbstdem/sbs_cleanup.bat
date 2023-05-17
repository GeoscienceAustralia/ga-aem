REM batch file for deleting output from different processors
REM DO NOT run this until you have run sbs_sort.bat
cd output
del inversion.output.*.asc
del inversion.output.*.hdr
del inversion.output.*.dfn
del inversion.output.*.csv
del inversion.output.*.i3
del inversion.output.*.log
pause
