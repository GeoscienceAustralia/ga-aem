REM batch file for sorting outputs from different processors into one file
cd output
copy inversion.output.0000.log inversion.output.log
copy inversion.output.0000.hdr inversion.output.hdr
copy inversion.output.0000.dfn inversion.output.dfn
copy inversion.output.0000.csv inversion.output.csv
copy inversion.output.0000.i3  inversion.output.i3
copy inversion.output.*.asc temp.txt
sort /REC 65535 temp.txt /O inversion.output.dat
del temp.txt
pause
