mkdir lines

copy inversion.output.0000.log lines\inversion.output.log
copy inversion.output.0000.hdr lines\inversion.output.hdr
copy inversion.output.0000.dfn lines\inversion.output.dfn
copy inversion.output.*.asc lines\temp.txt

cd lines
sort /REC 65535 temp.txt /O inversion.output.dat
del temp.txt

pause