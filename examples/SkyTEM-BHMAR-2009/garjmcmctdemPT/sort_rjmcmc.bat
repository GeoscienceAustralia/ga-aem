cd output

mkdir lines
copy rjmcmc.*.asc    lines\rjmcmc.dat
copy rjmcmc.0000.hdr lines\rjmcmc.hdr
copy rjmcmc.0000.dfn lines\rjmcmc.dfn
copy rjmcmc.0000.log lines\rjmcmc.log

cd lines
sort /REC 65535 rjmcmc.dat /O rjmcmc.dat

pause