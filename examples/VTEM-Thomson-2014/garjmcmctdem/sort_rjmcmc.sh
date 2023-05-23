#!/bin/tcsh

cd output
if( !(-d lines) ) then
   mkdir lines
endif

sort rjmcmc.*.asc >  lines/rjmcmc.dat
cp   rjmcmc.0000.hdr lines/rjmcmc.hdr
cp   rjmcmc.0000.dfn lines/rjmcmc.dfn
cp   rjmcmc.0000.log lines/rjmcmc.log

