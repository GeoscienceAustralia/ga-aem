#!/bin/tcsh

if( !(-d lines) ) then
   mkdir lines
endif


sort -n inversion.output.*.asc > lines/inversion.output.dat
cp inversion.output.0000.hdr     lines/inversion.output.hdr
cp inversion.output.0000.dfn     lines/inversion.output.dfn
cp inversion.output.0000.log     lines/inversion.output.log

