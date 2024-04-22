#!/bin/tcsh

cd output
sort -n inversion.output.*.asc > inversion.output.dat
cp inversion.output.0000.hdr     inversion.output.hdr
cp inversion.output.0000.dfn     inversion.output.dfn
cp inversion.output.0000.log     inversion.output.log

