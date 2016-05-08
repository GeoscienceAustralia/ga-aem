#!/bin/tcsh
clear

set SRC="./src"
set SRCFILES=""
set SRCFILES="$SRCFILES $SRC/tdaem1d_shrlib_tester.cpp"
set SRCFILES="$SRCFILES $SRC/tdaem1d_shrlib.cpp"
set SRCFILES="$SRCFILES $SRC/rb_utils.cpp"
set SRCFILES="$SRCFILES $SRC/directoriesandfiles.cpp"
set SRCFILES="$SRCFILES $SRC/blocklanguage.cpp"
set SRCFILES="$SRCFILES $SRC/geometry3d.cpp"
set SRCFILES="$SRCFILES $SRC/le.cpp"
set SRCFILES="$SRCFILES $SRC/tdemsystemclass.cpp"
set INCLUDES="-I$SRC -I./.."
set CFLAGS="-pg -std=c++0x -O3 -DTESTOUTSIDEMATLAB"
set LINK_FFTW="-L/opt/cluster/fftw-3.1.2/lib -lfftw3 -lm"

set OUTDIR="./"
set CMD="icpc $CFLAGS $INCLUDES $SRCFILES $LINK_FFTW -o $OUTDIR/tdaem1d_shrlib_tester_intel-pg.exe"
echo $CMD
$CMD

set OUTDIR="./"
set CMD="g++ $CFLAGS $INCLUDES $SRCFILES $LINK_FFTW -o $OUTDIR/tdaem1d_shrlib_tester_gnu-pg.exe"
echo $CMD
$CMD
