#!/bin/tcsh
clear

set SRC="./src"
set SRCFILES=""
set SRCFILES="$SRCFILES $SRC/tdaem1d_shrlib.cpp"
set SRCFILES="$SRCFILES $SRC/rb_utils.cpp"
set SRCFILES="$SRCFILES $SRC/directoriesandfiles.cpp"
set SRCFILES="$SRCFILES $SRC/blocklanguage.cpp"
set SRCFILES="$SRCFILES $SRC/geometry3d.cpp"
set SRCFILES="$SRCFILES $SRC/le.cpp"
set SRCFILES="$SRCFILES $SRC/tdemsystemclass.cpp"
set INCLUDES="-I$SRC -I./.."
set CFLAGS="-O3 -fPIC -shared"

set OUTDIR="../linux_intel"
set CMD="icpc $CFLAGS $INCLUDES $SRCFILES -lfftw3 -o $OUTDIR/tdaem1d_shrlib.so"
echo $CMD
$CMD

set OUTDIR="../linux_gnu"
set CMD="g++ $CFLAGS $INCLUDES $SRCFILES -lfftw3 -o $OUTDIR/tdaem1d_shrlib.so"
echo $CMD
$CMD
