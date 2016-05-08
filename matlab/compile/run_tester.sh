#!/bin/tcsh
clear

set nloops=100
set stm1="../examples/stmfiles/SkytemLM-SLV_newwave.stm"
set stm2="../examples/stmfiles/SkytemHM-SLV_newwave.stm"

echo $nloops
echo $stm1
echo $stm2
echo ""
echo "Intel xyz comps"
tdaem1d_shrlib_tester_intel.exe $nloops $stm1 $stm2
echo "GNU xyz comps"
tdaem1d_shrlib_tester_gnu.exe $nloops $stm1 $stm2


set stm1="../examples/stmfiles/SkytemLM-SLV_newwave_zonly.stm"
set stm2="../examples/stmfiles/SkytemHM-SLV_newwave_zonly.stm"
echo $nloops
echo $stm1
echo $stm2
echo ""
echo "Intel z comp only"
tdaem1d_shrlib_tester_intel.exe $nloops $stm1 $stm2
echo "GNU z comp only"
tdaem1d_shrlib_tester_gnu.exe $nloops $stm1 $stm2

