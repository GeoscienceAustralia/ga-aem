#!/bin/tcsh

set nloops=100
#set stm1="../examples/stmfiles/SkytemLM-SLV_newwave.stm"
#set stm2="../examples/stmfiles/SkytemHM-SLV_newwave.stm"
set stm1="../examples/stmfiles/SkytemLM-SLV_newwave_zonly.stm"
set stm2="../examples/stmfiles/SkytemHM-SLV_newwave_zonly.stm"

echo $nloops
echo $stm1
echo $stm2

rm gmon.out
tdaem1d_shrlib_tester_intel-pg.exe $nloops $stm1 $stm2
gprof -c tdaem1d_shrlib_tester_intel-pg.exe > profiling_results.txt


