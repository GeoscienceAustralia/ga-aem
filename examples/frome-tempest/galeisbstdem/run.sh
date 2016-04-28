#!/bin/tcsh

#Standalone
#galeisbstdem.exe galeisbstdem-do-not-solve-geometry.con
#galeisbstdem.exe galeisbstdem-solve-rxpitch.con
#galeisbstdem.exe galeisbstdem-solve-rxpitch-and-offsets.con

#Use 4 MPI processes
#mpirun -np 4 galeisbstdem.exe galeisbstdem-do-not-solve-geometry.con
#mpirun -np 4 galeisbstdem.exe galeisbstdem-solve-rxpitch.con
mpirun -np 4 galeisbstdem.exe galeisbstdem-solve-rxpitch-and-offsets.con

#Use 4 OpenMP threads
#galeisbstdem.exe galeisbstdem-do-not-solve-geometry.con 4
#galeisbstdem.exe galeisbstdem-solve-rxpitch.con 4
#galeisbstdem.exe galeisbstdem-solve-rxpitch-and-offsets.con 4

