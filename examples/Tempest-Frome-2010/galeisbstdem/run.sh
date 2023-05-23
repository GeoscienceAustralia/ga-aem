#!/bin/tcsh

#Run standalone
galeisbstdem.exe galeisbstdem-xzamplitude-solve-offsets.con
#galeisbstdem.exe galeisbstdem-xzcomponents-solve-offsets-and-rxpitch.con
#galeisbstdem.exe galeisbstdem-zcomponentonly-do-not-solve-geometry.con

#Run using 4 OpenMP threads
#galeisbstdem.exe galeisbstdem-xzamplitude-solve-offsets.con 4
#galeisbstdem.exe galeisbstdem-xzcomponents-solve-offsets-and-rxpitch.con 4
#galeisbstdem.exe galeisbstdem-zcomponentonly-do-not-solve-geometry.con 4

#Run using 4 MPI processes
#mpirun -np 4 galeisbstdem-xzamplitude-solve-offsets.con
#mpirun -np 4 galeisbstdem-xzcomponents-solve-offsets-and-rxpitch.con
#mpirun -np 4 galeisbstdem-zcomponentonly-do-not-solve-geometry.con

