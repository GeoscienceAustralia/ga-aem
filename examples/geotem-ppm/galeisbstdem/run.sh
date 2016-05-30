#!/bin/tcsh

#Standalone
galeisbstdem.exe galeisbstdem_z_nosolvegeometry.con
#galeisbstdem.exe galeisbstdem_xz_nosolvegeometry.con
#galeisbstdem.exe galeisbstdem_xz_solvegeometry.con

#Use 4 MPI processes
#mpiexec -np 4 galeisbstdem.exe galeisbstdem_z_nosolvegeometry.con
#mpiexec -np 4 galeisbstdem.exe galeisbstdem_xz_nosolvegeometry.con
#mpiexec -np 4 galeisbstdem.exe galeisbstdem_xz_solvegeometry.con

#Use 4 OpenMP threads
#galeisbstdem.exe galeisbstdem_z_nosolvegeometry.con 4
#galeisbstdem.exe galeisbstdem_xz_nosolvegeometry.con 4
#galeisbstdem.exe galeisbstdem_xz_solvegeometry.con 4

