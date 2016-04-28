#!/bin/tcsh

# Standalone
# galeisbstdem.exe galeisbstdem.con

# Use 4 MPI processes
mpiexec -np 4 galeisbstdem.exe galeisbstdem.con

# Use 4 OpenMP threads
# galeisbstdem.exe galeisbstdem.con 4



