#!/bin/tcsh

# Reconstruct the primary field from in input TFR geometry and then invert the amplitude of total-field data in the X&Z-component plane and solve for the Dx and Dz Tx-Rx offsets

# Run standalone
galeisbstdem.exe galeisbstdem-reconstruct-primaryfield-xzamplitude-solve-offsets.con

# Run using 4 OpenMP threads
# galeisbstdem.exe galeisbstdem-reconstruct-primaryfield-xzamplitude-solve-offsets.con 4

# Run using 4 MPI processes
# mpirun -np 4 galeisbstdem-reconstruct-primaryfield-xzamplitude-solve-offsets.con

