#!/bin/bash

#This is an example PBS shell script for running on gadi.nci.org.au

#PBS -P cr78
#PBS -q normal
#PBS -l ncpus=48
#PBS -l mem=48GB
#PBS -l walltime=01:00:00
#PBS -l wd
#PBS -N galeisbstdem
#PBS -o galeisbstdem.out
#PBS -e galeisbstdem.err
#PBS -j oe

module load openmpi/4.0.1
module load fftw3/3.3.8
module load eigen/3.3.7
module list

export PATH=/home/547/rcb547/code/repos/ga-aem-develop/bin/gadi/intel:$PATH
echo   $PATH
mpirun galeisbstdem.exe galeisbstdem.con


