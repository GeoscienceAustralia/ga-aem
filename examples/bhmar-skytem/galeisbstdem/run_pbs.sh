#!/bin/bash

#This is an example PBS shell script for running on raijin.nci.org.au

#PBS -P r17
#PBS -q normal
#PBS -l ncpus=16
#PBS -l mem=16GB
#PBS -l walltime=01:00:00
#PBS -l wd
#PBS -N galeisbstdem
#PBS -o galeisbstdem.out
#PBS -e galeisbstdem.err
#PBS -j oe

cd $PBS_O_WORKDIR
module load gaaem/1.0.beta
mpirun galeisbstdem.exe galeisbstdem.con

