#!/bin/bash

#This is an example PBS shell script for running on gadi.nci.org.au

#PBS -P cr78
#PBS -q express
#PBS -l ncpus=48
#PBS -l mem=48GB
#PBS -l walltime=01:00:00
#PBS -l storage=gdata/qi71+gdata/cr78
#PBS -l wd
#PBS -N galeisbstdem
#PBS -o galeisbstdem.out
#PBS -e galeisbstdem.err
#PBS -j oe

module load ga-aem/v2.0.0-Release20240424
module list

cat run_pbs.sh
cat galeisbstdem.con

echo Working directory
echo $PWD

which galeisbstdem.exe
mpirun galeisbstdem.exe galeisbstdem.con

