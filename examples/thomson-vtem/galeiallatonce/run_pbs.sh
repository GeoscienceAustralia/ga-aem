#!/bin/bash

#This is an example PBS shell script for running on raijin.nci.org.au
#PBS -P cr78
#PBS -q express
#PBS -l ncpus=16
#PBS -l mem=16GB
#PBS -l walltime=01:00:00
#PBS -l wd
#PBS -N galeiallatonce
#PBS -o galeiallatonce.out
#PBS -e galeiallatonce.err
#PBS -j oe

#module load ga-aem/dev-intel
#module list
#mpirun ../../../bin/raijin/intel/galeiallatonce.exe galeiallatonce.con

module load ga-aem/dev-gnu
module list
mpirun ../../../bin/raijin/gnu/galeiallatonce.exe galeiallatonce.con
