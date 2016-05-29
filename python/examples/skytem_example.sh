#!/bin/bash

module load gcc/5.2.0
#module load openmpi/1.6.3
#module load fftw3/3.3.3
module load python3/3.3.0
module load python3/3.3.0-matplotlib

export PYTHONPATH=/short/public/rcb547/apps/gaaem-1.0/python

python3 skytem_example.py

