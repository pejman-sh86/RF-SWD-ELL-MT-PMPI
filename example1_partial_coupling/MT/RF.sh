#!/bin/bash

#SBATCH --nodes=2
#SBATCH --ntasks-per-node=6   
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00



mpirun ../bin/prjmh_temper_rf
