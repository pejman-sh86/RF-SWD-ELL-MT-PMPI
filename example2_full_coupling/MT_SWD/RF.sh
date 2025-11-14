#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12  
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00

module load gcc/7.3.0   
module load nvhpc/nvhpc_21.7  

mpirun ../bin/prjmh_temper_rf
