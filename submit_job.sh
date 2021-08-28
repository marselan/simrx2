#!/bin/bash
#SBATCH --job-name=simrx
#SBATCH --partition=capability
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=18
#SBATCH --time 2-00:10
. /etc/profile
module load mpi/openmpi/1.8.3-gcc_4.8.3
srun simrx
