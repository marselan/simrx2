#!/bin/bash
#SBATCH --job-name=simrx
#SBATCH --partition=multi
#SBATCH --ntasks=20
#SBATCH --time 2-00:10
. /etc/profile
module load gcc
module load openmpi
srun simrx > output

