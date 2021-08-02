#!/bin/bash
#SBATCH --job-name=simrx
#SBATCH --partition=phi
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time 2-00:10
. /etc/profile
# 57/114/171/228, mas para hyperthreading
# 228 se queda sin memoria
mpiexec -np 2 ./simrx-phi
