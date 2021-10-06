#!/bin/bash
#SBATCH --mem=50GB
#SBATCH --time=1-0
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --partition=xeon24

module purge
module use /home/energy/modules/modules/all
module --ignore-cache load "binutils/2.31.1-GCCcore-8.2.0"
module load iomkl
module load QChem/5.2-multicore
export $QCLOCALSCR=/scratch/folder

bash $1
