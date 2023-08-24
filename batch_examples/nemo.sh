#!/bin/bash
#SBATCH --mem=100GB
#SBATCH --time=1-0
#SBATCH -N 1
#SBATCH -n 40
#SBATCH --partition=xeon40
module load QChem/5.2-multicore
export $QCLOCALSCR=/scratch/

#This is the key line:
bash $1


rm -rf /scratch/*
rm slurm*out 