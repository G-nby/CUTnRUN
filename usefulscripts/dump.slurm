#!/bin/bash 
#SBATCH -J dump
#SBATCH -p bio
#SBATCH -n 1 
#SBATCH --output=%j.out
#SBATCH --error=%j.err

srun bash "/Share/home/gby/CUTnRUN/usefulscripts/dump.sh" /Data/gby/antigen/mm/data /Data/gby/antigen/mm/output 
