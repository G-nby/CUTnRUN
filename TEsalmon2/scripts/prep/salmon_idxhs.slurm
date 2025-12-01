#!/bin/bash
#SBATCH -p Acluster
#SBATCH -J sidxm
#SBATCH -o smi.%j.out
#SBATCH -e smi.%j.err
#SBATCH -N 1
#SBATCH --cpus-per-task=20


srun salmon index -t "/Share/home/gby/CUTnRUN/TEsalmon2/REdiscoverTE/newref/genome.fasta" -i reference/hs --threads 20

