#!/bin/bash
#SBATCH -p Acluster
#SBATCH -J mkexpr
#SBATCH -o 1.%j.out
#SBATCH -e 1.%j.err
#SBATCH -N 1
#SBATCH --cpus-per-task=20


python /Share/home/gby/CUTnRUN/TEsalmon2/scripts/TEwithSalmon.py mk_expr --inpath=/Data/gby/zhao/250911RNAseq --outpath=/Data/gby/zhao/250911RNAseq --num_threads=8 --ref_name=hs

