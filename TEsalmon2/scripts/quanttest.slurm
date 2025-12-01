#!/bin/bash
#SBATCH -p Acluster
#SBATCH -J quanttest
#SBATCH -o 1.%j.out
#SBATCH -e 1.%j.err
#SBATCH -N 1
#SBATCH --cpus-per-task=8


python /Share/home/gby/CUTnRUN/TEsalmon2/scripts/TEwithSalmon.py quant \
	--reference=mouse1023 \
	--outpath=/Work/lht/lilab_users/25.GBY/BCandmm/1026mm/temouse1023out --num_threads=8 \
	"/Work/lht/lilab_users/25.GBY/BCandmm/1026mm/clean" --exprtype count

