#!/bin/bash
#SBATCH -p Acluster
#SBATCH -J quanttest
#SBATCH -o 1.%j.out
#SBATCH -e 1.%j.err
#SBATCH -N 1
#SBATCH --cpus-per-task=20


python TEwithSalmon.py quant --reference=hs --outpath=/Share/home/gby/datagby/TEsalmon2/hsout --num_threads=8 /Share/home/gby/datagby/TEsalmon2/hsdata

