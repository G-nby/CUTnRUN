#!/bin/bash
#SBATCH -p bio
#SBATCH -J rna_pca
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH -N 1

Rscript /Share/home/gby/CUTnRUN/usefulscripts/pcaVST.R EXPR.csv

