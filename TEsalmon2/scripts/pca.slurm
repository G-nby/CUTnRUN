#!/bin/bash
#SBATCH -p AFat
#SBATCH -J rna_pca
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH -N 1
#SBATCH --cpus-per-task=20

Rscript /Share/home/gby/CUTnRUN/TEsalmon2/scripts/pcatmp.R /Data/gby/zhao/250821RNAseqnotrim/out/EXPR.csv

