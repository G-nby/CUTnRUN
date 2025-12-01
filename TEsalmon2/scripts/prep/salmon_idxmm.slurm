#!/bin/bash
#SBATCH -p Acluster
#SBATCH -J sidxmm
#SBATCH -o mm.%j.out
#SBATCH -e mm.%j.err
#SBATCH -N 1
#SBATCH --cpus-per-task=20


python "/Share/home/gby/CUTnRUN/TEsalmon2/scripts/TEwithSalmon.py" index --ref_name=mm --input_fasta="/Share/home/gby/CUTnRUN/TEsalmon2/SalmonTE/SalmonTE-main/reference/mm/mm_upper.fa" --te_only

