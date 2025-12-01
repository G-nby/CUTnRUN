#!/bin/bash
#SBATCH -p Acluster
#SBATCH -J ttest
#SBATCH -o t.%j.out
#SBATCH -e t.%j.err
#SBATCH -N 1
#SBATCH --cpus-per-task=8


python /Share/home/gby/CUTnRUN/TEsalmon2/scripts/TEwithSalmon.py test_single --inpath=/Data/gby/test_hisatDE/ \
	--outpath=/Data/gby/test_hisatDE/corout1 \
	--table=tsv \
	--figtype=png \
	--samples=ctrl1,treat1 \
	--ref_name=hs

