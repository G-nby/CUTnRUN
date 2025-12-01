#!/bin/bash
#SBATCH -p Acluster
#SBATCH -J ttest
#SBATCH -o t.%j.out
#SBATCH -e t.%j.err
#SBATCH -N 1
#SBATCH --cpus-per-task=20


python /Share/home/gby/CUTnRUN/TEsalmon2/scripts/TEwithSalmon.py test --inpath=/Data/gby/TEsalmon2/hsout2 \
	--outpath=/Data/gby/TEsalmon2/hsout2/testout4 \
	--table=tsv \
	--figtype=png \
	--analysis_type=DE \
	--conditions=WT,KO

