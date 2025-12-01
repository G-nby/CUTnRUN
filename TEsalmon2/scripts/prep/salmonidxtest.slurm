#!/bin/bash
#SBATCH -p Acluster
#SBATCH -J sidxm
#SBATCH -o smi.%j.out
#SBATCH -e smi.%j.err
#SBATCH -N 1
#SBATCH --cpus-per-task=20


srun salmon index -t /Share/home/gby/0506antigen/data/mmref/mm10.fa -i /Share/home/gby/datagby/ZHAO0804/tmp/mmsalmonidx --threads 20

