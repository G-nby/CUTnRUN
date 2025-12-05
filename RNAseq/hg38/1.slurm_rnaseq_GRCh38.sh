#!/bin/bash 
#SBATCH -J RNA_GRCh38 
#SBATCH -p bio
#SBATCH -n 8 
#SBATCH --output=%j.out 
#SBATCH --error=%j.err 

DATA_DIR=$(pwd)

count=$(find "$DATA_DIR" -maxdepth 1 -type f \( -iname "*.fq" -o -iname "*.fastq" \) | wc -l)
#count=`ls *.f*q | wc -l`
if [ $# > 0 ]; then
        if[[ $1 == "trim" ]]; then
                echo "Data will be trimmed to remove adaptors."
        elif[[ $1 == "notrim" ]]; then
                echo "Data will not be trimmed."
        else
                echo "Wrong parameter for trimming."
                exit
        fi
else
        echo "Wrong parameter number."
        exit
fi

if [ $count == 1 ]; then
        i=`ls *.f*q`
        file_prefix=${i%.f*q}
        echo "File prefix=${file_prefix}"
        echo "Only one sequencing data file detected, single-end mode applied."
        sh /Share/home/gby/CUTnRUN/RNAseq/hg38/2.1.RNAseq_GRCh38_SE.sh ${file_prefix} $1
elif [ $count == 2 ]; then
        i=`ls *_1*f*q`
        file_prefix=${i%_1*f*q}
        echo "File prefix=${file_prefix}"
        echo "Two sequencing data files detected, paired-end mode applied."
        sh /Share/home/gby/CUTnRUN/RNAseq/hg38/2.2.RNAseq_GRCh38_PE.sh ${file_prefix} $1
else
        echo "Wrong file count, failed to submit the task."
fi

