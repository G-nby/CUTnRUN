#!/bin/bash
#SBATCH -p bio
#SBATCH -J cobed3loose
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH -N 1

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <filelist>"
    exit 1
fi

mkdir -p results

filelist=$1

mapfile -t files < "$filelist"

for ((i=0; i<${#files[@]}; i++)); do
    f="${files[$i]}"
    bn=$(basename "$f")

    echo "Processing $bn ..."

    tmpB="results/tmp_${bn}_others.bed"
    rm -f "$tmpB"

    for ((j=0; j<${#files[@]}; j++)); do
        if [ $i -ne $j ]; then
            cat "${files[$j]}" >> "$tmpB"
        fi
    done
    
    bedtools intersect \
        -a "$f" \
        -b "$tmpB" \
        -u > "results/${bn}_LooseOverlap.bed"

    bedtools intersect \
        -a "$f" \
        -b "$tmpB" \
        -v > "results/${bn}_Unique.bed"

    rm "$tmpB"
done



