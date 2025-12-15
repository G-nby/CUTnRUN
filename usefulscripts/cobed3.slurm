#!/bin/bash
#SBATCH -p bio
#SBATCH -J cobed3
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH -N 1

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <filelist>"
    exit 1
fi

mkdir -p tmp

filelist=$1
fileA=$(head -n 1 "$filelist")

while read -r f; do
    bn=$(basename "$f" )
    awk 'BEGIN{OFS="\t"}{c=int(($2+$3)/2); print $1, c-50, c+50, $4}' "$f" \
        > tmp/${bn}_center100.bed
done < "$filelist"

bedtools multiinter  -i tmp/*_center100.bed > tmp/multi_common.bed

num_files=$(wc -l < "$filelist")
echo "total file num : $num_files"

awk -v n="$num_files" '$4==n' tmp/multi_common.bed > tmp/common_intervals.bed

bedtools intersect -F 1 -a "$fileA" -b tmp/common_intervals.bed > common_all_peaks.bed

rm -r tmp
