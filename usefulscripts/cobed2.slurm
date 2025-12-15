#!/bin/bash
#SBATCH -p bio
#SBATCH -J cobed
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH -N 1

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <fileA> <fileB> <output_directory>"
    exit 1
fi

fileA=$1
fileB=$2
output_dir=$3

if [ ! -d "$output_dir" ]; then
    mkdir -p "$output_dir"
fi

awk 'BEGIN{OFS="\t"}{c=int(($2+$3)/2); print $1, c, c, $4}' "$fileA" > A_center.bed 
awk 'BEGIN{OFS="\t"}{c=int(($2+$3)/2); print $1, c, c, $4}' "$fileB" > B_center.bed
bedtools closest -d -a A_center.bed \
	-b B_center.bed > tmp.bed1
 awk '$NF <= 100 {print $4}' tmp.bed1 > "$output_dir/commonNames.txt"

grep -F -w -f "$output_dir/commonNames.txt" "$fileA" > "$output_dir/commonPeaks.bed"

rm A_center.bed
rm B_center.bed
rm tmp.bed1
rm "$output_dir/commonNames.txt"
#bedtools closest -d -a "$fileA" -b "$fileB" | awk '$NF <= 100' > "$output_dir/commonPeaks.bed"
#bedtools closest -d -a "$fileA" -b "$fileB" > "$output_dir/commonPeaks.bed"

#bedtools intersect \
#        -a "$fileA" \
#        -b "$fileB" \
#        -w 100 -r > "$output_dir/commonPeaks.bed"

bedtools intersect \
        -a "$fileA" \
        -b "$fileB" \
        -v > "$output_dir/uniqueA.bed"

bedtools intersect \
        -a "$fileB" \
        -b "$fileA" \
        -v > "$output_dir/uniqueB.bed"
        
bedtools intersect \
        -a "$fileA" \
        -b "$fileB" \
        -u > "$output_dir/LooseOverlapA.bed"

bedtools intersect \
        -a "$fileB" \
        -b "$fileA" \
        -u > "$output_dir/LooseOverlapB.bed"


