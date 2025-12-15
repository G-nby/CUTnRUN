#!/bin/bash

#SBATCH -p bio
#SBATCH -J homer
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH -N 1


if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_file> <ref>"
    exit 1
fi

input_file=$1
inref=$2
output="${input_file%.*}.${inref}.anno.xls"


tmp_peak=tmp_chr.peak
#awk 'BEGIN{OFS=FS="\t"} !/^chr/ {$1="chr"$1} {print $4, $1, $2, $3, $6}' "$input_file" > "$tmp_peak"
awk 'BEGIN{OFS=FS="\t"} !/^chr/ {$1="chr"$1} {print}' "$input_file" > "$tmp_peak"

if [[ $inref == TE* ]]; then
    ref=${inref#TE}
    echo "$ref"
    annotatePeaks.pl "$tmp_peak" "$ref" -ann /Share/home/gby/datagby/homer_ref/${inref}/${inref}_anno.final.txt > "$output"
else
    ref=$inref
    annotatePeaks.pl "$tmp_peak" "$ref" > "$output"
fi

rm -f "$tmp_peak"
