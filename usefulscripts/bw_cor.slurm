#!/bin/bash
#SBATCH -p bio
#SBATCH -J bw_cor
#SBATCH -o 1.%j.out
#SBATCH -e 1.%j.err
#SBATCH -N 1
#SBATCH --cpus-per-task=8

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <bw_dir> <outprefix> <binsize>"
    exit 1
fi

bw_dir=$1
prefix=$2
#binsize=$3
binsize=${3:-1000}

bw_files=$(ls ${bw_dir}/*.bw 2>/dev/null)

if [ -z "$bw_files" ]; then
    echo "no bw file found in $bw_dir"
    exit 1
fi

echo "mk rawcounts"
multiBigwigSummary bins \
    --bwfiles $bw_files \
    --binSize $binsize \
    --numberOfProcessors 8 \
    --outFileName ${prefix}_summary.npz \
    --outRawCounts ${prefix}_rawCounts.tab

echo "drawing... yes i will skip the trap haha"

plotCorrelation \
    --corData ${prefix}_summary.npz \
    --corMethod pearson \
    --skipZeros \
    --removeOutliers \
    --whatToPlot scatterplot \
    --colorMap RdYlBu \
    -o ${prefix}_pearson_correlation.png \
    --outFileCorMatrix ${prefix}_pearson_values.tab

plotCorrelation \
    --corData ${prefix}_summary.npz \
    --corMethod pearson \
    --skipZeros \
    --removeOutliers \
    --whatToPlot scatterplot \
    --colorMap RdYlBu \
    --log1p \
    -o ${prefix}_pearson_correlation_log1p.png \
    --outFileCorMatrix ${prefix}_pearson_values_log1p.tab

