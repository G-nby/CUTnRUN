#!/bin/bash
#SBATCH -p bio
#SBATCH -J bam_cor
#SBATCH -o 1.%j.out
#SBATCH -e 1.%j.err
#SBATCH -N 1
#SBATCH --cpus-per-task=8

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <bam_dir> <outprefix> [binsize]"
    exit 1
fi

bam_dir=$1
prefix=$2
binsize=${3:-1000}

bam_files=$(ls ${bam_dir}/*.bam 2>/dev/null)

if [ -z "$bam_files" ]; then
    echo "no bam file found in $bam_dir"
    exit 1
fi

for bam in $bam_files
do
    bai="${bam}.bai"
    if [ ! -f "$bai" ]; then
        samtools index "$bam"
    fi
done
echo "bai check done"

echo "mk rawcounts"
multiBamSummary bins \
    --bamfiles $bam_files \
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

