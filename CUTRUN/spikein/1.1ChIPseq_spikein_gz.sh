#!/bin/bash

file_prefix=$1
is_trim=$2
SEorPE=$3
#add_bowtie2_index=$4
star_index=$4
add_genome_index=$5
clean_string=$6
add_bw_script=$7
blacklist=$8
genomesize=$9
add_bowtie2_index=${10}
cpu_num=${11}

echo "9 is $9"
echo "10 is $10"
echo "11 is $11"
echo "$cpu_num is cpu_num"
echo "/2 = $((cpu_num / 2))"

fruitfly_star_index="/Data/lht/lilab_users/25.GBY/Ensembl_genome/BDGP6/STAR_index_BDGP6"
fruitfly_clean_string="Unmapped|mapped|2110000|rDNA|mitochondrion"
fruitfly_blacklist="/Data/lht/lilab_users/25.GBY/blacklist/dm6-blacklist.v2.bed"

mkdir -p STAR_mapped_${file_prefix}
mkdir -p STAR_spikein_${file_prefix}


if [ $is_trim == "trim" ] && [ ${SEorPE} == "PE" ]; then
        # Trimming
        echo 'Mode: Trimming + PE'
        echo 'Start trimming'
        trim_galore --paired ${file_prefix}*1.*f*q.gz ${file_prefix}*2.*f*q.gz  --dont_gzip
        echo "Trimming done"

        # Mapping
        echo "Start mapping"
        #bowtie2 -p 16 -t -q  -N 1 -L 25 -X 700 --no-unal -x ${add_bowtie2_index} -1 ${file_prefix}_1*_val_1.fq -2 ${file_prefix}_2*_val_2.fq > ${file_prefix}.sam 2>>${file_prefix}_log.txt
        STAR --genomeDir ${star_index} --runThreadN ${cpu_num} --readFilesIn ${file_prefix}_1*_val_1.fq ${file_prefix}_2*_val_2.fq --outFileNamePrefix STAR_mapped_${file_prefix}/${file_prefix}
        STAR --genomeDir ${fruitfly_star_index} --runThreadN ${cpu_num} --readFilesIn ${file_prefix}_1*_val_1.fq ${file_prefix}_2*_val_2.fq --outFileNamePrefix STAR_spikein_${file_prefix}/spikein_${file_prefix}
        wait
        rm ${file_prefix}*_val_*.fq
elif [ $is_trim == "notrim" ] && [ ${SEorPE} == "PE" ]; then
        echo 'Mode: No trimming + PE'
        #bowtie2 -p 16 -t -q  -N 1 -L 25 -X 700 --no-unal -x ${add_bowtie2_index} -1 ${file_prefix}_1*.f*q.gz -2 ${file_prefix}_2*.f*q.gz > ${file_prefix}.sam 2>>${file_prefix}_log.txt
        STAR --genomeDir ${star_index} --runThreadN ${cpu_num} --readFilesIn ${file_prefix}_1*.f*q.gz ${file_prefix}_2*.f*q.gz --outFileNamePrefix STAR_mapped_${file_prefix}/${file_prefix}
        STAR --genomeDir ${fruitfly_star_index} --runThreadN ${cpu_num} --readFilesIn ${file_prefix}_1*.f*q.gz ${file_prefix}_2*.f*q.gz --outFileNamePrefix STAR_spikein_${file_prefix}/spikein_${file_prefix}
        wait
elif [ $is_trim == "trim" ] && [ ${SEorPE} == "SE" ]; then
        # Trimming
        echo 'Mode: Trimming + SE'
        echo 'Start trimming'
        trim_galore ${file_prefix}.f*q.gz --dont_gzip
        echo "Trimming done"

        # Mapping
        #bowtie2 -p 16 -t -q  -N 1 -L 25 -X 700 --no-unal -x ${add_bowtie2_index} -U ${file_prefix}_trimmed.fq > ${file_prefix}.sam 2>>${file_prefix}_log.txt
        STAR --genomeDir ${star_index} --runThreadN ${cpu_num} --readFilesIn ${file_prefix}_trimmed.fq --outFileNamePrefix STAR_mapped_${file_prefix}/${file_prefix}
        STAR --genomeDir ${fruitfly_star_index} --runThreadN ${cpu_num} --readFilesIn ${file_prefix}_trimmed.fq --outFileNamePrefix STAR_spikein_${file_prefix}/spikein_${file_prefix}
        wait
elif [ $is_trim == "notrim" ] && [ ${SEorPE} == "SE" ]; then
        echo 'Mode: No trimming + SE'
        #bowtie2 -p 16 -t -q  -N 1 -L 25 -X 700 --no-unal -x ${add_bowtie2_index} -U ${file_prefix}.f*q.gz > ${file_prefix}.sam 2>>${file_prefix}_log.txt
        STAR --genomeDir ${star_index} --runThreadN ${cpu_num} --readFilesIn ${file_prefix}.f*q.gz --outFileNamePrefix STAR_mapped_${file_prefix}/${file_prefix}
        STAR --genomeDir ${fruitfly_star_index} --runThreadN ${cpu_num} --readFilesIn ${file_prefix}.f*q.gz --outFileNamePrefix STAR_spikein_${file_prefix}/spikein_${file_prefix}
        wait
fi
echo "Mapping done"

# Remove duplicates
echo "Start removing duplicates"
#samtools view -bt ${add_genome_index} ${file_prefix}.sam | samtools sort  >${file_prefix}_sorted.bam
samtools view -F 256 -O bam STAR_mapped_${file_prefix}/${file_prefix}Aligned.out.sam | samtools sort -@ $((cpu_num / 2)) > main_${file_prefix}_prim.bam
samtools view -F 256 -O bam STAR_spikein_${file_prefix}/spikein_${file_prefix}Aligned.out.sam   | samtools sort -@ $((cpu_num / 2)) > spikein_${file_prefix}_prim.bam
wait

#picard MarkDuplicates INPUT=${file_prefix}_sorted.bam OUTPUT=${file_prefix}_sorted_rmdup.bam METRICS_FILE=marked_dup_metrics.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true CREATE_INDEX=true
picard MarkDuplicates I=main_${file_prefix}_prim.bam O=main_${file_prefix}_prim_rmdup.bam M=main_${file_prefix}_dups.txt REMOVE_DUPLICATES=true
picard MarkDuplicates I=spikein_${file_prefix}_prim.bam  O=spikein_${file_prefix}_prim_rmdup.bam  M=spikein_${file_prefix}_dups.txt  REMOVE_DUPLICATES=true
wait
echo "Remvoe duplicates done"

#calculate spikein scalefactor
main_count=$(samtools view -@ 2 main_${file_prefix}_prim_rmdup.bam | wc -l)
spikein_count=$(samtools view -@ 2 spikein_${file_prefix}_prim_rmdup.bam | wc -l)

echo "$main_count" > ${file_prefix}_readcounts.txt
echo "$spikein_count" >> ${file_prefix}_readcounts.txt

# Remove blacklisted ones
echo "Start removing blacklist"
bedtools intersect -v -abam main_${file_prefix}_prim_rmdup.bam -b ${blacklist} > main_${file_prefix}_prim_rmdup_BLrm.bam
wait

samtools view -@ $((cpu_num / 2)) -h -b -q 30 -F 1804 -f 2 main_${file_prefix}_prim_rmdup_BLrm.bam > main_${file_prefix}_prim_rmdup_BLrm_Filtered.bam
wait
echo "Remove blacklist done"

# Clean bam file, remove reads with 2 or more alignment results
echo "Start bam cleanning"
samtools view main_${file_prefix}_prim_rmdup_BLrm_Filtered.bam | awk "\$3 !~ /${clean_string}/" - | grep -v 'XS:i:' | samtools view -bt ${add_genome_index} - | samtools sort - > ${file_prefix}_clean.bam
wait
echo "Bam cleanning done"

samtools index ${file_prefix}_clean.bam

#bash ${add_bw_script} ${file_prefix}_clean.bam ${SEorPE}
#wait
#echo "100bp bw done"

# Add the rpgc one from shuai's pipeline
#bamCoverage --binSize 10 -p 8 --normalizeUsing RPGC --effectiveGenomeSize ${genomesize} --ignoreForNormalization chrX --bam ${file_prefix}_clean.bam -o ${file_prefix}_clean_Filtered_10bp_rpgc.bw
#wait
#echo "10bp rpgc bw done"



# Callpeak
# /Share/app/anaconda3-20211117/bin/macs2 callpeak -t ${file_prefix}_clean.bam  -f BAM -g hs --nomodel --nolambda -q 0.05 -n ${file_prefix}_clean --outdir MACS2
mkdir -p ChIPseq_PEAKCALL
mkdir -p ChIPseq_PEAKCALL_broad

if [ ${SEorPE} == "PE" ] && echo "$add_bowtie2_index" | grep -qE 'hg..$'; then
        macs2 callpeak -t ${file_prefix}_clean.bam -f BAMPE -n ${file_prefix}_clean -g hs --keep-dup 1 --cutoff-analysis --outdir ChIPseq_PEAKCALL &
        macs2 callpeak -t ${file_prefix}_clean.bam -f BAMPE -n ${file_prefix}_clean -g hs --keep-dup 1 --cutoff-analysis --broad  --outdir ChIPseq_PEAKCALL_broad &
        wait
elif [ ${SEorPE} == "SE" ] && echo "$add_bowtie2_index" | grep -qE 'hg..$'; then
        macs2 callpeak -t ${file_prefix}_clean.bam -f BAM -n ${file_prefix}_clean -g hs --keep-dup 1 --cutoff-analysis --outdir ChIPseq_PEAKCALL &
        macs2 callpeak -t ${file_prefix}_clean.bam -f BAM -n ${file_prefix}_clean -g hs --keep-dup 1 --cutoff-analysis --broad  --outdir ChIPseq_PEAKCALL_broad &
        wait
elif [ ${SEorPE} == "PE" ] && echo "$add_bowtie2_index" | grep -qE 'GRCm..$'; then
        macs2 callpeak -t ${file_prefix}_clean.bam -f BAMPE -n ${file_prefix}_clean -g mm --keep-dup 1 --cutoff-analysis --outdir ChIPseq_PEAKCALL &
        macs2 callpeak -t ${file_prefix}_clean.bam -f BAMPE -n ${file_prefix}_clean -g mm --keep-dup 1 --cutoff-analysis --broad  --outdir ChIPseq_PEAKCALL_broad &
        wait
elif [ ${SEorPE} == "SE" ] && echo "$add_bowtie2_index" | grep -qE 'GRCm..$'; then
        macs2 callpeak -t ${file_prefix}_clean.bam -f BAM -n ${file_prefix}_clean -g mm --keep-dup 1 --cutoff-analysis --outdir ChIPseq_PEAKCALL &
        macs2 callpeak -t ${file_prefix}_clean.bam -f BAM -n ${file_prefix}_clean -g mm --keep-dup 1 --cutoff-analysis --broad  --outdir ChIPseq_PEAKCALL_broad &
        wait
fi
echo "peak calling done"

# cut MACS2/${file_prefix}_clean_peaks.narrowPeak -f 1-3 > MACS2/${file_prefix}_clean_peaks.narrowPeak.bed
cut ChIPseq_PEAKCALL/${file_prefix}_clean_peaks.narrowPeak -f 1-3 > ChIPseq_PEAKCALL/${file_prefix}_clean_peaks.narrowPeak.bed
cut ChIPseq_PEAKCALL_broad/${file_prefix}_clean_peaks.broadPeak -f 1-3 > ChIPseq_PEAKCALL_broad/${file_prefix}_clean_peaks.broadPeak.bed
cut ChIPseq_PEAKCALL_broad/${file_prefix}_clean_peaks.gappedPeak -f 1-3 > ChIPseq_PEAKCALL_broad/${file_prefix}_clean_peaks.gappedPeak.bed
wait

# Clean files
rm ${file_prefix}_trimmed.fq
rm ${file_prefix}.sam
rm ${file_prefix}_sorted_rmdup.bam
rm ${file_prefix}_sorted_rmdup_BLrm.bam
rm ${file_prefix}_sorted_rmdup_BLrm_Filtered.bam
rm ${file_prefix}_sorted.bam
rm ${file_prefix}*_val_*.fq
rm ${file_prefix}_sorted_rmdup.bai

