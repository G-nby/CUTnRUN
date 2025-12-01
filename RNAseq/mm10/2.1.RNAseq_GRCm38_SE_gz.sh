#!/bin/bash

file_prefix=$1

is_trim=$2

if [ "$is_trim" = "trim" ]; then
        # trim
        echo 'Start trimming'
        trim_galore ${file_prefix}.f*q.gz --dont_gzip
        echo "Trimming done"

        # mapping
        echo "Start mapping"
        hisat2 -p 8 --dta -x /Data/lht/lilab_users/25.GBY/Ensembl_genome/GRCm38/hisat2_index/GRCm38 -U ${file_prefix}_val*.fq -S ${file_prefix}.sam 2>>${file_prefix}_log.txt &
        wait
        echo "Mapping done"
	rm ${file_prefix}_val*.fq

elif [ "$is_trim" = "notrim" ]; then
        # mapping
        echo "Start mapping"
        hisat2 -p 8 --dta -x /Data/lht/lilab_users/25.GBY/Ensembl_genome/GRCm38/hisat2_index/GRCm38 -U ${file_prefix}.f*q.gz -S ${file_prefix}.sam 2>>${file_prefix}_log.txt &
        wait
        echo "Mapping done"
fi

# To bam
echo "Start bam conversion"
samtools view -uhS ${file_prefix}.sam |samtools sort -o ${file_prefix}_sorted.bam &
wait

samtools index ${file_prefix}_sorted.bam &
wait
echo "Bam conversion done"

# htseq-count
# Gene expression counting
echo "Gene expression counting start"
htseq-count -s no -r pos -t exon -n 16 -i gene_id -f bam ${file_prefix}_sorted.bam /Data/lht/lilab_users/25.GBY/Ensembl_genome/GRCm38/Mus_musculus.GRCm38.102.gtf   > ${file_prefix}_Mm_GRCm38.102.chr.counts &
wait
echo "Gene expression counting done"

# Clean bam
echo "Start bam cleanning"
samtools view ${file_prefix}_sorted.bam |awk '$3 !~ /GL|JH|MT/' - |grep -v 'XS:i:' | samtools view -bt /Data/lht/lilab_users/25.GBY/Ensembl_genome/GRCm38/Mus_musculus.GRCm38.dna.primary_assembly.fa.fai - |samtools sort - > ${file_prefix}_clean.bam
wait
echo "Bam cleanning done"

bash /Share/home/gby/CUTnRUN/RNAseq/mm10/3.bam_to_tsv_GRCm38_RNA.sh ${file_prefix}_clean.bam
wait

echo "Clean the sam file"
rm ${file_prefix}.sam
echo "File cleanning done"
