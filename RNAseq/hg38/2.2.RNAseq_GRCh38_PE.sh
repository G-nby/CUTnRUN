#!/bin/bash

file_prefix=$1
is_trim=$2

if [ "$is_trim" = "trim" ]; then

	# trim
	echo 'Start trimming'
	trim_galore --paired ${file_prefix}_1*f*q ${file_prefix}_2*f*q
	echo "Trimming done"

	# mapping
	echo "Start mapping"
	hisat2 -p 8 --dta -x /Share/home/lht/bioinformatic_tools/Ensembl_genome/GRCh38/hisat2_index/GRCh38 -1 ${file_prefix}_1_val*.fq -2 ${file_prefix}_2_val*.fq -S ${file_prefix}.sam 2>>${file_prefix}_log.txt &
	wait
	echo "Mapping done"
	#rm ${file_prefix}_1_val*.fq  ${file_prefix}_2_val*.fq

elif [ "$is_trim" = "notrim" ]; then
	# mapping
	echo "Start mapping"
  echo "	hisat2 -p 8 --dta -x /userdata/Ensembl_genome/GRCh38/hisat2_index/GRCh38 -1 ${file_prefix}_1*f*q -2 ${file_prefix}_2*f*q -S ${file_prefix}.sam 2>>${file_prefix}_log.txt &"
	hisat2 -p 8 --dta -x /Share/home/lht/bioinformatic_tools/Ensembl_genome/GRCh38/hisat2_index/GRCh38 -1 ${file_prefix}_1*f*q -2 ${file_prefix}_2*f*q -S ${file_prefix}.sam 2>>${file_prefix}_log.txt &
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
htseq-count -s no -r pos -t exon -n 16 -i gene_id -f bam ${file_prefix}_sorted.bam /Share/home/lht/bioinformatic_tools/Ensembl_genome/GRCh38/Homo_sapiens.GRCh38.109.chr.gtf  > ${file_prefix}_Hs_GRCh38.109.chr.counts &
wait
echo "Gene expression counting done"

# Clean bam
echo "Start bam cleanning"
samtools view ${file_prefix}_sorted.bam |awk '$3 !~ /GL|KI|MT/' - |grep -v 'XS:i:' | samtools view -bt /Share/home/lht/bioinformatic_tools/Ensembl_genome/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai - |samtools sort - > ${file_prefix}_clean.bam
wait
echo "Bam cleanning done"

bash /Share/home/gby/CUTnRUN/RNAseq/hg38/3.bam_to_tsv_GRCh38_RNA.sh ${file_prefix}_clean.bam
wait


echo "Clean the sam"
rm ${file_prefix}.sam 
echo "File cleanning done"
