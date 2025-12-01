echo "The bam must be sorted!!!!"
bam=$1
SEorPE=$2
samtools index $bam

if [ ${SEorPE} == "PE" ];then
        bamCoverage -b $bam -bs 100 -p 8 --normalizeUsing RPKM -o ${bam%.bam}_100bp_rpkm.bw
elif [ ${SEorPE} == "SE" ];then
        bamCoverage -b $bam -bs 100 -p 8 --normalizeUsing RPKM --extendReads 250 -o ${bam%.bam}_100bp_rpkm.bw
else
        echo wrong SEorPE
        echo ${SEorPE} iis not SE or PE
        exit
fi

#echo ${SEorPE} and we start bigwigtowig
bigWigToWig ${bam%.bam}_100bp_rpkm.bw ${bam%.bam}_100bp_rpkm.wig
#head ${bam%.bam}_100bp_rpkm.wig3
grep -v "#" ${bam%.bam}_100bp_rpkm.wig > ${bam%.bam}_100bp_rpkm.wig.tmp
mv ${bam%.bam}_100bp_rpkm.wig.tmp ${bam%.bam}_100bp_rpkm.wig
awk '{print $0"\t"$4}' ${bam%.bam}_100bp_rpkm.wig | sort-bed - | bedmap --mean --echo /Share/home/lht/bioinformatic_tools/Ensembl_genome/GRCh38/GRCh38_window_100_step_100.bed - | sed "s/\(^[^:]*\)|/\1\t/g" | awk '{print $0"\t"$1}' | cut -f2- |sed "s/\tNAN$/\t0/" | awk '{printf "%s\t%d\t%.2f\n",$1,($2+$3)/2,$4}' > ${bam%.bam}_100bp_rpkm.tsv
