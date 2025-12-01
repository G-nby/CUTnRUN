echo "The bam must be sorted!!!!"
bam=$1
SEorPE=$2

samtools index $bam

if [ ${SEorPE} == "PE" ]; then
        bamCoverage -b $bam -bs 100 -p 8 --normalizeUsing RPKM -o ${bam%.bam}_100bp_rpkm.bw
elif [ ${SEorPE} == "SE" ]; then
        bamCoverage -b $bam -bs 100 -p 8 --normalizeUsing RPKM --extendReads 250 -o ${bam%.bam}_100bp_rpkm.bw
else
        echo wrong SEorPE
        exit
fi

bigWigToWig ${bam%.bam}_100bp_rpkm.bw ${bam%.bam}_100bp_rpkm.wig
grep -v "#" ${bam%.bam}_100bp_rpkm.wig > ${bam%.bam}_100bp_rpkm.wig.tmp
mv ${bam%.bam}_100bp_rpkm.wig.tmp ${bam%.bam}_100bp_rpkm.wig
awk '{print $0"\t"$4}' ${bam%.bam}_100bp_rpkm.wig | sort-bed - | bedmap --mean --echo /Data/lht/lilab_users/25.GBY/Ensembl_genome/GRCm38/GRCm38_window_100_step_100.bed - | sed "s/\(^[^:]*\)|/\1\t/g" | awk '{print $0"\t"$1}' | cut -f2- |sed "s/\tNAN$/\t0/" | awk '{printf "%s\t%d\t%.2f\n",$1,($2+$3)/2,$4}' > ${bam%.bam}_100bp_rpkm.tsv