#!/bin/bash 
#SBATCH -J ChIP_common 
#SBATCH -p Acluster 
#SBATCH -n 16
#SBATCH --exclusive
#SBATCH --output=%j.out 
#SBATCH --error=%j.err 

count=`ls *.f*q | wc -l`

SEorPE=''
add_bowtie2_index=''
add_genome_index=''
clean_string=''
add_bw_script=''
blacklist=''
genomesize=0

# File count
if [ $count == 1 ]; then
	i=`ls *.f*q`
	file_prefix=${i%.f*q}
	echo "File prefix=${file_prefix}"
	echo "Only one sequencing data file detected, single-end mode applied."
        SEorPE='SE'
elif [ $count == 2 ]; then
	i=`ls *_1.f*q`
        file_prefix=${i%_1.f*q}
	echo "File prefix=${file_prefix}"
	echo "Two sequencing data files detected, paired-end mode applied."
        SEorPE='PE'
else
	echo "Abnormal gz file number."
        exit
fi

# Trim mode
if [ $1 != 'trim' ] && [ $1 != 'notrim' ]; then
        echo "The first parameter should either be trim or notrim"
        exit
fi

# Specieis
if [ $2 == 'ensembl_hg38' ]; then
        add_bowtie2_index='/Share/home/lht/bioinformatic_tools/Ensembl_genome/GRCh38/bowtie2_index/hg38'
        add_genome_index='/Share/home/lht/bioinformatic_tools/Ensembl_genome/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai'
        add_bw_script='/Share/home/gby/CUTnRUN/newscripts/bam_to_tsv_GRCh38.sh'
        clean_string='GL|KI|MT'
        blacklist='/Data/lht/lilab_users/25.GBY/blacklist/hg38-blacklist.v2.bed'
        genomesize=2913022398
elif [ $2 == 'ucsc_hg19' ]; then
        add_bowtie2_index='/Share/home/lht/bioinformatic_tools/UCSC_genome/hg19/bowtie2_index/hg19'
        add_genome_index='/Share/home/lht/bioinformatic_tools/UCSC_genome/hg19/hg19.fa.fai'
        add_bw_script='/Share/home/lht/bioinformatic_tools/bam_to_tsv_hg19.sh'
        clean_string='random|chrM|fix|chrMT|chrUn|_alt|hap'
        blacklist=''
        genomesize=3137161264
elif [ $2 == 'ensembl_mm39' ]; then
        add_bowtie2_index='/Share/home/lht/bioinformatic_tools/Ensembl_genome/GRCm39/bowtie2_index/GRCm39'
        add_genome_index='/Share/home/lht/bioinformatic_tools/Ensembl_genome/GRCm39/Mus_musculus.GRCm39.dna.primary_assembly.fa.fai '
        add_bw_script='/Share/home/lht/bioinformatic_tools/bam_to_tsv_GRCm39.sh'
        clean_string='GL|KI|JH|MU|MT'
        blacklist=''
        genomesize=2654621783
elif [ $2 == 'ensembl_mm10' ]; then
        add_bowtie2_index='/Data/lht/lilab_users/25.GBY/Ensembl_genome/GRCm38/bowtie2_index/GRCm38'
        add_genome_index='/Data/lht/lilab_users/25.GBY/Ensembl_genome/GRCm38/Mus_musculus.GRCm38.dna.primary_assembly.fa.fai'
        add_bw_script='/Share/home/gby/CUTnRUN/newscripts/bam_to_tsv_GRCm38.sh' 
        clean_string='GL|JH|MT'
        blacklist='/Data/lht/lilab_users/25.GBY/blacklist/mm10-blacklist.v2.bed'
        genomesize=2652783500
elif [ $2 == 'ensembl_rn7' ]; then
        add_bowtie2_index='/Share/home/lht/bioinformatic_tools/Ensembl_genome/mRatBN7.2/bowtie2_index/mRatBN7.2'
        add_genome_index='/Share/home/lht/bioinformatic_tools/Ensembl_genome/mRatBN7.2/Rattus_norvegicus.mRatBN7.2.dna_sm.toplevel.fa.fai'
        add_bw_script='/Share/home/lht/bioinformatic_tools/bam_to_tsv_mRatBN7.2.sh'
        clean_string='MU|JACYVU|MT'
        blacklist=''
elif [ $2 == 'ucsc_rn5' ]; then
        add_bowtie2_index='/Share/home/lht/bioinformatic_tools/UCSC_genome/rn5/bowtie2_index/rn5'
        add_genome_index='/Share/home/lht/bioinformatic_tools/UCSC_genome/rn5/rn5_assembly.fasta.fai'
        add_bw_script='/Share/home/lht/bioinformatic_tools/bam_to_tsv_rn5.sh'
        clean_string='MT'
        blacklist=''
elif [ $2 == 'ensembl_bdgp6' ]; then
        add_bowtie2_index='/Share/home/lht/bioinformatic_tools/Ensembl_genome/BDGP6.46/bowtie2_index/BDGP6.46'
        add_genome_index='/Share/home/lht/bioinformatic_tools/Ensembl_genome/BDGP6.46/Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa.fai'
        add_bw_script='/Share/home/lht/bioinformatic_tools/bam_to_tsv_BDGP6.46.sh'
        clean_string='Unmapped|mapped|2110000|rDNA|mitochondrion'
        blacklist=''
else
        echo wrong species. Only accepts: ensembl_hg38, ucsc_hg19, ensembl_mm39, ensembl_mm10, ensembl_rn7, ucsc_rn5, ensembl_bdgp6
        exit
fi

/Share/home/gby/CUTnRUN/newscripts/ecoli/2.ChIPseq_common.sh ${file_prefix} $1 ${SEorPE} ${add_bowtie2_index} ${add_genome_index} ${clean_string} ${add_bw_script} ${blacklist} ${genomesize}
