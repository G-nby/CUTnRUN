#!/bin/bash
#SBATCH --job-name=cutrun_ana
#SBATCH --output=%j.out 
#SBATCH --error=%j.err 
#SBATCH --cpus-per-task=8      
#SBATCH --partition=bio

# ======== the params you shoudf modify when use ======
reference=hg38                    
#bam_dir="/Share/home/gby/workgby/peak_cor_test/bamdir"
BAM_GROUPS=(
    "DIR2:/data/zwj/260323CUTRUN/bam/MTA1"
    "DIR3:/data/zwj/260323CUTRUN/bam/H3panAC"
) 
SPECIFIC_PLOT_PAIRS="GFP_MTA1p_DIR2,GFP_H3panAC_DIR3 MTA1_MTA1p_DIR2,MTA1_H3panAC_DIR3 Y13A_MTA1p_DIR2,Y13A_H3panAC_DIR3 singleP_MTA1p_DIR2,singleP_H3panAC_DIR3 F296P_MTA1p_DIR2,F296P_H3panAC_DIR3"
outputdir="/data/zwj/260323CUTRUN/bam/MTA1_H3panAC_new/"

TE_classes_of_interest="LTR,LINE,SINE" 
SAMPLES_FOR_FILTERING="MK9me3,MHDAC2,MH3panAC"
SAMPLES_FOR_GENE_BOXPLOT="GFP_MTA1p_DIR2,MTA1_MTA1p_DIR2,Y13A_MTA1p_DIR2,singleP_MTA1p_DIR2,F296P_MTA1p_DIR2,GFP_H3panAC_DIR3,MTA1_H3panAC_DIR3,Y13A_H3panAC_DIR3,singleP_H3panAC_DIR3,F296P_H3panAC_DIR3"

SAF_FROM_FILE="TRUE" #BED or PEAK or SAF, if more formats are wanted plz contact
FROM_FILE="/Work/lht/lilab_users/25.GBY/peak_cor_test/MGFP_yimin_broad/MGFP_IgG__Ctrl_clean_peaks.broadPeak" #BED or PEAK or SAF, if more formats are wanted plz contact

# ===== advanced modify params===
MIN_TOTAL_NORMALIZED_COUNTS=10
N_LABELS_GENES=20
N_LABELS_TE_FAMILY=20
N_LABELS_TE_REPNAME=40


SLURM_CPUS_PER_TASK=8

# ------0.prepareation---
mkdir -p "${outputdir}/featurecount"
mkdir -p "${outputdir}/results"
#mkdir -p "${outputdir}/merged_results"

#BASE=/Share/home/gby/CUTnRUN/CUTRUN/count_draw/CUTRUN_analysis   # 项目根目录
BASE=/opt/cutrun/gbyscripts/count_draw/CUTRUN_analysis 
if [[ "$reference" =~ ^(hg38|hs|human)$ ]]; then
    REPEAT_ANNO="${BASE}/anno/te_anno_hg38.tsv"
    GENE_ANNO="${BASE}/anno/gene_anno_ensembl_hg38.tsv"
    DENY="${BASE}/anno/hg38-blacklist.v2.bed"
    SAF="${BASE}/anno/gene_te_hg38.saf"
elif [[ "$reference" =~ ^(mm10|mm|mouse)$ ]]; then
    REPEAT_ANNO="${BASE}/anno/te_anno_mm10.tsv"
    GENE_ANNO="${BASE}/anno/gene_anno_ensembl_mm10.tsv"
    DENY="${BASE}/anno/mm10-blacklist.v2.bed"
    SAF="${BASE}/anno/gene_te_mm10.saf"
elif [[ "$reference" =~ ^(mm39)$ ]]; then
    REPEAT_ANNO="${BASE}/anno/te_anno_mm39.tsv"
    GENE_ANNO="${BASE}/anno/gene_anno_ensembl_mm39.tsv"
    DENY="${BASE}/anno/mm39-blacklist.v2.bed"
    SAF="${BASE}/anno/gene_te_mm39.saf"
else
    echo "ERROR: unknown reference '$reference'"; exit 1
fi

shopt -s nullglob
#bam_files=(${bam_dir}/*.bam)
#if [ ${#bam_files[@]} -eq 0 ]; then echo "ERROR: no bam in $bam_dir"; exit 1; fi


# -----1.featurecount----
min_frag=30
max_frag=1200

# 如果上面的SAF_FROM_FILE为TRUE，则覆盖SAF路径，根据文件名后缀判断格式，目前支持BED、PEAK、SAF三种格式，如果需要其他格式请联系我
if [[ "$SAF_FROM_FILE" == "TRUE" ]]; then
    if [[ "$FROM_FILE" =~ \.bed$ ]]; then
        SAF="${outputdir}/featurecount/saffrombed.saf"
        awk 'BEGIN{OFS="\t"} {if($1 !~ /^chr/) $1="chr"$1; print "ENSG"$4,$1,$2,$3,$6}' "$FROM_FILE" > "$SAF"
        GENE_ANNO="${outputdir}/featurecount/geneanno_frombed.tsv"
        echo -e "GeneID\tGeneName\tGeneBiotype" > "$GENE_ANNO"
        awk 'BEGIN{OFS="\t"} {print "ENSG"$4, "chr"$1":"$2"-"$3, "bed_sourced"}' "$FROM_FILE" >> "$GENE_ANNO"
    elif [[ "$FROM_FILE" =~ \Peak$ ]]; then
        SAF="${outputdir}/featurecount/saffrompeak.saf"
        awk 'BEGIN{OFS="\t"} {if($1 !~ /^chr/) $1="chr"$1; print "ENSG"$4,$1,$2,$3,$6}' "$FROM_FILE" > "$SAF"
        GENE_ANNO="${outputdir}/featurecount/geneanno_frompeak.tsv"
        echo -e "GeneID\tGeneName\tGeneBiotype" > "$GENE_ANNO"
        awk 'BEGIN{OFS="\t"} {print "ENSG"$4, "chr"$1":"$2"-"$3, "peak_sourced"}' "$FROM_FILE" >> "$GENE_ANNO"
    elif [[ "$FROM_FILE" =~ \.saf$ ]]; then
        SAF="$FROM_FILE"
        awk 'BEGIN{OFS="\t"} {if($2 !~ /^chr/) $2="chr"$2; print "ENSG"$1,$2,$3,$4,$5}' "$FROM_FILE" > "$SAF"
        GENE_ANNO="${outputdir}/featurecount/geneanno_from_saf.tsv"
        echo -e "GeneID\tGeneName\tGeneBiotype" > "$GENE_ANNO"
        awk 'BEGIN{OFS="\t"} {print "ENSG"$1, $2":"$3"-"$4, "saf_sourced"}' "$FROM_FILE" >> "$GENE_ANNO"
    else
        echo "ERROR: 无法识别的FROM_FILE格式，请使用.bed、.peak或.saf后缀的文件"; exit 1
    fi
    echo ">>> 使用来自文件的SAF: $SAF"
fi



for grp in "${BAM_GROUPS[@]}"; do
    gname="${grp%%:*}"
    gdir="${grp##*:}"
    
    shopt -s nullglob
    bam_files=(${gdir}/*.bam)
    if [ ${#bam_files[@]} -eq 0 ]; then echo "警告: ${gdir} 中没有bam文件，跳过"; continue; fi

    TARGET_FC="${outputdir}/featurecount/fc_${gname}.txt"
    if [ -s "${TARGET_FC}" ]; then
        echo ">>> [缓存命中] 发现已存在非空文件，跳过 featureCounts: ${TARGET_FC}"
    else
        echo ">>> 正在为组别 ${gname} 运行 featureCounts..."
        featureCounts -a "${SAF}" -o "${TARGET_FC}" \
            -F SAF --ignoreDup -p -d ${min_frag} -D ${max_frag} -T ${SLURM_CPUS_PER_TASK} -M --fraction "${bam_files[@]}"
    fi
done

# ---- 2. 生成临时 config_slurm.R ----
COUNTS_TXT="${outputdir}/featurecount/featurecounts.txt"
CONFIG_R="${outputdir}/featurecount/config_slurm.R"
#RUN_R="/Share/home/gby/CUTnRUN/CUTRUN/count_draw/run.R"
#RUN_MULTI_R="/Share/home/gby/CUTnRUN/CUTRUN/count_draw/run_multi.R"
RUN_R="/opt/cutrun/gbyscripts/count_draw/run.R"

TE_R=$(echo $TE_classes_of_interest | tr ',' ' ' | awk '{printf("c(\""); for(i=1;i<NF;i++) printf("%s\",\""),$i; printf("%s\")",$NF)}')
FILT_R=$(echo $SAMPLES_FOR_FILTERING | tr ',' ' ' | awk '{printf("c(\""); for(i=1;i<NF;i++) printf("%s\",\""),$i; printf("%s\")",$NF)}')
BOX_R=$(echo $SAMPLES_FOR_GENE_BOXPLOT | tr ',' ' ' | awk '{printf("c(\""); for(i=1;i<NF;i++) printf("%s\",\""),$i; printf("%s\")",$NF)}')

GROUPS_R="list("
for grp in "${BAM_GROUPS[@]}"; do
    gname="${grp%%:*}"
    GROUPS_R="${GROUPS_R}\"${gname}\"=\"${outputdir}/featurecount/fc_${gname}.txt\","
done
GROUPS_R="${GROUPS_R%,})" # 去掉结尾逗号

# 构建画图配对的 R 列表
PAIRS_R="list("
for pair in $SPECIFIC_PLOT_PAIRS; do
    p1="${pair%,*}"
    p2="${pair#*,}"
    PAIRS_R="${PAIRS_R}c(\"${p1}\",\"${p2}\"),"
done
PAIRS_R="${PAIRS_R%,})"

cat > "${CONFIG_R}" <<EOF
# ===== 自动生成configR =====
#input_files <- list(
#  repeat_annotations_file = "$REPEAT_ANNO",
#  gene_annotations_file   = "$GENE_ANNO",
#  counts_file             = "$COUNTS_TXT",
#  deny_list_file          = "$DENY",
#  saf_file_path           = "$SAF"
#)

repeat_annotations_file = "$REPEAT_ANNO"
gene_annotations_file   = "$GENE_ANNO"
counts_file             = "$COUNTS_TXT"
deny_list_file          = "$DENY"
saf_file_path           = "$SAF"

output_dir       <- "${outputdir}/results"
intermediate_dir <- file.path(output_dir, "intermediate_data")
figures_dir      <- file.path(output_dir, "figures")

MULTI_COUNTS_FILES  <- $GROUPS_R
SPECIFIC_PLOT_PAIRS <- $PAIRS_R

MIN_TOTAL_NORMALIZED_COUNTS <- $MIN_TOTAL_NORMALIZED_COUNTS
N_LABELS_GENES       <- $N_LABELS_GENES
N_LABELS_TE_FAMILY   <- $N_LABELS_TE_FAMILY
N_LABELS_TE_REPNAME  <- $N_LABELS_TE_REPNAME
TE_CLASSES_OF_INTEREST <- $TE_R

SAMPLES_FOR_FILTERING    <- $FILT_R
SAMPLES_FOR_GENE_BOXPLOT <- $BOX_R

EOF

# ---- 3. 运行 ----
#Rscript "${RUN_R}"  "${CONFIG_R}"
echo "start"
#RUN_MULTI_R="/Share/home/gby/CUTnRUN/CUTRUN/count_draw/run_multi.R"
RUN_MULTI_R="/opt/cutrun/gbyscripts/count_draw/run_multi.R"
Rscript "${RUN_MULTI_R}"  "${CONFIG_R}"


echo "======  Job finished at $(date)  ======"
