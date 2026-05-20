#!/bin/bash
#SBATCH --job-name=cutrun_ana
#SBATCH --output=%j.out 
#SBATCH --error=%j.err 
#SBATCH --cpus-per-task=8      
#SBATCH --partition=bio

# ======== the params you shoudf modify when use ======
reference=hg38                    
#bam_dir="/Share/home/gby/workgby/peak_cor_test/bamdir"
DATA_GROUPS=(
    "WT:/Share/home/gby/workgby/peak_cor_test/bamdir_wt"
    "GFP:/Share/home/gby/workgby/peak_cor_test/bamdir_gfp"
) 
SPECIFIC_PLOT_PAIRS="MHDAC2_WT,MH3panAC_WT MHDAC2_GFP,MH3panAC_GFP MHDAC2_WT,MHDAC2_GFP"
outputdir="/Share/home/gby/workgby/peak_cor_test/bamdirout"

TE_classes_of_interest="LTR,LINE,SINE" 
SAMPLES_FOR_FILTERING="MK9me3,MHDAC2,MH3panAC"
SAMPLES_FOR_GENE_BOXPLOT="MHDAC2,MGFP,MK9me3,MH3panAC"

# ===== advanced modify params===
MIN_TOTAL_NORMALIZED_COUNTS=10
N_LABELS_GENES=20
N_LABELS_TE_FAMILY=20
N_LABELS_TE_REPNAME=40


# ------0.prepareation---
mkdir -p "${outputdir}/featurecount"
mkdir -p "${outputdir}/results"
mkdir -p "${outputdir}/merged_results"

BASE=/Share/home/gby/CUTnRUN/CUTRUN/count_draw/CUTRUN_analysis   # 项目根目录
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
bam_files=(${bam_dir}/*.bam)
if [ ${#bam_files[@]} -eq 0 ]; then echo "ERROR: no bam in $bam_dir"; exit 1; fi


# -----1.featurecount----
min_frag=30
max_frag=1200


for grp in "${BAM_GROUPS[@]}"; do
    gname="${grp%%:*}"
    gdir="${grp##*:}"
    
    shopt -s nullglob
    bam_files=(${gdir}/*.bam)
    if [ ${#bam_files[@]} -eq 0 ]; then echo "警告: ${gdir} 中没有bam文件，跳过"; continue; fi
    
    echo ">>> 正在为组别 ${gname} 运行 featureCounts..."
    featureCounts -a "${SAF}" -o "${outputdir}/featurecount/fc_${gname}.txt" \
        -F SAF --ignoreDup -p -d ${min_frag} -D ${max_frag} -T ${SLURM_CPUS_PER_TASK} -M --fraction "${bam_files[@]}"
done

# ---- 2. 生成临时 config_slurm.R ----
COUNTS_TXT="${outputdir}/featurecount/featurecounts.txt"
CONFIG_R="${outputdir}/featurecount/config_slurm.R"
#RUN_R="/Share/home/gby/CUTnRUN/CUTRUN/count_draw/run.R"
RUN_MULTI_R="/Share/home/gby/CUTnRUN/CUTRUN/count_draw/run_multi.R"

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

cache_files <- list(
  clean_counts = file.path(intermediate_dir, "01_clean_counts.rds"),
  plot_data    = file.path(intermediate_dir, "02_plot_data.rds")
)
EOF

# ---- 3. 运行 ----
#Rscript "${RUN_R}"  "${CONFIG_R}"
RUN_MULTI_R="/Share/home/gby/CUTnRUN/CUTRUN/count_draw/run_multiigg.R"
Rscript "${RUN_MULTI_R}"  "${CONFIG_R}"


echo "======  Job finished at $(date)  ======"
