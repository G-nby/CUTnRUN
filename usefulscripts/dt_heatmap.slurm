#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH -p bio
#SBATCH --cpus-per-task=20
#SBATCH --output=DeepTools-%A.out
#SBATCH -e DeepTools-%A.err

computeMatrix reference-point \
--missingDataAsZero \
--sortUsing median \
--skipZeros -p 20 \
-R /Share/home/gby/CUTnRUN/assist/data/HEK293_TNRC18_common_peaks.bed \
-S /Data/gby/cutnrun/Data/SRR18761785/SRR18761785_clean_Filtered_10bp_rpgc.bw /Data/gby/cutnrun/Data/SRR18761784/SRR18761784_clean_Filtered_10bp_rpgc.bw \
-o ../Results/SETDB1GFP_Tcontrol3.gz \
-a 5000 -b 5000 \
--referencePoint center \
--samplesLabel HEK293_T18 HEK293_GFP_T18

plotProfile -m ../Results/SETDB1GFP_Tcontrol3.gz -out ../Results/SETDB1GFP_TcontrolpltProfile3.pdf --plotFileFormat 'pdf' --perGroup --colors black red --plotTitle SETDB1GFP_Tcontrol --refPointLabel center

plotHeatmap -m ../Results/SETDB1GFP_Tcontrol3.gz -out ../Results/SETDB1GFP_TcontrolpltHeatmap3.pdf --plotFileFormat 'pdf' --colorMap YlGnBu --plotTitle SETDB1GFP_Tcontrol --refPointLabel center

