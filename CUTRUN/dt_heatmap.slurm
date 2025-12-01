#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --partition=Acluster
#SBATCH --cpus-per-task=20
#SBATCH --output=DeepTools-%A.out

module load 

computeMatrix reference-point \
--missingDataAsZero \
--sortUsing median \
--skipZeros -p 20 \
-R HEK293_WT_T18_peaks.broadPeak \
-S HEK293_WT_T18_AlignedSortedChrMrmDuprmBLrmFiltered.bw HEK293_SETDB1KO_T18_AlignedSortedChrMrmDuprmBLrmFiltered.bw \
-o SETDB1KO_WTcontrol.gz \
-a 5000 -b 5000 \
--referencePoint center \
--samplesLabel HEK293WT_T18 HEK293_SETDB1KO_T18

plotProfile -m SETDB1KO_WTcontrol.gz -out SETDB1KO_WTcontrolpltProfile.pdf --plotFileFormat 'pdf' --perGroup --colors black red --plotTitle SETDB1KO_WTcontrol --refPointLabel PKC

plotHeatmap -m SETDB1KO_WTcontrol.gz -out SETDB1KO_WTcontrolpltHeatmap.pdf --plotFileFormat 'pdf' --colorMap YlGnBu --plotTitle SETDB1KO_WTcontrol --refPointLabel PKC

