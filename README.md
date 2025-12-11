# CUTnRUN
backup basic pipelines on acluster

## 目录
- [CUTRUN](#cutrun)
  - [CUTRUN 基本 pipeline 用法（common/ecoli）](#cutrun-基本pipeline用法commonecoli)
  - [CUTRUN spikein](#cutrun-spikein)
  - [CUTRUN count_draw](#cutrun-countdraw)
- [RNAseq](#rnaseq)
  - [使用 TEsalmon2 处理 TE 数据](#使用tesalmon2处理te数据)
  - [使用 TEsalmon2 处理 regular gene 数据](#使用tesalmon2处理regular-gene数据)
- [其他工具脚本](#其他工具脚本)
  - [dt_heatmap](#dt_heatmap)
  - [cobed](#cobed)
  - [bw_cor/bam_cor](#bw_corbam_cor)
  - [homer](#homer)
  - [pca](#pca)
  - [dump](#dump)


## CUTRUN
### CUTRUN 基本pipeline用法（common/ecoli）
运行以下命令激活环境
```bash
source /Share/samples/Acluster.sh
```
每个样本的数据位于一个文件夹中  
脚本位于`./CUTRUN` 文件夹中，复制 `1.xxx.sh` 到数据文件夹，根据实际情况选择参数，运行下方命令进行数据处理（以1.slurm_chipseq_common_gz.sh为例）
```bash
sbatch 1.slurm_chipseq_common_gz.sh trim/notrim ensembl_hg38/ensembl_mm10
```
common处理使用 `./CUTRUN` 文件夹脚本即可
评估ecoli污染使用 `./CUTRUN/ecoli` 文件夹中的脚本，用法相同
### CUTRUN spikein
脚本位于`./CUTRUN/spikein` 文件夹中，每个样本单独一个文件夹。所有样本文件夹位于一个总文件夹中  
先使用 `1.xxx.sh` 对每个样本进行数据处理，然后将 `2.spikeinscale` 复制到总文件夹中。运行如下命令
```bash
python 2.spikeinscale
```
根据提示选择基因组、选择windowsize、选择ctrl group（选多个的话用空格分割，其他file会分别和每个ctrl做一个比较），然后等待结果生成。结果bigwig文件在大文件夹下 `SPIKE_IN_bigwigFile` 文件夹里  
### CUTRUN count&draw
运行以下命令激活环境
```bash
conda activate TEsalmon2
source ~/Acluster.sh
```
脚本位于`./CUTRUN/count_draw` 文件夹中。将 `count_draw.slurm` 复制到自己需要的地方，cutrun处理过的bamfile放在一个文件夹中，填写到bamdir参数处。
根据自己的样本和其他需要修改 `count_draw.slurm` 中的参数，然后运行如下命令进行提交
```bash
sbatch my_count_draw.slurm
```
结果可在输出文件夹的 `result` 文件夹中查看。


## RNAseq
### 使用TEsalmon2处理TE数据
运行以下命令激活环境
```bash
source ~/Acluster.sh
conda activate TEsalmon2
```
- raw data
每个样本数据位于一个小文件夹中，使用 `./TEsalmon2/scripts/quanttest.slurm` 进行处理（复制到自己的文件夹中，修改参数）。  
文件最后填写输入数据文件夹（包含所有样本所在小文件夹的大文件夹），`--outpath` 处填写希望结果输出的路径，`--reference` 处根据样本种类填写hs或mm，然后使用 `sbatch` 命令提交即可
```bash
sbatch my_quanttest.slurm
```
- clean data
所有clean数据（fq或gz）放于同一个文件夹中
使用 `./TEsalmon2/scripts/quantcleantest.slurm`进行处理（复制到自己的文件夹中，修改参数）。  
文件最后（无参数标记的输入）填写输入数据文件夹（包含所有样本的文件夹），`--outpath` 处填写希望结果输出的路径，`--reference` 处根据样本种类填写hs或mm，然后使用 `sbatch` 命令提交即可
```bash
sbatch quantcleantest.slurm
```
处理完成后，在输出文件夹中找到 `conditions.csv` 文件，根据样本分组在`condition`列填写分组信息（名称自定义，同组样本填写相同字段）
- 有重复
使用 `./TEsalmon2/scripts/test1115.slurm` 进行处理（复制到自己的文件夹中，修改参数）。  
`--inpath` 处填写输入数据文件夹（上一步的输出文件夹，路径下应有`EXPR.csv`等文件），`--outpath` 处填写希望分析结果输出的路径，`--conditions` 处填写希望进行差异分析的两个group（填写condition列的字段），**ctrl组在前，处理组在后**，根据需要使用 `--log2fc_min` 和 `--pmax` 参数调整log2fc和pvalue的阈值，然后使用 `sbatch` 命令提交即可
```bash
sbatch my_test.slurm
```
- 无重复，组中单样本
使用 `./TEsalmon2/scripts/cor.slurm`进行处理（复制到自己的文件夹中，修改参数）。  
`--inpath` 处填写输入数据文件夹（上一步的输出文件夹，路径下应有`EXPR.csv`等文件），`--outpath` 处填写希望分析结果输出的路径，`--samples` 处填写希望进行差异分析的两个sample（填写sample列的字段），**ctrl在前，处理在后**，`--ref_name` 处根据样本种类填写hs或mm，根据需要使用 `--log2fc_min` 和 `--pmax` 参数调整log2fc和pvalue的阈值，然后使用 `sbatch` 命令提交即可
```bash
sbatch my_cor.slurm
```
### 使用TEsalmon2处理regular gene数据
#### 上游处理
运行以下命令激活环境
```bash
source /Share/samples/Acluster.sh
```
每个样本的数据位于一个小文件夹中  
脚本位于`./RNAseq` 文件夹中，根据实际情况选择基因组，复制 `1.xxx.sh` 到数据文件夹，运行下方命令进行数据处理（以1.slurm_rnaseq_GRCh38_gz.sh为例）
```bash
sbatch 1.slurm_rnaseq_GRCh38_gz.sh trim/notrim 
```
#### 下游处理
重开窗口，并运行以下命令激活环境
```bash
source ~/Acluster.sh
conda activate TEsalmon2
```
使用 `./TEsalmon2/scripts/mk_expr.slurm`（复制到自己的文件夹中，修改参数）进行处理。  
`--inpath` 处填写数据所处的大文件夹，`--outpath` 处填写希望数据输出的路径，`--ref_name` 处根据样本种类填写hs或mm，然后使用 `sbatch` 命令提交即可
```bash
sbatch my_mk_expr.slurm
```
处理完成后，在输出文件夹中找到 `conditions.csv` 文件，根据样本分组在`condition`列填写分组信息（名称自定义，同组样本填写相同字段）
- 有重复
使用 `./TEsalmon2/scripts/test1115.slurm` 进行处理（复制到自己的文件夹中，修改参数）。  
`--inpath` 处填写输入数据文件夹（上一步的输出文件夹，路径下应有`EXPR.csv`等文件），`--outpath` 处填写希望分析结果输出的路径，`--conditions` 处填写希望进行差异分析的两个group（填写condition列的字段），**ctrl组在前，处理组在后**，根据需要使用 `--log2fc_min` 和 `--pmax` 参数调整log2fc和pvalue的阈值，然后使用 `sbatch` 命令提交即可
```bash
sbatch my_test.slurm
```
- 无重复，组中单样本
使用 `./TEsalmon2/scripts/cor.slurm`进行处理（复制到自己的文件夹中，修改参数）。  
`--inpath` 处填写输入数据文件夹（上一步的输出文件夹，路径下应有`EXPR.csv`等文件），`--outpath` 处填写希望分析结果输出的路径，`--samples` 处填写希望进行差异分析的两个sample（填写sample列的字段），**ctrl在前，处理在后**，`--ref_name` 处根据样本种类填写hs或mm，根据需要使用 `--log2fc_min` 和 `--pmax` 参数调整log2fc和pvalue的阈值，然后使用 `sbatch` 命令提交即可
```bash
sbatch my_cor.slurm
```

## 其他工具脚本
位于 `./usefulscripts` 文件夹下，使用 `source /Share/samples/Acluster.sh` 激活的环境即可
### dt_heatmap
（复制到自己的文件夹中，修改参数）vim进入文件，  
修改 `computeMatrix` 模块的 `-R`（参考位置，可提交bed或peak文件等格式） `-S`（需要画图的样本的peak文件或者bw文件等） `-o`（数据文件位置，xx.gz） `--samplesLabel`（样本命名）参数，  
修改 `plotProfile` `plotHeatmap`模块的 `-m`（上一部分的数据文件，xx.gz） `-out`（图片输出） `--plotTitle`（图片标题）参数，  
然后使用 `sbatch` 命令提交即可
```bash
sbatch my_dt_heatmap.slurm
```
### cobed
使用如下命令提交
```bash
sbatch cobed2.slurm <fileA> <fileB> <output_directory>
```
如有多个样本，希望寻找共有overlap，请使用
```bash
sbatch cobed3.slurm <filelist>
```
其中， `filelist` 为一个文件列表，是想要进行overlap分析的多个peak文件的路径
### bw_cor/bam_cor
使用如下命令提交
```bash
sbatch bw_cor.slurm <bw_dir> <outprefix> <binsize>
sbatch bam_cor.slurm <bam_dir> <outprefix> <binsize>
```
`outprefix` 为输出文件命名前缀，binsize为自定义窗口大小，默认为1000bp
### homer
使用如下命令提交
```bash
sbatch homer.slurm <input_file> <ref>
```
`ref` 为参考物种基因组名称，hg38/mm10为默认注释。可使用TEhg38/TEmm10进行TE注释
### ~~pca~~
~~（复制到自己的文件夹中，修改参数）vim进入文件，修改语句最后提交文件为自己的 `EXPR.csv` 文件路径，然后使用 `sbatch` 命令提交即可。~~  
计划更新normlize版本并加入到TEsalmon处理中
### dump
使用如下命令提交
```bash
sbatch homer.slurm <SRR文件夹所在路径> <输出路径>
```
