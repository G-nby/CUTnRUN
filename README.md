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
- [zhaolab server usage](#zhaolab-server-usage)


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
~~**注意：目前仅有 `hg38` 和 `mm39` 参考基因组可使用， `mm10` 将随后更新**~~
参考基因组 `hg38` 和 `mm39` 和 `mm10` 参考基因组可使用。


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
复制 `dump.slurm` 到自己的文件夹并修改参数，使用如下命令提交
```bash
sbatch my_dump.slurm
```


# zhaolab server usage

## 目录 for zhaolab server
- [CUTRUN](#zhaolabcutrun)
  - [CUTRUN 基本 pipeline 用法（common/ecoli）](#zhaolabcutrun-基本pipeline用法commonecoli)
  - [CUTRUN spikein](#zhaolabcutrun-spikein)
  - [CUTRUN count_draw](#zhaolabcutrun-countdraw)
- [RNAseq](#zhaolabrnaseq)
  - [使用 TEsalmon2 处理 TE 数据](#zhaolab使用tesalmon2处理te数据)
  - [使用 TEsalmon2 处理 regular gene 数据](#zhaolab使用tesalmon2处理regular-gene数据)
- [其他工具脚本](#zhaolab其他工具脚本)
  - [dt_heatmap](#zhaolabdt_heatmap)
  - [cobed](#zhaolabcobed)
  - [bw_cor/bam_cor](#zhaolabbw_corbam_cor)
  - [homer](#zhaolabhomer)
  - [pca](#zhaolabpca)
  - [dump](#zhaolabdump)
 
## zhaolabCUTRUN
### zhaolabCUTRUN 基本pipeline用法（common/ecoli）
运行以下命令激活环境
```bash
conda activate cutrun
```
每个样本的数据位于一个文件夹中  
根据测序数据的文件种类选择命令 `CUTRUN`或`CUTRUNgz`，根据实际情况选择参数，运行下方命令进行数据处理（以CUTRUN为例）
```bash
CUTRUN trim/notrim ensembl_hg38/ensembl_mm10
```
common处理使用 `CUTRUN`或`CUTRUNgz` 命令即可
评估ecoli污染使用 `CUTRUN_ecoli`或`CUTRUNgz_ecoli` 命令，用法相同
### zhaolabCUTRUN spikein
每个样本单独一个文件夹。所有样本文件夹位于一个总文件夹中  
先使用 `SPIKEIN_1`或 `SPIKEINgz_1` 对每个样本进行数据处理（用法同上），然后在总文件夹中，运行如下命令
```bash
SPIKEIN_2
```
根据提示选择基因组、选择windowsize、选择ctrl group（选多个的话用空格分割，其他file会分别和每个ctrl做一个比较），然后等待结果生成。结果bigwig文件在大文件夹下 `SPIKE_IN_bigwigFile` 文件夹里  
### zhaolabCUTRUN count&draw
运行以下命令激活环境
```bash
conda activate TEsalmon2
```
由于参数较多，可使用如下命令复制脚本到目标文件夹使用。
```bash
cp /opt/cutrun/gbyscripts/count_draw/count_draw.slurm .
```
cutrun处理过的bamfile放在一个文件夹中，填写到bamdir参数处。根据自己的样本和其他需要修改 `count_draw.slurm` 中的参数，然后运行如下命令进行提交
```bash
bash my_count_draw.slurm
```
结果可在输出文件夹的 `result` 文件夹中查看。  
~~**注意：目前仅有 `hg38` 和 `mm39` 参考基因组可使用， `mm10` 将随后更新**~~
参考基因组 `hg38` 和 `mm39` 和 `mm10` 参考基因组可使用。


## zhaolabRNAseq
### zhaolab使用TEsalmon2处理TE数据
运行以下命令激活环境
```bash
conda activate TEsalmon
```
- raw data
每个样本数据位于一个小文件夹中，使用如下命令进行处理，根据自己需求修改参数。
```bash
TEsalmon --reference hs/mm --outpath "/path/to/outputdir" "path/to/datadir"
```
最后的参数填写输入数据文件夹（包含所有样本所在小文件夹的大文件夹），`--outpath` 处填写希望结果输出的路径，`--reference` 处根据样本种类填写hs或mm，然后提交即可
- clean data
所有clean数据（fq或gz）放于同一个文件夹中
使用如下命令进行处理，根据自己需求修改参数。  
文件最后（无参数标记的输入）填写输入数据文件夹（包含所有样本的文件夹），`--outpath` 处填写希望结果输出的路径，`--reference` 处根据样本种类填写hs或mm，然后使用 `sbatch` 命令提交即可
```bash
TEsalmon2 quant \
        --reference=hs/mm \
        --outpath=/path/to/outputdir --num_threads=8 \
        path/to/datadir --exprtype count
```
处理完成后，在输出文件夹中找到 `conditions.csv` 文件，根据样本分组在`condition`列填写分组信息（名称自定义，同组样本填写相同字段）
- 有重复
使用下方命令进行处理，根据自己需求修改参数。  
`--inpath` 处填写输入数据文件夹（上一步的输出文件夹，路径下应有`EXPR.csv`等文件），`--outpath` 处填写希望分析结果输出的路径，`--conditions` 处填写希望进行差异分析的两个group（填写condition列的字段），**ctrl组在前，处理组在后**，根据需要使用 `--log2fc_min` 和 `--pmax` 参数调整log2fc和pvalue的阈值，然后使用 `sbatch` 命令提交即可。可不填写 `--log2fc_min` 和 `--pmax` 参数，默认值分别为 `0.5` 和 `0.05` .
```bash
TEsalmon2 test --inpath=path/to/inputdir \
        --outpath=/path/to/outputdir \
        --table=tsv \
        --figtype=png \
        --analysis_type=DE \
        --conditions=ctrl,treat \
        --log2fc_min=0.5 --pmax=0.05
```
- 无重复，组中单样本
使用下方命令进行处理，根据自己需求修改参数。  
`--inpath` 处填写输入数据文件夹（上一步的输出文件夹，路径下应有`EXPR.csv`等文件），`--outpath` 处填写希望分析结果输出的路径，`--samples` 处填写希望进行差异分析的两个sample（填写sample列的字段），**ctrl在前，处理在后**，`--ref_name` 处根据样本种类填写hs或mm，根据需要使用 `--log2fc_min` 和 `--pmax` 参数调整log2fc和pvalue的阈值，可不填写 `--log2fc_min` 和 `--pmax` 参数，默认值分别为 `0.5` 和 `0.05`.然后命令提交即可
```bash
TEsalmon2 test_single --inpath=path/to/inputdir \
        --outpath=/path/to/outputdir \
        --table=tsv \
        --figtype=png \
        --samples=ctrl1,treat1 \
        --ref_name=hs/mm
        --log2fc_min=0.5 --pmax=0.05
```
### zhaolab使用TEsalmon2处理regular gene数据
#### zhaolab上游处理
运行以下命令激活环境
```bash
conda activate cutrun
```
每个样本的数据位于一个小文件夹中  
在数据文件夹，运行 `RNA_RE` 或 `RNA_REgz` 命令进行数据处理（以RNA_RE为例），根据实际情况选择参数
```bash
RNA_RE trim/notrim ensembl_hg38/ensembl_mm10
```
#### zhaolab下游处理
运行以下命令激活环境
```bash
conda activate TEsalmon
```
使用如下命令（修改参数）进行处理。  
`--inpath` 处填写数据所处的大文件夹，`--outpath` 处填写希望数据输出的路径，`--ref_name` 处根据样本种类填写hs或mm，然后使用 `sbatch` 命令提交即可
```bash
mk_expr --inpath=path/to/datadir --outpath=/path/to/outdir --num_threads=8 --ref_name=hs/mm
```
处理完成后，在输出文件夹中找到 `conditions.csv` 文件，根据样本分组在`condition`列填写分组信息（名称自定义，同组样本填写相同字段）
- 有重复
使用下方命令进行处理，根据自己需求修改参数。  
`--inpath` 处填写输入数据文件夹（上一步的输出文件夹，路径下应有`EXPR.csv`等文件），`--outpath` 处填写希望分析结果输出的路径，`--conditions` 处填写希望进行差异分析的两个group（填写condition列的字段），**ctrl组在前，处理组在后**，根据需要使用 `--log2fc_min` 和 `--pmax` 参数调整log2fc和pvalue的阈值，然后使用 `sbatch` 命令提交即可。可不填写 `--log2fc_min` 和 `--pmax` 参数，默认值分别为 `0.5` 和 `0.05` .
```bash
TEsalmon2 test --inpath=path/to/inputdir \
        --outpath=/path/to/outputdir \
        --table=tsv \
        --figtype=png \
        --analysis_type=DE \
        --conditions=ctrl,treat \
        --log2fc_min=0.5 --pmax=0.05
```
- 无重复，组中单样本
使用下方命令进行处理，根据自己需求修改参数。  
`--inpath` 处填写输入数据文件夹（上一步的输出文件夹，路径下应有`EXPR.csv`等文件），`--outpath` 处填写希望分析结果输出的路径，`--samples` 处填写希望进行差异分析的两个sample（填写sample列的字段），**ctrl在前，处理在后**，`--ref_name` 处根据样本种类填写hs或mm，根据需要使用 `--log2fc_min` 和 `--pmax` 参数调整log2fc和pvalue的阈值，可不填写 `--log2fc_min` 和 `--pmax` 参数，默认值分别为 `0.5` 和 `0.05`.然后命令提交即可
```bash
TEsalmon2 test_single --inpath=path/to/inputdir \
        --outpath=/path/to/outputdir \
        --table=tsv \
        --figtype=png \
        --samples=ctrl1,treat1 \
        --ref_name=hs/mm
        --log2fc_min=0.5 --pmax=0.05
```

## zhaolab其他工具脚本
位于 `/opt/gbyscripts` 文件夹下，使用 `conda activate cutrun` 激活的环境即可
### zhaolabdt_heatmap
（复制到自己的文件夹中，修改参数）vim进入文件，  
修改 `computeMatrix` 模块的 `-R`（参考位置，可提交bed或peak文件等格式） `-S`（需要画图的样本的peak文件或者bw文件等） `-o`（数据文件位置，xx.gz） `--samplesLabel`（样本命名）参数，  
修改 `plotProfile` `plotHeatmap`模块的 `-m`（上一部分的数据文件，xx.gz） `-out`（图片输出） `--plotTitle`（图片标题）参数，  
然后使用 `bash` 命令提交即可
```bash
bash my_dt_heatmap.slurm
```
### zhaolabcobed
使用如下命令提交
```bash
cobed2 <fileA> <fileB> <output_directory>
```
如有多个样本，希望寻找共有overlap，请使用
```bash
cobed3 <filelist>
```
其中， `filelist` 为一个文件列表，是想要进行overlap分析的多个peak文件的路径
### zhaolabbw_cor/bam_cor
使用如下命令提交
```bash
bw_cor <bw_dir> <outprefix> <binsize>
bam_cor <bam_dir> <outprefix> <binsize>
```
`outprefix` 为输出文件命名前缀，binsize为自定义窗口大小，默认为1000bp
### zhaolabhomer
使用如下命令提交
```bash
homer_anno <input_file> <ref>
```
`ref` 为参考物种基因组名称，hg38/mm10为默认注释。可使用TEhg38/TEmm10进行TE注释
### ~~zhaolabpca~~
~~（复制到自己的文件夹中，修改参数）vim进入文件，修改语句最后提交文件为自己的 `EXPR.csv` 文件路径，然后使用 `sbatch` 命令提交即可。~~  
~~使用如下命令提交~~
~~```bash~~
~~PCA <input_file>~~
~~```~~
计划更新normlize版本并加入到TEsalmon处理中
### zhaolabdump
使用如下命令提交
```bash
DUMP <SRR文件夹所在路径> <输出路径>
```

