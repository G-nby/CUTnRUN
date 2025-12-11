# CUT&RUN数据featureCounts（Gene+TE），下游绘图
## 请将"CUT&RUN_analysis"文件夹下载到本地（里面包括各种注释文件）

1.  **安装R和R studio**
2.  **准备文件**
    经过上游bam文件 -> featureCounts得到featurecounts.txt文件，将其下载到本地，放在"CUT&RUN_analysis"文件夹下。
    **目录结构如下**：
    CUTRUN_analysis
    ├── analysis/ 
    │   ├── scripts/
    │   └── config.R
    │   └── run.R
    ├── anno/ 
    ├── results/ 
    ├── featurecounts.txt/ 

3.  **修改参数和路径**
    打开config.R文件，修改output_file为输出路径。其它路径可以按需改变，也可以不改（如果文件名和相对路径是正确的）。
    ！！如果你的样本是小鼠样本，那么将所有“hg38”替换成“mm39”！！
    修改分析参数！！特别是箱线图：
    SAMPLES_FOR_FILTERING <- c("MMp", "MHDAC2", "MK20me2")，会筛选第一和第三高，但是第二个低的基因/te，绘制箱形图。
    SAMPLES_FOR_GENE_BOXPLOT <- c("MMp", "MK20me2", "MHDAC2", "MGFP", "MK9me3", "MH3pan")，改成你所有的样品，或者你想看的那些样品。
    **ctrl + s保存！！**
4. **运行**
    打开run.R文件。
    在终端使用setwd("/your path to/analysis")将analysis这个文件夹设置为工作目录。你也可以在右下角进入analysis文件夹，点击"More/Set as Working Directory"。找不到你的文件夹，就点击这个框右侧的省略号。
    
    **直接在Rstudio的最上栏选择：Code/Run region/Run All**