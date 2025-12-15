#!/bin/bash


if [ $# -ne 2 ]; then
    echo "用法: $0 <SRR文件夹所在路径> <输出路径>"
    echo "示例: $0 /Data/gby/antigen/mm/data /Data/gby/antigen/mm/output"
    exit 1
fi

SRR_BASE_DIR="$1"
OUTPUT_BASE_DIR="$2"

if [ ! -d "$SRR_BASE_DIR" ]; then
    echo "错误: 输入目录 $SRR_BASE_DIR 不存在!"
    exit 1
fi

mkdir -p "$OUTPUT_BASE_DIR" || { echo "错误: 无法创建输出目录 $OUTPUT_BASE_DIR"; exit 1; }

find "$SRR_BASE_DIR" -maxdepth 1 -type d -name "SRR*" | while read -r srr_dir; do
    srr_id=$(basename "$srr_dir")
    
    sra_file=$(find "$srr_dir" -maxdepth 1 -type f -name "*.sra" -o -name "*.sralite" | head -n 1)
    
    if [ -n "$sra_file" ] && [ -f "$sra_file" ]; then
        output_dir="$OUTPUT_BASE_DIR/$srr_id"
        mkdir -p "$output_dir"
        
        echo "fasterq-dump --split-files -e 1 --outdir '$output_dir' '$sra_file'"
    else
        echo "警告: 在 $srr_dir 中未找到.sra文件，跳过该目录" >&2
    fi
done | parallel -j 5 

echo "dump complete！"

