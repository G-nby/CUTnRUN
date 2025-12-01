#!/bin/bash

set -e

usage() {
  echo "Usage: $0 --trim [true|false] --reference [hs/mm] --outpath <output_dir> --datapath <data_dir> --num_threads <num>"
  exit 1
}

# default params
trim="false"
reference="hs"
outpath=""
datapath=""
num_threads=8

echo "${trim}"

# get params
while [[ $# -gt 0 ]]; do
  case $1 in
    --trim)
      trim="true"
      shift
      ;;
    --reference)
      reference="$2"
      shift 2
      ;;
    --outpath)
      outpath="$2"
      shift 2
      ;;
    --datapath)
      datapath="$2"
      shift 2
      ;;
    --num_threads)
      num_threads="$2"
      shift 2
      ;;
    --*)
      echo "Unknown parameter: $1"
      usage
      ;;
    *)
      # if datapath is not specified
      if [[ -z "$datapath" ]]; then
        datapath="$1"
        shift
      else
        echo "Unexpected positional argument: $1"
        usage
      fi
      ;;
  esac
done

if [[ -z "$datapath" ]]; then
  echo "Error: datapath must be specified either as --datapath or positional argument."
  usage
fi


script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "${trim}"

if [[ "$trim" == "notrim" || "$trim" == "false" ]]; then
  input_path="$datapath"
else
  # trim mode
  clean_dir="${datapath}/clean"
  echo
  mkdir -p "$clean_dir"

  echo "Starting trimming with trim_galore..."
  
  for sample_dir in "$datapath"/*; do
    if [[ -d "$sample_dir" ]]; then
      sample_name=$(basename "$sample_dir")
      echo "Processing sample: $sample_name"

      fq_files=()
      while IFS= read -r -d $'\0' file; do
        fq_files+=("$file")
      done < <(find "$sample_dir" -maxdepth 1 -type f \( -iname "*.fq" -o -iname "*.fastq" -o -iname "*.fq.gz" -o -iname "*.fastq.gz" \) -print0)

      if [[ ${#fq_files[@]} -eq 0 ]]; then
        echo "Warning: No fastq files found in $sample_dir, skipping..."
        continue
      fi

      if [[ ${#fq_files[@]} -eq 1 ]]; then
        # single-end
        file_prefix=$(basename "${fq_files[0]}")
        file_prefix=${file_prefix%%.f*q*}
        echo "File prefix=${file_prefix}"
        echo "Only one sequencing data file detected, single-end mode applied."
        trim_galore --cores 4 --output_dir "$clean_dir" "${fq_files[0]}" --dont_gzip
        echo "${file_prefix} trim done"
  
      elif [[ ${#fq_files[@]} -eq 2 ]]; then
        # paired-end
        read1=""
        read2=""
        for fq in "${fq_files[@]}"; do
          fname=$(basename "$fq")
          if [[ "$fname" =~ _1 ]]; then
            read1="$fq"
            file_prefix=${fname%%_1.f*q*}
          elif [[ "$fname" =~ _2 ]]; then
            read2="$fq"
            file_prefix=${fname%%_2.f*q*}
          fi
        done
  
        if [[ -z "$read1" || -z "$read2" ]]; then
          echo "Could not identify both _1 and _2 files in $sample_dir, skipping..."
          continue
        fi
  
        echo "File prefix=${file_prefix}"
        echo "Two sequencing data files detected, paired-end mode applied."
        trim_galore --cores 4 --paired --output_dir "$clean_dir" "$read1" "$read2" --dont_gzip
        echo "${file_prefix} trim done"
  
      else
        echo "Abnormal number of FASTQ files in $sample_dir (${#fq_files[@]} files), skipping..."
        continue
      fi
    fi
  done

  input_path="$clean_dir"
fi

echo "Running TEwithSalmon.py..."

echo "python $script_dir/TEwithSalmon.py quant --reference=$reference --outpath=$outpath --num_threads=$num_threads $input_path --exprtype count"
python "$script_dir/TEwithSalmon.py" quant --reference="$reference" --outpath="$outpath" --num_threads="$num_threads" "$input_path" --exprtype count 
