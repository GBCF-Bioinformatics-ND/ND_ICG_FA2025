#!/bin/bash

#$ -M <netid>@nd.edu
#$ -m abe
#$ -q debug
#$ -N read_qc

# load the bio/2.0 module from CRC or intall fastqc

module load bio/2.0

# --- Configuration ---
# Set default input and results directories.

input_dir="../input"
results_base_dir="../results"
fastqc_output_subdir="fastqc_raw_reads" # Specific subdirectory for FastQC results


# Define the full path for FastQC specific results
fastqc_results_dir="${results_base_dir}/${fastqc_output_subdir}"

# 1. Prepare output directory
mkdir -p "$fastqc_results_dir"


# 2. Run FastQC directly from the input directory, outputting to the designated results folder
cd "$input_dir" || { echo "Error: Could not change to input directory '$input_dir'. Exiting."; exit 1; }

# FastQC will place its .zip and .html files directly into fastqc_results_dir
fastqc -o "$fastqc_results_dir" *.fastq*


# 3. Unzip FastQC reports and consolidate summaries
cd "$fastqc_results_dir"

# Check if there are any zip files before looping
for filename in *.zip; do
    echo "  - Unzipping '$filename'"
    unzip -q "$filename" # -q for quiet output
done
