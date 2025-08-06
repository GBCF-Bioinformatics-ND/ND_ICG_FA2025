#!/bin/bash

#$ -M <netid>@nd.edu
#$ -m abe
#$ -q debug
#$ -N read_trimming

# load the bio/2.0 module from CRC or intall fastqc

module load bio/2.0

# --- Configuration ---
input_dir="../input"
results_base_dir="../results"
trimmomatic_output_subdir="trimmed_reads"     # Specific subdirectory for Trimmomatic results
TRIMMOMATIC_ADAPTERS="../input/NexteraPE-PE.fa" # <<< ADJUST THIS PATH!

trimmomatic_results_dir="${results_base_dir}/${trimmomatic_output_subdir}"

# 1. Prepare Trimmomatic output directory
mkdir -p "$trimmomatic_results_dir"

# 2. Navigate to the input directory for Trimmomatic

cd "$input_dir" || { echo "Error: Could not change to input directory '$input_dir'. Exiting."; exit 1; }

# 3. Run Trimmomatic by Looping through the forward reads (*_1.fastq.gz)
for infile in *_1.fastq.gz; do
    # Extract base name (e.g., 'sampleA' from 'sampleA_1.fastq.gz')
    base=$(basename "${infile}" _1.fastq.gz)

    # Define paths for output files in the trimmomatic_results_dir
    forward_read_in="${infile}"
    reverse_read_in="${base}_2.fastq.gz" # Assumes _2.fastq.gz naming for reverse reads

    # Check if the reverse read file exists
    if [ ! -f "$reverse_read_in" ]; then
        echo "Warning: Reverse read file '${reverse_read_in}' not found for '${forward_read_in}'. Skipping this pair."
        continue # Skip to the next iteration of the loop
    fi

    forward_paired_out="${trimmomatic_results_dir}/${base}_1.trim.fastq.gz"
    forward_unpaired_out="${trimmomatic_results_dir}/${base}_1un.trim.fastq.gz"
    reverse_paired_out="${trimmomatic_results_dir}/${base}_2.trim.fastq.gz"
    reverse_unpaired_out="${trimmomatic_results_dir}/${base}_2un.trim.fastq.gz"

    # Execute Trimmomatic command
    # Assuming 'trimmomatic' is in your PATH, or provide full path if not (e.g., java -jar /path/to/trimmomatic.jar PE ...)
    trimmomatic PE -phred33 \
                   "${forward_read_in}" "${reverse_read_in}" \
                   "${forward_paired_out}" "${forward_unpaired_out}" \
                   "${reverse_paired_out}" "${reverse_unpaired_out}" \
                   ILLUMINACLIP:"${TRIMMOMATIC_ADAPTERS}":2:40:15 \
                   SLIDINGWINDOW:4:20 \
                   MINLEN:25
    if [ $? -ne 0 ]; then
        echo "Error: Trimmomatic failed for sample ${base}. Check previous messages."
        # You might choose to exit here or continue to the next sample
    fi

done

echo -e "\n--- Trimmomatic Analysis Complete! ---"
echo "Trimmed reads are located in: '$trimmomatic_results_dir'"
