#!/bin/bash

#$ -M <netid>@nd.edu
#$ -m abe
#$ -q long
#$ -N read_subsample

# load the bio/2.0 module from CRC or intall conda module

conda activate /users/bmishra2/.conda/envs/seqkit

#module load bio/2.0

# --- Configuration ---
# Set default input and results directories.

PERCENTAGE=0.7
SEED=100
INPUT_DIR="../raw_full" # Change this to the directory containing your FASTQ files
OUTPUT_DIR="../input" # Output directory for subsampled files

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Find all fastq.gz files in the input directory
find "$INPUT_DIR" -maxdepth 1 -name "*.fastq.gz" -print0 | while IFS= read -r -d $'\0' input_file; do
    if [[ -f "$input_file" ]]; then
        # Extract the base filename (e.g., "sampleA_R1.fastq.gz" -> "sampleA_R1")
        filename=$(basename "$input_file" .fastq.gz)

        # Construct the output filename
        output_file="${OUTPUT_DIR}/sub_${filename}.fastq.gz"

        echo "Processing: $input_file"
        echo "Output to: $output_file"

        # Run seqkit sample
        seqkit sample -p "$PERCENTAGE" -s "$SEED" "$input_file" -o "$output_file"

        if [[ $? -eq 0 ]]; then
            echo "Successfully subsampled $input_file."
        else
            echo "Error subsampling $input_file. Please check the seqkit command and input file."
        fi
        echo "-------------------------------------"
    else
        echo "Warning: Skipped non-existent file path: $input_file"
    fi
done

echo "FASTQ subsampling process completed."
