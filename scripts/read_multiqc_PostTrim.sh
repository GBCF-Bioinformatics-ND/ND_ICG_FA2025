#!/bin/bash
#
#$ -M zzeng4@nd.edu
#$ -m abe
#$ -q long
#$ -N multiqc

# load the bio module from CRC or intall multiqc

module load bio/0724


# --- Configuration ---
input_dir="../results" # MultiQC often searches results directories
results_base_dir="../results" # Output base for MultiQC report

# MultiQC specific configurations
multiqc_output_subdir="multiqc_report_PostTrim"       # Subdirectory for MultiQC report
multiqc_report_name="multiqc_report_PostTrim.html"    # Name of the generated HTML report

# Define full path for MultiQC specific results directory
multiqc_results_dir="${results_base_dir}/${multiqc_output_subdir}"

# --- MultiQC Section ---
# 1. Prepare MultiQC output directory
mkdir -p "$multiqc_results_dir" || { echo "Error: Could not create MultiQC results directory '$multiqc_results_dir'. Exiting."; exit 1; }

# 2. Run MultiQC
# MultiQC will search the provided directories for compatible log files
multiqc "$input_dir" -o "$multiqc_results_dir" -n "$multiqc_report_name" --force

if [ $? -ne 0 ]; then
    echo "Error: MultiQC failed to generate report."
    echo "Ensure MultiQC is installed and accessible in your PATH."
    echo "Also check if there are any compatible log files in '$input_dir'."
    exit 1
fi

echo -e "\n--- Script Execution Complete! ---"
