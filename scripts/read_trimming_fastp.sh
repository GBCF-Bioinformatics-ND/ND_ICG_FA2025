#!/bin/bash
#
#$ -M bmishra2@nd.edu      # Email address for notifications (replace with your actual email)
#$ -m abe                  # Send mail on begin, abort, and end
#$ -q debug                # Job queue (e.g., debug, long, normal)
#$ -N fastp_trim           # Job name for easy identification
#$ -pe threads 4           # Request 4 cores/threads for parallel tasks (adjust as needed)
#$ -V                      # Export all environment variables to the job

# Load the bio/2.0 module (assuming fastp is available here)
module load bio/2.0

# Configuration (adjust paths as needed)
input_dir="../input"
results_base_dir="../results"

TRIMMED_READS_SUBDIR="fastp_trimmed"     # Subdirectory for trimmed reads (matching your example)
FASTP_REPORTS_SUBDIR="fastp_reports"     # Subdirectory for fastp HTML/JSON reports
FASTP_LOGS_SUBDIR="fastp_logs"           # Subdirectory for detailed fastp command logs

# Define full paths for all output directories
TRIMMED_READS_DIR="${results_base_dir}/${TRIMMED_READS_SUBDIR}"
FASTP_REPORTS_DIR="${results_base_dir}/${FASTP_REPORTS_SUBDIR}"
FASTP_LOGS_DIR="${results_base_dir}/${FASTP_LOGS_SUBDIR}"

# Main script log file (all echoes from this script will go here)
MAIN_LOG_FILE="${FASTP_LOGS_DIR}/pipeline_fastp_run_$(date +%Y%m%d_%H%M%S).log"

mkdir -p "$TRIMMED_READS_DIR" "$FASTP_REPORTS_DIR" "$FASTP_LOGS_DIR" # Create output directory if it doesn't exist

# Loop through all forward raw FASTQ files
for infile in ${input_dir}/*_1.fastq.gz; do
  # Extract base name (e.g., 'Sample1' from 'data/Sample1_1.fastq.gz')
  base=$(basename "${infile}" _1.fastq.gz)

  # Define input and output file paths for the current sample
  fq1_in="${infile}"
  fq2_in="${input_dir}/${base}_2.fastq.gz"
  fq1_out="${TRIMMED_READS_DIR}/${base}_1.trim.fastq.gz"
  fq2_out="${TRIMMED_READS_DIR}/${base}_2.trim.fastq.gz"
  html_report="${FASTP_REPORTS_DIR}/${base}.fastp.html"
  json_report="${FASTP_REPORTS_DIR}/${base}.fastp.json"
  FASTP_COMMAND_LOG="${FASTP_LOGS_DIR}/${base}_fastp.log"
  
  echo "Processing sample: ${base}"

  # Run fastp command
  fastp \
    -i "${fq1_in}" \
    -I "${fq2_in}" \
    -o "${fq1_out}" \
    -O "${fq2_out}" \
    --html "${html_report}" \
    --json "${json_report}" \
    --thread 4 \
    --detect_adapter_for_pe \
    --overrepresentation_analysis \
    --correction \
    --cut_right \
    &> "$FASTP_COMMAND_LOG" # Redirect stdout and stderr to the log file
  if [ $? -ne 0 ]; then
        echo "ERROR: fastp failed for sample ${base}. Check '$FASTP_COMMAND_LOG' for details. Exiting."
        exit 1
    fi
    echo "  fastp trimming complete for ${base}."

done # End of per-sample loop