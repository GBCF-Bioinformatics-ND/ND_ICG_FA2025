#!/bin/bash
#
#$ -M <netid>@nd.edu      # Email address for notifications
#$ -m abe                  # Send mail on begin, abort, and end
#$ -q long                # Job queue (e.g., long)
#$ -N variant_auto         # Job name for easy identification
#$ -pe smp 4           # Request 4 cores/threads for parallel tasks (adjust as needed)
#$ -V                      # Export all environment variables to the job

# Load the required bioinformatics module
# This command is specific to environments like CRC
module load bio/2.0

# --- Configuration ---
# Genome file name (must be unzipped .fasta in reference_dir)
cd ../results
genome=../references/ecoli_rel606.fasta

bwa index $genome

mkdir -p sam bam bcf vcf

for fq1 in ../results/fastp_trimmed/*_1.trim.fastq.gz
    do
    echo "working with file $fq1"

    base=$(basename $fq1 _1.trim.fastq.gz)
    echo "base name is $base"

    fq1=../results/fastp_trimmed/${base}_1.trim.fastq.gz
    fq2=../results/fastp_trimmed/${base}_2.trim.fastq.gz
    sam=../results/sam/${base}.aligned.sam
    bam=../results/bam/${base}.aligned.bam
    sorted_bam=../results/bam/${base}.aligned.sorted.bam
    raw_bcf=../results/bcf/${base}_raw.bcf
    variants=../results/vcf/${base}_variants.vcf
    final_variants=../results/vcf/${base}_final_variants.vcf

    bwa mem $genome $fq1 $fq2 > $sam
    samtools view -S -b $sam > $bam
    samtools sort -o $sorted_bam $bam
    samtools index $sorted_bam
    bcftools mpileup -O b -o $raw_bcf -f $genome $sorted_bam
    bcftools call --ploidy 1 -m -v -o $variants $raw_bcf
    vcfutils.pl varFilter $variants > $final_variants

    done
