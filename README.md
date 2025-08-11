# ND_ICG_FA2025
Course content for Introduction to Computational Genomics (BIOS 40132-01 and BIOS 60132-01) at the University of Notre Dame, Notre Dame, Indiana, USA.


```markdown
# Quick‑Start Guide – Module 01: Genomics QC  
*(University HPC – “CRC” – `/scratch365/$USER`)*

> **Tip** – Replace `$USER` with your actual HPC username.

---

## 1. Log In & Prepare Working Directory  

```bash
# 1. Log in to the CRC HPC (replace <hpc-host> with the real host name)
ssh <hpc-host>

# 2. Go to your scratch space
cd /scratch365/$USER
```

---

## 2. Clone the GitHub Repository  

```bash
git clone https://github.com/GBCF-Bioinformatics-ND/ND_ICG_FA2025.git
cd ND_ICG_FA2025
```

The repo contains all helper scripts (`read_fastqc.sh`, `read_trimming_fastp.sh`, …) and this very README.

---

## 3. Download Input & Reference Data  

The data are hosted on Box.  
1. Open the link in a browser:  
   <https://notredame.app.box.com/folder/333412290959>  
2. Log in with your Not‑ReDame credentials.  
3. Download the following items to your local machine (or to the HPC via `curl`/`wget` if you prefer).  

| Item | Description |
|------|-------------|
| `paired_end_fastq.zip` | Raw paired‑end FASTQ files (e.g., `sample_R1.fastq.gz`, `sample_R2.fastq.gz`) |
| `reference_genome.zip` | Reference genome (FASTA, index files, etc.) |

> **Note**: The scripts assume a naming convention:  
> *read 1* → `*_R1_*.fastq.gz` (or `*_1.fastq.gz`)  
> *read 2* → `*_R2_*.fastq.gz` (or `*_2.fastq.gz`)  
> If your files use a different pattern, adjust the scripts or rename the files accordingly.

---

## 4. Transfer the Files to HPC with Globus  

1. **Create a Globus transfer** from the Box folder to `/scratch365/$USER/ND_ICG_FA2025`.  
2. In the Globus UI, point the *source* to the folder you just downloaded and *destination* to:  

```
/scratch365/$USER/ND_ICG_FA2025/
```

3. Start the transfer and wait until it completes.  
4. Verify on the HPC:

```bash
cd /scratch365/$USER/ND_ICG_FA2025
ls -R input/
ls -R reference/
```

You should see all FASTQ files and the reference genome file inside those directories.

---

## 5. Run the Analysis Pipeline – *Module 01: Genomics QC*  

All scripts are located in the `scripts` folder of the repo.  
Before running, make sure the required tools are available.  
On CRC you’ll usually load modules; check the list below or ask staff.

```bash
# Example module loads – adapt to your system
module load bio/2.0
module load bio/0724
```

### 5.1. FASTQ Quality Control (Pre‑Trim)

```bash
cd /scratch365/$USER/ND_ICG_FA2025/scripts
bash read_fastqc.sh
```

*Script output*: `fastqc/` directory containing FastQC reports for each raw FASTQ file.

### 5.2. Adapter/Quality Trimming

Choose one of the two trimming scripts:

| Script                | Tool |
|-----------------------|------|
| `read_Trimming.sh`    | Trimmomatic |
| `read_trimming_fastp.sh` | fastp |

#### Using fastp (recommended)

```bash
qsub read_trimming_fastp.sh
```

*Output*: Trimmed FASTQ files (e.g., `sample_R1_trimmed.fastq.gz`).

### 5.3. Post‑Trim FASTQ QC & MultiQC Summary

```bash
qsub read_fastqc_postTrim.sh
```

This script runs FastQC on the trimmed files.  
After it finishes, run MultiQC to aggregate all reports:

```bash
qsub read_multiqc_PostTrim.sh
```

> **Remember**: `multiqc` will only run if the previous `fastqc` step succeeded; otherwise it will abort.

### 5.4. Alignment & Variant Calling

```bash
qsub read_variants.sh
```

*What it does*  
1. Aligns trimmed reads to the reference genome using BWA.  
2. Sorts with Samtools.  
3. Calls variants with BCFtools.  
4. Generates VCF and filter variants.

All intermediate files (BAMs, VCFs, logs) are stored in the `results/` directory.

---

## 6. Inspect the Results

```bash
# FastQC reports
open fastqc/*.html

# MultiQC summary
open multiqc_report.html

# BAM/VCF files
ls results/
```

You can use `samtools view` or `bcftools view` to explore the alignment and variant data.

---

## 7. Troubleshooting & Tips

| Symptom | Likely Cause | Fix |
|---------|--------------|-----|
| `module not found` | Wrong module name | Use `module avail` to list available tools |
| FastQC fails on gzipped files | Corrupt download | Re‑download from Box |
| Trimming script errors | File paths wrong | Verify that `data/` contains the correct FASTQ names |
| `multiqc` aborts | No FastQC output | Re‑run `read_fastqc_postTrim.sh` |
| Alignment takes too long | Large genome | Use `bwa-mem2` or a faster index |

For any other questions, reach out to the course TA or post on the class Slack channel.

---

## 8. Credits

- **Module 1: Genomics QC** 
  - **Repo**: <https://datacarpentry.github.io/wrangling-genomics/index.html>  
  - **Data**: <https://www.pnas.org/doi/10.1073/pnas.0803151105> Blount et al. 2008: Historical contingency and the evolution of a key innovation in an experimental population of Escherichia coli.  

### Fell Free to contact me or TAs if you have any questions.

- Instructor: Dr. Bharat Mishra bmishra2@nd.edu; 3028C McCourtney Hall East.
- Co-Instructor and TA 1: Elizabeth Brooks ebrooks5@nd.edu
- TA 2: Ziyu Zeng zzeng4@nd.edu
- TA 3: Nirjhar Bhattacharyya nbhattac@nd.edu
- Undergraduate TA: Colleen Farrell cfarre23@nd.edu

#### *Happy analyzing!* 