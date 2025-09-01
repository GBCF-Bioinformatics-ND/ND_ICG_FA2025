# ND_ICG_FA2025
Course content for Introduction to Computational Genomics (BIOS 40132-01 and BIOS 60132-01) at the University of Notre Dame, Notre Dame, Indiana, USA.


```markdown
# Quick‑Start Guide – Module 01: Genomics QC  
*(University HPC – “CRC” – `/scratch365/$USER`)*

> **Tip** – Replace `$USER` with your actual HPC username.
```
---

## 1. Log In & Prepare Working Directory  

To log in, you can use Terminal (Mac), PowerShell/Command Prompt (Win), or 3rd party SSH clients like Termius or MobaXTerm.
To use CRC, you must be on campus and connected to eduroam. If you are off campus, consider using bastion (see below) or ND VPN.

```bash
# 1. Log in to the CRC front end
# replace netid with your real netid, e.g. jsmith
ssh netid@crcfe01.crc.nd.edu
# or
ssh netid@crcfe02.crc.nd.edu
# if you're off campus
ssh netid@bastion.crc.nd.edu

# When prompted, type in your ND password and hit enter when finished.
# You won't be able to see what you typed.

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
2. Log in with your Notre Dame credentials.  
3. Download the following items to your local machine.

  - `input` folder with `fastq.gz` files and adapter seqenece file
  - `references` folder with `.fasta` file

> **Note**: The scripts assume a naming convention:  
> *read 1* → `*_R1_*.fastq.gz` (or `*_1.fastq.gz`)  
> *read 2* → `*_R2_*.fastq.gz` (or `*_2.fastq.gz`)  
> If your files use a different pattern, adjust the scripts or rename the files accordingly.

---

## 4. Transfer the Files to HPC with Globus  

### **Set up a local endpoint**

1. Go to [https://www.globus.org/get-started](https://www.globus.org/get-started) and select “log in”  
2. Search for “University of Notre Dame”, and authenticate using ND credentials.
3. Go to “Collections” \> click on “Get globus personal” on the upper right and download to your local machine. Alternatively, use [this link](https://www.globus.org/globus-connect-personal) to download.
4. Open the application on your computer. A Set-up window will pop up for first-time users. Follow the on-screen instructions to log in and allow access to files.  
5. Then, enter a collection name and description of your choice, and click “Save” \> “Exit Setup”.  
6. To enable the local endpoint, keep the globus connect app running on your local machine.

### **Transfer files to CRC**

1. Now log in to the web app ([https://app.globus.org/](https://app.globus.org/)), and on the “File Transfer” tab, select the your local collection on one side. Navigate or enter the path to the file or folder you want to transfer. 
2. On the other side, search for and select “ND Center for Research Computing Collection” (Remember to check “Search All Collections”) .
3. For CRC, you can only transfer files to/from scratch space. On the CRC side, enter /scratch365/netid in the “Path” field. Replace netid with your real netid.   
   Note: First-time users will be prompted to consent to data access for CRC. Click “Continue” and authenticate using your ND email (1st option). Do **NOT** click “Link a new identity” (2nd option).
4. Click “Start” to schedule a transfer. You can track the status under the “Activity” tab. When the transfer is complete, you will also receive an email notification.
5. Verify on the HPC:

```bash
cd /scratch365/$USER/ND_ICG_FA2025
ls -R input/
ls -R references/
```

You should see all FASTQ files and the reference genome file inside those directories. Additional documentation and troubleshooting steps can be found at [CRC Documetation](https://docs.crc.nd.edu/resources/globus.html) and [Globus Tutorial](https://docs.globus.org/guides/tutorials/manage-files/transfer-files/).

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

Please also change the email in the script headers, for example:

```bash
nano read_fastqc.sh
# change the line to your real netid
#$ -M netid@nd.edu
# save and exit
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

### Feel Free to contact me or TAs if you have any questions.

- Instructor: Dr. Bharat Mishra bmishra2@nd.edu; 3028C McCourtney Hall East.
- Co-Instructor and TA 1: Elizabeth Brooks ebrooks5@nd.edu
- TA 2: Ziyu Zeng zzeng4@nd.edu
- TA 3: Nirjhar Bhattacharyya nbhattac@nd.edu
- Undergraduate TA: Colleen Farrell cfarre23@nd.edu

#### *Happy analyzing!* 
