# Polygenic Risk Scores for Psychiatric Disorders in ABCD

This repository documents the pipeline used to compute polygenic risk scores (PRS) for several psychiatric disorders in the ABCD genetic dataset using PRS-CS.

The current workflow uses PRS-CS with the UK Biobank European LD reference panel and ABCD Release 3 imputed genotype data. GWAS summary statistics are harmonized to the PRS-CS format and then converted into individual-level PRS using PLINK.

---

## Traits and current status

The pipeline is set up to compute PRS for the following psychiatric disorders:

- Anorexia nervosa (AN)
- Obsessive-compulsive disorder (OCD)
- Attention-deficit/hyperactivity disorder (ADHD)
- Post-traumatic stress disorder (PTSD)
- Tourette syndrome (TS)
- Opioid dependence (OD) – currently problematic (see note below)
- Anxiety disorder (AD) – not yet successfully computed (see note below)

Successfully computed PRS in the current version:

- AN: `AN_prs_score.*`
- OCD: `OCD_prs_score.*`
- ADHD: `adhd2022_PRS.*`
- PTSD: `ptsd_prs_score.*`
- TS: `ts_prs_score.*`

OD and AD are still under troubleshooting and are not included as final usable PRS in the current pipeline.

---

## Input data

### Target data (ABCD)

- ABCD Release 3 imputed genotype data
- Binary PLINK files:
  - `abcd_r3_impute_rsq03_maf001_rsid.bed`
  - `abcd_r3_impute_rsq03_maf001_rsid.bim`
  - `abcd_r3_impute_rsq03_maf001_rsid.fam`

Example path (adjust to your system):

- `/data/neuromark2/Data/ABCD/Data_genetic/abcd_r3_impute_rsq03_maf001_rsid`

### LD reference

- PRS-CS UK Biobank European LD panels:
  - Directory: `ldblk_ukbb_eur`
  - Files: `ldblk_1.hdf5`, ..., `ldblk_22.hdf5` (depending on PRS-CS distribution)

Example path:

- `/data/users1/<user>/ldblk_ukbb_eur`

### GWAS summary statistics

Each trait uses publicly available GWAS summary statistics, converted to the PRS-CS input format. Two basic formats are used:

1. Odds ratio (OR) + standard error (SE):

   Required columns (PRS-CS convention):

   - `SNP` (rsID)
   - `A1` (effect allele)
   - `A2` (other allele)
   - `OR`
   - `SE`

   Example file names:
   - `pgcAN2.ORSE.clean.tsv`
   - `OCD.ORSE.clean.tsv`
   - `TS_Oct2018.SNP_A1_A2_OR_SE.txt`
   - `ADHD_meta_Jan2022_iPSYCH1_iPSYCH2_deCODE_PGC.meta` → converted to `adhd2022_SNP_A1_A2_OR_SE.txt`

2. Beta + P-value:

   Required columns:

   - `SNP`
   - `A1`
   - `A2`
   - `BETA`
   - `P`

   Example file name:
   - `ptsd.prscs.tsv`
   - `OD.prscs.tsv` (note: this one currently leads to empty PRS after harmonization; see limitation note below)

All GWAS files are stored under a working directory, e.g.:

- `/data/users1/<user>/`

---

## Software and dependencies

- Python 3 (via conda)
- PRS-CS (Python script `PRScs.py`)
- PLINK 1.9 or 2.0
- SLURM (for running jobs on an HPC cluster)
- Optional: `awk` or other command-line tools for preprocessing GWAS files

Example environment:

- Conda environment: `/data/users1/<user>/envs/prscs`
- PRS-CS script: `/data/users1/<user>/PRScs/PRScs.py`
- PLINK binaries: `/data/users1/<user>/bin/plink`

---

## Directory structure (example)

A typical layout on the cluster:

- `/data/users1/<user>/`
  - `PRScs/` – PRS-CS code (`PRScs.py`, etc.)
  - `ldblk_ukbb_eur/` – LD reference panels
  - `prs_output/` – folder to store PRS-CS and PLINK output
    - `AN_pst_eff_*`
    - `AN_allchr_pst_eff.txt`
    - `AN_prs_score.*`
    - `OCD_pst_eff_*`
    - `OCD_allchr_pst_eff.txt`
    - `OCD_prs_score.*`
    - `adhd2022_pst_eff_*`
    - `adhd2022_allchr.txt`
    - `adhd2022_PRS.*`
    - `ptsd_pst_eff_*`
    - `ptsd_allchr_pst_eff.txt`
    - `ptsd_prs_score.*`
    - `ts_pst_eff_*`
    - `ts_allchr_pst_eff.txt`
    - `ts_prs_score.*`
    - `od_pst_eff_*` (currently leads to empty PRS)
    - `od_allchr_pst_eff.txt` (header only / empty)
  - GWAS summary files:
    - `pgcAN2.ORSE.clean.tsv`
    - `OCD.ORSE.clean.tsv`
    - `ptsd.prscs.tsv`
    - `TS_Oct2018.SNP_A1_A2_OR_SE.txt`
    - `ADHD_meta_Jan2022_iPSYCH1_iPSYCH2_deCODE_PGC.meta` (and converted file)
    - `OD.prscs.tsv`

---

## General PRS-CS workflow

The workflow is the same for each trait:

1. Prepare GWAS summary statistics
2. Run PRS-CS per chromosome (1–22)
3. Concatenate per-chromosome posterior effect sizes
4. Use PLINK `--score` to compute individual-level PRS

### 1. Prepare GWAS summary statistics

Example: converting a multi-column GWAS file (ADHD) into PRS-CS format:

```bash
awk 'NR==1{print "SNP\tA1\tA2\tOR\tSE"; next} {print $2"\t"$4"\t"$5"\t"$10"\t"$11}' \
  ADHD_meta_Jan2022_iPSYCH1_iPSYCH2_deCODE_PGC.meta \
  > adhd2022_SNP_A1_A2_OR_SE.txt
