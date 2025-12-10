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

- `AN_prs_score.*`  
- `OCD_prs_score.*`  
- `adhd2022_PRS.*`  
- `ptsd_prs_score.*`  
- `ts_prs_score.*`  

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

   Example file names:  
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
- PLINK binary: `/data/users1/<user>/bin/plink`

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

~~~bash
awk 'NR==1{print "SNP\tA1\tA2\tOR\tSE"; next} {print $2"\t"$4"\t"$5"\t"$10"\t"$11}' \
  ADHD_meta_Jan2022_iPSYCH1_iPSYCH2_deCODE_PGC.meta \
  > adhd2022_SNP_A1_A2_OR_SE.txt
~~~

For traits already in the required format (e.g., `pgcAN2.ORSE.clean.tsv`, `OCD.ORSE.clean.tsv`, `ptsd.prscs.tsv`, `OD.prscs.tsv`), only minor cleaning or renaming is needed.

### 2. Run PRS-CS per chromosome

A generic SLURM script skeleton for running PRS-CS:

~~~bash
#!/bin/bash -l
#SBATCH -J prscs_<TRAIT>
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 64:00:00
#SBATCH -p qTRD
#SBATCH -A <ACCOUNT>
#SBATCH -o /data/users1/<user>/prs_output/prscs_%j.out
#SBATCH -e /data/users1/<user>/prs_output/prscs_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=<USER_EMAIL>

set -euxo pipefail
export MKL_NUM_THREADS=8 OMP_NUM_THREADS=8 NUMEXPR_NUM_THREADS=8

# Paths
OUTDIR=/data/users1/<user>/prs_output
REFDIR=/data/users1/<user>/ldblk_ukbb_eur
BIM_PREFIX=/data/neuromark2/Data/ABCD/Data_genetic/abcd_r3_impute_rsq03_maf001_rsid
PRSCSPY=/data/users1/<user>/PRScs/PRScs.py
PLINK_BIN_DIR=/data/users1/<user>/bin

mkdir -p "$OUTDIR"
export PATH="$PLINK_BIN_DIR:$PATH"

module load miniconda3/4.12.0
eval "$(conda shell.bash hook)"
conda activate /data/users1/<user>/envs/prscs
export CONDA_PKGS_DIRS=/data/users1/<user>/conda_pkgs

# Trait-specific settings
TRAIT=<TRAIT>                          # e.g. AN, OCD, adhd2022, ptsd, ts
GWAS_SST=/data/users1/<user>/<SST>.tsv # GWAS sumstats in PRS-CS format
N_GWAS=<N_GWAS>                        # effective sample size
OUTPREFIX="$OUTDIR/${TRAIT}_pst_eff"

# Run PRS-CS for chromosomes 1–22
for CHR in {1..22}; do
  python "$PRSCSPY" \
    --ref_dir="$REFDIR" \
    --bim_prefix="$BIM_PREFIX" \
    --sst_file="$GWAS_SST" \
    --n_gwas="$N_GWAS" \
    --out_dir="$OUTPREFIX" \
    --phi=1e-2 \
    --chrom=${CHR}
done
~~~

Here `<TRAIT>`, `<SST>`, `<N_GWAS>`, `<ACCOUNT>`, `<USER_EMAIL>`, and `<user>` should be replaced with trait-specific or user-specific values.

### 3. Concatenate effect sizes and compute PRS

After PRS-CS finishes, per-chromosome posterior effect sizes are combined and fed to PLINK.

~~~bash
# Concatenate PRS-CS weight files
cat ${OUTPREFIX}_pst_eff_a1_b0.5_phi1e-02_chr*.txt > "$OUTDIR/${TRAIT}_allchr_pst_eff.txt"

# Compute PRS with PLINK
plink \
  --bfile "$BIM_PREFIX" \
  --score "$OUTDIR/${TRAIT}_allchr_pst_eff.txt" 2 4 6 header sum \
  --out "$OUTDIR/${TRAIT}_prs_score"
~~~

The column indices `2 4 6` correspond to `SNP`, `A1`, and `BETA_post` in the PRS-CS output (with `header` indicating that the first row contains column names).

---

## Outputs

For each successfully computed trait, the pipeline produces:

- Per-chromosome PRS-CS weight files:  
  - `<TRAIT>_pst_eff_a1_b0.5_phi1e-02_chr<CHR>.txt`  
- Merged weight file:  
  - `<TRAIT>_allchr_pst_eff.txt`  
- PLINK PRS files:  
  - `<TRAIT>_prs_score.profile`  
  - `<TRAIT>_prs_score.log`  

Example final PRS prefixes:

- `AN_prs_score`  
- `OCD_prs_score`  
- `adhd2022_PRS`  
- `ptsd_prs_score`  
- `ts_prs_score`  

For OD, the current pipeline generates `od_allchr_pst_eff.txt` with a header only (no SNPs), which leads to an empty PRS in the PLINK `--score` step.

---

## Current limitations and troubleshooting notes

- OD (`OD.prscs.tsv`):  
  - After harmonization and matching to the ABCD genotype set, there is effectively no overlap between GWAS SNPs and target SNPs passing QC.  
  - As a result, PRS-CS produces no valid SNP weights, and the merged file `od_allchr_pst_eff.txt` is empty (header only).  
  - PLINK then reports `Error: Empty --score file.` when attempting to compute `od_prs_score`.  

- AD (anxiety disorder):  
  - At the time of this documentation, AD PRS has not yet been successfully computed.  
  - Possible causes include issues with the source GWAS file format or missing required columns for PRS-CS.  

These traits are therefore not included as usable PRS in the current version of the pipeline.

---

## Reproducibility and usage

To reproduce PRS for a given trait:

1. Place or generate the GWAS summary statistics in PRS-CS format (`SNP`, `A1`, `A2`, and either `OR`+`SE` or `BETA`+`P`).  
2. Update the paths, `TRAIT`, `GWAS_SST`, and `N_GWAS` in the SLURM script.  
3. Submit the job to the cluster via `sbatch`.  
4. After completion, use the `.profile` file from the `--out` prefix (for example, `AN_prs_score.profile`) as the trait’s PRS in downstream analyses.

---

## References

- PRS-CS: Ge T, Chen CY, Ni Y, Feng YA, Smoller JW (2019). Polygenic prediction via Bayesian regression and continuous shrinkage priors.  
- ABCD Study: For details on the ABCD genetic data, see the official ABCD documentation.  
- GWAS sources: Trait-specific GWAS summary statistics are taken from large consortia (e.g., PGC and related studies) as indicated by the file names above.
