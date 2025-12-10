Polygenic Risk Scores for Psychiatric Disorders in ABCD
This repository documents the pipeline used to compute polygenic risk scores (PRS) for several psychiatric disorders in the ABCD genetic dataset using PRS-CS.
The current workflow uses PRS-CS with the UK Biobank European LD reference panel and ABCD Release 3 imputed genotype data. GWAS summary statistics are harmonized to the PRS-CS format and then converted into PRS at the individual level using PLINK.

Data
Target data
ABCD Release 3 imputed genotype data
Binary PLINK files:
abcd_r3_impute_rsq03_maf001_rsid.bed
abcd_r3_impute_rsq03_maf001_rsid.bim
abcd_r3_impute_rsq03_maf001_rsid.fam
LD reference panel
PRS-CS LD reference:
ldblk_ukbb_eur/ (UK Biobank European ancestry HapMap3 SNPs)
GWAS summary statistics (PRS-CS formatted or converted)
Anorexia nervosa: pgcAN2.ORSE.clean.tsv
Obsessive-compulsive disorder: OCD.ORSE.clean.tsv
ADHD: adhd2022_SNP_A1_A2_OR_SE.txt
PTSD: ptsd.prscs.tsv
Tourette syndrome: TS_Oct2018.SNP_A1_A2_OR_SE.txt
Opioid dependence (OD): OD.prscs.tsv
Anxiety disorder (AD): source file not yet finalized (see notes below)
Methods
Prepare GWAS summary statistics
Columns are harmonized to the PRS-CS expected format (for example, SNP, A1, A2, effect size, standard error or odds ratio/SE depending on the GWAS).
SNP IDs are rsIDs and aligned to the ABCD genotype data and the UKBB LD reference where possible.
Run PRS-CS per chromosome
Script template (example):
python PRScs.py \
  --ref_dir=/data/users1/csun10/ldblk_ukbb_eur \
  --bim_prefix=/data/neuromark2/Data/ABCD/Data_genetic/abcd_r3_impute_rsq03_maf001_rsid \
  --sst_file=<GWAS_SUMMARY_TSV> \
  --n_gwas=<N_eff> \
  --out_dir=<OUTPREFIX> \
  --phi=1e-2 \
  --chrom=<1-22>
PRS-CS performs Bayesian shrinkage using the LD structure from the UKBB reference and outputs posterior SNP effect sizes for each chromosome.
Concatenate per-chromosome outputs
All 22 chromosomes are merged into a single file containing SNP-level weights (posterior betas):
cat <OUTPREFIX>_pst_eff_a1_b0.5_phi1e-02_chr*.txt \
  > <TRAIT>_allchr_pst_eff.txt
Compute individual-level PRS using PLINK
PLINK --score is used to project SNP weights onto ABCD genotypes:
plink \
  --bfile /data/neuromark2/Data/ABCD/Data_genetic/abcd_r3_impute_rsq03_maf001_rsid \
  --score <TRAIT>_allchr_pst_eff.txt 2 4 6 header sum \
  --out <TRAIT>_prs_score
Column indices correspond to:
2: SNP ID
4: effect allele (A1)
6: posterior effect size (BETA_post)
Traits and PRS availability
Currently:
Anorexia nervosa (AN)
Status: PRS successfully computed
Example output files:
AN_allchr_pst_eff.txt
AN_prs_score.profile
AN_prs_score.log
Obsessive-compulsive disorder (OCD)
Status: PRS successfully computed
Example output: OCD_prs_score.*
ADHD
Status: PRS successfully computed
Example output: adhd2022_PRS.profile (naming may vary slightly)
PTSD
Status: PRS successfully computed
Example output: ptsd_prs_score.*
Tourette syndrome (TS)
Status: PRS successfully computed
Example output: ts_prs_score.*
Anxiety disorder (AD)
Status: not yet finalized
Reason: the initial GWAS summary file format for AD requires additional cleaning and harmonization before it can be used in PRS-CS. This PRS is planned but not yet included in the current set.
Opioid dependence (OD)
Status: PRS could not be computed (see details below)
Example attempted output prefix: od_prs_score.*
PLINK error:
Error: Empty --score file.
Why OD PRS could not be computed
For PRS-CS to produce usable SNP weights, each SNP must be present in three places simultaneously:
The LD reference panel (ldblk_ukbb_eur),
The GWAS summary statistics for the trait (for example, OD.prscs.tsv),
The target genotype data (ABCD; abcd_r3_impute_rsq03_maf001_rsid.bim).
PRS-CS internally intersects these three SNP sets. Only SNPs that are present in all three datasets are retained and assigned posterior effect sizes.
For AN, OCD, ADHD, PTSD, and TS:

The underlying GWAS are large, standard consortium GWAS that rely mostly on common HapMap3 SNPs.
These SNP sets overlap well with both the UKBB LD reference and the ABCD genotype data.
As a result, the per-chromosome PRS-CS outputs contain many SNPs, and the merged <TRAIT>_allchr_pst_eff.txt files have a large number of rows.
PLINK can then successfully compute PRS.
For OD:
The OD GWAS (OD.prscs.tsv) appears to use a SNP set that has very poor overlap with:
The HapMap3 SNPs in the UKBB LD reference, and/or
The SNPs present in the ABCD imputed genotype data.
When PRS-CS intersects the three SNP sets (LD reference, OD GWAS, ABCD genotypes), almost no SNPs remain in common.
PRS-CS still runs to completion, but each per-chromosome output file contains only a header line or effectively zero SNP rows.
After concatenation, od_allchr_pst_eff.txt is essentially an empty weight file (no SNPs with valid posterior betas).
PLINK then fails with:
Error: Empty --score file.
This behavior indicates that the pipeline itself is working as expected, but the OD GWAS SNP set is not compatible with the current LD reference plus target genotype combination. In other words, there are not enough overlapping SNPs across all three datasets to construct a meaningful PRS for OD with this setup.
Known limitations and potential next steps
Anxiety disorder (AD)
The AD PRS is not yet included due to remaining format and harmonization issues with the AD GWAS summary statistics. Once a clean, PRS-CS compatible file is prepared, the same pipeline can be applied.
Opioid dependence (OD)
Potential strategies to obtain an OD PRS include:
Using a different LD reference that better matches the OD GWAS SNP set.
Re-mapping or liftover of SNP IDs if there are build inconsistencies.
Filtering or reformatting the OD GWAS to a SNP set that overlaps more strongly with HapMap3 and the ABCD genotype data.
At the moment, this repository reports OD as “PRS not available due to insufficient SNP overlap between OD GWAS, LD reference panel, and ABCD genotypes.”
Reproducibility
All PRS-CS and PLINK commands are executed on an HPC cluster with SLURM. Example SLURM script structure:
#!/bin/bash -l
#SBATCH -J prscs_job
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 64:00:00
#SBATCH -p qTRD
#SBATCH -A trends53c17
#SBATCH -o /data/users1/csun10/prs_output/prscs_%j.out
#SBATCH -e /data/users1/csun10/prs_output/prscs_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=csun10@student.gsu.edu

set -euxo pipefail
export MKL_NUM_THREADS=8 OMP_NUM_THREADS=8 NUMEXPR_NUM_THREADS=8

OUTDIR=/data/users1/csun10/prs_output
GWAS_SST=<GWAS_SUMMARY_TSV>
REFDIR=/data/users1/csun10/ldblk_ukbb_eur
BIM_PREFIX=/data/neuromark2/Data/ABCD/Data_genetic/abcd_r3_impute_rsq03_maf001_rsid
PRSCSPY=/data/users1/csun10/PRScs/PRScs.py
PLINK_BIN_DIR=/data/users1/csun10/bin

mkdir -p "$OUTDIR"
export PATH="$PLINK_BIN_DIR:$PATH"

module load miniconda3/4.12.0
eval "$(conda shell.bash hook)"
conda activate /data/users1/csun10/envs/prscs
export CONDA_PKGS_DIRS=/data/users1/csun10/conda_pkgs

OUTPREFIX="$OUTDIR/<TRAIT>_pst_eff"
N_GWAS=<N_eff>

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

cat ${OUTPREFIX}_pst_eff_a1_b0.5_phi1e-02_chr*.txt > "$OUTDIR/<TRAIT>_allchr_pst_eff.txt"

plink \
  --bfile "$BIM_PREFIX" \
  --score "$OUTDIR/<TRAIT>_allchr_pst_eff.txt" 2 4 6 header sum \
  --out "$OUTDIR/<TRAIT>_prs_score"
Paths and filenames in this README are examples based on the current working setup and can be adapted as needed for other environments.



Thinking



ChatGPT can make mistakes. Check import

