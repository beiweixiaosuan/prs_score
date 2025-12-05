RS pipeline README (PGC GWAS: AD, OCD, ADHD, PTSD, TS) — PRS-CS + PLINK
Overview
This project computes polygenic risk scores (PRS) for five psychiatric traits using PGC GWAS summary statistics:
AD
OCD
ADHD
PTSD
TS
We use PRS-CS to infer posterior SNP effect sizes under an external LD reference panel, then compute individual-level PRS in the target cohort using PLINK --score.
This README documents: inputs, software environment, processing steps, outputs, and common QC/interpretation notes.
Methods summary
Start from cleaned PGC GWAS summary statistics (*.tsv) for each trait.
Run PRS-CS per chromosome (1–22) using:
EUR LD reference blocks (UKBB EUR)
Target cohort BIM for variant matching
GWAS effective sample size (median Neff)
A fixed global shrinkage parameter (phi=1e-2)
Concatenate per-chromosome PRS-CS outputs into a single weight file.
Score the target cohort genotypes using PLINK --score to produce a .profile file containing raw PRS (SCORESUM) per individual.
Software and compute environment
Scheduler: Slurm
Node request (example):
1 node, 1 task, 4 CPU cores, 32GB RAM, 64h walltime
Conda environment:
Miniconda module: miniconda3/4.12.0
Conda env: /data/users1/csun10/envs/prscs
PRS-CS:
Script: /data/users1/csun10/PRScs/PRScs.py
PLINK:
Binary directory added to PATH: /data/users1/csun10/bin
Note on threads:
The Slurm script requests -c 4, but the environment variables in the example set MKL/OMP/NUMEXPR threads to 8. For cluster compliance and efficiency, set these thread variables to match the Slurm CPU request (e.g., 4). This does not change results, but avoids oversubscribing CPUs.
Inputs (paths as used in the example script)
1) Target cohort genotypes (PLINK bed/bim/fam)
Prefix:
/data/neuromark2/Data/ABCD/Data_genetic/abcd_r3_impute_rsq03_maf001_rsid
Required files:
${BIM_PREFIX}.bed
${BIM_PREFIX}.bim
${BIM_PREFIX}.fam
2) LD reference panel (EUR)
Directory:
/data/users1/csun10/ldblk_ukbb_eur
This contains per-chromosome LD blocks in the format expected by PRS-CS.
3) GWAS summary statistics (PGC)
Example (trait-specific):
/data/users1/csun10/pgcAN2.ORSE.clean.tsv
For production runs, each trait should have its own cleaned summary statistics file, e.g.:
pgcAD.clean.tsv
pgcOCD.clean.tsv
pgcADHD.clean.tsv
pgcPTSD.clean.tsv
pgcTS.clean.tsv
Note: In the provided workflow, PRS-CS was run on uncompressed .tsv files (not .gz) due to tool/version compatibility.
Outputs (per trait)
Assuming a trait name <TRAIT> and output directory <OUTDIR>:
1) PRS-CS per-chromosome weight files
Produced by PRS-CS using:
--out_dir <OUTPREFIX>
Expected naming pattern:
<OUTPREFIX>*_chr1.txt, ..., <OUTPREFIX>*_chr22.txt
Each output line typically includes (no header in some cases):
CHR SNP(rsID) BP A1 A2 BETA_post
Example line:
1 rs4475691 846808 T C -7.138128e-04
2) Concatenated all-chromosome weight file
<OUTDIR>/<TRAIT>_allchr_pst_eff.txt
3) PLINK scoring output
<OUTDIR>/<TRAIT>_prs_score.profile (main PRS output)
<OUTDIR>/<TRAIT>_prs_score.log (QC / run log)
Other PLINK side outputs (e.g., .nosex) may also be produced.
Step-by-step pipeline (example commands)
This reflects the structure of the Slurm script used for one trait.
Step A. Run PRS-CS (chromosomes 1–22)
Variables (example):
OUTDIR=/data/users1/csun10/prs_output
GWAS_SST=/data/users1/csun10/pgc<TRAIT>.clean.tsv
REFDIR=/data/users1/csun10/ldblk_ukbb_eur
BIM_PREFIX=/data/neuromark2/Data/ABCD/Data_genetic/abcd_r3_impute_rsq03_maf001_rsid
PRSCSPY=/data/users1/csun10/PRScs/PRScs.py
N_GWAS=<median_Neff_for_trait>
OUTPREFIX=$OUTDIR/<TRAIT>_pst_eff
Run:
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
Step B. Concatenate per-chromosome outputs
cat ${OUTPREFIX}*_chr*.txt > "$OUTDIR/<TRAIT>_allchr_pst_eff.txt"
Critical note:
If the glob pattern does not match any files, the output redirection (>) will still create/truncate the target file, resulting in an empty <TRAIT>_allchr_pst_eff.txt. This will later cause PLINK to fail with “Error: Empty --score file.” Always verify the file is non-empty after concatenation.
Recommended checks:
ls -lh "$OUTDIR/<TRAIT>_allchr_pst_eff.txt"
wc -l  "$OUTDIR/<TRAIT>_allchr_pst_eff.txt"
ls -lh ${OUTPREFIX}*_chr*.txt | wc -l
Header note:
If the concatenated weight file has no header line (i.e., line 1 is a data row), do not use the PLINK header modifier, otherwise PLINK will skip the first SNP. If you prefer using header, add a single header line to the top of the file and ensure no duplicate headers exist after concatenation.
Step C. Compute PRS with PLINK
The PRS-CS concatenated file is structured as:
col2 = SNP (rsID)
col4 = A1 (scoring/effect allele)
col6 = BETA_post (posterior effect size)
PLINK scoring:
plink \
  --bfile "$BIM_PREFIX" \
  --score "$OUTDIR/<TRAIT>_allchr_pst_eff.txt" 2 4 6 sum \
  --out "$OUTDIR/<TRAIT>_prs_score"
If and only if your weight file has a header line:
plink \
  --bfile "$BIM_PREFIX" \
  --score "$OUTDIR/<TRAIT>_allchr_pst_eff.txt" 2 4 6 header sum \
  --out "$OUTDIR/<TRAIT>_prs_score"
Interpreting the PLINK .profile output
The main output for downstream analyses is:
<TRAIT>_prs_score.profile
Typical columns:
FID, IID: family and individual identifiers
PHENO: phenotype field echoed from the 6th column of the target .fam
CNT, CNT2: counting fields produced by PLINK for the scoring run
SCORESUM: the raw PRS (sum of dosage × weight across variants)
Downstream analyses generally use IID + SCORESUM (often after standardization).
Important note about the PHENO column (and why it may contain “3”)
PHENO in <TRAIT>_prs_score.profile is not used to compute PRS in the plink --score step. It is simply carried over from the 6th column of the target genotype dataset’s .fam file.
Therefore:
PHENO values such as “3” do not indicate an error in PRS computation.
They reflect how the phenotype field was encoded in the genotype .fam file.
The PLINK log may report “Phenotype data is quantitative” if the .fam phenotype field is not in standard case/control encoding.
For downstream association analyses, collaborators should not rely on the .fam PHENO field unless it has been explicitly curated. Instead, analyses should use a separate phenotype/covariate file provided by the analysis team (e.g., UConn collaborators).
QC and troubleshooting
1) PLINK error: “Error: Empty --score file.”
Most common causes:
The concatenated weight file <TRAIT>_allchr_pst_eff.txt is empty (e.g., cat/glob mismatch, missing PRS-CS outputs).
The weight file contains only non-parseable lines (format issues).
Checks:
wc -l <TRAIT>_allchr_pst_eff.txt
head <TRAIT>_allchr_pst_eff.txt
ls -lh <OUTPREFIX>*_chr*.txt | wc -l  # should be ~22
2) Variant matching and allele alignment
PRS correctness depends on proper variant and allele alignment between:
PGC summary stats
PRS-CS output weights
Target cohort genotypes
Always retain:
PRS-CS logs/output
PLINK .log files
because they provide counts of matched vs. skipped variants.
3) Duplicate headers after concatenation
If per-chromosome outputs include header lines, concatenation may introduce duplicate headers, which PLINK may treat as invalid records. Prefer concatenating without duplicate headers (keep one header or none). Verify by searching for repeated header tokens if your PRS-CS outputs include them.
Recommendations for collaborators (UConn) using the PRS
Use SCORESUM as the raw PRS.
Standardize PRS (z-score) before regression for interpretability and comparability.
Include appropriate covariates (ancestry PCs, sex, age, site/batch if relevant).
Use a separate curated phenotype file instead of the .fam phenotype field unless explicitly validated.
Trait list and naming convention
For reproducibility and clarity, use consistent naming:
Input GWAS: pgc<TRAIT>.clean.tsv
PRS-CS output prefix: <TRAIT>_pst_eff
Concatenated weights: <TRAIT>_allchr_pst_eff.txt
PLINK output: <TRAIT>_prs_score.profile, <TRAIT>_prs_score.log
Traits in this project:
AD, OCD, ADHD, PTSD, TS
