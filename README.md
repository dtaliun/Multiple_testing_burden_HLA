# Evaluating multiple testing burden in HLA association studies
***


## 1. Contents

### 1.1 Data generation

- `scripts/regression` - linear and Firth logistic regression models implemented in R for simulation-derived estimates of the effective number of tests.

- `scripts/phenotype_simulations` - Python 3 scripts to simulate binary and continuous phenotypes from VCF and covariate (optional) files.

- `scripts/reads_extraction` - Python 3 script that extracts reads from short-read BAM/CRAM files that are either mapped to HLA‑relevant contigs or unmapped. It is used to reduce BAM/CRAM file size prior to HLA typing in order to speed up downstream computations, as recommended by major HLA typing tools. The implementation is inspired by similar routines from `HLA-HD` and `SpecHLA`.

- `scripts/parsing` - Python 3 scripts that parse outputs from `HLA-HD`, `SpecHLA`, and `SpecImmune` and convert them to VCF files for downstream analyses, with built‑in data sanity and quality checks.

- `pipelines/typing` - Nextflow pipelines with correspinding configuration files used to run `HLA-HD` and `SpecHLA` typing tools.

- `pipelines/gwas_simulation` - Nextflow pipelines and configuration files used to run HLA association analyses on simulated phenotypes (30,000 simulations in our case).

- `docker` - Definition files in Docker format used to build images for the `HLA-HD` and `SpecHLA` typing tools. These files specify the exact tool versions to ensure reproducibility. Images can be built with either `Docker` or `Apptainer` (we used `Apptainer`).

### 1.2. Data

- `data/1000G/Covariates` - Covariate files used to simulate phenotypes for individuals with short‑read sequencing data from the 1000 Genomes Project Phase 3. These files include the Sex variable from the 1000 Genomes metadata and the top 10 principal components computed from common LD‑pruned SNVs.

### 1.3. Analyses

## 2. Citation

You are free to re-use data and code in this repository. If you do so, please cite: TBA.
