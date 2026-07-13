# Evaluating multiple testing burden in HLA association studies

We are currently reorganizing and refining our code and expect to fully populate this repository by August 1, 2026.

## 1. Contents
- `scripts/regression` - linear and Firth logistic regression models implemented in R for simulation-derived estimates of the effective number of tests.

- `scripts/phenotype_simulations` - Python 3 scripts to simulate binary and continuous phenotypes from VCF and covariate (optional) files.

- `scripts/reads_extraction` - Python 3 script that extracts reads from short-read BAM/CRAM files that are either mapped to HLA‑relevant contigs or unmapped. It is used to reduce BAM/CRAM file size prior to HLA typing in order to speed up downstream computations, as recommended by major HLA typing tools. The implementation is inspired by similar routines from `HLA-HD` and `SpecHLA`.

- `scripts/parsing` - Python 3 scripts that parse outputs from `HLA-HD`, `SpecHLA`, and `SpecImmune` and convert them to VCF files for downstream analyses, with built‑in data sanity and quality checks.

## 2. Citation

You are free to re-use data and code in this repository. If you do so, please cite: TBA.
