# Evaluating multiple testing burden in HLA association studies
***


## 1. Contents

### 1.1 Data generation

- `scripts/regression` - linear and Firth logistic regression models implemented in R for simulation-derived estimates of the effective number of tests.

- `scripts/phenotype_simulations` - Python 3 scripts to simulate binary and continuous phenotypes from VCF and covariate (optional) files.

- `scripts/reads_extraction` - Python 3 script that extracts reads from short-read BAM/CRAM files that are either mapped to HLA‑relevant contigs or unmapped. It is used to reduce BAM/CRAM file size prior to HLA typing in order to speed up downstream computations, as recommended by major HLA typing tools. The implementation is inspired by similar routines from `HLA-HD` and `SpecHLA`.

- `scripts/parsing` - Python 3 scripts that parse outputs from `HLA-HD`, `SpecHLA`, and `SpecImmune` and convert them to VCF files for downstream analyses, with built‑in data sanity and quality checks.

- `scripts/plink2_permutations` - R and Bash scripts used to obtain an estimate of the effective number of tests by permuting the true phenotype, re-running `PLINK2 --glm` with the same genotypes and covariates, and recording the minimal P-value in each permutation.

- `pipelines/typing` - Nextflow pipelines with correspinding configuration files used to run `HLA-HD` and `SpecHLA` typing tools.

- `pipelines/gwas_simulation` - Nextflow pipelines and configuration files used to run HLA association analyses on simulated phenotypes (30,000 simulations in our case).

- `docker` - Definition files in Docker format used to build images for the `HLA-HD` and `SpecHLA` typing tools. These files specify the exact tool versions to ensure reproducibility. Images can be built with either `Docker` or `Apptainer` (we used `Apptainer`).

### 1.2. Data

- `data/1000G/Covariates` - Covariate files used to simulate phenotypes for individuals with short‑read sequencing data from the 1000 Genomes Project Phase 3. These files include the Sex variable from the 1000 Genomes metadata and the top 10 principal components computed from common LD‑pruned SNVs.

> [!CAUTION]
> If you want to replicate the phenotypes we simulated, use the following starting random seeds in the `pipelines/gwas_simulation` pipelines. For binary phenotypes: ALL – 20, AFR – 21, AMR – 22, EAS – 23, EUR – 24, SAS – 25. For continuous phenotypes: ALL – 20, AFR – 21, AMR – 25, EAS – 22, EUR – 23, SAS – 24.

- `data/1000G/HLA_typed` - 2-field HLA alleles in VCF format, typed from high-depth short-read sequencing data in the 1000 Genomes Project Phase 3 using `HLA-HD`. These files were used in our main analyses of 2-field HLA alleles.

- `data/1000G_ONT/Covariates` - Covariate files used to simulate phenotypes for individuals with Oxford Nanopore Technologies (ONT) long‑read sequencing data from the 1000 Genomes Project Phase 3. These files include the Sex variable from the 1000 Genomes metadata and the top 10 principal components computed from common LD‑pruned SNVs derived from short‑read sequencing data.

> [!CAUTION]
> If you want to replicate the phenotypes we simulated for HLA alleles derived from long-read sequencing data, then use the following starting random seeds in the `pipelines/gwas_simulation` pipelines: 2-field - 421, 4-field - 421.

- `data/1000G_ONT/HLA_typed` - 2-field and 4-field HLA alleles in original CSV and converted VCF formats, typed from Oxford Nanopore Technologies (ONT) long-read sequencing data for 1000 Genomes Project Phase 3 using `SpecImmune`. These files were used in our secondary analyses to compare the multiple testing burden between 2-field and 4-field HLA alleles.

> [!IMPORTANT]
> The HLA typing of long‑read sequencing data was performed by Wang et al. (2026) and shared in the official `SpecImmune` GitHub repository [https://github.com/deepomicslab/SpecImmune](https://github.com/deepomicslab/SpecImmune). We gratefully acknowledge the work of Wang et al. (2026), which sped up our analyses. If you use these data, please be sure to cite the following:
>
> Wang, S., Wang, X., Wang, M., Zhou, Q., Wang, L., & Li, S. C. (2026). A scalable framework for comprehensive typing of polymorphic immune genes from long-read data. Advanced Science (Weinheim, Baden-Wurttemberg, Germany), e21531, e21531.

- `data/1000G/Simulated_*_phenotypes` - simulated continuous and binary phenotypes used to perform association tests with 2-field HLA alleles derived from short-read sequencing data.

- `data/1000G_ONT/Simulated_continuous_phenotypes` - simulated continuous phenotypes used to perform association tests with 2-field and 4-field HLA alleles derived from long-read sequencing data.

### 1.3. Analyses

- `notebooks` - Jupyter Notebook files containing the code used for data analyses and for generating the figures and tables included in the corresponding manuscript.

## 2. Citation
You are free to re-use data and code in this repository. If you do so, please cite the following:

- If you used our scripts or simulation results cite: TBD

- If you used only the raw data in `data/1000G_ONT/HLA_typed` you must cite the original [SpecImmune](https://github.com/deepomicslab/SpecImmune) repository and citation:

  Wang, S., Wang, X., Wang, M., Zhou, Q., Wang, L., & Li, S. C. (2026). A scalable framework for comprehensive typing of polymorphic immune genes from long-read data. Advanced Science (Weinheim, Baden-Wurttemberg, Germany), e21531, e21531.


