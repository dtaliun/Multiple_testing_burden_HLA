#!/bin/bash

n_permutations=1000
phenotype_file="<tab-delimited phenotype file where the first two columns are FID and IID, respectively>"
phenotype_name="<column name of phenotype of interest in phenotype_file>"
covar_file="<tab-delimited covariate file where the first two columns are FID and IID, respectively>"
bfile="<plink bfile prefix>"
plink_exec="plink2"

output_file="minimal_pvalues.txt"

tail -n+2 ${phenotype_file} | cut -f1,2 > temp_sample_ids.txt 

for i in `seq 1 ${n_permutations}`; do
	# First: permute phenotype
	head -n1 ${phenotype_file} > permuted_pheno_${i}.txt
	tail -n+2 ${phenotype_file} | cut -f3 | shuf > temp_values_${i}.txt
	paste temp_sample_ids.txt temp_values_${i}.txt >> permuted_pheno_${i}.txt
	
	# Cleanup
	rm temp_values_${i}.txt

	# Run PLINK2
	${plink_exec} --bfile ${bfile} --glm hide-covar --mac 5 --covar ${covar_file} --covar-variance-standardize --pheno permuted_pheno_${i}.txt --pheno-name ${phenotype_name} --out gwas_${i}
	
	# Take the smallest p-value
	cut -f15 gwas_${i}.${phenotype_name}.glm.linear | tail -n+2 | LC_ALL=C sort -g | head -n1 >> ${output_file}

	# Cleanup
	rm gwas_${i}.log
	rm gwas_${i}.${phenotype_name}.glm.linear
	rm permuted_pheno_${i}.txt
done
