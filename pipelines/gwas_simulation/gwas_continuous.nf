process create_phenotypes {
	cache "lenient"
	scratch true

	cpus 1
        memory "4 GB"
        time "1h"

	input:
	each seed
	val n_pheno
	path covars
	path vcf

	output:
	path("phenos_${seed}_${n_pheno}_cov.tsv"), emit: cov
	path("phenos_${seed}_${n_pheno}_cov.diagnostic.tsv"), emit: cov_diagnostics
	path("phenos_${seed}_${n_pheno}_nocov.tsv"), emit: nocov
	path("phenos_${seed}_${n_pheno}_nocov.diagnostic.tsv"), emit: nocov_diagnostics

	publishDir "${params.output_dir}", pattern: "phenos_${seed}_${n_pheno}*tsv", mode: "copy"

	"""
	${params.simulate_pheno_exec} -i ${vcf} -c ${covars} -s ${seed} -n ${n_pheno} -o phenos_${seed}_${n_pheno}_cov.tsv -d phenos_${seed}_${n_pheno}_cov.diagnostic.tsv
	${params.simulate_pheno_exec} -i ${vcf} -s ${seed} -n ${n_pheno} -o phenos_${seed}_${n_pheno}_nocov.tsv -d phenos_${seed}_${n_pheno}_nocov.diagnostic.tsv
	"""
}


def combine_gwas_output(base_name) {
	"""
	find . -name "*.linear.tsv" > files.txt
	head -n1 `head -n1 files.txt` | gzip -c > pvalues.tsv.gz
	while read -r f; do tail -n+2 \${f}; done < files.txt | gzip -c >> pvalues.tsv.gz
	mv pvalues.tsv.gz ${base_name}.pvalues.tsv.gz
	"""
}


process run_gwas_with_covariates {
	cache "lenient"
	scratch true

	cpus 1
	memory "2 GB"
	time "1d"

	input:
	each path(phenos)
	path vcf
	path covars

	output:
	path("${phenos.baseName}.pvalues.tsv.gz")

	"""
	${params.regression_exec} -v ${vcf} -p ${phenos} -c ${covars} -o ${phenos.baseName}
	${combine_gwas_output(phenos.baseName)}
	"""
}


process run_gwas_without_covariates {
	cache "lenient"
	scratch true

	cpus 1
	memory "2 GB"
	time "1d"

	input:
	each path(phenos)
	path vcf

	output:
	path("${phenos.baseName}.pvalues.tsv.gz")

	"""
	${params.regression_exec} -v ${vcf} -p ${phenos} -o ${phenos.baseName}
	${combine_gwas_output(phenos.baseName)}
	"""
}


process combine {
	cache "lenient"
	scratch false

	cpus 1
	memory "4 GB"
	time "1h"

	input:
	tuple val(prefix), path(gwas_files)

	output:
	path("${prefix}.tsv.gz")

	publishDir "${params.output_dir}", pattern: "${prefix}.tsv.gz", mode: "copy"

	"""
	find . -name "*.pvalues.tsv.gz" > files.txt
	gzip -dc `head -n1 files.txt` | head -n1 | gzip -c > pvalues.tsv.gz
	while read -r f; do gzip -dc \${f} | tail -n+2; done < files.txt | gzip -c >> pvalues.tsv.gz 
	
	mv pvalues.tsv.gz ${prefix}.tsv.gz
	"""
}


workflow {
	// Generate unique random integers for random seeds in each independent simulation job
	def rng = new Random(params.seed)
	def uniqueSimSeeds = new LinkedHashSet()
	while (uniqueSimSeeds.size() < params.n_sim_jobs) {
		uniqueSimSeeds.add(rng.nextInt().abs())
	}
	seeds = Channel.fromList(uniqueSimSeeds.asList())
	//seeds.view()

	// Input VCF with HLA alleles. We expect a single file.
	vcf = Channel.fromPath(params.vcf).first()
	//vcf.view()

	// Input covariates TSV file in PLINK-compatible format. We expect a single file.
	covariates = Channel.fromPath(params.covariates_tsv).first()
	//covariates.view()

	// Simulate phenotype files
	phenos = create_phenotypes(seeds, params.n_pheno_per_sim_job, covariates, vcf)
	//phenos.cov.view()
	//phenos.nocov.view()

	// Run GWASs
	gwas_cov = run_gwas_with_covariates(phenos.cov, vcf, covariates)
	gwas_nocov = run_gwas_without_covariates(phenos.nocov, vcf)
	//gwas_cov.view()
	//gwas_nocov.view()

	combine(gwas_cov.collect().map { files -> ["seed_${params.seed}_cov", files] }.mix(gwas_nocov.collect().map { files -> ["seed_${params.seed}_nocov", files] }))
}
