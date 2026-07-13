process create_phenotypes {
	cache "lenient"
	scratch true

	cpus 1
        memory "4 GB"
        time "1h"

	input:
	each seed
	val case_fraction
	val n_pheno
	path covars
	path vcf

	output:
	tuple val(case_fraction), path("phenos_s${seed}_n${n_pheno}_f${case_fraction}_cov.tsv"), path("phenos_s${seed}_n${n_pheno}_f${case_fraction}_nocov.tsv"), emit: data
	tuple path("phenos_s${seed}_n${n_pheno}_f${case_fraction}_cov.diagnostic.tsv"), path("phenos_s${seed}_n${n_pheno}_f${case_fraction}_nocov.diagnostic.tsv"), emit: diagnostics

	publishDir "${params.output_dir}", pattern: "phenos_s${seed}_n${n_pheno}_f${case_fraction}_*tsv", mode: "copy"

	"""
	${params.simulate_pheno_exec} -i ${vcf} -f ${case_fraction} -c ${covars} -s ${seed} -n ${n_pheno} -o phenos_s${seed}_n${n_pheno}_f${case_fraction}_cov.tsv -d phenos_s${seed}_n${n_pheno}_f${case_fraction}_cov.diagnostic.tsv
	${params.simulate_pheno_exec} -i ${vcf} -f ${case_fraction} -s ${seed} -n ${n_pheno} -o phenos_s${seed}_n${n_pheno}_f${case_fraction}_nocov.tsv -d phenos_s${seed}_n${n_pheno}_f${case_fraction}_nocov.diagnostic.tsv
	"""
}


def combine_gwas_output(base_name) {
	"""
	find . -name "*.firth.tsv" > files.txt
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
	time "2d"

	input:
	tuple val(case_fraction), path(phenos)
	path vcf
	path covars

	output:
	tuple val(case_fraction), path("${phenos.baseName}.pvalues.tsv.gz")

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
	tuple val(case_fraction), path(phenos)
	path vcf

	output:
	tuple val(case_fraction), path("${phenos.baseName}.pvalues.tsv.gz")

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
	phenos = create_phenotypes(seeds, Channel.fromList([0.1, 0.2, 0.3]), params.n_pheno_per_sim_job, covariates, vcf)
	//phenos.data.view()

	phenos_split = phenos.data.multiMap { it ->
		cov: [ it[0], it[1] ]
		nocov: [ it[0], it[2] ]
    	}

	//phenos_split.cov.view()
	//phenos_split.nocov.view()

	// Run GWASs
	gwas_cov = run_gwas_with_covariates(phenos_split.cov, vcf, covariates)
	gwas_nocov = run_gwas_without_covariates(phenos_split.nocov, vcf)
	//gwas_cov.view()
	//gwas_nocov.view()

	combine(gwas_cov.groupTuple().map { it -> ["seed_${params.seed}_ratio_${it[0]}_cov", it[1]] }.mix(
		gwas_nocov.groupTuple().map { it -> ["seed_${params.seed}_ratio_${it[0]}_nocov", it[1]] }))
}
