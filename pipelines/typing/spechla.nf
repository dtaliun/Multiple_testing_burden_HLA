process type_sample {
	cache "lenient"
	scratch true

	errorStrategy "retry"
	maxRetries 2
	maxErrors 50

	cpus params.spechla_ncpus
	memory "4 GB"
	time "1h"

	container "${params.spechla_container}"

	input:
	tuple val(sample_name), path(bam_file)

	output:
	path "${sample_name}.spechla.tar.gz", emit: out

	publishDir "${params.output_dir}", pattern: "${sample_name}.spechla.tar.gz", mode: "copy"


	"""
	# Let's make sure we fail entire pipeline if at lease one command fails in case we ever use piped commands
	set -o pipefail

	bash /opt/SpecHLA/script/ExtractHLAread.sh -b ${bam_file} -r hg38 -T ${params.reference_fasta} -s ${sample_name} -o fastq_files
	bash /opt/SpecHLA/script/whole/SpecHLA.sh -j ${params.spechla_ncpus} -n ${sample_name} -1 fastq_files/${sample_name}_extract_1.fq.gz -2 fastq_files/${sample_name}_extract_2.fq.gz -u 0 -p nonuse -w 0.3 -y 0.2 -l 0 -o output

	cd output/
	tar -czvf ${sample_name}.spechla.tar.gz -C ${sample_name} .
	cd ..

	mv output/${sample_name}.spechla.tar.gz .
	"""

}


workflow {
	// Input BAM files. We assume that the sample name is encoded as a prefix in the file name i.e. SAMPLE_XYZ.my.bam.
	bams = Channel.fromPath(params.bams).map(f -> [f.getSimpleName(), f])

	// Run
	type_sample(bams)
}
