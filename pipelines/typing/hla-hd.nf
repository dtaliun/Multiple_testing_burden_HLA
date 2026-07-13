process type_sample {
	cache "lenient"
	scratch true

	errorStrategy "retry"
	maxRetries 2
	maxErrors 50

	cpus params.hlahd_ncpus
	memory "16 GB"
	time "1h"

	container "${params.hlahd_container}"

	input:
	tuple val(sample_name), path(bam_file)

	output:
	path "${sample_name}.hla-hd.tar.gz", emit: out

	publishDir "${params.output_dir}", pattern: "${sample_name}.hla-hd.tar.gz", mode: "copy"


	"""
	# Let's make sure we fail entire pipeline if at lease one command fails in case we ever use piped commands
	set -o pipefail

	# Optional: Sort BAM by read name, so that when translating back to FASTQ files, the reads follow same order. 
	samtools sort -n ${bam_file} -o tmp.name_sorted.bam

	# Translate back to FASTQ files
	gatk SamToFastq I=tmp.name_sorted.bam F=tmp.R1.fastq F2=tmp.R2.fastq
	rm tmp.name_sorted.bam

	# Change read name suffixes from "/1" and "/2" to " 1" and " 2". Here we follow HLA-HD documentation.
	awk 'NR%4==1{gsub("/1"," 1")}1' tmp.R1.fastq > R1.fastq
	awk 'NR%4==1{gsub("/2"," 2")}1' tmp.R2.fastq > R2.fastq
	rm tmp.R*.fastq

	# Run HLA-HD
	hlahd.sh -t ${task.cpus} -m 100 -c 0.95 -f /opt/hlahd/freq_data R1.fastq R2.fastq /opt/hlahd/HLA_gene.split.txt /opt/hlahd/dictionary ${sample_name} . > ${sample_name}_runtime.log 2>&1


   	if grep -iq "killed" "${sample_name}_runtime.log"; then
        	exit 1
    	fi
   	
	if grep -iq "error" "${sample_name}_runtime.log"; then
        	exit 1
    	fi

	# Bundle HLA-HD's results
	mkdir results
	mv ${sample_name}_runtime.log results
	mv ${sample_name}/log/*.log results
	mv ${sample_name}/result/*.txt results
	tar -czvf ${sample_name}.hla-hd.tar.gz -C results .
	"""

}


workflow {
	// Input BAM files. We assume that the sample name is encoded as a prefix in the file name i.e. SAMPLE_XYZ.my.bam.
	bams = Channel.fromPath(params.bams).map(f -> [f.getSimpleName(), f])

	// Run
	type_sample(bams)
}
