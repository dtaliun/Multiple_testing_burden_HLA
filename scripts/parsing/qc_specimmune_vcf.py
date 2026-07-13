import argparse
from collections import OrderedDict


argparser = argparse.ArgumentParser(description = 'Quality checks the VCF file with HLA alleles typed using SpecImmune.  Run it on the output from `specimmune_to_vcf.py`.')
argparser.add_argument('-i', '--input', metavar = 'file', dest = 'in_vcf', type = str, required = True, help = 'Input VCF, generated using `specimmune_to_vcf.py`.')
argparser.add_argument('-o', '--output', metavar = 'file', dest = 'out_vcf', type = str, required = True, help = 'Output VCF.')
argparser.add_argument('-s', '--fmiss-sample', metavar = 'number', dest = 'fmiss_sample', type = float, default = 0.1, help = 'Fraction of missing genotypes per sample.')
argparser.add_argument('-g', '--fmiss-gene', metavar = 'number', dest = 'fmiss_gene', type = float, default = 0.1, help = 'Fraction of missing genotypes per gene.')
argparser.add_argument('-r', '--min-read-count', metavar = 'number', dest = 'min_read_count', type = int, default = 10, help = 'Minimal read count supporting the allele.')


HG38_HLA_GENE_APPROX_POSITION = dict({
    
    'F': 29678000,
    'L': 29880000,
    'V': 29900000,
    'K': 29920000,
    'H': 29940000,
    'A': 29970000,
    'G': 29800000,
    'J': 29910000,
    'E': 30475000,
    'C': 31240000,
    'B': 31310000,
    'DRA': 32450000,
    'DRB9': 32500000,
    'DRB6': 32505000,
    'DRB7': 32510000,
    'DRB8': 32515000,
    'DRB3': 32540000,
    'DRB4': 32545000,
    'DRB5': 32550000,
    'DRB1': 32580000,
    'DRB2': 32520000,
    'DQA1': 32645000,
    'DQB1': 32665000,
    'DPA1': 33040000,
    'DPB1': 33080000,
    'DQA1': 32645000,
    'DQB1': 32665000,
    'DPA1': 33040000,
    'DPB1': 33080000,
    'DMA': 32980000,
    'DMB': 32910000,
    'DOA': 33030000,
    'DOB': 32820000
})

FORMAT_FIELD = ['GT', 'RC']

GENE_COVERAGE_RATIO_THRESHOLD = 0.10

if __name__ == '__main__':
    args = argparser.parse_args()

    SAMPLE_MISSINGNESS_THRESHOLD = args.fmiss_sample
    GENE_MISSINGNESS_THRESHOLD = args.fmiss_gene
    COVERAGE_THRESHOLD = args.min_read_count

    header = None
    samples = []

    samples_high_missingness = set()
    genes_high_missingness = set()
    genes_low_coverage = set()
    
    with open(args.in_vcf, 'rt') as ifile:
        for line in ifile:
            if line.startswith('#CHROM'):
                header = line.rstrip().split('\t')
                # At least one sample
                assert len(header) > 9
                samples = header[9:]
                break

    assert header is not None
    assert len(samples) > 0

    # First pass: missingness per sample
    missingness_by_sample = dict()
    with open(args.in_vcf, 'rt') as ifile:
        for line in ifile:
            if line.startswith('#'):
                continue
            fields =  dict(zip(header, line.rstrip().split('\t')))
            
            allele = fields['ID']
            gene = allele.split('*')[0]
            if gene.removeprefix('HLA-') not in {'A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRA', 'DRB1'}:
                continue

            for sample in samples:
                if sample not in missingness_by_sample:
                    missingness_by_sample[sample] = {'n_missing': 0, 'n_total': 0}
                sample_info = dict(zip(FORMAT_FIELD, fields[sample].split(':')))
                missingness_by_sample[sample]['n_total'] += 1
                if sample_info['GT'] == './.':
                    missingness_by_sample[sample]['n_missing'] += 1

    for sample, counts in missingness_by_sample.items():
        missingness = counts['n_missing'] / counts['n_total']
        if missingness > SAMPLE_MISSINGNESS_THRESHOLD:
            samples_high_missingness.add(sample)
    print(f'{len(samples_high_missingness)} samples have >{SAMPLE_MISSINGNESS_THRESHOLD * 100:.2f}% missing genotypes in classical HLA genes and will be excluded.')


    updated_samples = []
    for sample in samples:
        if sample not in samples_high_missingness:
            updated_samples.append(sample)
    samples = updated_samples
    print(f'{len(samples)} samples remaining.')

    # Second pass: check missigness per HLA gene
    missingness_by_gene = dict()
    with open(args.in_vcf, 'rt') as ifile:
        for line in ifile:
            if line.startswith('#'):
                continue
            fields = dict(zip(header, line.rstrip().split('\t')))

            allele = fields['ID']
            gene = allele.split('*')[0]
  
            # FORMAT field should be as expected
            assert all(x == y for x, y in zip(FORMAT_FIELD,  fields['FORMAT'].split(':')))

            n_missing = 0
            for sample in samples:
                sample_info = dict(zip(FORMAT_FIELD, fields[sample].split(':')))
                if sample_info['GT'] == './.':
                    n_missing += 1
            
            gene_missingness = missingness_by_gene.setdefault(gene, {'n_missing': 0, 'n_total': 0})
            gene_missingness['n_missing'] += n_missing
            gene_missingness['n_total'] += len(samples)


    print("Ratio of missing alleles per HLA gene:")
    for gene, gene_missingness in missingness_by_gene.items():
        missingness = gene_missingness['n_missing'] / gene_missingness['n_total']
        if missingness > GENE_MISSINGNESS_THRESHOLD:
            genes_high_missingness.add(gene)
        print(f'\t{gene} has {missingness * 100:.2f}% of alleles (GT field) missing')
    print(f"{len(genes_high_missingness)} genes will be removed due to high missigness (>{GENE_MISSINGNESS_THRESHOLD * 100:.0f}%):", ",".join(genes_high_missingness))


    # Third pass: set GT to missing for alleles with minimal read count below the threshold
    coverage_by_gene = dict()
    with open(args.in_vcf, 'rt') as ifile:
        for line in ifile:
            if line.startswith('#'):
                continue
            fields = dict(zip(header, line.rstrip().split('\t')))

            allele = fields['ID']
            gene = allele.split('*')[0]
            if gene in genes_high_missingness:
                continue

            # FORMAT field should be as expected
            assert all(x == y for x, y in zip(FORMAT_FIELD,  fields['FORMAT'].split(':')))

            n_total = 0
            n_low_coverage = 0
            for sample in samples:
                sample_info = dict(zip(FORMAT_FIELD, fields[sample].split(':')))
                if sample_info['GT'] == './.':
                    continue
                n_total += 1
                rc = min(int(x) for x in sample_info['RC'].split(","))
                if rc < COVERAGE_THRESHOLD:
                    n_low_coverage += 1
            gene_coverage = coverage_by_gene.setdefault(gene, {'n_low_coverage': 0, 'n_total': 0})
            gene_coverage['n_low_coverage'] += n_low_coverage
            gene_coverage['n_total'] += n_total

    print("Ratio of alleles with low coverage per HLA gene:")
    for gene, gene_coverage in coverage_by_gene.items():
        low_coverage_ratio = gene_coverage['n_low_coverage'] / gene_coverage['n_total']
        if low_coverage_ratio > GENE_COVERAGE_RATIO_THRESHOLD:
            genes_low_coverage.add(gene)
        print(f'\t{gene} has {low_coverage_ratio * 100:.2f}% of alleles with low coverage (<{COVERAGE_THRESHOLD})')
    print(f'{len(genes_low_coverage)} genes will be removed due to high ratio of alleles with low coverage (>{GENE_COVERAGE_RATIO_THRESHOLD * 100:.0f}%):', ",".join(genes_low_coverage))


    with open(args.in_vcf, 'rt') as ifile, open(args.out_vcf, 'wt') as ofile:
        for line in ifile:
            if line.startswith('#'):
                if line.startswith('#CHROM'):
                    for i, x in enumerate(header, 1):
                        if i == 1:
                            ofile.write(x)
                        elif i <= 9:
                            ofile.write(f'\t{x}')
                        else:
                            if x not in samples:
                                continue
                            ofile.write(f'\t{x}')
                    ofile.write('\n')
                else:
                    ofile.write(line)
                continue
            fields = dict(zip(header, line.rstrip().split('\t')))
            
            allele = fields['ID']
            gene = allele.split('*')[0]
            if gene in genes_high_missingness or gene in genes_low_coverage:
                continue

            for i, x in enumerate(header, 1):
                if i == 1:
                    ofile.write(f'{fields[x]}')
                elif i <= 9:
                    ofile.write(f'\t{fields[x]}')
                else:
                    if x not in samples:
                        continue
                    ofile.write(f'\t')
                    sample_info = dict(zip(FORMAT_FIELD, fields[x].split(':')))
                    if sample_info['GT'] != './.':
                        rc = min(int(x) for x in sample_info['RC'].split(","))
                        if rc < COVERAGE_THRESHOLD:
                            sample_info['GT'] = './.'
                    ofile.write(':'.join(sample_info[f] for f in FORMAT_FIELD))

            ofile.write('\n')


