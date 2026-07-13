import argparse
from collections import OrderedDict


argparser = argparse.ArgumentParser(description = 'Quality checks the VCF file with HLA alleles typed using HLA-HD.  Run it on the output from `parsed_hlahd_to_vcf.py`.')
argparser.add_argument('-i', '--input', metavar = 'file', dest = 'in_vcf', type = str, required = True, help = 'Input VCF, generated using `parsed_hlahd_to_vcf.py`.')
argparser.add_argument('-o', '--output', metavar = 'file', dest = 'out_vcf', type = str, required = True, help = 'Output VCF.')


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

FORMAT_FIELD = ['GT', 'MC', 'AC', 'CC', 'CB', 'NCP', 'NP']
GENE_MISSIGNESS_THRESHOLD = 0.01
GENE_COVERAGE_RATIO_THRESHOLD = 0.01
COVERAGE_THRESHOLD = 10
ALLELE_BALANCE_THRESHOLD = 0.5

if __name__ == '__main__':
    args = argparser.parse_args()

    header = None
    samples = []

    genes_high_missigness = set()
    genes_low_coverage = set()
    genes_umbiguous_alleles = set()

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

    # First pass: check missigness per HLA gene
    missigness_by_gene = dict()

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
            
            gene_missigness = missigness_by_gene.setdefault(gene, {'n_missing': 0, 'n_total': 0})
            gene_missigness['n_missing'] += n_missing
            gene_missigness['n_total'] += len(samples)

    print("Ratio of missing alleles per HLA gene:")
    for gene, gene_missigness in missigness_by_gene.items():
        missigness = gene_missigness['n_missing'] / gene_missigness['n_total']
        if missigness > GENE_MISSIGNESS_THRESHOLD:
            genes_high_missigness.add(gene)
        print(f'\t{gene} has {missigness * 100:.2f}% of alleles (GT field) missing')

    print(f"{len(genes_high_missigness)} genes will be removed due to high missigness (>{GENE_MISSIGNESS_THRESHOLD * 100:.0f}%):", ",".join(genes_high_missigness))


    # Second pass: set GT to missing for alleles with minimal exon coverage below the threshold or incomplete exon coverage
    coverage_by_gene = dict()

    with open(args.in_vcf, 'rt') as ifile:
        for line in ifile:
            if line.startswith('#'):
                continue
            fields = dict(zip(header, line.rstrip().split('\t')))

            allele = fields['ID']
            gene = allele.split('*')[0]
            if gene in genes_high_missigness:
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
                mc = min(int(x) for x in sample_info['MC'].split(","))
                cc = min(int(x) for x in sample_info['CC'].split(","))
                if mc < COVERAGE_THRESHOLD or cc == 0:
                    n_low_coverage += 1
                
            gene_coverage = coverage_by_gene.setdefault(gene, {'n_low_coverage': 0, 'n_total': 0})
            gene_coverage['n_low_coverage'] += n_low_coverage
            gene_coverage['n_total'] += n_total

    print("Ratio of alleles with low coverage per HLA gene:")
    for gene, gene_coverage in coverage_by_gene.items():
        low_coverage_ratio = gene_coverage['n_low_coverage'] / gene_coverage['n_total']
        if low_coverage_ratio > GENE_COVERAGE_RATIO_THRESHOLD:
            genes_low_coverage.add(gene)
        print(f'\t{gene} has {low_coverage_ratio * 100:.2f}% of alleles with low coverage (<{COVERAGE_THRESHOLD} or incomplete exon coverage)')
    print(f'{len(genes_low_coverage)} genes will be removed due to high ratio of alleles with low coverage (>{GENE_COVERAGE_RATIO_THRESHOLD * 100:.0f}%):', ",".join(genes_low_coverage))


    # Third pass: set GT to missing for alleles with high inbalance
    umbiguous_alleles_by_gene = dict()

    with open(args.in_vcf, 'rt') as ifile:
        for line in ifile:
            if line.startswith('#'):
                continue
            fields = dict(zip(header, line.rstrip().split('\t')))

            allele = fields['ID']
            gene = allele.split('*')[0]
            if gene in genes_high_missigness or gene in genes_low_coverage:
                continue

            # FORMAT field should be as expected
            assert all(x == y for x, y in zip(FORMAT_FIELD,  fields['FORMAT'].split(':')))

            n_total = 0
            n_alleles = 0
            for sample in samples:
                sample_info = dict(zip(FORMAT_FIELD, fields[sample].split(':')))
                if sample_info['GT'] == './.':
                    continue
                mc = min(int(x) for x in sample_info['MC'].split(","))
                cc = min(int(x) for x in sample_info['CC'].split(","))
                if mc < COVERAGE_THRESHOLD or cc == 0:
                    continue
                n_total += 1
                cb = float(sample_info['CB'])
                ncp = int(sample_info['NCP'])
                if cb < ALLELE_BALANCE_THRESHOLD or ncp > 0:
                    n_alleles += 1                
                
            gene_umbiguous = umbiguous_alleles_by_gene.setdefault(gene, {'n_alleles': 0, 'n_total': 0})
            gene_umbiguous['n_alleles'] += n_alleles
            gene_umbiguous['n_total'] += n_total

    print("Ratio of alleles with inbalance or alternative pairs:")
    for gene, gene_umbiguous in umbiguous_alleles_by_gene.items():
        umbiguous_ratio = gene_umbiguous['n_alleles'] / gene_umbiguous['n_total']
        if umbiguous_ratio > 0.01:
            genes_umbiguous_alleles.add(gene)
        print(f'\t{gene} has {umbiguous_ratio * 100:.2f}% of alleles with inbalance (<{ALLELE_BALANCE_THRESHOLD}) or alternative pairs')
    print(f'{len(genes_umbiguous_alleles)} genes will be removed:', ','.join(genes_umbiguous_alleles))

    with open(args.in_vcf, 'rt') as ifile, open(args.out_vcf, 'wt') as ofile:
        for line in ifile:
            if line.startswith('#'):
                ofile.write(line)
                continue
            fields = dict(zip(header, line.rstrip().split('\t')))
            
            allele = fields['ID']
            gene = allele.split('*')[0]
            if gene in genes_high_missigness or gene in genes_low_coverage or gene in genes_umbiguous_alleles:
                continue

            for i, x in enumerate(header, 1):
                if i == 1:
                    ofile.write(f'{fields[x]}')
                elif i <= 9:
                    ofile.write(f'\t{fields[x]}')
                else:
                    ofile.write(f'\t')
                    sample_info = dict(zip(FORMAT_FIELD, fields[x].split(':')))
                    if sample_info['GT'] != './.':
                        mc = min(int(x) for x in sample_info['MC'].split(","))
                        cc = min(int(x) for x in sample_info['CC'].split(","))
                        cb = float(sample_info['CB'])
                        ncp = int(sample_info['NCP'])
                        if mc < COVERAGE_THRESHOLD or cc == 0:
                            sample_info['GT'] = './.'
                        elif cb < ALLELE_BALANCE_THRESHOLD or ncp > 0:
                            sample_info['GT'] = './.'
                    ofile.write(':'.join(sample_info[f] for f in FORMAT_FIELD))

            ofile.write('\n')


