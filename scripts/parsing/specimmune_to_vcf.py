import argparse
from collections import OrderedDict


argparser = argparse.ArgumentParser(description = 'This script creates a VCF file from SpecImmune output shared in https://doi.org/10.1002/advs.202521531.')
argparser.add_argument('-i', '--input', metavar = 'file', dest = 'in_csv', type = str, required = True, help = 'Input CSV file.')
argparser.add_argument('-o', '--output', metavar = 'file', dest = 'out_vcf', type = str, required = True, help = 'Output VCF.')
argparser.add_argument('-d', '--digits', metavar = 'number', dest = 'in_digits', type = int, default = 8, help = 'Resolution in digits: 2, 4, 6, 8.')

HEADER = ['Locus', 'Chromosome', 'Genotype', 'Match_info', 'Reads_num', 'Step1_type', 'One_guess', 'Sample']

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

VCF_HEADER = [
    '##fileformat=VCFv4.2',
    '##contig=<ID=chr6,length=170805979>',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=RC,Number=2,Type=Integer,Description="Read count for the allele represented by this record and the other allele at the same HLA gene. For homozygous genotypes (0/0 or 1/1), the assignment of values to alleles is arbitrary. For heterozygous genotypes (0/1), the first value corresponds to the allele represented by this record and the second value corresponds to the other allele at the same gene.">',
]


def get_resolution(allele, n_digits):
    gene, digits = allele.split('*')

    assert n_digits in {2, 4, 6, 8}

    fields = digits.split(':')
    if len(fields) * 2 < n_digits:
        return None

    # Drop a suffix (e.g. N, L, S, C, Q, A, P, G) if present
    if not fields[1][-1].isdigit():
        fields[1] = fields[1][:-1]

    new_allele = f'{gene}*' + ':'.join(f'{fields[i]}' for i in range(0, n_digits // 2))
    return new_allele


if __name__ == '__main__':
    args = argparser.parse_args()

    n_digits = args.in_digits

    samples = OrderedDict()
    gene_names = OrderedDict()
    raw_alleles = set()
    alleles = dict()

    # Collect all sample and allele names
    with open(args.in_csv, 'rt') as ifile:
        for i, line in enumerate(ifile, 0):
            if i == 0:
                fields = line.rstrip().split(',')
                assert all(x == y for x, y in zip(fields, HEADER))
                continue

            fields = line.rstrip().split(',')
            assert len(fields) == len(HEADER)
            fields = dict(zip(HEADER, fields))

            samples[fields['Sample']] = dict()
            gene_name = fields['Locus'].removeprefix('HLA-')
            if gene_name in HG38_HLA_GENE_APPROX_POSITION: # we will skip loci like MICA, MICB, TAP1, TAP2
                gene_names[gene_name] = None
                gene_pos = HG38_HLA_GENE_APPROX_POSITION[gene_name]
            else:
                continue
 
            raw_allele = fields['Genotype'].split(';')
            assert len(raw_allele) > 0
            if len(raw_allele) == 1: # only one candidate allele
                raw_allele = raw_allele[0]
            else: # multiple candidate alleles, so we take best guess field
                raw_allele = fields['One_guess']
                
            if raw_allele == '': # skip missing genotypes
                continue

            assert fields['Chromosome'] in {'1', '2'}, fields

            raw_alleles.add(raw_allele)
            allele = get_resolution(raw_allele, n_digits)
            if allele is None: # did not type at the desired resolution
                continue
            alleles[allele] = { 'POS': gene_pos, 'ID': allele }
 

    print(f'No. of samples: {len(samples)}')
    print(f'No. of HLA genes: {len(gene_names)}')
    print(f'No. of all HLA alleles: {len(raw_alleles)}')
    print(f'No. of {n_digits // 2}-field HLA alleles: {len(alleles)}')

    # Order alleles by apporximate HLA gene position and assign a unique pseudo-position withn the gene. 
    alleles = OrderedDict(sorted(alleles.items(), key = lambda item: (item[1]['POS'], item[1]['ID']) ))
    current_gene_name = None
    for allele, record in  alleles.items():
        gene_name = allele.split('*')[0]
        if current_gene_name == None or current_gene_name != gene_name:
            current_gene_name = gene_name
            i = 0
        else:
            i += 1
        record['POS'] += i
    
    # Collect HLA genotypes:
    with open(args.in_csv, 'rt') as ifile:
        # skip header line because we already checked it
        line = ifile.readline()
        for line in ifile:
            fields = line.rstrip().split(',')
            fields = dict(zip(HEADER, fields))

            gene_name = fields['Locus'].removeprefix('HLA-')
            sample_name = fields['Sample']
            chr_copy = fields['Chromosome']

            if gene_name not in gene_names:
                continue

            if gene_name not in samples[sample_name]:
                samples[sample_name][gene_name] = {}

            raw_allele = fields['Genotype'].split(';')
            assert len(raw_allele) > 0
            if len(raw_allele) == 1: # only one candidate allele
                raw_allele = raw_allele[0]
            else: # multiple candidate alleles, so we take best guess field
                raw_allele = fields['One_guess']

            if raw_allele == '':
                samples[sample_name][gene_name][f'Allele{chr_copy}'] = None
                samples[sample_name][gene_name][f'Reads{chr_copy}'] = None
            else:
                samples[sample_name][gene_name][f'Allele{chr_copy}'] = get_resolution(raw_allele, n_digits)
                samples[sample_name][gene_name][f'Reads{chr_copy}'] = int(fields['Reads_num'])

    # Sanity check: number of genes should be  the same across all individuals
    assert all(len(gene_names) == len(sample_genes) for sample_name, sample_genes in samples.items())

    # Transform from per-sample to per-allele genotypes
    for sample_name, sample_genes in samples.items():
        for gene_name, sample_alleles in sample_genes.items():
            allele1 = sample_alleles['Allele1']
            allele2 = sample_alleles['Allele2']
            reads1 = sample_alleles['Reads1']
            reads2 = sample_alleles['Reads2']

            if allele1 is None or allele2 is None:
                # If sample had no alleles typed for this gene, then we should set missing GT for all correspinding alleles in the data.
                for allele, record in alleles.items():
                    if gene_name == allele.split('*')[0].removeprefix('HLA-'):
                        record[sample_name] = './.:.'
            elif allele1 == allele2:
                for allele, record in alleles.items():
                    if gene_name == allele.split('*')[0].removeprefix('HLA-'):
                        if allele1 == allele: # Assign 1/1 genotype
                            record[sample_name] = f'1/1:{reads1},{reads2}'
                        else: # Assign 0/0 genotype to all other alleles from this HLA gene in this sample
                            record[sample_name] = f'0/0:{reads1},{reads2}'
            else: 
                for allele, record in alleles.items():
                    if gene_name == allele.split('*')[0].removeprefix('HLA-'):
                        if allele1 == allele: # Assign '0/1' to allele1
                            record[sample_name] = f'0/1:{reads1},{reads2}'
                        elif allele2 == allele: # Assigne '0/1' to allele2
                            record[sample_name] = f'0/1:{reads2},{reads1}'
                        else: # Assign 0/0 genotype to all other alleles from this HLA gene in this sample
                            record[sample_name] = f'0/0:{reads1},{reads2}'

    # Sanity check: according to our logic, all sample names must be present inside each record
    for allele, record in alleles.items():
        assert all(x in record for x in samples.keys())

    with open(args.out_vcf, 'wt') as ofile:
        for line in VCF_HEADER:
            ofile.write(f'{line}\n')
        ofile.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t')
        ofile.write('\t'.join(samples.keys()))
        ofile.write('\n')
        for allele, record in  alleles.items():
            ofile.write(f'chr6\t{record["POS"]}\t{record["ID"]}\tA\tC\t.\t.\t.\tGT:RC\t')
            ofile.write('\t'.join(record[x] for  x in samples.keys()))
            ofile.write('\n')
