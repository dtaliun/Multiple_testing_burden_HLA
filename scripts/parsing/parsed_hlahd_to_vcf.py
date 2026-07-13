import argparse
from collections import OrderedDict


argparser = argparse.ArgumentParser(description = 'This script creates a VCF file from HLA-HD typing results. It transforms all HLA alleles to 2-field resolution.  Run it on the output from `parse_hlahd_results.py`.')
argparser.add_argument('-i', '--input', metavar = 'file', dest = 'in_tsv', type = str, required = True, help = 'Input TSV file with HLA alleles typed using HLA-HD and parsed using `parse_hlahd_results.py`.')
argparser.add_argument('-o', '--output', metavar = 'file', dest = 'out_vcf', type = str, required = True, help = 'Output VCF.')

HEADER = ['#sample', 'gene', 'allele1', 'allele2', 'n_candidate_pairs', 'pair_count', 'allele1_min_cov', 'allele1_mean_cov', 'allele1_cov_complete', 'allele2_min_cov', 'allele2_mean_cov', 'allele2_cov_complete', 'alleles_cov_balance']

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
    '##FORMAT=<ID=MC,Number=2,Type=Integer,Description="Minimal exon coverage for the allele represented by this record and the other allele at the same HLA gene. For homozygous genotypes (0/0 or 1/1), the assignment of values to alleles is arbitrary. For heterozygous genotypes (0/1), the first value corresponds to the allele represented by this record and the second value corresponds to the other allele at the same gene.">',
    '##FORMAT=<ID=AC,Number=2,Type=Float,Description="Average exon coverage for the allele represented by this record and the other allele at the same HLA gene. For homozygous genotypes (0/0 or 1/1), the assignment of values to alleles is arbitrary. For heterozygous genotypes (0/1), the first value corresponds to the allele represented by this record and the second value corresponds to the other allele at the same gene.">',
    '##FORMAT=<ID=CC,Number=2,Type=Integer,Description="Complete exon coverage for the allele represented by this record and the other allele at the same HLA gene. For homozygous genotypes (0/0 or 1/1), the assignment of values to alleles is arbitrary. For heterozygous genotypes (0/1), the first value corresponds to the allele represented by this record and the second value corresponds to the other allele at the same gene.">',
    '##FORMAT=<ID=CB,Number=1,Type=Float,Description="Coverage balance between the two alleles at the gene, defined as the ratio of their average exon coverages. Values close to 1 indicate balanced coverage.">',
    '##FORMAT=<ID=NCP,Number=1,Type=Integer,Description="Number of candidate allele pairs reported in the final HLA-HD typing result file. Reflects inference ambiguity and may be used as a QC indicator.">',
    '##FORMAT=<ID=NP,Number=1,Type=Integer,Description="Number of allele pair entries reported in the HLA-HD .est.txt file. Represents intermediate pair enumeration from the HLA-HD typing procedure and may be used as a QC indicator.">'
]


def get_4digit(allele):
    gene, digits = allele.split('*')

    fields = digits.split(':')

    # We expect that HLA-HD was able to type to at least 2 field (4-digit) resolution
    assert len(fields) >= 2

    # Drop a suffix (e.g. N, L, S, C, Q, A, P, G) if present
    if not fields[1][-1].isdigit():
        fields[1] = fields[1][:-1]

    return f'{gene}*{fields[0]}:{fields[1]}'


if __name__ == '__main__':
    args = argparser.parse_args()



    sample_names = OrderedDict()
    gene_names = OrderedDict()
    raw_alleles = set()
    alleles = dict()

    # Collect all sample and allele names
    with open(args.in_tsv, 'rt') as ifile:
        for i, line in enumerate(ifile, 0):
            if i == 0:
                fields = line.rstrip().split('\t')
                assert all(x == y for x, y in zip(fields, HEADER))
                continue

            fields = line.rstrip().split('\t')
            assert len(fields) == len(HEADER)
            fields = dict(zip(HEADER, fields))

            sample_names[fields['#sample']] = None
            gene_names[fields['gene']] = None
            gene_pos = HG38_HLA_GENE_APPROX_POSITION[fields['gene']]

            if fields['allele1'] != 'NA':
                raw_allele = fields['allele1']
                raw_alleles.add(raw_allele)
                allele = get_4digit(raw_allele)
                alleles[allele] = { 'POS': gene_pos, 'ID': allele }
            
            if fields['allele2'] != 'NA':
                raw_allele = fields['allele2']
                raw_alleles.add(raw_allele)
                allele = get_4digit(raw_allele)
                alleles[allele] = { 'POS': gene_pos, 'ID': allele } 

    print(f'No. of samples: {len(sample_names)}')
    print(f'No. of HLA genes: {len(gene_names)}')
    print(f'No. of all HLA alleles: {len(raw_alleles)}')
    print(f'No. of 2-field HLA alleles: {len(alleles)}')

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
    with open(args.in_tsv, 'rt') as ifile:
        # skip header line because we already checked it
        line = ifile.readline()
        for line in ifile:
            fields = line.rstrip().split('\t')
            fields = dict(zip(HEADER, fields))

            sample_name = fields['#sample']
            #if sample_name != 'HG00096': continue
            gene_name = fields['gene']
            allele1 = fields['allele1']
            allele2 = fields['allele2']

            #print(sample_name, allele1, allele2)

            if allele1 == 'NA' or allele2 == 'NA':
                # Can't have only one allele typed (i.e. only on one chromosome)
                assert allele1 == 'NA' and allele2 == 'NA'

                # If sample had no alleles typed for this gene, then we should set missing GT for all correspinding alleles in the data.
                for allele, record in alleles.items():
                    if gene_name == allele.split('*')[0].removeprefix('HLA-'):
                        record[sample_name] = './.:.:.:.:.:.:.' 
            elif allele1 == allele2:
                allele1 = get_4digit(allele1)
                allele2 = get_4digit(allele2)
                for allele, record in alleles.items():
                    if gene_name == allele.split('*')[0].removeprefix('HLA-'):
                        mc = f'{float(fields["allele1_min_cov"]):.0f},{float(fields["allele2_min_cov"]):.0f}'
                        ac = f'{float(fields["allele1_mean_cov"]):.0f},{float(fields["allele2_mean_cov"]):.0f}'
                        cc = f'{1 if fields["allele1_cov_complete"] == "True" else 0},{1 if fields["allele2_cov_complete"] == "True" else 0}'
                        cb = f'{float(fields["alleles_cov_balance"]):.2f}'
                        ncp = fields['n_candidate_pairs']
                        np = fields['pair_count']
                        if allele1 == allele: # Assign 1/1 genotype
                            record[sample_name] = f'1/1:{mc}:{ac}:{cc}:{cb}:{ncp}:{np}'
                        else: # Assign 0/0 genotype to all other alleles from this HLA gene in this sample
                            record[sample_name] = f'0/0:{mc}:{ac}:{cc}:{cb}:{ncp}:{np}'
            else:
                allele1 = get_4digit(allele1)
                allele2 = get_4digit(allele2)
                for allele, record in alleles.items():
                    if gene_name == allele.split('*')[0].removeprefix('HLA-'):
                        if allele1 == allele: # Assign '0/1' to allele1
                            mc = f'{float(fields["allele1_min_cov"]):.0f},{float(fields["allele2_min_cov"]):.0f}'
                            ac = f'{float(fields["allele1_mean_cov"]):.0f},{float(fields["allele2_mean_cov"]):.0f}'
                            cc = f'{1 if fields["allele1_cov_complete"] == "True" else 0},{1 if fields["allele2_cov_complete"] == "True" else 0}'
                            cb = f'{float(fields["alleles_cov_balance"]):.2f}'
                            ncp = fields['n_candidate_pairs']
                            np = fields['pair_count']
                            record[sample_name] = f'0/1:{mc}:{ac}:{cc}:{cb}:{ncp}:{np}'
                        elif allele2 == allele: # Assign '0/1' to allele2
                            mc = f'{float(fields["allele2_min_cov"]):.0f},{float(fields["allele1_min_cov"]):.0f}'
                            ac = f'{float(fields["allele2_mean_cov"]):.0f},{float(fields["allele1_mean_cov"]):.0f}'
                            cc = f'{1 if fields["allele2_cov_complete"] == "True" else 0},{1 if fields["allele1_cov_complete"] == "True" else 0}'
                            cb = f'{float(fields["alleles_cov_balance"]):.2f}'
                            ncp = fields['n_candidate_pairs']
                            np = fields['pair_count']
                            record[sample_name] = f'0/1:{mc}:{ac}:{cc}:{cb}:{ncp}:{np}'
                        else: # Assign 0/0 genotype to all other alleles from this HLA gene in this sample
                            mc = f'{float(fields["allele1_min_cov"]):.0f},{float(fields["allele2_min_cov"]):.0f}'
                            ac = f'{float(fields["allele1_mean_cov"]):.0f},{float(fields["allele2_mean_cov"]):.0f}'
                            cc = f'{1 if fields["allele1_cov_complete"] == "True" else 0},{1 if fields["allele2_cov_complete"] == "True" else 0}'
                            cb = f'{float(fields["alleles_cov_balance"]):.2f}'
                            ncp = fields['n_candidate_pairs']
                            np = fields['pair_count']
                            record[sample_name] = f'0/0:{mc}:{ac}:{cc}:{cb}:{ncp}:{np}'
            #print(record[allele1])
            #print(record[allele2])


    # Sanity check: according to our logic, all sample names must be present inside each record
    for allele, record in alleles.items():
        assert all(x in record for x in sample_names.keys())

    with open(args.out_vcf, 'wt') as ofile:
        for line in VCF_HEADER:
            ofile.write(f'{line}\n')
        ofile.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t')
        ofile.write('\t'.join(sample_names.keys()))
        ofile.write('\n')
        for allele, record in  alleles.items():
            ofile.write(f'chr6\t{record["POS"]}\t{record["ID"]}\tA\tC\t.\t.\t.\tGT:MC:AC:CC:CB:NCP:NP\t')
            ofile.write('\t'.join(record[x] for  x in sample_names.keys()))
            ofile.write('\n')
