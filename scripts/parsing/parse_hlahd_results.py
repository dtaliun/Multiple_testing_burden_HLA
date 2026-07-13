import argparse
import tarfile
from collections import OrderedDict

argparser = argparse.ArgumentParser(description = 'This script parses results from HLA-HD and combines them into a single tab-delimited file.')
argparser.add_argument('-i', '--input', metavar = 'file', dest = 'in_tar', type = str, required = True, help = 'Input tar.gz file with HLA-HD results for a sample.')
argparser.add_argument('-s', '--sample', metavar = 'name', dest = 'sample', type = str, required = True, help = 'Sample name used when running HLA-HD.')
argparser.add_argument('-o', '--output', metavar = 'file', dest = 'out_tsv', type = str, required = True, help = 'Output TSV file name.')

def parse_allele_coverage_string(allele_coverage):
    coverages = []
    complete_coverage = True
    for c in  allele_coverage.split(','):
        exon, exon_coverage, status = c.split(':')
        coverages.append(float(exon_coverage))
        # status can be only 'comp.' or 'incomp.'
        assert status.startswith('comp.') or status.startswith('incomp.')
        if status.startswith('incomp.'):
            complete_coverage = False
    return coverages, complete_coverage


if __name__ == '__main__':
    args = argparser.parse_args()

    genes = OrderedDict()

    with tarfile.open(args.in_tar, 'r:gz') as tar:
        main_result_file = f'./{args.sample}_final.result.txt' 
        with tar.extractfile(main_result_file) as ifile:
            for line in ifile:
                line = line.decode('utf-8').strip()
                
                fields = line.split('\t')
 
                # Should be at least 3 columns
                assert len(fields) >= 3

                # If there are more than 3 columns, then the difference should be multiple by 2
                assert (len(fields) - 3) % 2 == 0

                gene_name = fields[0]
                allele1 = fields[1]
                allele2 = fields[2]

                # No duplicated genes
                assert gene_name not in genes

                # No partially typed genotypes
                assert (allele1 == 'Not typed' and allele2 == 'Not typed') or (allele1 != 'Not typed' and allele2 != 'Not typed')

                # Allele1 can't be '-'
                assert allele1 != '-'

                # Replace 'Not typed' with 'NA'
                if allele1 == 'Not typed' or allele2 == 'Not typed':
                    allele1 = 'NA'
                    allele2 = 'NA'

                genes[gene_name] = { 
                    "allele1": allele1,
                    "allele2": allele2 if allele2 != "-" else allele1,
                    "n_candidate_pairs": (len(fields) - 3) % 2
                }


        for gene_name, gene_result in genes.items():
            gene_result_file = f'./{args.sample}_{gene_name}.est.txt'
            with tar.extractfile(gene_result_file) as ifile:
                line = ifile.readline().decode('utf-8').strip()
                if line.startswith('No candidate.'):
                    genes[gene_name]['pair_count'] = 'NA'
                    genes[gene_name]['allele1_min_cov'] = 'NA'
                    genes[gene_name]['allele1_mean_cov'] = 'NA'
                    genes[gene_name]['allele1_cov_complete'] = 'NA'
                    genes[gene_name]['allele2_min_cov'] = 'NA'
                    genes[gene_name]['allele2_mean_cov'] = 'NA'
                    genes[gene_name]['allele2_cov_complete'] = 'NA'
                    genes[gene_name]['alleles_cov_balance'] = 'NA'
                    continue
                assert line.startswith('#Pair count')
                pair_count = int(line.split()[-1])

                line = ifile.readline().decode('utf-8').strip()
                assert line.startswith('#Best allele pair')
                best_pair = int(line.split()[-1])

                allele1_cov_field = None
                allele2_cov_field = None
                for i in range(0, min(best_pair, pair_count)):
                    line = ifile.readline().decode('utf-8').strip()
                    fields = line.split('\t')
                    assert len(fields) == 4, gene_result_file
                    if i == best_pair - 1:
                        allele1 = fields[0]
                        allele2 = fields[1]
                        allele1_cov_field = fields[2]
                        allele2_cov_field = fields[3]

                        # Allele1 can't be '-'
                        assert (allele1 != '-' and allele1_cov_field != '-')

                        # If allele2 is '-', the exon coverages are '-' as well
                        assert (allele2 == '-' and allele2_cov_field == '-') or (allele2 != '-' and allele2_cov_field != '-')

                # No missing coverages
                assert allele1_cov_field is not None
                assert allele2_cov_field is not None 

                allele1_cov, allele1_cov_status = parse_allele_coverage_string(allele1_cov_field)
                if allele2_cov_field == '-':
                    allele2_cov = allele1_cov
                    allele2_cov_status = allele1_cov_status
                else:
                    allele2_cov, allele2_cov_status = parse_allele_coverage_string(allele2_cov_field) 

                genes[gene_name]['pair_count'] = pair_count
                genes[gene_name]['allele1_min_cov'] = min(allele1_cov)
                genes[gene_name]['allele1_mean_cov'] = sum(allele1_cov) / len(allele1_cov)
                genes[gene_name]['allele1_cov_complete'] = allele1_cov_status
                genes[gene_name]['allele2_min_cov'] = min(allele2_cov)
                genes[gene_name]['allele2_mean_cov'] = sum(allele2_cov) / len(allele2_cov)
                genes[gene_name]['allele2_cov_complete'] = allele2_cov_status
                genes[gene_name]['alleles_cov_balance'] = min(sum(allele1_cov) / len(allele1_cov), sum(allele2_cov) / len(allele2_cov)) / max(sum(allele1_cov) / len(allele1_cov), sum(allele2_cov) / len(allele2_cov))
               
    header = ['#sample', 'gene', 'allele1', 'allele2', 'n_candidate_pairs', 'pair_count', 'allele1_min_cov', 'allele1_mean_cov', 'allele1_cov_complete', 'allele2_min_cov', 'allele2_mean_cov', 'allele2_cov_complete', 'alleles_cov_balance'] 
    with open(args.out_tsv, 'wt') as ofile:
        ofile.write('\t'.join(header))
        ofile.write('\n')
        for gene_name, values in genes.items():
            ofile.write(f'{args.sample}\t{gene_name}\t')
            ofile.write('\t'.join( str(genes[gene_name][x]) for x in header[2:] ))
            ofile.write('\n')
