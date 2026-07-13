import argparse
import tarfile
import pandas as pd
from collections import OrderedDict

argparser = argparse.ArgumentParser(description = 'This script parses results from SpecHLA and combines them into a single tab-delimited file.')
argparser.add_argument('-i', '--input', metavar = 'file', dest = 'in_tar', type = str, required = True, help = 'Input tar.gz file with HLA-HD results for a sample.')
argparser.add_argument('-s', '--sample', metavar = 'name', dest = 'sample', type = str, required = True, help = 'Sample name used when running HLA-HD.')
argparser.add_argument('-o', '--output', metavar = 'file', dest = 'out_tsv', type = str, required = True, help = 'Output TSV file name.')


def get_blast_bit_score(tar, filename, allele):
    with tar.extractfile(filename) as ifile:
        blast_df = pd.read_csv(ifile, sep = '\t', header = None, comment = '#',
            names = ["query", "subject", "percent_identity", "alignment_length", "n_mismatches", "n_gaps", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score"])
        
        blast_df = blast_df[blast_df.subject == allele]

        assert len(blast_df) >= 1, f'{filename} {allele}'
        
        return blast_df["bit_score"].sum()
        

if __name__ == '__main__':
    args = argparser.parse_args()

    genes = OrderedDict()

    with tarfile.open(args.in_tar, 'r:gz') as tar:
        main_result_file = f'./hla.result.details.txt' 
        with tar.extractfile(main_result_file) as ifile:
            header = None
            for line in ifile:
                line = line.decode('utf-8').strip()

                if line.startswith('#'):
                    continue
                elif line.startswith('Gene'):
                    header = line.split('\t')
                    assert header[0] == 'Gene'
                    assert header[1] == 'G_best'
                    assert header[2] == 'allele'
                    assert header[3].startswith('details:allele;Score;')
                    continue

                fields = line.split('\t')

                # Should be consistent with header
                assert len(fields) >= len(header)
                
                prefix, gene, allele_idx = fields[0].split('_')
                assert prefix == 'HLA'
                assert gene in {'A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1'} 
                assert allele_idx in {'1', '2'}

                best_allele = fields[1].rstrip(";").split(';')[0] # In many cases it is only one allele. If there are multiple, then we take the first and the rest will be in anyway checked inside `all_candidate_alleles`.
                all_candidate_alleles = fields[2].rstrip(';').split(';')
                
                alleles_2field = []
                alleles_2field_map_scores = []

                details = fields[3:]
                for detail in details:
                    detail = detail.split(';')
                    assert len(detail) >= 2
                    alleles_2field.append(detail[0])
                    alleles_2field_map_scores.append(float(detail[1]))
                
                blast_filename = f'./tmp/{fields[0]}.blast.out1'
 
                best_allele_bit_score = get_blast_bit_score(tar, blast_filename, best_allele)

                best_allele_tied_candidates = 0
                if len(all_candidate_alleles) == 1:
                    assert best_allele == all_candidate_alleles[0]
                else:
                    bit_scores = []
                    for allele in all_candidate_alleles:
                        if allele == best_allele: # first best allele will re-appear here, so we skip it
                            continue
                        allele_bit_score = get_blast_bit_score(tar, blast_filename, allele)
                        if allele_bit_score >= best_allele_bit_score:
                            best_allele_tied_candidates += 1

                best_allele_2field = None
                best_allele_2field_map_score = None
                best_allele_2field_tied_candidates = 0
                for allele_2field, allele_2field_map_score in zip(alleles_2field, alleles_2field_map_scores):
                    if best_allele_2field is None:
                        best_allele_2field = allele_2field
                        best_allele_2field_map_score = allele_2field_map_score
                    elif best_allele_2field_map_score < allele_2field_map_score:
                        best_allele_2field = allele_2field
                        best_allele_2field_map_score = allele_2field_map_score
                        best_allele_2field_tied_candidates = 0
                    elif best_allele_2field_map_score == allele_2field_map_score:
                        best_allele_2field_tied_candidates += 1
                
                gene_record = genes.setdefault(gene, {})

                gene_record[f'allele{allele_idx}'] = best_allele
                gene_record[f'allele{allele_idx}_bit_score'] = best_allele_bit_score
                gene_record[f'allele{allele_idx}_tied_candidates'] = best_allele_tied_candidates
                gene_record[f'allele{allele_idx}_2field'] = best_allele_2field
                gene_record[f'allele{allele_idx}_2field_map_score'] = best_allele_2field_map_score
                gene_record[f'allele{allele_idx}_2field_tied_candidates'] = best_allele_2field_tied_candidates



    header = ['#sample', 'gene', 'allele1', 'allele2', 'allele1_2field', 'allele2_2field',
        'allele1_bit_score', 'allele1_tied_candidates', 'allele1_2field_map_score', 'allele1_2field_tied_candidates',
        'allele2_bit_score', 'allele2_tied_candidates', 'allele2_2field_map_score', 'allele2_2field_tied_candidates'] 
    with open(args.out_tsv, 'wt') as ofile:
        ofile.write('\t'.join(header))
        ofile.write('\n')
        for gene_name, gene_record in genes.items():
            ofile.write(f'{args.sample}\t{gene_name}\t')
            ofile.write('\t'.join( str(gene_record[x]) for x in header[2:] ))
            ofile.write('\n')
            
