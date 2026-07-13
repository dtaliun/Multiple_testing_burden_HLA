import argparse
import os
import pysam
import time

argparser = argparse.ArgumentParser(description = 'Extract all reads aligned to the region of interest and any unmapped read from BAM/CRAM.')
argparser.add_argument('-i', '--input', metavar = 'file', dest = 'in_bam', type = str, required = True, help = 'Input BAM/CRAM file.')
argparser.add_argument('-o', '--output', metavar = 'file', dest = 'out_bam', type = str, required = True, help = 'Output BAM file.')
argparser.add_argument('-r', '--reference', metavar = 'file', dest = 'in_fasta', type = str, default = "", help = 'Reference FASTA file if CRAM input is used.')
argparser.add_argument('-t', '--threads', metavar = 'number', dest = 'n_threads', type = int, default = 1, help = 'Number of threads (default: %(default)s).')

HG38_HLA_BEGIN_BP = 28510120
HG38_HLA_END_BP = 33480577

pysam.set_verbosity(4)

if __name__ == '__main__':
    args = argparser.parse_args()

    if args.in_fasta:
        if not os.path.isfile(args.in_fasta):
            raise Exception(f'Can\'t access {args.in_fasta}')

    selected_query_names = set()

    print(f'Threads specified: {args.n_threads}')
    
    with pysam.AlignmentFile(args.in_bam, mode = "r", reference_filename = args.in_fasta, threads = args.n_threads) as ifile:
        hla_decoy_ids = set()
        for reference in ifile.references:
            if 'HLA-' in reference:
                hla_decoy_ids.add(ifile.get_tid(reference))

        print(f'HLA decoy contigs found: {len(hla_decoy_ids)}')
        
        chrom6_id = ifile.get_tid('6')
        if chrom6_id < 0:
            chrom6_id = ifile.get_tid('chr6')
            if chrom6_id < 0:
                raise Exception(f'Chromosome 6 was not found in {args.in_bam}')
   
        print(f'Looking for candidate segments... ')
        start_time = time.perf_counter()
        for alignment in ifile:
            if alignment.is_unmapped:
                selected_query_names.add(alignment.query_name)
            elif alignment.reference_id == chrom6_id:
                if alignment.reference_start <= HG38_HLA_END_BP and alignment.reference_end >= HG38_HLA_BEGIN_BP:
                    selected_query_names.add(alignment.query_name)
            elif alignment.reference_id in hla_decoy_ids:
                selected_query_names.add(alignment.query_name)
        end_time = time.perf_counter()
        elapsed_time = end_time - start_time
        print(f'Done ({elapsed_time:.3f} seconds)')

        print(f'Unique segments selected: {len(selected_query_names)}')
        
    with pysam.AlignmentFile(args.in_bam, mode = "r", reference_filename = args.in_fasta, threads = args.n_threads) as ifile, pysam.AlignmentFile(args.out_bam, mode = 'wb', template = ifile) as ofile:
        print(f'Exporting reads... ')
        start_time = time.perf_counter()
        for alignment in ifile:
            if alignment.query_name in selected_query_names:
                ofile.write(alignment)
        end_time = time.perf_counter()
        elapsed_time = end_time - start_time
        print(f'Done ({elapsed_time:.3f} seconds)')


