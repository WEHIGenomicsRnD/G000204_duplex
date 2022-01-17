'''
Script for generating stats on duplex reads.
Counts the number of families, their min, mean,
media, and max, as well as the proportion of
families with reads > 1, the number of paired
families (+ and - strand) and number of paired
families where each family has reads > 1.
'''
import pysam
import numpy as np
import pandas as pd
import seaborn as sns
import os
import sys
from argparse import ArgumentParser
from glob import iglob

def parse_args(args):
    description = 'Generate stats on reads from duplex sequencing data.'
    parser = ArgumentParser(description=description)
    parser.add_argument('bam_path',
                        help='''BAM path. May contain wildcard to match multiple files,
                                but in this case you MUST enclose the argument in quotes.''')

    return parser.parse_args(args)

def get_read_info(bam):
    rinfo = []
    for read in bam:
        if read.pos < 0:
            # read not mapped
            continue
        if read.mate_is_unmapped:
            continue
            
        strand = ''
        if read.is_read1 and read.pos < read.next_reference_start:
            # read 1 is left-most read
            if not read.is_reverse and read.mate_is_reverse:
                # check that orentation is proper
                strand = '+'
        elif read.is_read2 and read.pos > read.next_reference_start:
            # read 2 is right-most read
            if read.is_reverse and not read.mate_is_reverse:
                # same as above, we've just selected read 2
                strand = '+'
        elif read.is_read1 and read.pos > read.next_reference_start:
            # read 1 is right-most read
            if read.is_reverse and not read.mate_is_reverse:
                strand = '-'
        elif read.is_read2 and read.pos < read.next_reference_start:
            # read 2 is left-most read
            if not read.is_reverse and read.mate_is_reverse:
                # same as above, we've just selected read 2
                strand = '-'

        if strand not in ['-', '+']:
            continue
        
        pos = min(read.pos, read.next_reference_start)
        umi = read.get_tag('ZB') + '-' + read.get_tag('ZA') if strand == '-' else read.get_tag('RX')
        r = (read.query_name,
             read.reference_name,
             pos,
             strand,
             umi)
        rinfo.append(r)
        if len(rinfo) % 1000000 == 0:
            print('Processed %d reads' % len(rinfo))
    rinfo = pd.DataFrame(rinfo, columns=['rname', 'chrom', 'pos', 'strand', 'umi']).drop_duplicates()
    return rinfo

def process_reads(bam_path):
    stats = []
    for bam_file in iglob(bam_path):
        bam = pysam.AlignmentFile(bam_file)
        
        sample = os.path.basename(bam_file)
        sample = '_'.join(sample.split('_')[:2])
        
        print('Processing', sample, '...')
        
        rinfo = get_read_info(bam)
        rinfo_summary = rinfo.groupby(['chrom', 'pos', 'strand', 'umi']).size().reset_index()
        #rinfo_summary.to_csv('%s_rinfo.txt' % sample, sep='\t', index=False)    

        family_size = rinfo_summary[0].values    
        by_umi = rinfo_summary.groupby(['pos', 'umi']).size().reset_index()
        by_umi_gt1 = rinfo_summary[family_size > 1].groupby(['pos', 'umi']).size().reset_index()

        stats.append({'sample': sample,
                      'len': len(family_size),
                      'mean': np.mean(family_size),
                      'median': np.median(family_size),
                      'max': max(family_size),
                      'min': min(family_size),
                      'frac_fam_gt1': sum(rinfo_summary[0] > 1) / len(rinfo_summary),
                      'frac_fam_paired': sum(by_umi[0] > 1) / len(by_umi),
                      'frac_fam_paired_gt1': sum(by_umi_gt1[0] > 1) / len(by_umi)
                    }
        )
        return pd.DataFrame(stats)

def main():
    args = parse_args(sys.argv[1:])
    bam_path = args.bam_path

    stats = process_reads(bam_path)
    stats.to_csv('family_sizes.txt', sep='\t', index=False)

if __name__ == '__main__':
    main()
