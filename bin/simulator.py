import pandas as pd
import numpy as np
import argparse
import textwrap
from tqdm import tqdm
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from random import randint

import helper

def parse_arg():
    parser = argparse.ArgumentParser(
        description=textwrap.dedent(
        '''
        This script can be used to simulate reads (perfect reads with no errors)
        from 10X single-cell Nanopore RNA sequencing data. 
        '''),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # Required positional argument
    parser.add_argument('-i', '--umi-count-file', type=str, default='../data/SR_bc.csv',
                        help='Filename of the BC UMI count csv file. ')
    parser.add_argument('--amp-rate', type=int, default=None,
                        help=textwrap.dedent('''
                        The number of reads (X) a unique molecule will be 1 if --amp-rate 
                        is not specified. Otherwise, X will be randomly drawn from poisson 
                        distribution with specified rate.
                        '''
                        ))
    parser.add_argument('-o','--output-filename',  type=str, default='template.fa',
                        help=textwrap.dedent('''
                        Filename of simulated reads.
                        '''
                        ))
    
    args = parser.parse_args()

    if args.amp_rate  and args.amp_rate <= 0:
        helper.err_msg(textwrap.dedent(
                f'''
                Invalid number of --amp-rate. It takes positive value only...
                '''))

    # check file 
    helper.check_exist([args.umi_count_file])
    return args


def polyT(length):
    return 'T'*length
    # do I need to randomize the length?

def random_seq(length):
    '''
    Generate random DNA seq with specified length
    '''
    return ''.join(np.random.choice(['A','G','C', 'T'],length))

def artificial_template(BC_list, repeats = None,
                        adaptor='ACTAAAGGCCATTACGGCCTACACGACGCTCTTCCGATCT',
                        TSO='TGTACTCTGCGTTGATACCACTGCTT',
                        output_fn='template.fa',
                        amp_rate = None):
    '''
    Generate a barcode 
        BC_list: list of unique BC
        repeats: list of count each BC is simulated
        adaptor/TSO: adapter/TSO sequence to attach before/after the barcodes, 
                            (use 10X 3' solution V3 adaptor by default)
        output_fn: output filename (append to the end of the file if exist)
        num_read: if specified, each UMI will be repected x time where x~pois(amp_rate)
    output:
        FASTA:
            seqname: barcode sequence
            seq: adapter + BC + 12nt random seq + 15nt polyT + 500 random sequence + TSO
    ''' 
    if repeats is None:
        repeats = [1]*len(BC_list)
        
        
    # v3 adaptor
    with open(output_fn, 'a') as f:
        random_frag = get_random_frag(
                    '/home/ubuntu/vol_data/project/SC_analysis/data/Genome_n_Anno/gencode.v31.transcripts.fa',200)
        read_index = 0
        pbar = tqdm(total=sum(repeats), desc='Total UMI count')
        for BC,count in tqdm(zip(BC_list, repeats), total = len(BC_list), desc='BCs'):
            for i in tqdm(range(count), leave=False, desc = f'Reads for BC: {BC}'):
                seq = adaptor+BC+random_seq(12) + polyT(15) + next(random_frag) +TSO
                
                if amp_rate is None:
                    # seq (randomly reversed) 
                    read_index += 1
                    read_id = BC + '_' +  str(read_index)
                    f.write(">" + read_id + '\n')
                    if np.random.choice([False,True]):
                        f.write(seq + '\n')
                    else:
                        f.write(helper.reverse_complement(seq) + '\n')

                elif amp_rate<=0:
                    print('invalide number of amp_rate, please use amp_rate > 0')
                    sys.exit()
                else:
                    # random repeat each UMI x times (x~poi(amp_rate))
                    for j in range(np.random.poisson(amp_rate)):
                        read_index += 1
                        read_id = BC + '_' +  str(read_index)
                        f.write(">" + read_id + '\n')
                        if np.random.choice([False,True]):
                            f.write(seq + '\n')
                        else:
                            f.write(helper.reverse_complement(seq) + '\n')
                pbar.update(1)
                    
def get_random_frag(ref, frag_len):                    
    # There should be one and only one record, the entire genome:
    records = list(SeqIO.parse(ref, "fasta"))
    num_records = len(records)
    limits = [len(r.seq) for r in records]
    while True:
        i = randint(0,num_records-1)
        record = records[i]
        limit = limits[i]
        if frag_len > limit:
            frag = record.seq
            if 'N' in frag:
                continue
            yield helper.reverse_complement(frag)
        else:
            start = randint(0, limit - frag_len)
            end = start + frag_len
            frag = record.seq[start:end]
            if 'N' in frag:
                continue
            yield helper.reverse_complement(frag)
                            
def main():
    args = parse_arg()
    BC_df = pd.read_csv(args.umi_count_file)
    artificial_template(BC_df.BC, 
                        BC_df.counts, 
                        output_fn=args.output_filename,
                        amp_rate=args.amp_rate)
   
if __name__=='__main__':
    main()
