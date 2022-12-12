import sys
import argparse
import Bio.SeqIO
import textwrap
from tqdm import tqdm
from badread.simulate import sequence_fragment, ErrorModel, QScoreModel, Identities
import os
import multiprocessing as mp
from six.moves import xrange
import random
import numpy as np
import math

BASES = ['A', 'T', 'C', 'G']

import  nanosim_mixed_model as mm
import helper
import config

def parse_arg():
    parser = argparse.ArgumentParser(
        description=textwrap.dedent(
        '''
        Descript come later
        '''),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # Required positional argument
    parser.add_argument('-t', '--template-fasta', type=str,required=True,
                        help='Template fastq file.')
    parser.add_argument('--nanosim-model-prefix', type=str, required=True,
                        help=textwrap.dedent('''
                        Prefix to the nanosim model prefix
                        '''
                        ))
    # Required positional argument    
    parser.add_argument('-o','--output-filename',  type=str, default='sim_read.fq',
                        help=textwrap.dedent('''
                        Filename of simulated reads with errors.
                        '''
                        ))
    parser.add_argument('--batch-size',  type=int, default=500,
                        help=textwrap.dedent('''
                        Batch size
                        '''
                        ))
    parser.add_argument('--thread',  type=int, default = mp.cpu_count()-1,
                        help=textwrap.dedent('''
                        Thread
                        '''
                        ))

    args = parser.parse_args()
    # if args.amp_rate  and args.amp_rate <= 0:
    #     helper.err_msg(textwrap.dedent(
    #             f'''
    #             Invalid number of --amp-rate. It takes positive value only...
    #             '''))
    # check file 
    helper.check_exist([args.template_fasta])
    return args


def sim_read(perfect_read, args):
    """Simulate error into perfect read using Badread error and qscore model

        read (str): perfect reads
    """
    match_ht_list, error_par, trans_error_pr, match_markov_model = \
                        read_nanosim_profile(args.nanosim_model_prefix)
    middle_read, middle_ref, error_dict, error_count = error_list(len(perfect_read), match_markov_model,
                                                            match_ht_list, error_par, trans_error_pr,
                                                            fastq = True)
    print('func:sim_read')
    print(error_dict)
    print('---------------')
    print(error_count)

    seq, quals = mutate_read(perfect_read, error_dict, error_count, 
                                            basecaller = 'guppy', 
                                            read_type = 'linear', 
                                            fastq=True, read_name=None)
    

    return seq, quals

def read_batch_generator(fasta_fns, batch_size):
    """Generator of barches of reads from list of fasta files
    Args:
        fasta_fns (list): fasta filenames
        batch_size (int, optional):  Defaults to 100.
    """
    for fn in fasta_fns:
        if str(fn).endswith('.gz'):
            with gzip.open(fn, "rt") as handle:
                fasta = Bio.SeqIO.parse(handle, "fasta")
                read_batch = helper.batch_iterator(fasta, batch_size=batch_size)
                for batch in read_batch:
                    yield batch
        else:
            fasta = Bio.SeqIO.parse(fn, "fasta")
            read_batch = helper.batch_iterator(fasta, batch_size=batch_size)
            for batch in read_batch:
                yield batch

def sim_read_batch(read_batch, args):
    fastq_lines = []
    for read in read_batch:
        failed = True
        #while failed:
        seq, qscore = sim_read(str(read.seq), args)
        print('sim_read_batch')
        print(len(seq))
        print(len(qscore))
        qscore = ''.join([chr(x+33) for x in qscore])
        failed = False
            # except:
            #     failed = True
            #     pass
        # write read.id, se quals\
        fastq_lines.extend([f'@{read.id}', seq, '+' , qscore])
    return fastq_lines

def read_nanosim_profile(model_prefix):
    """Write error profile from nanosim model
    """
    def read_ecdf(profile):
        # We need to count the number of zeros. If it's over 10 zeros, l_len/l_ratio need to be changed to higher.
        # Because it's almost impossible that the ratio is much lower than the lowest historical value.
        header = profile.readline()
        header_info = header.strip().split()
        ecdf_dict = {}
        lanes = len(header_info[1:])

        for i in header_info[1:]:
            boundaries = i.split('-')
            ecdf_dict[(int(boundaries[0])), int(boundaries[1])] = {}

        ecdf_key = sorted(ecdf_dict.keys())
        l_prob = [0.0] * lanes
        l_ratio = [0.0] * lanes

        for line in profile:
            new = line.strip().split('\t')
            ratio = [float(x) for x in new[0].split('-')]
            prob = [float(x) for x in new[1:]]
            for i in xrange(lanes):
                if prob[i] == l_prob[i]:
                    continue
                else:
                    if l_prob[i] != 0:
                        ecdf_dict[ecdf_key[i]][(l_prob[i], prob[i])] = (l_ratio[i], ratio[1])
                    else:
                        ecdf_dict[ecdf_key[i]][(l_prob[i], prob[i])] \
                            = (max(l_ratio[i], ratio[1] - 10 * (ratio[1] - ratio[0])), ratio[1])
                    l_ratio[i] = ratio[1]
                    l_prob[i] = prob[i]

        for i in xrange(0, len(ecdf_key)):
            last_key = sorted(ecdf_dict[ecdf_key[i]].keys())[-1]
            last_value = ecdf_dict[ecdf_key[i]][last_key]
            ecdf_dict[ecdf_key[i]][last_key] = (last_value[0], ratio[1])

        return ecdf_dict
    error_par = {}
    model_profile = model_prefix + "_model_profile"
    with open(model_profile, 'r') as mod_profile:
        mod_profile.readline()
        for line in mod_profile:
            new_line = line.strip().split("\t")
            if "mismatch" in line:
                error_par["mis"] = [float(x) for x in new_line[1:]]
            elif "insertion" in line:
                error_par["ins"] = [float(x) for x in new_line[1:]]
            else:
                error_par["del"] = [float(x) for x in new_line[1:]]

    trans_error_pr = {}
    with open(model_prefix + "_error_markov_model", "r") as error_markov:
        error_markov.readline()
        for line in error_markov:
            info = line.strip().split()
            k = info[0]
            trans_error_pr[k] = {}
            trans_error_pr[k][(0, float(info[1]))] = "mis"
            trans_error_pr[k][(float(info[1]), float(info[1]) + float(info[2]))] = "ins"
            trans_error_pr[k][(1 - float(info[3]), 1)] = "del"

    with open(model_prefix + "_first_match.hist", 'r') as fm_profile:
        match_ht_list = read_ecdf(fm_profile)

    with open(model_prefix + "_match_markov_model", 'r') as mm_profile:
        match_markov_model = read_ecdf(mm_profile)

    return match_ht_list, error_par, trans_error_pr, match_markov_model

def error_list(m_ref, m_model, m_ht_list, error_p, trans_p, fastq):
    # l_old is the original length, and l_new is used to control the new length after introducing errors
    l_new = m_ref
    pos = 0
    e_dict = {}
    middle_ref = m_ref
    prev_error = "start"
    e_count = {"mis": 0, "ins": 0, "match": 0}

    # The first match come from m_ht_list
    p = random.random()
    k1 = list(m_ht_list.keys())[0]
    for k2, v2 in m_ht_list[k1].items():
        if k2[0] < p <= k2[1]:
            prev_match = int(np.floor((p - k2[0]) / (k2[1] - k2[0]) * (v2[1] - v2[0]) + v2[0]))
            if prev_match < 2:
                prev_match = 2
    pos += prev_match
    if fastq:
        if prev_match > middle_ref: 
            e_count["match"] += middle_ref
        else:
            e_count["match"] += prev_match

    # Select an error, then the step size, and then a match and so on so forth.
    while pos < middle_ref:
        # pick the error based on Markov chain
        p = random.random()
        for k in trans_p[prev_error].keys():
            if k[0] <= p < k[1]:
                error = trans_p[prev_error][k]
                break

        if error == "mis":
            step = mm.pois_geom(error_p[error][0], error_p[error][2], error_p[error][3])
        elif error == "ins":
            step = mm.wei_geom(error_p[error][0], error_p[error][1], error_p[error][2], error_p[error][3])
            l_new += step
        else:
            step = mm.wei_geom(error_p[error][0], error_p[error][1], error_p[error][2], error_p[error][3])
            l_new -= step

        if error != "ins":
            e_dict[pos] = [error, step]
            pos += step
            if pos >= middle_ref:
                l_new += pos - middle_ref
                middle_ref = pos
        else:
            e_dict[pos - 0.5] = [error, step]

        prev_error = error

        if fastq:
            if error == "mis" or error == "ins":
                e_count[error] += step

        # Randomly select a match length
        for k1 in m_model.keys():
            if k1[0] <= prev_match < k1[1]:
                break
        p = random.random()
        for k2, v2 in m_model[k1].items():
            if k2[0] < p <= k2[1]:
                step = int(np.floor((p - k2[0]) / (k2[1] - k2[0]) * (v2[1] - v2[0]) + v2[0]))
                break
        # there are no two 0 base matches together
        if prev_match == 0 and step == 0:
            step = 1

        prev_match = step

        if fastq:
            e_count["match"] += step

        if pos + prev_match > middle_ref:
            l_new += pos + prev_match - middle_ref
            middle_ref = pos + prev_match

        pos += prev_match
        if prev_match == 0:
            prev_error += "0"

    return l_new, middle_ref, e_dict, e_count

def mutate_read(read, e_dict, e_count, basecaller = 'guppy', read_type = 'linear', 
                        fastq=True, error_log = None, read_name=None):
    new_e_dict = e_dict

    if fastq:  # Sample base qualities for mis/ins/match
        mis_quals = mm.trunc_lognorm_rvs("mis", read_type, basecaller, e_count["mis"]).tolist()
        ins_quals = mm.trunc_lognorm_rvs("ins", read_type, basecaller, e_count["ins"]).tolist()
        match_quals = mm.trunc_lognorm_rvs("match", read_type, basecaller, e_count["match"]).tolist()

    # Mutate read
    quals = []
    prev = len(read)
    for key in sorted(new_e_dict.keys(), reverse=True):
        val = new_e_dict[key]
        key = math.ceil(key)  # Ceil instead of round for consistent match calculations during base qual sim
        err_quals = []

        if val[0] == "mis":
            ref_base = read[key: key + val[1]]
            new_bases = ""
            for i in xrange(val[1]):
                tmp_bases = list(BASES)
                tmp_bases.remove(read[key + i])
                # tmp_bases.remove(read[key]) ## Edited this part for testing
                new_base = random.choice(tmp_bases)
                new_bases += new_base
                if fastq:
                    err_quals.append(mis_quals.pop())

            new_read = read[:key] + new_bases + read[key + val[1]:]
            err_end = key + val[1]

        elif val[0] == "del":
            new_bases = val[1] * "-"
            ref_base = read[key: key + val[1]]
            new_read = read[: key] + read[key + val[1]:]
            err_end = key + val[1]

        elif val[0] == "ins":
            ref_base = val[1] * "-"
            new_bases = ""
            for i in xrange(val[1]):
                new_base = random.choice(BASES)
                new_bases += new_base
                if fastq:
                    err_quals.append(ins_quals.pop())
            new_read = read[:key] + new_bases + read[key:]
            err_end = key

        if fastq:
            if err_end != prev:  # Match after error
                for j in xrange(prev - err_end):
                    quals.append(match_quals.pop())
            quals += err_quals

        read = new_read
        prev = key

        # if val[0] != "match" and error_log:
        #     error_log.write(read_name + "\t" + str(key) + "\t" + val[0] + "\t" + str(val[1]) +
        #                     "\t" + ref_base + "\t" + new_bases + "\n")

    if fastq:  # Add first match quals
        while len(match_quals) > 0:
            quals.append(match_quals.pop())

    quals.reverse()
    return read, quals

def main(output=sys.stderr):
    args = parse_arg()

    # input template file
    if os.path.isdir(args.template_fasta):
        fns = helper.get_files(fastq_dir, ['*.fastq', '*.fq', '*.fastq.gz', '*.fq.gz',
                                    '*.fasta', '*.fa', '*.fasta.gz', '*.fa.gz'])
    elif os.path.isfile(args.template_fasta):
        fns = [args.template_fasta]
    else:
        helper.err_msg("Invalid input of template fasta file")
        sys.exit(1)

    read_batchs = read_batch_generator(fns, args.batch_size)
    

    # single threading mode (good for debuging)
    if args.thread == 1:
        pbar = tqdm(unit = 'read', desc='Processed')
        for batch in read_batchs:
            fastq_records = sim_read_batch(batch, args)
            
            with open(args.output_filename, 'w') as f:
                f.write('\n'.join(fastq_records) + '\n')
            pbar.update(args.batch_size)
        
        helper.green_msg(f'Simulation finished!')
        return None
    else:
        # multiprocessing mode
        pass
    
    rst_futures = helper.multiprocessing_submit(sim_read_batch,
                    read_batchs, n_process=args.thread, 
                    pbar_update=args.batch_size,
                    args=args)
    
    for idx, f in enumerate(rst_futures):
        fastq_records = f.result() #write rst_df out
        if idx == 0:
            with open(args.output_filename, 'w') as f:
                f.write('\n'.join(fastq_records)+'\n')
        else:
            with open(args.output_filename, 'a') as f:
                f.write('\n'.join(fastq_records)+'\n')

    helper.green_msg(f'Simulation finished!')
    return None    


if __name__ == '__main__':
    main()


# test chunk
# sys.argv = ['sim_error.py', '-t', '../test/test_50_reads.fa']
# args = parse_arg()
# fasta = Bio.SeqIO.parse(args.template_fasta, "fasta")
