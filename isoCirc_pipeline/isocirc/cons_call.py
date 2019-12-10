import argparse
import sys
import os
import mappy as mp
import itertools

import isocirc.utils as ut

cons_header_ele = ['repeat_start', 'repeat_end', 'period_size', 'copy_number', 'consensus_length', 'match_frac',
                   'indel_frac', 'align_score', 'A', 'C', 'G', 'T', 'entropy', 'cons_seq', 'repeat_seq',
                   'left_flank_seq', 'right_flank_seq']
cons_header_idx = {cons_header_ele[i]: i for i in range(len(cons_header_ele))}

trf = 'trf409.legacylinux64'
tidehunter = 'TideHunter'

fxtools = 'fxtools' #dir_path + '/bin/fxtools'
cons_min_len = 30
cons_min_copy = 2.0
cons_min_frac = 0.0
threads = 8

# TRF parameters
match = 2
mismatch = 7
indel = 7
match_frac = 80
indel_frac = 10
min_score = 100
max_period = 2000

# TideHunter parameters
w=1
k=8
s=1
l='' #'-l'

def get_relation(seq1, seq2):
    len1, len2 = len(seq1), len(seq2)
    if 0.9 * len1 > len2 or 1.1 * len1 < len2: return '', 0, 0
    res = ''
    min_len = min(len1, len2)
    # seq1 = seq1 + seq1
    # seq2 = seq2 + seq2
    a = mp.Aligner(seq=seq1)
    f_iden = r_iden = 0
    for h in a.map(seq2):
        # only use primary alignment
        if not h.is_primary: continue
        f_iden = h.mlen / (min_len+ 0.0)
        f_strand = h.strand
        break
    rseq2 = seq2[::-1]
    for h in a.map(rseq2):
        # only use primary alignment
        if not h.is_primary: continue
        r_iden = h.mlen / (min_len + 0.0)
        r_strand = h.strand
        break
    if max(f_iden, r_iden) < 0.8: return 'NA', f_iden, r_iden
    if f_iden > r_iden:
        if f_strand == 1:
            res = 'ID'  # 'identical'
        else:
            res = 'RC'
    else:
        if r_strand == 1:
            res = 'R'
        else:
            res = 'C'
    return res

def get_read_len(read_fn, out_dir, fxtools=fxtools):
    read_len = {}
    read_len_fn = out_dir + '/' + os.path.basename(os.path.abspath(read_fn)) + '.len'
    if not os.path.exists(read_len_fn):
        ut.exec_cmd(sys.stderr, 'fxtools', '{} lp {} > {} 2> /dev/null'.format(fxtools, read_fn, read_len_fn))
    with open(read_len_fn, 'r') as len_fp:
        for line in len_fp:
            read_len[line.rsplit()[0]] = int(line.rsplit()[1])
    return read_len

# regular chimeric, inverted, identical(non-chimeric)
def get_chimeric_info(seqs):
    relations = []
    for i in range(len(seqs)):
        relations.append(str(i))
    for i, j in itertools.combinations(range(len(seqs)), 2):
        s1, s2 = seqs[i], seqs[j]
        res = get_relation(s1, s2)
        if res == 'RC':
            relations[i] += 'RC{}'.format(j)
            relations[j] += 'RC{}'.format(i)
        elif res == 'ID':
            relations[i] += 'ID{}'.format(j)
            relations[j] += 'ID{}'.format(i)
    return relations

def overlap(s1, e1, s2, e2):
    if e1 < s2 or e2 < s1:
        return 0.0
    else:
        return min(e2-s1+1, e1-s2+1)

def is_embedded(c, l, s, e, cons_seqs, span, lens):
    for i in range(len(cons_seqs)):
        ovlp_len = overlap(s, e, span[i][0], span[i][1])
        if ovlp_len / (e-s+1.0) >= 0.8 and lens[i] < l:
            return True # igore l,s,e
        elif ovlp_len / (span[i][1]-span[i][0]+1.0) >= 0.8 and lens[i] > l:
            cons_seqs[i] = c
            lens[i] = l
            span[i] = (s, e)
            return True
    return False


def extract_tidehunter_cons(th_out, cons_fa, cons_info, cons_min_len=cons_min_len, cons_min_copy=cons_min_copy, cons_min_frac=cons_min_frac):
    last_name, last_len, read_name, cons_seqs, cons_lens, copy_nums, spans, fracs = '', 0, '', [], [], [], [], []
    with open(th_out, 'r') as in_fp, open(cons_info, 'w') as info_fp, open(cons_fa, 'w') as cons_fp:
        cons_i = 0
        for line in in_fp:
            # print line
            [read_name, cons_id, read_len, cons_start, cons_end, cons_len, copy_num, is_full_len, cons_seq] = line.rsplit()
            read_len, cons_start, cons_end, cons_len, copy_num = int(read_len), int(cons_start), int(cons_end), int(cons_len), float(copy_num)
            frac = (cons_end-cons_start+1.0) / read_len
            if cons_len < cons_min_len or frac < cons_min_frac or copy_num < cons_min_copy:
                continue
            if read_name == last_name and cons_i > 0 and is_embedded(cons_seq, cons_len, cons_start, cons_end, cons_seqs, spans, cons_lens): continue
            if read_name != last_name:
                if cons_i > 0:
                    chimeric_info = ['NA'] if cons_i == 1 else get_chimeric_info(cons_seqs)
                    for i in range(len(cons_lens)):
                        c, cl, copy_n, fr = cons_seqs[i], cons_lens[i], copy_nums[i], fracs[i]
                        # print c, cl, copy_n, sp, chimeric_info[i]
                        info_fp.write('{}_cons{}\t{}\t{}\t{}\t{}\t{}\n'.format(last_name, i, last_len, cl, copy_n, fr, chimeric_info[i]))
                        cons_fp.write('>{}_cons{}\n{}{}\n'.format(last_name, i, c, c))  # 2 copies of consensus sequence
                last_name, last_len, cons_i, cons_seqs, cons_lens, copy_nums, spans, fracs = read_name, read_len, 1, [cons_seq], [cons_len], [copy_num], [(cons_start, cons_end)], [frac]
            else:
                cons_i += 1
                cons_seqs.append(cons_seq)
                cons_lens.append(cons_len)
                copy_nums.append(copy_num)
                spans.append((cons_start, cons_end))
                fracs.append(frac)
        if cons_i > 0:
            chimeric_info = ['NA'] if cons_i == 1 else get_chimeric_info(cons_seqs)
            for i in range(len(cons_lens)):
                c, cl, copy_n, fr = cons_seqs[i], cons_lens[i], copy_nums[i], fracs[i]
                # print c, cl, copy_n, sp, chimeric_info[i]
                info_fp.write('{}_cons{}\t{}\t{}\t{}\t{}\t{}\n'.format(last_name, i, last_len, cl, copy_n, fr, chimeric_info[i]))
                cons_fp.write('>{}_cons{}\n{}{}\n'.format(last_name, i, c, c))  # 2 copies of consensus sequence



def extract_multi_cons(cons_fa, cons_info, trf_out, read_len, cons_min_len=cons_min_len, cons_min_copy=cons_min_copy, cons_min_frac=cons_min_frac):
    read_name, cons_seqs, cons_lens, copy_number, span = '', [], [], [], []
    with open(trf_out, 'r') as trf_fp, open(cons_fa, 'w') as cons_fp, open(cons_info, 'w') as info_fp:
        cons_i = 0
        for line in trf_fp:
            ele = line.rsplit()
            if line.startswith('@'):
                if cons_i > 0:
                    chimeric_info = ['NA'] if cons_i == 1 else get_chimeric_info(cons_seqs)
                    for i in range(len(cons_lens)):
                        c, cl, copy_n, sp = cons_seqs[i], cons_lens[i], copy_number[i], span[i]
                        # print c, cl, copy_n, sp, chimeric_info[i]
                        info_fp.write('{}_cons{}\t{}\t{}\t{}\t{}\t{}\n'.format(read_name, i, read_len[read_name], cl, copy_n, (sp[1]-sp[0]+1.0)/read_len[read_name], chimeric_info[i]))
                        cons_fp.write('>{}_cons{}\n{}{}\n'.format(read_name, i, c, c))  # 2 copies of consensus sequence
                read_name = line[1:].rsplit()[0]
                cons_i, cons_seqs, cons_lens, copy_number, span = 0, [], [], [], []
            else:  # repeat
                cons_len = int(ele[cons_header_idx['consensus_length']])
                if cons_len < cons_min_len: continue
                cons_end = int(ele[cons_header_idx['repeat_end']])
                cons_start = int(ele[cons_header_idx['repeat_start']])
                if (cons_end - cons_start + 1.0) / read_len[read_name] < cons_min_frac: continue
                copy_n = float(ele[cons_header_idx['copy_number']])
                if copy_n < cons_min_copy: continue
                cons_seq = ele[cons_header_idx['cons_seq']]
                if is_embedded(cons_seq, cons_len, cons_start, cons_end, cons_seqs, span, cons_lens): continue
                cons_seqs.append(cons_seq), cons_lens.append(cons_len), copy_number.append(copy_n), span.append((cons_start, cons_end))
                cons_i += 1

        if cons_i > 0:
            chimeric_info = ['NA'] if cons_i == 1 else get_chimeric_info(cons_seqs)
            for i in range(len(cons_lens)):
                c, cl, copy_n, sp = cons_seqs[i], cons_lens[i], copy_number[i], span[i]
                # print c, cl, copy_n, sp, chimeric_info[i]
                info_fp.write('{}_cons{}\t{}\t{}\t{}\t{}\t{}\n'.format(read_name, i, read_len[read_name], cl, copy_n, (sp[1]-sp[0]+1.0)/read_len[read_name], chimeric_info[i]))
                cons_fp.write('>{}_cons{}\n{}{}\n'.format(read_name, i, c, c))  # 2 copies of consensus sequence

def extract_cons(cons_fa, cons_info, trf_out, read_len, cons_min_len=cons_min_len, cons_min_copy=cons_min_copy, cons_min_frac=cons_min_frac):
    read_name, cons_seq, max_copy_number = '', '', 0
    all_cons_len = 0
    with open(trf_out, 'r') as trf_fp, open(cons_fa, 'w') as cons_fp, open(cons_info, 'w') as info_fp:
        for line in trf_fp:
            ele = line.rsplit()
            if line.startswith('@'):
                if read_name and min_cons_len >= cons_min_len and min_cons_len != 10000:
                    # 2 copies of consensus sequence
                    info_fp.write('{}\t{}\t{}\t{}\t{}\n'.format(read_name, read_len[read_name], len(cons_seq), max_copy_number, all_cons_len / read_len[read_name]))
                    cons_fp.write('>{}\n{}{}\n'.format(read_name, cons_seq, cons_seq))
                    # cons_fp.write('>{}\n{}\n'.format(read_name, cons_seq))
                read_name = line[1:].rsplit()[0]
                min_cons_len = 10000
                cons_seq = ''
            else:  # repeat
                if int(ele[cons_header_idx['consensus_length']]) < min_cons_len \
                        and int(ele[cons_header_idx['consensus_length']]) >= cons_min_len \
                        and read_len[read_name] * cons_min_frac <= (1.0 + float(
                    ele[cons_header_idx['repeat_end']]) - float(ele[cons_header_idx['repeat_start']])) \
                        and float(ele[cons_header_idx['copy_number']]) >= cons_min_copy:
                    all_cons_len = (1.0 + float(ele[cons_header_idx['repeat_end']]) - float(ele[cons_header_idx['repeat_start']]))
                    max_copy_number = float(ele[cons_header_idx['copy_number']])
                    cons_seq = ele[cons_header_idx['cons_seq']]
                    min_cons_len = int(ele[cons_header_idx['consensus_length']])
        if read_name and min_cons_len >= cons_min_len and min_cons_len != 10000:
            info_fp.write('{}\t{}\t{}\t{}\t{}\n'.format(read_name, read_len[read_name], len(cons_seq), max_copy_number, all_cons_len / read_len[read_name]))
            cons_fp.write('>{}\n{}{}\n'.format(read_name, cons_seq, cons_seq))


def run_tidehunter(in_long, cons_fa, cons_info, th_out, tidehunter=tidehunter, cons_min_len=cons_min_len, cons_min_copy=cons_min_copy, cons_min_frac=cons_min_frac, threads=threads):
    ut.exec_cmd(sys.stderr, 'TideHunter', '{} {} -t {} -f 2 -w {} -k {} -s {} {} -p {} -c {} -o {}'.format(tidehunter, in_long, threads, w, k, s, l, cons_min_len, cons_min_copy, th_out))
    extract_tidehunter_cons(th_out, cons_fa, cons_info, cons_min_len, cons_min_copy, cons_min_frac)

def run_trf(in_long, cons_fa, cons_info, trf_out, trf=trf, match=match, mismatch=mismatch, indel=indel, match_frac=match_frac, indel_frac=indel_frac, min_score=min_score, max_period=max_period, cons_min_len=cons_min_len, cons_min_copy=cons_min_copy, cons_min_frac=cons_min_frac):
    ut.exec_cmd(sys.stderr, 'Tandem Repeats Finder', '{} {} {} {} {} {} {} {} {} -h -ngs > {}'.format(trf, in_long, match, mismatch, indel, match_frac,
                                                     indel_frac, min_score, max_period, trf_out))
    read_len = get_read_len(in_long, os.path.dirname(os.path.abspath(cons_fa)))
    # extract_cons(cons_fa, cons_info, trf_out, read_len, cons_min_len, cons_min_copy, cons_min_frac)
    extract_multi_cons(cons_fa, cons_info, trf_out, read_len, cons_min_len, cons_min_copy, cons_min_frac)

    return

def run_trf_parall(in_long, cons_fa, cons_info, trf_out, trf=trf, match=match, mismatch=mismatch, indel=indel, match_frac=match_frac, indel_frac=indel_frac, min_score=min_score, max_period=max_period, cons_min_len=cons_min_len, cons_min_copy=cons_min_copy, cons_min_frac=cons_min_frac, threads=threads):
    out_dir = os.path.dirname(os.path.abspath(trf_out)) + '/'
    ut.exec_cmd(sys.stderr, 'fxtools', '{} sx {} {} {}'.format(fxtools, in_long, threads, out_dir))
    fp_list = [sys.stderr] * threads
    header_list = ['Tandem Repeats Finder'] * threads
    tmp_in_long = [out_dir+os.path.basename(in_long) + '.' + str(i) for i in range(1, threads + 1)]
    tmp_trf_out = [trf_out + '.' + str(i) for i in range(1, threads + 1)]

    cmd_list = ['{} {} {} {} {} {} {} {} {} -h -ngs > {}; rm {}'.format(trf, in_tmp_long, match, mismatch, indel, match_frac, indel_frac, min_score, max_period, trf_tmp_out, in_tmp_long)
                for in_tmp_long, trf_tmp_out in zip(tmp_in_long, tmp_trf_out)]
    ut.exec_cmd_parall(fp_list, header_list, cmd_list)
    if os.path.exists(trf_out):
        ut.exec_cmd(sys.stderr, 'Tandem Repeats Finder', 'rm {}'.format(trf_out))
    for tmp in tmp_trf_out:
        ut.exec_cmd(sys.stderr, 'Tandem Repeats Finder', 'cat {} >> {}; rm {}'.format(tmp, trf_out, tmp))
    read_len = get_read_len(in_long, os.path.dirname(os.path.abspath(cons_fa)))
    # extract_cons(cons_fa, cons_info, trf_out, read_len, cons_min_len, cons_min_copy, cons_min_frac)
    extract_multi_cons(cons_fa, cons_info, trf_out, read_len, cons_min_len, cons_min_copy, cons_min_frac)


def run_rf_core(args):
    # if args.use_tidehunter:
        # run_tidehunter(args.in_long_read, args.cons_fa, args.cons_info, args.tr_out, args.tidehunter, args.min_len, args.min_copy, args.min_frac, args.threads)
    # else:
    if args.threads > 1:
        run_trf_parall(args.in_long_read, args.cons_fa, args.cons_info, args.tr_out, args.trf, args.match, args.mismatch, args.indel, args.match_frac, args.indel_frac, args.min_score, args.max_period, args.min_len, args.min_copy, args.min_frac, args.threads)
    else:
        run_trf(args.in_long_read, args.cons_fa, args.cons_info, args.tr_out, args.trf, args.match, args.mismatch, args.indel, args.match_frac, args.indel_frac, args.min_score, args.max_period, args.min_len, args.min_copy, args.min_frac)


def parser_argv():
    # parse command line arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="Consensus calling with Tandem Repeats Finder(TRF)")
    parser.add_argument("in_long_read", metavar='long.fa', type=str, help='Long read sequencing data generated with isoCirc.')
    parser.add_argument("tr_out", metavar='tr.out', type=str, help='Output file.')
    parser.add_argument("cons_fa", metavar='cons.fa', type=str, help='Consensus sequence file.')
    parser.add_argument("cons_info", metavar='cons.info', type=str, help='Consensus information file.')

    # parser.add_argument('-T', '--use-tidehunter', default=False, help='Use TideHunter as the tandem repeats detection tool.', action='store_true')

    parser.add_argument('--trf', type=str, help='Path to trf program.', default=trf)
    # parser.add_argument('--tidehunter', type=str, help='Path to TideHunter.', default=tidehunter)

    parser.add_argument('--threads', type=int, help='Number of thread to use.', default=threads)
    parser.add_argument('-m', '--match', type=int, help='Match score.', default=match)
    parser.add_argument('-x', '--mismatch', type=int, help='Mismatch penalty.', default=mismatch)
    parser.add_argument('-i', '--indel', type=int, help='Indel penalty.', default=indel)
    parser.add_argument('-M', '--match-frac', type=int, help='Match probability.', default=match_frac)
    parser.add_argument('-I', '--indel-frac', type=int, help='Indel probability.', default=indel_frac)
    parser.add_argument('-s', '--min-score', type=int, help='Minimum alignment score to report.', default=min_score)
    parser.add_argument('-p', '--max-period', type=int, help='Maximum period size to report.', default=max_period)
    parser.add_argument('-l', '--min-len', type=int, help='Minimum consensus length to keep.', default=cons_min_len)
    parser.add_argument('-c', '--min-copy', type=float, help='Minimum copy number of consensus to keep.', default=cons_min_copy)
    parser.add_argument('-f', '--min-frac', type=float, help='Minimum fraction of original long read to keep.', default=cons_min_frac)


    return parser.parse_args()


if __name__ == '__main__':
    args = parser_argv()
    run_rf_core(args)
    # read_len = get_read_len(sys.argv[4], './')
    # extract_multi_cons(sys.argv[1], sys.argv[2], sys.argv[3], read_len)