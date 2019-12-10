#!/usr/bin/env python
import argparse
from collections import defaultdict as dd
import itertools

import isocirc.utils as ut
import isocirc.hcBSJ_fullIso as hf
from isocirc.__init__ import __program__
from isocirc.__init__ import __version__
from isocirc.__init__ import whole_output_header
from isocirc.__init__ import whole_output_header_idx
from isocirc.__init__ import isoform_output_header
from isocirc.__init__ import isoform_output_header_idx

whole_output_header = hf.whole_output_header
whole_output_header_idx = hf.whole_output_header_idx
isoform_output_header = hf.isoform_output_header
isoform_output_header_idx = hf.isoform_output_header_idx
'''
Compare PARRIS results between different long-read datasets, or long-read and short-read datasets
'''

ciri_header = ['circRNA_ID','chr','circRNA_start','circRNA_end', '#junction_reads', 'SM_MS_SMS','#non_junction_reads',
               'junction_reads_ratio', 'circRNA_type', 'gene_id', 'strand', 'junction_reads_ID']
ciri_header_idx = {h: i for i, h in enumerate(ciri_header)}

# whole_output_header = ['#readID', 'chrom', 'startCoor0based', 'endCoor', 'mapStrand', 'geneStrand', 'geneID', 'geneName', # 'transID', 'transName',
#                        'blockCount', 'blockSize', 'blockStarts', 'refMapLen',  # mapping information
#                        'blockType', 'blockAnno',  # #evaluation with whole annotation for each block
#                        'readLen', 'consLen', 'consMapLen', 'copyNum', 'consFrac',  # original read and consensus sequence information
#                        'novelFlag', 'circID', 'circLen',  # annotated circRNA information
#                        'isKnownBSJ', 'isCanoBSJ', 'disToCanoBSJ', 'canoBSJMotif', 'alignAroundCanoBSJ', # back-splicing junction
#                        'isKnownCircSS', 'isKnownCircSJ', 'isKnownCircExon', 'disToKnownCircSS', # evaluation with circRNA annotation
#                        'isKnownSS', 'isKnownSJ', 'isKnownExon', 'disToKnownSS',
#                        'isCanoSJ', 'disToCanoSJ', 'canoSJStrand', # evaluation with whole annotation
#                        'CDS', 'UTR', 'lincRNA', 'antisense', # $3; gene_type/biotype
#                        'rRNA', 'Alu', 'allRepeat', 'upFlankAlu', 'downFlankAlu']  # repeat element
# whole_output_header_idx = {h: i for i, h in enumerate(whole_output_header)}

# isoform_output_header = ['#isoformID'] + whole_output_header[1:] + ['readCount', 'readIDs']
# isoform_output_header_idx = {h: i for i, h in enumerate(isoform_output_header)}


def print_ovlp(out_fp, ovlp_entry, line1, line2):
    if ovlp_entry == 'A':
        out_fp.write(line1)
    elif ovlp_entry == 'B':
        out_fp.write(line2)
    elif ovlp_entry == 'AB':
        out_fp.write(line1[:-1]+'\t'+line2)


def is_bsj_in_dict(coor_tuple, d, bsj_dis):
    S_range = [0]
    for i in range(1, bsj_dis + 1):
        S_range.append(i)
        S_range.append(-i)

    ret = []
    for s1 in S_range:
        for s2 in S_range:
            if (coor_tuple[0], coor_tuple[1]+s1, coor_tuple[2]+s2) in d:
                ret.append((coor_tuple[0], coor_tuple[1]+s1, coor_tuple[2]+s2))
    return ret


def get_overlap(a_dict, b_dict, bsj_dis, ovlp_entry, ovlp_fn):
    with open(ovlp_fn, 'w') as out_fp:
        for a in a_dict:
            bs = is_bsj_in_dict(a, b_dict, bsj_dis)
            if bs:
                if ovlp_entry == 'A':
                    out_fp.write(''.join(a_dict[a]))
                elif ovlp_entry == 'B':
                    out_fp.write(''.join(b_dict[a]))
                elif ovlp_entry == 'AB':
                    for b in bs:
                        for r in itertools.product(a_dict[a], b_dict[b]):
                            out_fp.write(r[0][:-1]+'\t'+r[1])


def is_iden_coor(coor1, coor2, site_dis):
    if len(coor1) != len(coor2): return False
    for c1, c2 in zip(coor1[1:-1], coor2[1:-1]):
        if abs(c1-c2) > site_dis: return False
    return True


def get_detailed_overlap(a_dict, b_dict, bsj_dis, site_dis, ovlp_entry, ovlp_fn):
    with open(ovlp_fn, 'w') as out_fp:
        for a in a_dict:
            bsjs = is_bsj_in_dict(a, b_dict, bsj_dis)
            for bsj in bsjs:
                for r in itertools.product(a_dict[a], b_dict[bsj]):
                    if is_iden_coor(r[0], r[1], site_dis):
                        print_ovlp(out_fp, ovlp_entry, a_dict[a][r[0]], b_dict[bsj][r[1]])


def get_only(a_dict, b_dict, bsj_dis, only_fn):
    with open(only_fn, 'w') as out_fp:
        for a in a_dict:
            if not is_bsj_in_dict(a, b_dict, bsj_dis):
                out_fp.write(''.join(a_dict[a]))

def get_short_ciri_input(in_fn):
    short_dict = dd(lambda : []) # coordinate : BSJ read count
    with open(in_fn) as in_fp:
        for line in in_fp:
            if line.startswith('#') or line.startswith('circRNA_ID'): continue
            ele = line.split()
            short_dict[(ele[ciri_header_idx['chr']], int(ele[ciri_header_idx['circRNA_start']]), int(ele[ciri_header_idx['circRNA_end']]))].append(line) # int(ele[ciri_header_idx['#junction_reads']])
    return short_dict

def get_short_bed_input(in_fn):
    short_dict = dd(lambda : [])
    with open(in_fn) as in_fp:
        for line in in_fp:
            if line.startswith('#'): continue
            ele = line.split()
            short_dict[(ele[0], int(ele[1]), int(ele[2]))].append(line) # int(ele[3])
    return short_dict

def get_detailed_long_whole_input(in_fn):
    long_dict = dd(lambda : dd(lambda: ''))
    with open(in_fn) as in_fp:
        for line in in_fp:
            if line.startswith('#'): continue
            ele = line.split()
            chrom = ele[whole_output_header_idx['chrom']]
            start = int(ele[whole_output_header_idx['startCoor0based']])
            end = int(ele[whole_output_header_idx['endCoor']])
            coors = []
            for s, l in zip(map(int,ele[whole_output_header_idx['blockStarts']].rsplit(',')), map(int,ele[whole_output_header_idx['blockSize']].rsplit(','))):
                coors.append(start + s)
                coors.append(start + s + l)
            long_dict[(chrom, start, end)][tuple(coors)] = line
    return long_dict

def get_long_whole_input(in_fn):
    long_dict = dd(lambda : [])
    with open(in_fn) as in_fp:
        for line in in_fp:
            if line.startswith('#'): continue
            ele = line.split()
            chrom = ele[whole_output_header_idx['chrom']]
            start = int(ele[whole_output_header_idx['startCoor0based']])
            end = int(ele[whole_output_header_idx['endCoor']])
            long_dict[(chrom, start, end)].append(line)
    return long_dict

def get_detailed_long_isoform_input(in_fn):
    long_dict = dd(lambda: dd(lambda: ''))
    with open(in_fn) as in_fp:
        for line in in_fp:
            if line.startswith('#'): continue
            ele = line.split()
            chrom = ele[whole_output_header_idx['chrom']]
            start = int(ele[whole_output_header_idx['startCoor0based']])
            end = int(ele[whole_output_header_idx['endCoor']])
            coors = []
            for s, l in zip(map(int,ele[isoform_output_header_idx['blockStarts']].rsplit(',')), map(int,ele[isoform_output_header_idx['blockSize']].rsplit(','))):
                coors.append(start + s)
                coors.append(start + s + l)
            long_dict[(chrom, start, end)][tuple(coors)] = line
    return long_dict

def get_long_isoform_input(in_fn):
    long_dict = dd(lambda: [])
    with open(in_fn) as in_fp:
        for line in in_fp:
            if line.startswith('#'): continue
            ele = line.split()
            chrom = ele[whole_output_header_idx['chrom']]
            start = int(ele[whole_output_header_idx['startCoor0based']])
            end = int(ele[whole_output_header_idx['endCoor']])
            long_dict[(chrom, start, end)].append(line)
    return long_dict

def get_input(in_fn, in_type, detailed):
    if detailed and (in_type == 'ciri' or in_type == 'bed'):
        ut.fatal_format_time('get_input', 'ciri does not support detailed alignment information.')
    if in_type == 'whole':
        return get_detailed_long_whole_input(in_fn) if detailed else get_long_whole_input(in_fn)
    elif in_type == 'isoform':
        return get_detailed_long_isoform_input(in_fn) if detailed else get_long_isoform_input(in_fn)
    elif in_type == 'ciri':
        return get_short_ciri_input(in_fn)
    elif in_type == 'bed':
        return get_short_bed_input(in_fn)

# type: 'whole'/'isoform'/'ciri'/'bed'
def isocirc_comp_core(args):
    a_fn = args.a_input
    a_type = args.a_type
    a_dict = get_input(a_fn, a_type, args.detailed)
    b_fn = args.b_input
    b_type = args.b_type
    b_dict = get_input(b_fn, b_type, args.detailed)
    if args.overlap:
        ut.err_format_time('get_overlap', '{} {} ... '.format(a_fn, b_fn))
        if args.detailed:
            get_detailed_overlap(a_dict, b_dict, args.back_dis, args.inter_dis, args.overlap_entry, args.overlap)
        else:
            get_overlap(a_dict, b_dict, args.back_dis, args.overlap_entry, args.overlap)
        ut.err_format_time('get_overlap', '{} {} done!'.format(a_fn, b_fn))
    if args.a_only:
        ut.err_format_time('get_A_only', '{} {} ... '.format(a_fn, b_fn))
        get_only(a_dict, b_dict, args.back_dis, args.a_only)
        ut.err_format_time('get_A_only', '{} {} done!'.format(a_fn, b_fn))
    if args.b_only:
        ut.err_format_time('get_B_only', '{} {} ... '.format(a_fn, b_fn))
        get_only(b_dict, a_dict, args.back_dis, args.b_only)
        ut.err_format_time('get_B_only', '{} {} done!'.format(a_fn, b_fn))
    return

# parse command line arguments
def parser_argv():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="{}-comp: compare {} result with other long-read and/or short-read dataset".format(__program__, __program__))
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    input_par = parser.add_argument_group('Input options')
    input_par.add_argument('-a', '--a-input', required=True, type=str, help='{} long-read result file or short-read result file(Default: {} whole read-wise output).'.format(__program__, __program__))
    input_par.add_argument('-b', '--b-input', required=True, type=str, help='{} long-read result file or short-read result file(Default: {} whole read-wise output).'.format(__program__, __program__))

    input_par.add_argument('--a-type', default='isoform', choices=['whole', 'isoform', 'ciri', 'bed'], help='Input file type of A.')
    input_par.add_argument('--b-type', default='isoform', choices=['whole', 'isoform', 'ciri', 'bed'], help='Input file type of B.')

    cust_par = parser.add_argument_group('Customized options')
    cust_par.add_argument('-s', '--inter-dis', type=int, default=0, help='Maximum allowed distance between two internal-splice-sites.')
    cust_par.add_argument('-S', '--back-dis', type=int, default=10, help='Maximum allowed distance between two back-spice-sites.')

    output_par = parser.add_argument_group('Output options')

    output_par.add_argument('-O', '--overlap', type=str, help='Output file name of shared records.')
    output_par.add_argument('-d', '--detailed', default=False, action='store_true')
    output_par.add_argument('-w', '--overlap-entry', type=str, default='A', choices=['A', 'B', 'AB'], help='Write the entry of A/B/AB in overlap output file.')

    output_par.add_argument('-A', '--a-only', type=str, help='Output file name of records that only show up in A')
    output_par.add_argument('-B', '--b-only', type=str, help='Output file name of records that only show up in B')

    return parser.parse_args()

def main():
    args = parser_argv()
    isocirc_comp_core(args)

if __name__ == '__main__':
    main()