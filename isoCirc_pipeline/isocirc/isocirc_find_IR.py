#!/usr/bin/env python
import argparse
from collections import defaultdict as dd
import itertools
import re

import isocirc.utils as ut
import isocirc.hcBSJ_fullIso as hf
from isocirc.__init__ import __version__
from isocirc.__init__ import __program__

'''
Compare PARRIS results between different long-read datasets, or long-read and short-read datasets
'''

ciri_header = ['circRNA_ID','chr','circRNA_start','circRNA_end', '#junction_reads', 'SM_MS_SMS','#non_junction_reads',
               'junction_reads_ratio', 'circRNA_type', 'gene_id', 'strand', 'junction_reads_ID']
ciri_header_idx = {h: i for i, h in enumerate(ciri_header)}


whole_output_header = hf.whole_output_header
whole_output_header_idx = hf.whole_output_header_idx
isoform_output_header = hf.isoform_output_header
isoform_output_header_idx = hf.isoform_output_header_idx

def isocirc_find_IR_core(in_fn, min_exon_len, idx):
    with open(in_fn) as in_fp:
        for line in in_fp:
            if line.startswith('#'): continue
            ele = line.rsplit()
            if int(ele[idx['blockCount']]) != 1: continue
            anno = ele[idx['blockAnno']]
            if anno == 'NA': continue
            # find intron retaion
            trans = anno.rsplit(';')
            for t in trans:
                # print t
                [tname, exon_intron] = t.rsplit(':')
                eis = exon_intron.rsplit('+')
                if len(eis) != 3: continue
                eids, elens, ilens = [], [], []
                for ei in eis:
                    if ei[0] == 'E': # exon
                        [eid, elen] = list(map(int, re.split('\(|\)', ei[1:])[:-1]))
                        eids.append(eid)
                        elens.append(elen)
                    else: # intron
                        ilens.append(int(re.split('\(|\)', ei)[1]))
                if len(eids) == 2 and abs(eids[0]-eids[1]) == 1 and min(elens) >= min_exon_len:
                    print(line)

# type: 'whole'/'isoform'/'ciri'/'bed'
def isocirc_find_IR(args):
    in_fn = args.in_out
    min_len = args.min_len

    # idx = isoform_output_header_idx
    idx = whole_output_header_idx
    isocirc_find_IR_core(in_fn, min_len, idx)
    return

# parse command line arguments
def parser_argv():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="{}-IR: find intron-retation from {} result".format(__program__, __program__))
    parser.add_argument('in_out', metavar="in.out", type=str, help='{} long-read result file.'.format(__program__))
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)
    parser.add_argument('-m', '--min-len', type=int, default=10, help='Minimum length of each exon block.')

    return parser.parse_args()

def main():
    args = parser_argv()
    isocirc_find_IR(args)

if __name__ == '__main__':
    main()