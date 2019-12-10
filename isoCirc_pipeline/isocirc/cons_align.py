import argparse
import os
import sys
from pyfaidx import Fasta
import pysam as ps

import isocirc.utils as ut

minimap2 = 'minimap2' #dir_path + '/bin/minimap2'
threads = 8

# def rotated_cons_based_on_sam(cons_fa, rotated_cons, sam_fn):
    # with Fasta(cons_fa) as in_cons, open(rotated_cons, 'w') as out_cons, ps.AlignmentFile(sam_fn) as in_sam:


# cons_fa: 2 copies of consensus sequence
def align_cons(cons_fa, ref_fa, cons_all_sam, minimap2=minimap2, threads=threads):
    ut.exec_cmd(sys.stderr, 'Mapping', '{} -ax splice -ub --MD --eqx {} {} -t {} > {}'.format(minimap2, ref_fa, cons_fa, threads, cons_all_sam))

# cons_fa: 2 copies of consensus sequence
def new_align_cons(cons_fa, ref_fa, cons_all_sam, minimap2=minimap2, threads=threads):
    # 1st round alignment
    tmp_all_sam = cons_all_sam + '.tmp'
    ut.exec_cmd(sys.stderr, 'Mapping', '{} -ax splice -ub --MD --eqx {} {} -t {} > {}'.format(minimap2, ref_fa, cons_fa, threads, tmp_all_sam))
    # rotate cons based on the alignment
    rotated_cons = cons_fa + '.rotated'
    rotate_cons_based_on_sam(cons_fa, rotated_cons, tmp_all_sam)

    # 2nd round alignment
    ut.exec_cmd(sys.stderr, 'Mapping', '{} -ax splice -ub --MD --eqx {} {} -t {} > {}'.format(minimap2, ref_fa, rotated_cons, threads, cons_all_sam))

def cons_align_core(args):
    if args.two_rounds:
        new_align_cons(args.cons_fa, args.ref_fa, args.cons_all_sam, args.minimap2, args.threads)
    else:
        align_cons(args.cons_fa, args.ref_fa, args.cons_all_sam, args.minimap2, args.threads)


def parser_argv():
    # parse command line arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Mapping of consensus sequence (minimap2)")
    parser.add_argument("ref_fa", type=str, metavar='ref.fa', help='Reference genome sequence file.')
    parser.add_argument("cons_fa", metavar='cons.fa', type=str, help='Consensus sequence file.')
    parser.add_argument("cons_all_sam", metavar='cons.fa.sam', type=str, help='Consensus SAM alignment file.')

    parser.add_argument('--minimap2', help='Path to minimap2.', default=minimap2)
    parser.add_argument('--two-rounds', help='Do two-rounds of alignment.', default=False, action='store_true')
    parser.add_argument('-t', '--threads', type=int, help='Number of threads.', default=threads)

    return parser.parse_args()


if __name__ == '__main__':
    args = parser_argv()

    cons_align_core(args)
