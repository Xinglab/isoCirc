import argparse
import sys
import os

import isocirc.utils as ut

lordec = 'lordec-correct' #dir_path + '/deps/lordec-src_0.9/lordec-correct'
kmer=21
solid=3
threads=8

def error_corr(in_long, in_short, out_long, lordec=lordec, kmer=kmer, solid=solid, threads=threads):
    ut.exec_cmd(sys.stderr, 'LoRDEC', '{} -2 {} -i {} -o {} -k {} -s {} -T {}'.format(lordec, in_short, in_long, out_long, kmer, solid, threads))
    return out_long


def error_corr_core(args):
    error_corr(args.in_long_read, args.in_short_read, args.out_long_read, args.lordec, args.kmer, args.solid, args.threads)


def parser_argv():
    # parse command line arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Hybrid error-correction for long-read data")
    parser.add_argument("in_long_read", metavar='long.fa', type=argparse.FileType('r'), help='Raw long-read sequencing data.')
    parser.add_argument("in_short_read", metavar='short.fa', type=argparse.FileType('r'), help='Short-read data. Use \',\' to connect multiple or paired-end short read data.')
    parser.add_argument("out_long_read", metavar='correct_long.fa', type=argparse.FileType('w'), help='Corrected long-read data.')

    parser.add_argument('-l', '--lordec', type=str, help='Path to LoRDEC-correct.', default=lordec)
    parser.add_argument('-k', '--kmer', type=int, help='k-mer size.', default=kmer)
    parser.add_argument('-s', '--solid', type=int, help='Solid k-mer abundance threshold.', default=solid)
    parser.add_argument('-T', '--threads', type=int, help='Number of threads.', default=threads)

    return parser.parse_args()


if __name__ == '__main__':
    args = parser_argv()

    error_corr_core(args)
