#!/usr/bin/python
import sys, re, os, argparse
from collections import defaultdict as dd
header = ['#isoformID', 'chrom', 'startCoor0base', 'endCoor', # 1-4
          'geneStrand', 'geneID', 'geneName', # 5-7
          'blockCount', 'blockSize', 'blockStarts', 'refMapLen', # 8-11
          'blockType', 'blockAnno', # 12-13
          'isKnownSS', 'isKnownSJ', 'isCanoSJ', 'isHighSJ', 'isKnownExon', # 14-18
          'isKnownBSJ', 'isCanoBSJ', 'canoBSJMotif', # 19-21
          'isFullLength', 'BSJCate', 'interIsoCate', # FSM: full splice match, NIC: novel in catelog, NNC, novel and not in catelog
          'CDS', 'UTR', 'lincRNA', 'antisense',  # 22-25 gene_type/biotype
          'rRNA', 'Alu', 'allRepeat', 'upFlankAlu', 'downFlankAlu', # 26-30 repeat element
          'readCount', 'readIDs'] # 31-32

header_idx = {h: i for i, h in enumerate(header)}

ciri_sum_header = ['circRNA_ID', 'XRJ_circ_S1', 'XRJ_circ_S3',
                   'XRJ_circ_R3', 'XRJ_circ_R2', 'XRJ_circ_S2',
                   'XRJ_circ_R1']
ciri_sum_idx = {h: i for i,h in enumerate(ciri_sum_header)}

# ciri_header = ['circRNA_ID', 'chr', 'circRNA_start', 'circRNA_end',
#                '#junction_reads',   'SM_MS_SMS', '#non_junction_reads', 
#                'junction_reads_ratio', 'circRNA_type', 'gene_id', 'strand', 
#                'junction_reads_ID']
# ciri_idx = {h: i for i,h in enumerate(ciri_header)}

def parser_argv():
    # parse command line arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="")
    parser.add_argument("isocirc_out", metavar='IsoCirc.out', type=str, help='IsoCirc output file.')
    parser.add_argument('ciri_out', metavar='ciri.out', type=str, help='CIRI output file.')
    parser.add_argument('min_cnt', metavar='min_cnt', type=int, help='Minimum read count for each BSJ')
    parser.add_argument('scatter_png', metavar='circShortLong.png', type=str, help='Scatter plot file.')
    return parser.parse_args()

def circ_spec_core(args):
    isocirc_out, ciri_sum_out, min_cnt, fig = args.isocirc_out, args.ciri_out, int(args.min_cnt), args.scatter_png 
    tmp_dir=os.path.dirname(os.path.abspath(fig))
    tmp_dat=tmp_dir+'/short_long.dat'

    bsj_cnt = dd(lambda : dd(lambda:0)) # shortCnt, longCnt
    with open(tmp_dat, 'w') as out_fp, open(isocirc_out) as isocirc_fp, open(ciri_sum_out) as ciri_fp:
        out_fp.write('BSJ\tShortCnt\tLongCnt\n')
        # 1. get ciri bsjs, cmp with ciri bsjs, output in_ciri.bsj
        for line in ciri_fp:
            ele = line.rsplit()
            if len(ele) != len(ciri_sum_header): continue
            if not ele[ciri_sum_idx['XRJ_circ_R1']].isdigit(): continue
            cnt = int(ele[ciri_sum_idx['XRJ_circ_R1']]) + int(ele[ciri_sum_idx['XRJ_circ_R2']]) + int(ele[ciri_sum_idx['XRJ_circ_R3']])
            if cnt == 0: continue
            circ_id = ele[ciri_sum_idx['circRNA_ID']]
            _bsj = re.split(':|\|', circ_id)
            bsj = (_bsj[0], str(int(_bsj[1])-1), _bsj[2])
            bsj_cnt[bsj]['short'] += cnt
        # 2. get all bsjs, output all bsjs
        for line in isocirc_fp:
            if line.startswith('#'): continue
            ele = line.rsplit()
            bsj = (ele[header_idx['chrom']], ele[header_idx['startCoor0base']], ele[header_idx['endCoor']])
            bsj_cnt[bsj]['long'] += int(ele[header_idx['readCount']])
        for bsj, cnt_dict in bsj_cnt.items():
            if (cnt_dict['short'] >= min_cnt and cnt_dict['long'] >= min_cnt):
                out_fp.write('{}\t{}\t{}\n'.format('_'.join(bsj), cnt_dict['short'], cnt_dict['long']))

    # 5. R plot venn diagram to fig
    cmd='Rscript /home/gaoy1/program/circ_plot/circShortLong.R {} {}'.format(tmp_dat, fig)
    print(cmd)
    os.system(cmd)

if __name__ == '__main__':
    args = parser_argv()
    circ_spec_core(args)