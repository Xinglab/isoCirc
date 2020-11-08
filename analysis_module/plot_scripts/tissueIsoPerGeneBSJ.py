#!/usr/bin/python
import os, sys, re, argparse
import math
from collections import defaultdict as dd

header = [
        'isoformID', 'chrom', 'startCoor0based', 'endCoor', # 1-4
        'geneStrand', 'geneID', 'geneName', # 5-7
        'blockCount', 'blockSize', 'blockStarts', 'refMapLen', # 8-11
        'blockType', 'blockAnno', # 12-13
        'isKnownSS', 'isKnownSJ', 'isCanoSJ', 'isHighSJ', 'isKnownExon', # 14-18
        'isKnownBSJ', 'isCanoBSJ', 'canoBSJMotif', # 19-21
        'isFullLength', 'BSJCate', 'FSJCate', # 22-24
        'CDS', 'UTR', 'lincRNA', 'antisense',  # 25-28 gene_type/biotype
        'rRNA', 'Alu', 'allRepeat', 'upFlankAlu', 'downFlankAlu', 
        'Adipose_1', 'Adrenal_1', 'Blood_1', 'Brain_1', 'Heart_1', 'Kidney_1',
        'Liver_1', 'Lung_1', 'Prostate_1', 'SkeletalMuscle_1', 'SmoothMuscle_1', 'Testis_1'
]

tissues = ['Adipose_1', 'Adrenal_1', 'Blood_1', 'Brain_1', 'Heart_1', 'Kidney_1',
        'Liver_1', 'Lung_1', 'Prostate_1', 'SkeletalMuscle_1', 'SmoothMuscle_1', 'Testis_1']

idx = {h: i for i, h in enumerate(header)}


def parser_argv():
    # parse command line arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="")
    parser.add_argument("isocirc_out", metavar='IsoCirc.out', type=str, help='IsoCirc output file.')
    parser.add_argument('genecdfPlot', metavar='geneIsoNum.pdf', type=str, help='CDF plot file.')
    parser.add_argument('BSJcdfPlot', metavar='BSJIsoNum.pdf', type=str, help='CDF plot file.')

    return parser.parse_args()


def circ_iso_num_core(args):
    isocirc_out, gene_fig, bsj_fig = args.isocirc_out, args.genecdfPlot, args.BSJcdfPlot

    tmp_dir = os.path.dirname(os.path.abspath(gene_fig))
    tmp_dat = tmp_dir + '/tissue_iso_num_per_gene_bsj.dat'

    rep_id = 0
    iso_per_gene_dict, iso_per_bsj_dict = dd(lambda: dd(lambda:0)), dd(lambda: dd(lambda:0))
    bsj_rep_dict, gene_rep_dict = dd(lambda:dict()), dd(lambda:dict())

    with open(tmp_dat, 'w') as out_fp, open(isocirc_out) as isocirc_fp:
        out_fp.write('Type\tRepCnt\tIsoCnt\tReadCnt\n')
        first_line = True
        for line in isocirc_fp:
            if first_line:
                first_line = False
                continue
            ele = line.rsplit()
            read_cnt = 0
            for t in tissues:
                read_cnt += int(ele[idx[t]])
            bsj = (ele[idx['chrom']], ele[idx['startCoor0based']], ele[idx['endCoor']])
            iso = (bsj[0], bsj[1], bsj[2], ele[idx['blockCount']], ele[idx['blockSize']], ele[idx['blockStarts']])
            iso_per_bsj_dict[bsj][iso] += read_cnt
            bsj_rep_dict[bsj][rep_id] = 1
            if ele[idx['geneID']] != 'NA':
                for _id in ele[idx['geneID']].rsplit(','):
                    iso_per_gene_dict[_id][iso] += read_cnt
                    gene_rep_dict[_id][rep_id] = 1
        for gene_id, iso_cnt in iso_per_gene_dict.items():
            iso_n, read_n = 0, 0
            for iso, cnt in iso_cnt.items():
                iso_n += 1
                read_n += cnt
            if iso_n > 0:
                out_fp.write("Gene\t{}\t{}\t{}\n".format(len(gene_rep_dict[gene_id]), iso_n, read_n))
        for bsj, iso_cnt in iso_per_bsj_dict.items():
            iso_n, read_n = 0, 0
            for iso, cnt in iso_cnt.items():
                iso_n += 1
                read_n += cnt
            if iso_n > 0:
                out_fp.write("BSJ\t{}\t{}\t{}\n".format(len(bsj_rep_dict[bsj]), iso_n, read_n))

    # 5. R plot venn diagram to fig
    # cmd = 'Rscript /home/gaoy1/program/circ_plot/circIsoNumScatter.R {} {}'.format(tmp_dat, fig)
    # print(cmd)
    cmd = 'Rscript /home/gaoy1/program/circ_plot/circIsoNumPerGeneCDF.R {} {}'.format(tmp_dat, gene_fig)
    print(cmd)
    os.system(cmd)
    cmd = 'Rscript /home/gaoy1/program/circ_plot/circIsoNumPerBSJCDF.R {} {}'.format(tmp_dat, bsj_fig)
    print(cmd)
    os.system(cmd)


def main():
    args = parser_argv()
    circ_iso_num_core(args)


if __name__ == '__main__':
    main()
