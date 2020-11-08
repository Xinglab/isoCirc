import sys, os, re
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


if len(sys.argv) != 5:
    print('Usage:')
    print('{} min_cnt in.out cateTableCnt2Fig cateTableCnt3Fig'.format(sys.argv[0]))
    sys.exit(1)

min_read_cnt = int(sys.argv[1])
in_fn = sys.argv[2]
fig1 = sys.argv[3]
fig2 = sys.argv[4]
tmp_dir=os.path.dirname(os.path.abspath(fig1))
tmp_out = tmp_dir + '/tissue_bsj_iso_cate.out'

iso_cate = dd(lambda:0)
cate_list = ['FSM', 'NIC', 'NNC']

with open(tmp_out, 'w') as out_fp, open(in_fn) as in_fp:
    out_fp.write('IsoCnt\tReadCnt\tBSJ\tISO\n')
    first_line = True
    for line in in_fp:
        if first_line:
            first_line = False
            continue
        ele = line.rsplit()
        iso_cate[(ele[idx['BSJCate']], ele[idx['FSJCate']])] += 1
    for k, v in iso_cate.items():
        out_fp.write('{}\t2\t{}\t{}\n'.format(v, k[0], k[1]))

cmd='Rscript /home/gaoy1/program/circ_plot/circCateTable_tissue.R {} {}'.format(tmp_out, fig1)
print(cmd)
os.system(cmd)

# cmd='Rscript /home/gaoy1/program/circ_plot/circCateTableCnt3.R {} {}'.format(tmp_out, fig2)
# print(cmd)
# os.system(cmd)
