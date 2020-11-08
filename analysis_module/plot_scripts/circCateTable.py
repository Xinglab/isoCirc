import sys, os, re
from collections import defaultdict as dd

header = ['#isoformID', 'chrom', 'startCoor0base', 'endCoor', # 1-4
          'geneStrand', 'geneID', 'geneName', # 5-7
          'blockCount', 'blockSize', 'blockStarts', 'refMapLen', # 8-11
          'blockType', 'blockAnno', # 12-13
          'isKnownSS', 'isKnownSJ', 'isCanoSJ', 'isHighSJ', 'isKnownExon', # 14-18
          'isKnownBSJ', 'isCanoBSJ', 'canoBSJMotif', # 19-21
          'isFullLength', 'BSJCate', 'FSJCate', # FSM: full splice match, NIC: novel in catelog, NNC, novel and not in catelog
          'CDS', 'UTR', 'lincRNA', 'antisense',  # 22-25 gene_type/biotype
          'rRNA', 'Alu', 'allRepeat', 'upFlankAlu', 'downFlankAlu', # 26-30 repeat element
          'readCount', 'readIDs'] # 31-32
idx = {h: i for i, h in enumerate(header)}


if len(sys.argv) != 4:
    print('Usage:')
    print('{} min_cnt in.out cateHeatTableFig'.format(sys.argv[0]))
    sys.exit(1)

min_read_cnt = int(sys.argv[1])
in_fn = sys.argv[2]
fig = sys.argv[3]
tmp_dir=os.path.dirname(os.path.abspath(fig))
tmp_out = tmp_dir + '/bsj_iso_cate.out'

iso_read_cnt = dd(lambda:0)
iso_cate = dict()
full_iso = dict()
cate_list = ['FSM', 'NIC', 'NNC']
with open(in_fn) as in_fp:
    for line in in_fp:
        if line.startswith('#'): continue
        ele = line.rsplit()
        iso = (ele[idx['chrom']], ele[idx['startCoor0base']], ele[idx['endCoor']], ele[idx['blockCount']], ele[idx['blockSize']], ele[idx['blockStarts']])
        if ele[idx['isFullLength']] == 'True': full_iso[iso] = 1
with open(tmp_out, 'w') as out_fp, open(in_fn) as in_fp:
    out_fp.write('IsoCnt\tReadCnt\tBSJ\tISO\n')
    for line in in_fp:
        if line.startswith('#'): continue
        ele = line.rsplit()
        iso = (ele[idx['chrom']], ele[idx['startCoor0base']], ele[idx['endCoor']], ele[idx['blockCount']], ele[idx['blockSize']], ele[idx['blockStarts']])
        if iso not in full_iso: continue
        read_cnt = int(ele[idx['readCount']])
        iso_read_cnt[iso] += read_cnt
        iso_cate[iso] = (ele[idx['BSJCate']], ele[idx['FSJCate']])
    for i in range(1, min_read_cnt+1):
        cate_cnt = dd(lambda:0)
        for c1 in cate_list:
            for c2 in cate_list:
                cate_cnt[(c1, c2)] = 0
        for iso, cate in iso_cate.items():
            if iso_read_cnt[iso] >= i:
                cate_cnt[cate] += 1
        for cate, cnt in cate_cnt.items():
            out_fp.write('{}\t{}\t{}\n'.format(cnt, i, '\t'.join(cate)))
# cmd='Rscript /home/gaoy1/program/circ_plot/circCateTable.R {} {}'.format(tmp_out, fig)
# print(cmd)
# os.system(cmd)
