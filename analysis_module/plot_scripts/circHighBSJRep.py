import sys, os, re
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

idx = {h: i for i, h in enumerate(header)}

if len(sys.argv) < 6:
    print('Usage:')
    print('{} min_cnt in1.out in2.out ... stackFig cdfFig pairFig'.format(sys.argv[0]))
    sys.exit(1)

min_read_cnt = int(sys.argv[1])
in_fn_list = sys.argv[2:-3]
[fig1, fig2, fig3] = sys.argv[-3:]
names= []
fn_n = len(in_fn_list)
tmp_dir=os.path.dirname(os.path.abspath(fig1))
bsj_out = tmp_dir + '/cnt{}_bsj.out'.format(min_read_cnt)
anno_out = tmp_dir + '/cnt{}_bsj_anno.out'.format(min_read_cnt)
cnt_out  = tmp_dir + '/cnt{}_bsj_read_cnt.out'.format(min_read_cnt)
pairs_out = tmp_dir + '/cnt{}_bsj_pair_cnt.out'.format(min_read_cnt)

in_fp_list = []
for fn in in_fn_list:
    fp = open(fn)
    names.append(fn.rsplit('/')[-2])
    in_fp_list.append(fp)

bsj_dict_list = []
all_bsj_dict = dict()
with open(bsj_out, 'w') as out_fp, open(anno_out, 'w') as anno_fp, open(cnt_out, 'w') as cnt_fp, open(pairs_out, 'w') as pair_fp:
    out_fp.write('Type\tValue\n')
    for fp, name in zip(in_fp_list,names):
        bsj_dict = dd(lambda: 0)
        for line in fp:
            if line.startswith('#'): continue
            ele = line.rsplit()
            bsj = (ele[idx['chrom']], ele[idx['startCoor0base']], ele[idx['endCoor']])
            cnt = int(ele[idx['readCount']])
            bsj_dict[bsj] += cnt
            if bsj not in all_bsj_dict:
                is_known = ele[idx['isKnownBSJ']]
                if ',' not in is_known:
                    print('Two circRNA annotations needed.')
                    sys.exit(1)
                if 'True' not in is_known:
                    all_bsj_dict[bsj] = 'Novel'
                elif 'False' not in is_known:
                    all_bsj_dict[bsj] = 'Both'
                elif is_known.startswith('True'):
                    all_bsj_dict[bsj] = 'circBase'
                else:
                    all_bsj_dict[bsj] = 'MiOncoCirc'
        bsj_dict_list.append(bsj_dict)
        for bsj, cnt in bsj_dict.items():
            if cnt >= min_read_cnt:
                out_fp.write('{}\t{}\n'.format(name, '_'.join(bsj)))
    for bsj_dict in bsj_dict_list:
        bsjs = list(bsj_dict.keys())
        for bsj in bsjs:
            if bsj_dict[bsj] < min_read_cnt:
                del bsj_dict[bsj]
    anno_fp.write('Cate\tAnno\tCnt\tFraction\n')
    anno_cnt_dict = dd(lambda : dd(lambda :0))
    cnt_dict = dict()
    for bsj, anno in all_bsj_dict.items():
        cnt = 0
        for bsj_dict in bsj_dict_list:
            if bsj in bsj_dict:
                cnt += 1
        if cnt == 0:
            continue
            print('Unexpected: 0')
            sys.exit(1)
        cnt_dict[bsj] = cnt
        anno_cnt_dict[cnt][anno]+=1
    for c,ad in anno_cnt_dict.items():
        s = 0.0
        for a,n in ad.items():
            s+=n
        for a,n in ad.items():
            anno_fp.write('{}\t{}\t{}\t{}\n'.format(c, a, n, n/s))

        # anno_fp.write('{}\t{}\n'.format(cnt, anno))
    cnt_fp.write('RepAnno\tCount\n')
    for bsj_dict in bsj_dict_list:
        for bsj, cnt in bsj_dict.items():
            anno = 'Novel' if all_bsj_dict[bsj] == 'Novel' else 'Known'
            cnt_fp.write('{}{}\t{}\n'.format(cnt_dict[bsj], anno, cnt))
    pair_fp.write('Sample1\tSample2\tCount1\tCount2\n')
    for i in range(len(names)):
        pair_fp.write('{}\t{}\tNA\tNA\n'.format(names[i], names[i]))
        for j in range(i+1, len(names)):
            name1, name2 = names[i], names[j]
            samp = '_'.join([name1, name2])
            bsj_dict1, bsj_dict2 = bsj_dict_list[i], bsj_dict_list[j]
            for bsj in all_bsj_dict:
                cnt1, cnt2 = 0, 0
                if bsj in bsj_dict1:
                    cnt1 = bsj_dict1[bsj]
                if bsj in bsj_dict2:
                    cnt2 = bsj_dict2[bsj]
                if cnt1 or cnt2:
                    pair_fp.write('{}\t{}\t{}\t{}\n'.format(name1, name2, cnt1, cnt2))
                    pair_fp.write('{}\t{}\t{}\t{}\n'.format(name2, name1, cnt2, cnt1))


cmd='Rscript /home/gaoy1/program/circ_plot/circHighBSJRepAnno.R {} {}'.format(anno_out, fig1)
print(cmd)
os.system(cmd)

cmd='Rscript /home/gaoy1/program/circ_plot/circHighBSJRepCnt.R {} {}'.format(cnt_out, fig2)
print(cmd)
os.system(cmd)

cmd='Rscript /home/gaoy1/program/circ_plot/circHighBSJRepPair.R {} {}'.format(pairs_out, fig3)
print(cmd)
os.system(cmd)

for fp in in_fp_list:
    fp.close()
