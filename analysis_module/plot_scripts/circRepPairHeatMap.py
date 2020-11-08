import sys, os, re
from collections import defaultdict as dd
from itertools import permutations as perm

header = ['#isoformID', 'chrom', 'startCoor0base', 'endCoor', # 1-4
         'geneStrand', 'geneID', 'geneName', # 5-7
         'blockCount', 'blockSize', 'blockStarts', 'refMapLen', # 8-11
         'blockType', 'blockAnno', # 12-13
         'isKnownSS', 'isKnownSJ', 'canoSJMotif', 'isHighSJ', 'isKnownExon', # 14-18
         'isKnownBSJ', 'isCanoBSJ', 'canoBSJMotif', # 19-21
         'isFullLength', 'BSJCate', 'FSJCate', # 22-24 FSM: full splice match, NIC: novel in catelog, NNC, novel and not in catelog
         'CDS', 'UTR', 'lincRNA', 'antisense',  # 25-28 gene_type/biotype 
         'rRNA', 'Alu', 'allRepeat', 'upFlankAlu', 'downFlankAlu', # 29-33 repeat element
         'readCount', 'readIDs'] # 34-35
idx = {h: i for i, h in enumerate(header)}

ciri_sum_header = ['circRNA_ID', 'XRJ_circ_S1', 'XRJ_circ_S3',
                   'XRJ_circ_R3', 'XRJ_circ_R2', 'XRJ_circ_S2',
                   'XRJ_circ_R1']
ciri_sum_idx = {h: i for i,h in enumerate(ciri_sum_header)}

def get_jaccard_index(dict1, dict2):
    itsct_cnt, union_cnt = 0, len(dict2)
    for i in dict1:
        if i in dict2:
            itsct_cnt += 1
        else:
            union_cnt += 1
    return itsct_cnt, union_cnt, itsct_cnt / (union_cnt+0.0)

# min_read_cnt = 3

# if len(sys.argv) < 8:
if len(sys.argv) < 5:
    print('Usage:')
    print('{} cnt_cutoff in1.out in2.out ... BSJHeatMap.fig IsoHeatMap.fig'.format(sys.argv[0]))
    sys.exit(1)

read_cutoff = int(sys.argv[1])
in_fn_list = sys.argv[2:-2]
[fig1, fig2] = sys.argv[-2:]
iso_names = []
bsj_names = []
fn_n = len(in_fn_list)
tmp_dir=os.path.dirname(os.path.abspath(fig1))
bsj_pair_out = tmp_dir + '/bsj_cnt.out'
iso_pair_out = tmp_dir + '/iso_cnt.out'

all_bsj_dict = dict()
all_iso_dict = dict()

full_iso = dict()
for fn in in_fn_list:
    name = fn.rsplit('/')[-2]
    if 'CIRI' in name: continue
    with open(fn) as fp:
        for line in fp:
            if line.startswith('#'): continue
            ele = line.rsplit()
            iso = (ele[idx['chrom']], ele[idx['startCoor0base']], ele[idx['endCoor']], ele[idx['blockCount']], ele[idx['blockSize']], ele[idx['blockStarts']])
            if ele[idx['isFullLength']] == 'True':
                full_iso[iso] = 1

with open(bsj_pair_out, 'w') as bsj_out_fp, open(iso_pair_out, 'w') as iso_out_fp: 
    bsj_out_fp.write('Cnt\tVar1\tVar2\tValue\n')
    iso_out_fp.write('Cnt\tVar1\tVar2\tValue\n')
    # bsj_out_fp.write('Cnt\tVar1\tVar2\tTot1\tTot2\tItstN\tUnionN\tValue\n')
    # iso_out_fp.write('Cnt\tVar1\tVar2\tTot1\tTot2\tItstN\tUnionN\tValue\n')
    for fn in in_fn_list:
        name = fn.rsplit('/')[-2]
        if 'CIRI' in name:
            r1_bsj_dict, r2_bsj_dict, r3_bsj_dict = dd(lambda:0), dd(lambda:0), dd(lambda:0)
            with open(fn) as fp:
                for line in fp:
                    ele = line.rsplit()
                    if len(ele) != len(ciri_sum_header): continue
                    if not ele[ciri_sum_idx['XRJ_circ_R1']].isdigit(): continue
                    circ_id = ele[ciri_sum_idx['circRNA_ID']]
                    _bsj = re.split(':|\|', circ_id)
                    bsj = (_bsj[0], str(int(_bsj[1])-1), _bsj[2])
                    r1_bsj_dict[bsj] += int(ele[ciri_sum_idx['XRJ_circ_R1']])
                    r2_bsj_dict[bsj] += int(ele[ciri_sum_idx['XRJ_circ_R2']])
                    r3_bsj_dict[bsj] += int(ele[ciri_sum_idx['XRJ_circ_R3']])
            all_bsj_dict['Illumina_1'] = r1_bsj_dict
            all_bsj_dict['Illumina_2'] = r2_bsj_dict
            all_bsj_dict['Illumina_3'] = r3_bsj_dict
            bsj_names.extend(['Illumina_1', 'Illumina_2', 'Illumina_3'])
        else:
            bsj_names.append(name)
            iso_names.append(name)
            bsj_dict, iso_dict = dd(lambda:0), dd(lambda:0)
            with open(fn) as fp:
                for line in fp:
                    if line.startswith('#'): continue
                    ele = line.rsplit()
                    cnt = int(ele[idx['readCount']])
                    bsj = (ele[idx['chrom']], ele[idx['startCoor0base']], ele[idx['endCoor']])
                    bsj_dict[bsj] += cnt
                    iso = (ele[idx['chrom']], ele[idx['startCoor0base']], ele[idx['endCoor']], ele[idx['blockCount']], ele[idx['blockSize']], ele[idx['blockStarts']])
                    if iso not in full_iso: continue
                    iso_dict[iso] += cnt
            all_bsj_dict[name] = bsj_dict
            all_iso_dict[name] = iso_dict
    # 
    for names, dicts, out_fp in zip([bsj_names, iso_names], [all_bsj_dict, all_iso_dict], [bsj_out_fp, iso_out_fp]):
        for (name1, name2) in list(perm(names, 2)):
            for min_cnt in range(1, read_cutoff+1):
                _dict1, _dict2 = dict(), dict()
                for i, cnt in dicts[name1].items():
                    if cnt >= min_cnt:
                        _dict1[i] = 1
                for i, cnt in dicts[name2].items():
                    if cnt >= min_cnt:
                        _dict2[i] = 1
                # out_fp.write('{}\t{}\t{}\t{:.3f}\n'.format(min_cnt, name1, name2, get_jaccard_index(_dict1, _dict2)))
                itsct_cnt, union_cnt, ratio = get_jaccard_index(_dict1, _dict2)
                # out_fp.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3f}\n'.format(min_cnt, name1, name2, len(_dict1), len(_dict2), itsct_cnt, union_cnt, ratio))
                out_fp.write('{}\t{}\t{}\t{:.3f}\n'.format(min_cnt, name1, name2, ratio))
        for name in names:
            for min_cnt in range(1, read_cutoff+1):
                out_fp.write('{}\t{}\t{}\t1\n'.format(min_cnt, name, name))
                out_fp.write('{}\t{}\t{}\t1\n'.format(min_cnt, name, name))

# cmd='Rscript /home/gaoy1/program/circ_plot/pairCorrHeatMap.R {} {}'.format(bsj_pair_out, fig1)
cmd='Rscript /home/gaoy1/program/circ_plot/pairHeatMap.R {} {}'.format(bsj_pair_out, fig1)
print(cmd)
os.system(cmd)

# cmd='Rscript /home/gaoy1/program/circ_plot/pairCorrHeatMap.R {} {}'.format(iso_pair_out, fig2)
cmd='Rscript /home/gaoy1/program/circ_plot/pairHeatMap.R {} {}'.format(iso_pair_out, fig2)
print(cmd)
os.system(cmd)