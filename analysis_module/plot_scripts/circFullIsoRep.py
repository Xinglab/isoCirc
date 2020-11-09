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

def get_count_bin(cnt):
    if cnt == 0: return 0
    if cnt == 1:
        return '1'
    elif cnt == 2:
        return '2'
    elif cnt >= 3 and cnt <= 5:
        return '3~5'
    elif cnt >=6 and cnt <= 10:
        return '6~10'
    elif cnt >=11 and cnt <=50:
        return '11~50'
    else: return '>50'

if len(sys.argv) < 8:
    print('Usage:')
    print('{} min_cnt in1.out in2.out ... pairFig'.format(sys.argv[0]))
    sys.exit(1)

min_read_cnt = int(sys.argv[1])
in_fn_list = sys.argv[2:-1]
fig = sys.argv[-1]
names= []
fn_n = len(in_fn_list)
tmp_dir=os.path.dirname(os.path.abspath(fig))
iso_out = tmp_dir + '/cnt{}_iso.out'.format(min_read_cnt)
anno_out = tmp_dir + '/cnt{}_iso_anno.out'.format(min_read_cnt)
cnt_out  = tmp_dir + '/cnt{}_iso_read_cnt.out'.format(min_read_cnt)
pairs_out = tmp_dir + '/cnt{}_iso_pair_cnt.out'.format(min_read_cnt)
len_out = tmp_dir + '/cnt{}_iso_len.out'.format(min_read_cnt)

for fn in in_fn_list:
    fp = open(fn)
    names.append(fn.rsplit('/')[-2])

iso_dict_list = []
all_iso_dict = dict()
with open(iso_out, 'w') as out_fp, open(anno_out, 'w') as anno_fp, open(cnt_out, 'w') as cnt_fp, open(pairs_out, 'w') as pair_fp, open(len_out, 'w') as len_fp:
    out_fp.write('Type\tValue\n')
    iso_cate = dd(lambda:'')
    full_iso = dict()
    for fn, name in zip(in_fn_list,names):
        with open(fn) as fp:
            for line in fp:
                if line.startswith('#'): continue
                ele = line.rsplit()
                iso = (ele[idx['chrom']], ele[idx['startCoor0base']], ele[idx['endCoor']], ele[idx['blockCount']], ele[idx['blockSize']], ele[idx['blockStarts']])
                if ele[idx['isFullLength']] == 'True':
                    full_iso[iso] = 1
    for fn, name in zip(in_fn_list,names):
        iso_dict = dd(lambda : dd(lambda :0))
        with open(fn) as fp:
            for line in fp:
                if line.startswith('#'): continue
                ele = line.rsplit()
                iso = (ele[idx['chrom']], ele[idx['startCoor0base']], ele[idx['endCoor']], ele[idx['blockCount']], ele[idx['blockSize']], ele[idx['blockStarts']])
                if iso not in full_iso: continue
                cnt = int(ele[idx['readCount']])
                tot_len = int(ele[idx['refMapLen']])
                block_cnt = int(ele[idx['blockCount']])
                iso_dict[iso]['read_cnt'] += cnt
                iso_dict[iso]['tot_len'] = tot_len
                iso_dict[iso]['block_cnt'] = block_cnt
                # iso_cate[iso] = ('FSM' if ele[idx['BSJCate']] == 'FSM' else 'NIC/NNC', 'FSM' if ele[idx['FSJCate']] == 'FSM' else 'NIC/NNC')
                iso_cate[iso] = ('NNC' if ele[idx['BSJCate']] == 'NNC' else 'FSM/NIC', 'NNC' if ele[idx['FSJCate']] == 'NNC' else 'FSM/NIC')
                if iso not in all_iso_dict:
                    known = []
                    is_known_bsj = ele[idx['isKnownBSJ']]
                    is_known_sj = ele[idx['isKnownSJ']]
                    known.append('KnownBSJ') if 'True' in is_known_bsj else known.append('NovelBSJ')
                    known.append('NovelFSJ') if 'False' in is_known_sj else known.append('KnownFSJ')
                    all_iso_dict[iso] = '_'.join(known)
        iso_dict_list.append(iso_dict)
    for iso_dict, name in zip(iso_dict_list, names):
        isos = list(iso_dict.keys())
        for iso in isos:
            cnt = iso_dict[iso]['read_cnt']
            if cnt >= min_read_cnt:
                out_fp.write('{}\t{}\n'.format(name, '_'.join(iso)))
            else: 
                del iso_dict[iso]
                del all_iso_dict[iso]
    anno_fp.write('Cate\tAnno\n')
    anno_dict = dd(lambda : dd(lambda :0))
    rep_cnt_dict = dict()
    for iso, anno in all_iso_dict.items():
        cnt = 0
        for iso_dict in iso_dict_list:
            if iso in iso_dict:
                cnt += 1
        if cnt == 0:
            print('Unexpected: 0')
            sys.exit(1)
        rep_cnt_dict[iso] = cnt
        anno_fp.write('{}\t{}\n'.format(cnt, anno))
    cnt_fp.write('Sample\tRep\tCount\n')
    len_fp.write('Sample\tRep\tCount\tCountBin\tCate\tLen\tBlockCount\n')
    for iso_dict, name in zip(iso_dict_list, names):
        for iso in iso_dict:
            cnt_fp.write('{}\t{}\t{}\n'.format(name, rep_cnt_dict[iso], iso_dict[iso]['read_cnt']))
            len_fp.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(name, rep_cnt_dict[iso], iso_dict[iso]['read_cnt'], get_count_bin(iso_dict[iso]['read_cnt']), '-'.join(iso_cate[iso]), iso_dict[iso]['tot_len'], iso_dict[iso]['block_cnt']))

            len_fp.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(name, rep_cnt_dict[iso], iso_dict[iso]['read_cnt'], get_count_bin(iso_dict[iso]['read_cnt']), 'All-All', iso_dict[iso]['tot_len'], iso_dict[iso]['block_cnt']))

    pair_fp.write('Sample1\tSample2\tCount1\tCount2\n')
    for i in range(len(names)):
        pair_fp.write('{}\t{}\tNA\tNA\n'.format(names[i], names[i]))
        for j in range(i+1, len(names)):
            name1, name2 = names[i], names[j]
            iso_dict1, iso_dict2 = iso_dict_list[i], iso_dict_list[j]
            for iso in all_iso_dict:
                cnt1, cnt2 = 0, 0
                if iso in iso_dict1:
                    cnt1 = iso_dict1[iso]['read_cnt']
                if iso in iso_dict2:
                    cnt2 = iso_dict2[iso]['read_cnt']
                if cnt1 or cnt2:
                    pair_fp.write('{}\t{}\t{}\t{}\n'.format(name1, name2, cnt1, cnt2))
                    pair_fp.write('{}\t{}\t{}\t{}\n'.format(name2, name1, cnt2, cnt1))



cmd='Rscript /home/gaoy1/program/circ_plot/circFullIsoRepPair.R {} {}'.format(pairs_out, fig)
print(cmd)
os.system(cmd)