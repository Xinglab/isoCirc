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

if len(sys.argv) !=4:
    print('Usage:')
    print('{} comb.out exonNumFig lenFig'.format(sys.argv[0]))
    sys.exit(1)

in_fn = sys.argv[1]
[fig1, fig2] = sys.argv[2:]
names= []
tmp_dir=os.path.dirname(os.path.abspath(fig1))
len_out = tmp_dir + '/tissue_exon_len.out'

iso_dict_list = []
all_iso_dict = dict()
with open(len_out, 'w') as len_fp, open(in_fn) as in_fp:
    len_fp.write('Cate\tLen\tExonNum\n')
    first_line = True
    for line in in_fp:
        if first_line:
            first_line = False
            continue
        ele = line.rsplit()
        tot_len = int(ele[idx['refMapLen']])
        exon_num = int(ele[idx['blockCount']])
        iso_cate = ('NNC' if ele[idx['BSJCate']] == 'NNC' else 'FSM/NIC', 'NNC' if ele[idx['FSJCate']] == 'NNC' else 'FSM/NIC')
        len_fp.write('{}\t{}\t{}\n'.format('-'.join(iso_cate), tot_len, exon_num))
        len_fp.write('{}\t{}\t{}\n'.format('All-All', tot_len, exon_num))



cmd='Rscript /home/gaoy1/program/circ_plot/circFullIsoRepBlockCnt.R {} {}'.format(len_out, fig1)
print(cmd)
os.system(cmd)
    
cmd='Rscript /home/gaoy1/program/circ_plot/circFullIsoRepLen.R {} {}'.format(len_out, fig2)
print(cmd)
os.system(cmd)
