import sys, os, re
import collections
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

raw_names=('Adipose_1', 'Adrenal_1', 'Blood_1', 'Brain_1', 'Heart_1', 'Kidney_1',
        'Liver_1', 'Lung_1', 'Prostate_1', 'SkeletalMuscle_1', 'SmoothMuscle_1', 'Testis_1')
names=('Adipose', 'Adrenal', 'Blood', 'Brain', 'Heart', 'Kidney', 'Liver', 'Lung', 'Prostate', 'SkeletalMuscle', 'SmoothMuscle', 'Testis')
tissue_name_dict = {
    x:y for x,y in zip(raw_names, names)
}

def is_exon_intron_exon1(anno1):
    if re.search(r'E[0-9]+\([0-9]+\)\+I\([0-9]+\)\+E[0-9]+\(', anno1):
        return True
    else:
        return False

def is_exon_intron_exon(anno):
    annos = anno.rsplit(',')
    # for each block
    for anno1 in annos:
        # for each transcipt's annotation
        trans_annos = anno1.rsplit(';')
        for trans_anno1 in trans_annos:
            if is_exon_intron_exon1(trans_anno1):
                return True
    return False

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print('Usage:')
        print('{} cnt_cutoff comb.out isoNumBar.png'.format(sys.argv[0]))
        sys.exit(1)

    cnt_cutoff = int(sys.argv[1])
    in_fn = sys.argv[2]
    fig = sys.argv[3]

    tmp_dir = os.path.dirname(os.path.abspath(fig))
    tmp_out = tmp_dir + '/bsj_iso_num.out'

    full_iso_num = dd(lambda:dd(lambda:0))

    with open(in_fn) as in_fp:
        first_line = True
        for line in in_fp:
            if first_line:
                first_line = False
                continue
            ele = line.rsplit()
            for rname, name in zip(raw_names, names):
                cnt = int(ele[idx[rname]])
                if cnt == 0: continue
                if cnt < cnt_cutoff:
                    full_iso_num[name][cnt] += 1
                else:
                    full_iso_num[name][cnt_cutoff] += 1
                full_iso_num[name]['all'] += 1

    sorted_iso_num = sorted(full_iso_num.items(), key=lambda kv:kv[1]['all'])
    sorted_iso_num.reverse()
    sorted_iso_num_dict = collections.OrderedDict(sorted_iso_num)
    with open(tmp_out, 'w') as fp:
        fp.write('Sample\tType\tNumber\n')
        for name, num_dict in sorted_iso_num_dict.items():
            for i in range(1, cnt_cutoff+1):
                # fp.write('{}\tBSJNum\t{}\n'.format(name, num_dict['bsj1']))
                fp.write('{}\tISONum{}\t{}\n'.format(name, i, num_dict[i]))
                # fp.write('{}\tEIEISONum\t{}\n'.format(name, num_dict['eie1']))


    cmd='Rscript /home/gaoy1/program/circ_plot/isonumBarPlot.R {} {}'.format(tmp_out, fig)
    print(cmd)
    os.system(cmd)