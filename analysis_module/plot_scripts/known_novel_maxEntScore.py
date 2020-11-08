import sys, os, re
from collections import defaultdict as dd
import mappy as mp
from pyfaidx import Fasta

header = ['#isoformID', 'chrom', 'startCoor0base', 'endCoor',  # 1-4
          'geneStrand', 'geneID', 'geneName',  # 5-7
          'blockCount', 'blockSize', 'blockStarts', 'refMapLen',  # 8-11
          'blockType', 'blockAnno',  # 12-13
          'isKnownSS', 'isKnownSJ', 'isCanoSJ', 'isHighSJ', 'isKnownExon',  # 14-18
          'isKnownBSJ', 'isCanoBSJ', 'canoBSJMotif',  # 19-21
          'isFullLength', 'BSJCate', 'interIsoCate',  # 22-24, FSM: full splice match, NIC: novel in catelog, NNC, novel and not in catelog
          'CDS', 'UTR', 'lincRNA', 'antisense',  # 25-28 gene_type/biotype
          'rRNA', 'Alu', 'allRepeat', 'upFlankAlu', 'downFlankAlu',  # 29-33 repeat element
          'readCount', 'readIDs']  # 34-35
idx = {h: i for i, h in enumerate(header)}

_5_perl='/home/gaoy1/software/maxent/score5.pl'
_3_perl='/home/gaoy1/software/maxent/score3.pl'
# 5_intron_len = 6
# 5_exon_len = 3

# 3_intron_len = 20
# 3_exon_len = 3

def extract_fa(fa, rep1, rep2, known_rep1_5_fa, known_rep2_5_fa, novel_rep1_5_fa, novel_rep2_5_fa, known_rep1_3_fa, known_rep2_3_fa, novel_rep1_3_fa, novel_rep2_3_fa):
    bsj_dict = dd(lambda: dd(lambda:0))

    rep = 0
    for fn in [rep1, rep2]:
        rep += 1
        with open(fn) as fp:
            for line in fp:
                if line.startswith('#'): continue
                ele = line.rsplit()
                bsj = (ele[idx['chrom']], ele[idx['startCoor0base']], ele[idx['endCoor']], ele[idx['canoBSJMotif']][0])
                is_known = ele[idx['isKnownBSJ']]
                if 'True' not in is_known:
                    bsj_dict[bsj]['Cate'] = 'Novel'
                else:
                    bsj_dict[bsj]['Cate'] = 'Known'
                bsj_dict[bsj][rep] = 1

    with open(known_rep1_5_fa, 'w') as k1_5, open(known_rep2_5_fa, 'w') as k2_5, open(novel_rep1_5_fa, 'w') as n1_5, open(novel_rep2_5_fa, 'w') as n2_5, open(known_rep1_3_fa, 'w') as k1_3, open(known_rep2_3_fa, 'w') as k2_3, open(novel_rep1_3_fa, 'w') as n1_3, open(novel_rep2_3_fa, 'w') as n2_3:
        for bsj in bsj_dict:
            (chrom, start, end, strand) = bsj
            start, end = int(start), int(end)
            ref_seq = fa[chrom]
            # 5'
            if strand == '+':
                seq = ref_seq[end-3:end+6].seq.upper()
                _5_str = '>{}:{}-{} {}\n{}\n'.format(chrom, end-2, end+6, strand, seq)
                seq = ref_seq[start-20:start+3].seq.upper()
                _3_str = '>{}:{}-{} {}\n{}\n'.format(chrom, start-19, start+3, strand, seq)
            elif strand == '-':
                seq = ref_seq[start-6:start+3].seq.upper()
                _5_str = '>{}:{}-{} {}\n{}\n'.format(chrom, start-5, start+3, strand, mp.revcomp(seq))
                seq = ref_seq[end-3:end+20].seq.upper()
                _3_str = '>{}:{}-{} {}\n{}\n'.format(chrom, end-2, end+20, strand, mp.revcomp(seq))
            else:
                print('Unexpected strand: {}'.format(strand))

            cate = bsj_dict[bsj]['Cate']
            rep_cnt = bsj_dict[bsj][1] + bsj_dict[bsj][2]
            if cate == 'Known':
                if rep_cnt == 1:
                    k1_5.write(_5_str)
                    k1_3.write(_3_str)
                elif rep_cnt == 2:
                    k2_5.write(_5_str)
                    k2_3.write(_3_str)
                else:
                    print('Unexpected rep_cnt: {}'.format(rep_cnt))
            elif cate == 'Novel':
                if rep_cnt == 1:
                    n1_5.write(_5_str)
                    n1_3.write(_3_str)
                elif rep_cnt == 2:
                    n2_5.write(_5_str)
                    n2_3.write(_3_str)
                else:
                    print('Unexpected rep_cnt: {}'.format(rep_cnt))
            else:
                print('Unexpected cate: {}'.format(cate))

if __name__ == '__main__':
    if len(sys.argv) != 5:
        print('{} ref.fa rep1.isocirc.out rep2.isocirc.out score.out'.format(sys.argv[0]))
        sys.exit(1)
    
    fa = Fasta(sys.argv[1])
    rep1, rep2 = sys.argv[2], sys.argv[3]
    score_out = sys.argv[4]
    
    known_rep1_5_fa, known_rep2_5_fa, novel_rep1_5_fa, novel_rep2_5_fa = score_out + 'known1.5sites.fa', score_out + 'known2.5sites.fa', score_out + 'novel1.5sites.fa', score_out + 'novel2.5sites.fa' 
    known_rep1_3_fa, known_rep2_3_fa, novel_rep1_3_fa, novel_rep2_3_fa = score_out + 'known1.3sites.fa', score_out + 'known2.3sites.fa', score_out + 'novel1.3sites.fa', score_out + 'novel2.3sites.fa' 
    
    for _5_fa in [known_rep1_5_fa, known_rep2_5_fa, novel_rep1_5_fa, novel_rep2_5_fa]:
        cmd = 'perl {} {} > {}.score'.format(_5_perl, _5_fa, _5_fa)
        print(cmd)
    
    for _3_fa in [known_rep1_3_fa, known_rep2_3_fa, novel_rep1_3_fa, novel_rep2_3_fa]:
        cmd = 'perl {} {} > {}.score'.format(_3_perl, _3_fa, _3_fa)
        print(cmd)

