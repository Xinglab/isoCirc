import os,sys
import random
from datetime import datetime
from collections import defaultdict as dd
bedtools = 'bedtools'

header = ['#isoformID', 'chrom', 'startCoor0base', 'endCoor', # 1-4
          'geneStrand', 'geneID', 'geneName', # 5-7
          'blockCount', 'blockSize', 'blockStarts', 'refMapLen', # 8-11
          'blockType', 'blockAnno', # 12-13
          'isKnownSS', 'isKnownSJ', 'canoFSJMotif', 'isHighSJ', 'isKnownExon', # 14-18
          'isKnownBSJ', 'isCanoBSJ', 'canoBSJMotif', # 19-21
          'isFullLength', 'BSJCate', 'FSJCate', # 22-24 FSM: full splice match, NIC: novel in catelog, NNC, novel and not in catelog
          'CDS', 'UTR', 'lincRNA', 'antisense',  # 25-28 gene_type/biotype
          'rRNA', 'Alu', 'allRepeat', 'upFlankAlu', 'downFlankAlu', # 29-33 repeat element
          'readCount', 'readIDs'] # 34-35
idx = {h: i for i, h in enumerate(header)} 

genepred_header = ['transName', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'score', 'geneName', 'cdsStartStatus', 'cdsEndStatus', 'exonFrames'] # similar to bed12
genepred_idx = {h: i for i, h in enumerate(genepred_header)}


def add_neg_control(neg_control_bsj, bsj_rep_cnt_dict):
    for bsj in neg_control_bsj:
        if bsj in bsj_rep_cnt_dict:
            print('Error in neg_control_bsj')
            sys.exit(1)
        else:
            bsj_rep_cnt_dict[bsj] = 0
    return bsj_rep_cnt_dict
# cate_out_dict = dd(lambda: dd(lambda:0)) #  {RepAnno : {Cate: num}}
# for negative control: {Neg: {Cate:num}}
def get_neg_control_non_bsj(all_non_bsj_sites, tot_n):
    seed = datetime.now() # 1024
    random.seed(seed)
    skip_name = dict()
    N = tot_n
    gene_names = list(all_non_bsj_sites.keys())
    gene_n = len(gene_names)
    neg_control_sites = dd(lambda:dict())
    neg_control_bsj = dict()
    while N > 0:
        name = gene_names[random.randint(0, gene_n-1)]
        if name in skip_name:
            continue
        up_sites, down_sites = list(all_non_bsj_sites[name]['up'].keys()), list(all_non_bsj_sites[name]['down'].keys())
        up_n, down_n = len(up_sites), len(down_sites)
        if up_n == 0 or down_n == 0:
            skip_name[name] = 1
            continue
        # print(name, up_n, down_n)
        up_site, down_site = up_sites[random.randint(0, up_n-1)], down_sites[random.randint(0, down_n-1)]
        if up_site not in neg_control_sites['up'] and down_site not in neg_control_sites['down'] and int(up_site[1]) < int(down_site[1]) and up_site[2] == down_site[2]:
            neg_control_sites['up'][up_site] = 1
            neg_control_sites['down'][down_site] = 1
            neg_control_bsj[(up_site[0], up_site[1], down_site[1])] = 1
            N -= 1
    return neg_control_bsj

def get_circ_sites(all_bsj_sites, circbase, mionco):
    with open(circbase) as circ_fp, open(mionco) as mionco_fp:
        for fp in [circ_fp, mionco_fp]:
            for line in fp:
                ele = line.rsplit()
                chrom, start, end = ele[0], ele[1], ele[2]
                all_bsj_sites['up'][(chrom, start)] = 1
                all_bsj_sites['down'][(chrom, end)] = 1

def get_all_sites(gene_pred, bsj_genes, all_bsj_sites):
    all_non_bsj_sites = dd(lambda : dd(lambda: dict())) # 'up': , 'down':
    with open(gene_pred) as fp:
        for line in fp:
            ele = line.rsplit()
            geneName = ele[genepred_idx['geneName']]
            if geneName not in bsj_genes:
                continue
            chrom, strand, starts, ends = ele[genepred_idx['chrom']], ele[genepred_idx['strand']], ele[genepred_idx['exonStarts']].rsplit(','), ele[genepred_idx['exonEnds']].rsplit(',')
            # collect all splice sites
            up_sites, down_sites = [], []
            for i in range(int(ele[genepred_idx['exonCount']])):
                up_sites.append((chrom, starts[i]))
                down_sites.append((chrom, ends[i]))
            up_sites = up_sites[1:]
            down_sites = down_sites[:-1]
            for sites, label in zip([up_sites, down_sites], ['up', 'down']):
                for s in sites:
                    if s not in all_bsj_sites[label]:
                        all_non_bsj_sites[geneName][label][(s[0], s[1], strand)] = 1
    return all_non_bsj_sites


def get_bsj_dict(min_cnt, rep1_out, rep2_out, all_bsj_sites):
    bsj_rep_cnt_dict, bsj_known_novel_dict, bsj_gene = dd(lambda:0), dict(), dict()
    for fn in [rep1_out, rep2_out]:
        bsj_dict = dd(lambda: 0)
        with open(fn) as fp:
            for line in fp:
                if line.startswith('#'): continue
                ele = line.rsplit()
                bsj = (ele[idx['chrom']], ele[idx['startCoor0base']], ele[idx['endCoor']])
                all_bsj_sites['up'][(ele[idx['chrom']], ele[idx['startCoor0base']])] = 1
                all_bsj_sites['down'][(ele[idx['chrom']], ele[idx['endCoor']])] = 1
                for g in ele[idx['geneID']].rsplit(','):
                    bsj_gene[g] = 1
                cnt = int(ele[idx['readCount']])
                bsj_dict[bsj] += cnt
                if bsj not in bsj_known_novel_dict:
                    is_known = ele[idx['isKnownBSJ']]
                    if ',' not in is_known:
                        print('Two circRNA annotations needed.')
                        sys.exit(1)
                    if 'True' not in is_known:
                        bsj_known_novel_dict[bsj] = 'Novel'
                    else:
                        bsj_known_novel_dict[bsj] = 'Known'
        for bsj, cnt in bsj_dict.items():
            if cnt >= min_cnt:
                bsj_rep_cnt_dict[bsj] += 1

    return bsj_rep_cnt_dict, bsj_known_novel_dict, bsj_gene

# up +, down - => convergent
# up -, down + => divergent
def get_flank_bed(bsj_dict, flank_len, up_flank_bed, down_flank_bed):
    up_dict, down_dict = dict(), dict()
    with open(up_flank_bed, 'w') as up_fp, open(down_flank_bed, 'w') as down_fp:
        for (chrom, start, end) in bsj_dict:
            start, end = int(start), int(end)
            up_start, up_end = start - flank_len, start
            down_start, down_end = end, end + flank_len
            if up_start >= 0:
                if '{}_{}'.format(chrom, start) not in up_dict:
                    up_fp.write('{}\t{}\t{}\t{}_{}\n'.format(chrom, up_start, up_end, chrom, start))
                    up_dict['{}_{}'.format(chrom, start)] = 1
            if '{}_{}'.format(chrom, end) not in down_dict:
                down_fp.write('{}\t{}\t{}\t{}_{}\n'.format(chrom, down_start, down_end, chrom, end))
                down_dict['{}_{}'.format(chrom, end)] = 1


def get_flank_Alu(flank_out, flank_dict):
    with open(flank_out, 'r') as in_fp:
        for line in in_fp:
            if line.startswith('#'): continue
            ele = line.rsplit()
            # 3:id, 7: name, 9: strand
            flank_dict[ele[9]][ele[5]] += 1



def itst_alu_bed(bsj_dict, alu_bed, flank_len, up_flank_bed, down_flank_bed, up_flank_out, down_flank_out):
    get_flank_bed(bsj_dict, flank_len, up_flank_bed, down_flank_bed)
    
    cmd = '{} intersect -a {} -b {} -f 1.0 -wa -wb > {}'.format(bedtools, alu_bed, up_flank_bed, up_flank_out)
    print(cmd)
    os.system(cmd)
    cmd = '{} intersect -a {} -b {} -f 1.0 -wa -wb > {}'.format(bedtools, alu_bed, down_flank_bed, down_flank_out)
    print(cmd)
    os.system(cmd)
    up_flank_Alu_dict, down_flank_Alu_dict = dd(lambda: dd(lambda:0)), dd(lambda: dd(lambda:0)) # site -> '-':#, '+':#
    get_flank_Alu(up_flank_out, up_flank_Alu_dict)
    get_flank_Alu(down_flank_out, down_flank_Alu_dict)
    return up_flank_Alu_dict, down_flank_Alu_dict


def get_alu_pair_cate(bsj_dict, up_flank_Alu_dict, down_flank_Alu_dict):
    bsj_alu_pair_cate = dd(lambda:'None')
    for (chrom, start, end) in bsj_dict:
        
        n_up_for = up_flank_Alu_dict['{}_{}'.format(chrom, start)]['+']
        n_up_rev = up_flank_Alu_dict['{}_{}'.format(chrom, start)]['-']
        n_down_for = down_flank_Alu_dict['{}_{}'.format(chrom, end)]['+']
        n_down_rev = down_flank_Alu_dict['{}_{}'.format(chrom, end)]['-']

        if n_up_for > 0 and n_up_rev > 0 and n_down_for > 0 and n_down_rev > 0:
            cate = 'Both'
        elif n_up_for > 0 and n_down_rev > 0:
            cate = 'Convergent'
        elif n_up_rev > 0 and n_down_for > 0:
            cate = 'Divergent'
        else:
            cate = 'None'
        bsj_alu_pair_cate[(chrom, start, end)] = cate

    return bsj_alu_pair_cate

if __name__ == '__main__':
    if len(sys.argv) != 10:
        print('{} min_cnt flank_len rep1.out rep2.out alu.bed gene.pred circBase.bed mionco.bed fig.png'.format(sys.argv[0]))
        sys.exit(1)
    min_cnt = int(sys.argv[1])
    flank_len = int(sys.argv[2])
    rep1_out, rep2_out = sys.argv[3], sys.argv[4]
    alu_bed = sys.argv[5]
    gene_pred = sys.argv[6]
    circbase = sys.argv[7]
    mionco = sys.argv[8]
    fig = sys.argv[-1]
    tmp_dir=os.path.dirname(os.path.abspath(fig))
    bsj_alu_out = tmp_dir + '/cnt{}_flank{}_bsj_alu_pair.out'.format(min_cnt, flank_len)

    all_bsj_sites = dd(lambda : dict())
    bsj_rep_cnt_dict, bsj_known_novel_dict, bsj_genes = get_bsj_dict(min_cnt, rep1_out, rep2_out, all_bsj_sites)
    get_circ_sites(all_bsj_sites, circbase, mionco)
    all_non_bsj_sites = get_all_sites(gene_pred, bsj_genes, all_bsj_sites)
    # negative control
    neg_control_bsj = get_neg_control_non_bsj(all_non_bsj_sites, 10000)
    bsj_rep_cnt_dict = add_neg_control(neg_control_bsj, bsj_rep_cnt_dict)
    # do intersect for up and for down
    up_flank_bed, down_flank_bed = tmp_dir + '/up.bed', tmp_dir + '/down.bed'
    up_flank_out, down_flank_out = tmp_dir + '/up_flank_Alu.out', tmp_dir + '/down_flank_Alu.out'
    up_flank_Alu_dict, down_flank_Alu_dict = itst_alu_bed(bsj_rep_cnt_dict, alu_bed, flank_len, up_flank_bed, down_flank_bed, up_flank_out, down_flank_out)

    bsj_alu_pair_cate = get_alu_pair_cate(bsj_rep_cnt_dict, up_flank_Alu_dict, down_flank_Alu_dict)

    # output
    cate_out_dict = dd(lambda: dd(lambda:0)) #  {RepAnno : {Cate: num}}
    for bsj, rep_n in bsj_rep_cnt_dict.items():
        label = 'Negative' if rep_n == 0 else '{}{}'.format(rep_n, bsj_known_novel_dict[bsj])
        cate_out_dict[label][bsj_alu_pair_cate[bsj]] += 1

    with open(bsj_alu_out, 'w') as out_fp:
        tot_cnt_dict = dd(lambda:0.0)
        out_fp.write('RepAnno\tCate\tCount\tFrac\tMinCount\tFlankLen\n') # Cate: Convergent/Divergent/Both/Other
        for repanno, cate_cnt in cate_out_dict.items():
            for cate, cnt in cate_cnt.items():
                tot_cnt_dict[repanno] += cnt
        for repanno, cate_cnt in cate_out_dict.items():
            for cate, cnt in cate_cnt.items():
                out_fp.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(repanno, cate, cnt, cnt/ tot_cnt_dict[repanno], min_cnt, flank_len))

