import sys, os, re
from collections import defaultdict as dd
from copy import copy
import parse_gff as pg

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
        'Liver_1', 'Lung_1', 'Prostate_1', 'SkeletalMuscle_1', 'SmoothMuscle_1', 'Testis_1',
        'Hek293_1',    'Hek293_2',    'Hek293_3',    'Hek293_4',    'Hek293_5',    'Hek293_6'
]

idx = {h: i for i, h in enumerate(header)}

def get_exon(ele):
    exons = []
    start = int(ele[idx['startCoor0based']])
    starts = list(map(int, ele[idx['blockStarts']].rsplit(',')))
    lens = list(map(int, ele[idx['blockSize']].rsplit(',')))
    n = int(ele[idx['blockCount']])
    for i in range(n):
        exons.append((start+starts[i]+1, start+starts[i]+lens[i]))
    return n, exons

def get_inter_site(ele):
    sites = []
    start = int(ele[idx['startCoor0based']])
    starts = list(map(int, ele[idx['blockStarts']].rsplit(',')))
    lens = list(map(int, ele[idx['blockSize']].rsplit(',')))
    n = int(ele[idx['blockCount']])
    for i in range(n):
        sites.append(start + starts[i] + 1)
        sites.append(start + starts[i] + lens[i])
    sites.remove(sites[0])
    sites.remove(sites[-1])
    return sites

def get_read_cnt(ele):
    N = 18
    cnt = 0
    for i in range(N):
        cnt += int(ele[idx['downFlankAlu'] + i + 1])
    return cnt

def get_es_coors(inc_exons, exc_exons):
    for i, me in enumerate(inc_exons):
        if me != exc_exons[i]: # e1 is skipped exon
            e1 = exc_exons[i-1]
            e2 = exc_exons[i]
            return [e1[1], me[0], me[1], e2[0]], [e1[0], e2[1]], i


# exons1 = [exon1, exon2, exon3]
# exons2 = [exon1, exon2]
def is_exon_skip(exons1, exons2):
    if exons1[0][0] == exons2[0][0] and exons1[0][1] == exons2[0][1] and exons1[-1][0] == exons2[-1][0] and exons1[-1][1] == exons2[-1][1]:
        return True#[exons1[0][1], exons1[1][0], exons1[1][1], exons1[2][0]]  # exons1 is longer isoform
    else:
        return None


# exons1 = exon
# exons2 = [exon1, exon2]
def is_intron_retention(exons1, exons2):
    if exons1[0] == exons2[0][0] and exons1[1] == exons2[1][1]:
        return True
    else:
        return None


# exons1 = [exon1, exon2]
# exons2 = [exon1, exon2]
def is_alt53_splice(exons1, exons2):
    if exons1[0][1] > exons2[0][1] and exons1[0][0] == exons2[0][0] and exons1[1][0] == exons2[1][0] and exons1[1][1] == exons2[1][1]:
        return 'A5S'
    elif exons1[1][0] < exons2[1][0] and exons1[0][0] == exons2[0][0] and exons1[0][1] == exons2[0][1] and exons1[1][1] == exons2[1][1]:
        return 'A3S'
    else:
        return None


def get_bsj_dis(major, minor, bsj_dis_fp):
    strand = major[idx['canoBSJMotif']][0]

    _3_dis = abs(int(major[idx['startCoor0based']]) - int(minor[idx['startCoor0based']]))
    _5_dis = abs(int(major[idx['endCoor']]) - int(minor[idx['endCoor']]))
    if strand == '-':
        tmp = _3_dis
        _3_dis = _5_dis
        _5_dis = tmp
    if _3_dis >= 0:
        bsj_dis_fp.write('3prime\t{}\n'.format(_3_dis))
    if _5_dis >= 0:
        bsj_dis_fp.write('5prime\t{}\n'.format(_5_dis))
    if _3_dis + _5_dis >= 0:
        bsj_dis_fp.write('53prime\t{}\n'.format(_3_dis + _5_dis))


def get_fsj_cate(all_trans, gene_ids, exons, known_ss):
    if len(exons) == 1: return 'FSM'
    iso_coor = []
    for (start, end) in exons:
        iso_coor.append(start)
        iso_coor.append(end)
    int_iso_coor_str = '_'.join(map(str, iso_coor[1:-1]))
    iso_cate = None
    for gene_id in gene_ids:
        for trans in all_trans[gene_id]:
            if int_iso_coor_str in trans:
                idx = trans.index(int_iso_coor_str)
                cnt = trans[:idx].count('_')
                if cnt % 2 == 1:
                    iso_cate = 'FSM'
                    break
        if iso_cate: break
    if not iso_cate:
        iso_cate = 'NIC' if 'False' not in known_ss else 'NNC'
    return iso_cate


# ES, IR, A3S, A5S, others
# output:
#   AS_cate, gene, strand, chrom, keyCoors, majorID, minorID, majorCnt, minorCnt, major/minorInc, novelSites, ucscCoors
def get_inter_cate(all_trans, major, minor, alt_fp):
    # if not ((major[idx['FSJCate']] == 'FSM' or major[idx['FSJCate']] == 'NIC') and \
        # (minor[idx['FSJCate']] == 'FSM' or minor[idx['FSJCate']] == 'NIC')):
    # if major[idx['FSJCate']] != 'FSM' or minor[idx['FSJCate']] != 'FSM':
        # return
    cate = 'others'
    major_exon_n, major_exon = get_exon(major)
    minor_exon_n, minor_exon = get_exon(minor)
    # major_fsj, minor_fsj = major[idx['FSJCate']], minor[idx['FSJCate']]
    major_read_cnt, minor_read_cnt = get_read_cnt(major), get_read_cnt(minor)
    major_ID, minor_ID = major[0], minor[0]
    gene = major_ID.rsplit('.circRNA.')[0]
    geneIDs = major[idx['geneID']].rsplit(',')
    strand = major[idx['canoBSJMotif']][0]
    novel_flag = 'NA'
    inc = ''
    for inc, exc, inc_exon, inc_n, exc_exon, exc_n in zip([major, minor], [minor, major], [major_exon, minor_exon], [major_exon_n, minor_exon_n], [minor_exon, major_exon], [minor_exon_n, major_exon_n]):
        # check if ES
        for i in range(inc_n-2):
            inc_exons = [inc_exon[i], inc_exon[i+1], inc_exon[i+2]]
            for j in range(exc_n-1):
                exc_exons = [exc_exon[j], exc_exon[j+1]]
                if is_exon_skip(inc_exons, exc_exons):
                    cate = 'ES'
                    inc_isoform = 'Predom' if inc_exon == major_exon else 'Minor'
                    keycoors = [inc_exons[0][1], inc_exons[1][0], inc_exons[1][1], inc_exons[2][0]]
                    ucsc_coors = [inc_exons[0][0], inc_exons[2][1]]
                    novel_n = 1 if inc[idx['isKnownSS']].rsplit(',')[(i+1)*2] == 'False' else 0
                    if inc[idx['isKnownSS']].rsplit(',')[(i+1)*2+1] == 'False':
                        novel_n += 1
                    novel_flag = 'novel_{}sites'.format(novel_n) if novel_n else 'NA'
                    inc_known_ss = inc[idx['isKnownSS']].rsplit(',')[i*2+1 : (i+2)*2+1]
                    exc_known_ss = exc[idx['isKnownSS']].rsplit(',')[j*2+1 : (j+1)*2+1]
                    (major_fsj, minor_fsj) = (get_fsj_cate(all_trans, geneIDs, inc_exons, inc_known_ss), get_fsj_cate(all_trans, geneIDs, exc_exons, exc_known_ss)) if inc_isoform == 'Predom' else (get_fsj_cate(all_trans, geneIDs, exc_exons, exc_known_ss), get_fsj_cate(all_trans, geneIDs, inc_exons, inc_known_ss))
                    alt_fp.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}:{}-{}\n'.format(cate, gene, strand, major[idx['chrom']], ','.join(map(str, keycoors)), major_ID, minor_ID, major_read_cnt, minor_read_cnt, major_fsj, minor_fsj, inc_isoform, novel_flag, major[idx['chrom']], ucsc_coors[0], ucsc_coors[1]))
        # check if IR
        for i in range(inc_n):
            inc_exons = inc_exon[i]
            for j in range(exc_n-1):
                exc_exons = [exc_exon[j], exc_exon[j+1]]
                if is_intron_retention(inc_exons, exc_exons):
                    cate = 'IR'
                    inc_isoform = 'Predom' if inc_exon == major_exon else 'Minor'
                    keycoors = [exc_exons[0][1], exc_exons[1][0]]
                    ucsc_coors = [exc_exons[0][0], exc_exons[1][1]]
                    novel_n = 1 if exc[idx['isKnownSS']].rsplit(',')[j*2+1] == 'False' else 0
                    if exc[idx['isKnownSS']].rsplit(',')[(j+1)*2] == 'False':
                        novel_n += 1
                    novel_flag = 'novel_{}sites'.format(novel_n) if novel_n else 'NA'
                    exc_known_ss = exc[idx['isKnownSS']].rsplit(',')[j*2+1 : (j+1)*2+1]
                    (major_fsj, minor_fsj) = ('FSM', get_fsj_cate(all_trans, geneIDs, exc_exons, exc_known_ss)) if inc_isoform == 'Predom' else (get_fsj_cate(all_trans, geneIDs, exc_exons, exc_known_ss), 'FSM')
                    alt_fp.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}:{}-{}\n'.format(cate, gene, strand, major[idx['chrom']], ','.join(map(str, keycoors)), major_ID, minor_ID, major_read_cnt, minor_read_cnt, major_fsj, minor_fsj,inc_isoform, novel_flag, major[idx['chrom']], ucsc_coors[0], ucsc_coors[1]))
        
        # check if +Alt 5' or -Alt 3'
        for i in range(inc_n-1):
            inc_exons = [inc_exon[i], inc_exon[i+1]]
            for j in range(exc_n-1):
                exc_exons = [exc_exon[j], exc_exon[j+1]]
                cate = is_alt53_splice(inc_exons, exc_exons)
                if cate:
                    inc_isoform = 'Predom' if inc_exon == major_exon else 'Minor'
                    if cate == 'A5S':
                        keycoors = [exc_exons[0][1], inc_exons[0][1], inc_exons[1][0]]
                        novel_n = 1 if exc[idx['isKnownSS']].rsplit(',')[j*2+1] == 'False' else 0
                        if inc[idx['isKnownSS']].rsplit(',')[i*2+1] == 'False':
                            novel_n += 1
                    else:
                        keycoors = [inc_exons[0][1], inc_exons[1][0], exc_exons[1][0]]
                        novel_n = 1 if inc[idx['isKnownSS']].rsplit(',')[(i+1)*2] == 'False' else 0
                        if exc[idx['isKnownSS']].rsplit(',')[(j+1)*2] == 'False':
                            novel_n += 1
                    ucsc_coors = [inc_exons[0][0], inc_exons[1][1]]
                    novel_flag = 'novel_{}sites'.format(novel_n) if novel_n else 'NA'
                    inc_known_ss = inc[idx['isKnownSS']].rsplit(',')[i*2+1 : (i+1)*2+1]
                    exc_known_ss = exc[idx['isKnownSS']].rsplit(',')[j*2+1 : (j+1)*2+1]
                    if strand == '-':
                        cate = 'A3S' if cate == 'A5S' else 'A5S'
                    (major_fsj, minor_fsj) = (get_fsj_cate(all_trans, geneIDs, inc_exons, inc_known_ss), get_fsj_cate(all_trans, geneIDs, exc_exons, exc_known_ss)) if inc_isoform == 'Predom' else (get_fsj_cate(all_trans, geneIDs, exc_exons, exc_known_ss), get_fsj_cate(all_trans, geneIDs, inc_exons, inc_known_ss))
                    alt_fp.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}:{}-{}\n'.format(cate, gene, strand, major[idx['chrom']], ','.join(map(str, keycoors)), major_ID, minor_ID, major_read_cnt, minor_read_cnt, major_fsj, minor_fsj, inc_isoform, novel_flag, major[idx['chrom']], ucsc_coors[0], ucsc_coors[1]))


def cal_diff_cate(lines, all_trans, alt_fp, bsj_dis_fp):
    diff_dict = dd(lambda: 0)
    if len(lines) < 2: return diff_dict

    major = None
    for ele in lines:
        # search for major isoform
        if ele[idx['isoformID']].endswith('.1'):
            major = ele
            break
    if not major:
        print('No major in {}'.format(ele[idx['isoformID']]))
        return diff_dict
    major_inter_site = get_inter_site(major)
    for ele in lines:
        if ele[idx['isoformID']].endswith('.1'): continue
        if major[idx['startCoor0based']] == ele[idx['startCoor0based']] and major[idx['endCoor']] == ele[idx['endCoor']]:
            diff_dict['diff_inter'] += 1
            get_inter_cate(all_trans, major, ele, alt_fp)
        else:
            minor_inter_site = get_inter_site(ele)
            # all internal same
            if major_inter_site == minor_inter_site:
                diff_dict['diff_bsj'] += 1
            else:
                diff_dict['diff_both'] += 1
                get_inter_cate(all_trans, major, ele, alt_fp)
            get_bsj_dis(major, ele, bsj_dis_fp)
    return diff_dict


def alt_splice_analysis(comb_out, all_trans, alt_splice_out, bsj_dis_out):
    with open(comb_out) as isocirc_fp, open(alt_splice_out, 'w') as alt_fp, open(bsj_dis_out, 'w') as bsj_dis_fp:
        last_gene = ''
        diff_cate_dict = dd(lambda: 0)
        lines = []

        alt_fp.write('Cate\tGene\tStrand\tChrom\tKeyCoors\tPredomID\tMinorID\tPredomCnt\tMinorCnt\tPredomFSJ\tMinorFSJ\tIncIsoform\tNovelSites\tUCSCCoors\n')
        bsj_dis_fp.write('Type\tDistance\n')  # Type: 5prime/3prime/53pirme

        for line in isocirc_fp:
            if line.startswith('NA.') or line.startswith('isoformID'): # or line.startswith('#'):
                continue
            ele = line.rsplit()
            # if ',' in ele[0]: continue # multip genes
            gene = ele[idx['isoformID']].rsplit('.circRNA.')[0]
            if gene == last_gene:
                lines.append(ele)
            else:
                diff_dict = cal_diff_cate(lines, all_trans, alt_fp, bsj_dis_fp)
                for c, n in diff_dict.items():
                    diff_cate_dict[c] += n
                lines = [ele]
                last_gene = gene

        diff_dict = cal_diff_cate(lines, all_trans, alt_fp, bsj_dis_fp)
        for c, n in diff_dict.items():
            diff_cate_dict[c] += n

        print('\nInternal: {}\tBSJ: {}\tBoth: {}\n'.format(diff_cate_dict['diff_inter'], diff_cate_dict['diff_bsj'], diff_cate_dict['diff_both']))

if __name__ == '__main__':
    if len(sys.argv) != 5:
        print('Usage:')
        print('{} comb.out gene.pred alt.out bsj_dis.out'.format(sys.argv[0]))
        sys.exit(1)

    comb_out = sys.argv[1]
    all_trans = pg.get_transcript_from_gene_pred(sys.argv[2])
    alt_out = sys.argv[3]
    bsj_dis = sys.argv[4]
    alt_splice_analysis(comb_out, all_trans, alt_out, bsj_dis)

# 5. R plot venn diagram to fig
# cmd='Rscript /home/gaoy1/program/circ_plot/RI_lenBox.R {} {}'.format(ri_dat, fig)
# print(cmd)
# os.system(cmd)