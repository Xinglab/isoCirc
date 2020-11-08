import sys, os, re
import math
from itertools import combinations as cmb
from collections import defaultdict as dd
from scipy.stats import chi2_contingency as cc
from scipy.stats import fisher_exact as fe
import statsmodels.stats.multitest as smt
from rpy2 import robjects as rb
from rpy2.robjects import r
from rpy2.robjects.packages import importr
stats = importr('stats')

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
raw_names = ['Adipose_1', 'Adrenal_1', 'Blood_1', 'Brain_1', 'Heart_1', 'Kidney_1',
         'Liver_1', 'Lung_1', 'Prostate_1', 'SkeletalMuscle_1', 'SmoothMuscle_1', 'Testis_1']
names = ['Adipose', 'Adrenal', 'Blood', 'Brain', 'Heart', 'Kidney',
         'Liver', 'Lung', 'Prostate', 'SkeletalMuscle', 'SmoothMuscle', 'Testis']

name_dict = {r:n for r,n in zip(raw_names, names)}

fdr_thres=0.05
p_thres = 0.05
delta_ratio = 0.05

def get_chisq_test_p_value(cnt=[]):
    cnt1, cnt2 = cnt[0], cnt[1]
    r1 = rb.IntVector(cnt1)
    r2 = rb.IntVector(cnt2)
    df = r.rbind(r1, r2)
    out = stats.chisq_test(df)
    return list(out[2])[0]

def get_fisher_test_p_value(cnt=[]):
    cnt1, cnt2 = cnt[0], cnt[1]
    r1 = rb.IntVector(cnt1)
    r2 = rb.IntVector(cnt2)
    df = r.rbind(r1, r2)
    out = stats.fisher_test(df)
    return list(out[0])[0]

def get_fdr(ps=[]):
    rp = rb.FloatVector(ps)
    return list(stats.p_adjust(rp, method='fdr'))

def specific_isoform(min_read_cnt=2, min_iso_cnt=2, in_fn='', names=[], gene_fdr_out='', iso_p_out='', iso_cnt_fn='', ratio_fp=None):
    iso_cnt_in_each_gene = dd(lambda : dd(lambda: dd(lambda: 0)))
    gene_to_isos, iso_to_genes = dd(lambda:dict()), dd(lambda:dict())
    gene_to_max_ratio = dd(lambda:0.0)
    gene_to_fdrs = dd(lambda:dd(lambda:1)) # p:1, fdr:1

    with open(in_fn) as in_fp:
        first_line = True
        for line in in_fp:
            if first_line:
                first_line = False
                continue
            if line.startswith('NA.'): continue
            ele = line.rsplit()
            iso = (ele[idx['chrom']], ele[idx['startCoor0based']], ele[idx['endCoor']], ele[idx['blockCount']], ele[idx['blockSize']], ele[idx['blockStarts']])
            gene = ele[idx['geneName']]
            for g in gene.rsplit(','):
                g = g.rsplit('.')[0]
                gene_to_isos[g][iso] = 0
                iso_to_genes[iso][g] = 0
                for name in names:
                    iso_cnt_in_each_gene[name][g][iso] += int(ele[idx[name]])

    genes = list(gene_to_isos.keys())
    for gene in genes:
        iso_dict = gene_to_isos[gene]
        isos = list(iso_dict.keys())
        if len(isos) < min_iso_cnt: 
            del gene_to_isos[gene]
            continue
        for iso in isos:
            filtered = 1
            for name in names:
                if iso_cnt_in_each_gene[name][gene][iso] >= min_read_cnt:
                    filtered = 0
            if filtered: del iso_dict[iso] # filtered out iso
    gene_list, p_list = [], []
    genes = list(gene_to_isos.keys())
    for gene in genes:
        iso_dict = gene_to_isos[gene]
        isos = list(iso_dict.keys())
        if len(isos) < min_iso_cnt:
            del gene_to_isos[gene]
            continue
        cnts = []
        for name in names:
            cnt = []
            for iso in isos:
                cnt.append(iso_cnt_in_each_gene[name][gene][iso])
            cnts.append(cnt)
        # delta_psi
        rat1 = [x / sum(cnts[0]) if sum(cnts[0]) > 0 else 0 for x in cnts[0]]
        rat2 = [x / sum(cnts[1]) if sum(cnts[1]) > 0 else 0 for x in cnts[1]]
        gene_to_max_ratio[gene] = max([abs(x-y) for x,y in zip(rat1, rat2)])
        p = get_chisq_test_p_value(cnts)
        gene_to_fdrs[gene]['p'] = p
        gene_list.append(gene)
        p_list.append(p)
    gene_fdrs = get_fdr(p_list)

    diff_gene_n = 0
    iso_ps, iso_list = [], []
    for gene, fdr in zip(gene_list, gene_fdrs):
        gene_to_fdrs[gene]['fdr'] = fdr
        if fdr > fdr_thres or gene_to_max_ratio[gene] < delta_ratio or math.isnan(fdr): continue
        # if fdr > fdr_thres or math.isnan(fdr): continue
        diff_gene_n += 1
        tot_cnt = dd(lambda:0.0)
        for name in names:
            for iso in gene_to_isos[gene]:
                tot_cnt[name] += iso_cnt_in_each_gene[name][gene][iso]
        for iso in gene_to_isos[gene]:
            cnts = []
            ratios = []
            for name in names:
                cnt = [iso_cnt_in_each_gene[name][gene][iso],tot_cnt[name] - iso_cnt_in_each_gene[name][gene][iso]]
                cnts.append(cnt)
                if tot_cnt[name] == 0: ratios.append(0)
                else: ratios.append(iso_cnt_in_each_gene[name][gene][iso]/tot_cnt[name])
            if abs(ratios[0]-ratios[1]) < delta_ratio: continue
            p = get_fisher_test_p_value(cnts)
            if math.isnan(p): continue
            iso_ps.append(p)
            iso_list.append(iso)
            iso_to_genes[iso][gene] = 1

    diff_iso_n, diff_iso_list, diff_iso_genes = 0, dict(), dict()
    cnt_fp = open(iso_cnt_fn, 'w') if iso_cnt_fn else None
    if cnt_fp: cnt_fp.write('Isoform\tP-value\tGene\tTissue1\tTissue2\tReadCnt\tIsoCnt\tCount1\tRatio1\tCount2\tRatio2\n')
    for iso_p, iso in zip(iso_ps, iso_list):
        if iso_p > p_thres or iso in diff_iso_list: continue
        diff_iso_n += 1
        diff_iso_list[iso] = 1
        for g in iso_to_genes[iso]:
            if iso_to_genes[iso][g] == 1:
                diff_iso_genes[g] = 1
                if cnt_fp:
                    cnt1, cnt2 = iso_cnt_in_each_gene[names[0]][g][iso], iso_cnt_in_each_gene[names[1]][g][iso]
                    tot_cnt1, tot_cnt2 = 0.0, 0.0
                    for iso1 in gene_to_isos[g]:
                        tot_cnt1 += iso_cnt_in_each_gene[names[0]][g][iso1]
                        tot_cnt2 += iso_cnt_in_each_gene[names[1]][g][iso1]
                    ratio1 = 0 if tot_cnt1 == 0 else cnt1 / tot_cnt1
                    ratio2 = 0 if tot_cnt2 == 0 else cnt2 / tot_cnt2
                    if abs(ratio1-ratio2) < delta_ratio:
                        print(names, g, iso)
                        sys.exit(1)
                    cnt_fp.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('_'.join(iso), iso_p, g, name_dict[names[0]], name_dict[names[1]],min_read_cnt, min_iso_cnt, cnt1, ratio1, cnt2, ratio2))
                    if ratio_fp:
                        ratio_fp.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(name_dict[names[0]], name_dict[names[1]],min_read_cnt, min_iso_cnt, cnt1, ratio1, cnt2, ratio2))
                        ratio_fp.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(name_dict[names[1]], name_dict[names[0]],min_read_cnt, min_iso_cnt, cnt2, ratio2, cnt1, ratio1))
    if cnt_fp: cnt_fp.close()

    # XXX for now, all the genes and isos are output 
    # if gene_fdr_out and iso_p_out:
    #     with open(gene_fdr_out, 'w') as gene_fp, open(iso_p_out, 'w') as iso_fp :
    #         gene_fp.write('Gene\tFDR\tP-value\tCnt1\tCnt2\n')
    #         iso_fp.write('Isoform\tP-value\tCnt1\tCnt2\tGene\tGene-FDR\tGene-P\n')
    #         for gene, isos in gene_to_isos.items():
    #             gene_fp.write('{}\t{}\t{}'.format(gene, gene_to_fdrs[gene]['fdr'], gene_to_fdrs[gene]['p']))
    #             for name in names:
    #                 cnt = []
    #                 for iso in isos:
    #                     cnt.append(iso_cnt_in_each_gene[name][gene][iso])
    #                 gene_fp.write('\t{}'.format('_'.join(map(str, cnt))))
    #             gene_fp.write('\n')
    #         for iso_p, iso in zip(iso_ps, iso_list):
    #             iso_fp.write('{}\t{}'.format('_'.join(iso), iso_p))
    #             g = list(iso_to_genes[iso].keys())[0]
    #             for name in names:
    #                 iso_fp.write('\t{}'.format(iso_cnt_in_each_gene[name][g][iso]))
    #             iso_fp.write('\t{}\t{}\t{}\n'.format(','.join(iso_to_genes[iso]), gene_to_fdrs[g]['fdr'], gene_to_fdrs[g]['p']))

    print('\n=={}==v.s.=={}=='.format(name_dict[names[0]], name_dict[names[1]]))
    print('FDR <= {}, read_count >= {}, isoform count >= {}'.format(fdr_thres, min_read_cnt, min_iso_cnt))
    print('Total diff circRNA genes: {}'.format(diff_gene_n))   
    print('P-value <= {}'.format(p_thres))
    print('Total diff circRNA isoforms: {}, from {} genes'.format(diff_iso_n, len(diff_iso_genes)))   
    return diff_gene_n, diff_iso_n

def pair_specific_circRNA(in_fn='', min_read_cnt=2, min_iso_cnt=2, spe_out='', all_iso_ratio_out='', tmp_dir=''):
    raw_name_cmb = list(cmb(raw_names, 2))
    with open(spe_out, 'w') as out_fp, open(all_iso_ratio_out, 'w') as ratio_fp:
        out_fp.write('ReadCnt\tIsoCnt\tVar1\tVar2\tCate\tValue\n')

        ratio_fp.write('Tissue1\tTissue2\tReadCnt\tIsoCnt\tCount1\tRatio1\tCount2\tRatio2\n')
        for name in names:
            ratio_fp.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(name, name, min_read_cnt, min_iso_cnt, 'NA', 'NA', 'NA', 'NA'))
        for raw_name in raw_name_cmb:
            name= [name_dict[raw_name[0]], name_dict[raw_name[1]]]
            tissue_spe_gene_fn = tmp_dir + '/{}_{}.spe.gene'.format(name[0], name[1])
            tissue_spe_iso_fn = tmp_dir + '/{}_{}.spe.iso'.format(name[0], name[1])
            iso_cnt_fn = tmp_dir + '/{}_{}.spe_iso.cnt'.format(name[0], name[1])
            spe_gene_n, spe_iso_n = specific_isoform(min_read_cnt, min_iso_cnt, in_fn, raw_name, tissue_spe_gene_fn, tissue_spe_iso_fn, iso_cnt_fn, ratio_fp)
            out_fp.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(min_read_cnt, min_iso_cnt, name[0], name[1], 'Gene', spe_gene_n))
            out_fp.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(min_read_cnt, min_iso_cnt, name[1], name[0], 'Gene', spe_gene_n))
            out_fp.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(min_read_cnt, min_iso_cnt, name[0], name[1], 'Isoform', spe_iso_n))
            out_fp.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(min_read_cnt, min_iso_cnt, name[1], name[0], 'Isoform', spe_iso_n))
        for name in names:
            out_fp.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(min_read_cnt, min_iso_cnt, name, name, 'Gene', 0))
            out_fp.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(min_read_cnt, min_iso_cnt, name, name, 'Isoform', 0))


if __name__ == '__main__':
    if len(sys.argv) != 6:
        print('Usage:')
        print('{} min_read_cnt min_iso_cnt_per_gene comb.out circSpeGeneIso.png circSpePairwise.png'.format(sys.argv[0]))
        sys.exit(1)

    min_read_cnt = int(sys.argv[1])
    min_iso_cnt = int(sys.argv[2])
    in_fn = sys.argv[3]
    fig1, fig2 = sys.argv[-2], sys.argv[-1]
    tmp_dir = os.path.dirname(os.path.abspath(fig1))
    gene_fdr_out = tmp_dir + '/gene_fdr.out'
    iso_p_out = tmp_dir + '/iso_p.out'
    spe_out = tmp_dir + '/spe.out'
    iso_ratio_out = tmp_dir + '/iso_ratio.out'


    pair_specific_circRNA(in_fn, min_read_cnt, min_iso_cnt, spe_out, iso_ratio_out, tmp_dir)

    # specific_isoform(min_read_cnt, min_iso_cnt, in_fp_list, names, gene_fdr_out, iso_p_out)
    cmd='Rscript /home/gaoy1/program/circ_plot/circSpecificIso.R {} {}'.format(spe_out, fig1)
    print(cmd)
    os.system(cmd)
    cmd='Rscript /home/gaoy1/program/circ_plot/circTissueFullIsoRepPair.R {} {}'.format(iso_ratio_out, fig2)
    print(cmd)
    os.system(cmd)