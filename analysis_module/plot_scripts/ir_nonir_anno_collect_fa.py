import sys, os, re
from collections import defaultdict as dd
from pyfaidx import Fasta
from copy import copy
import mappy as mp

fxtools = 'fxtools'

ir_header = ['Cate','Gene', 'Strand', 'Chrom','KeyCoors','MajorID','MinorID','MajorCount','MinorCount','Inclusion','NovelSite','UCSCCoors']
ir_idx = {h: i for i, h in enumerate(ir_header)}

intron_header = ['chrom', 'start0based', 'end', 'intronName', 'color', 'strand']
intron_idx = {h: i for i, h in enumerate(intron_header)}

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

_5_perl='/home/gaoy1/software/maxent/score5.pl'
_3_perl='/home/gaoy1/software/maxent/score3.pl'

def get_anno_intron(intron_bed):
    anno_3_sites = dict()
    anno_5_sites = dict()
    anno_intron = dict()
    with open(intron_bed) as fp:
        for line in fp:
            ele = line.rsplit()
            chrom = ele[intron_idx['chrom']]
            strand = ele[intron_idx['strand']]
            if ele[intron_idx['strand']] == '+':
                anno_5_sites[(chrom, strand, int(ele[intron_idx['start0based']]))] = 1
                anno_3_sites[(chrom, strand, int(ele[intron_idx['end']])+1)] = 1
            else:
                anno_5_sites[(chrom, strand, int(ele[intron_idx['end']])+1)] = 1
                anno_3_sites[(chrom, strand, int(ele[intron_idx['start0based']]))] = 1
            anno_intron[(chrom, strand, int(ele[intron_idx['start0based']]), int(ele[intron_idx['end']])+1)] = 1
    return anno_intron # anno_5_sites, anno_3_sites


# def get_ir_gene_and_sites(ir_list, anno_5_sites, anno_3_sites):
def get_ir_gene_and_sites(ir_list, anno_intron):
    ir_gene_dict = dict()
    # ir_5_sites = dict()
    # ir_3_sites = dict()
    ir = dict()
    with open(ir_list) as fp:
        for line in fp:
            ele = line.rsplit()
            if ele[0] != 'IR': continue
            gene = ele[ir_idx['Gene']]
            ir_gene_dict[gene] = 1
            chrom, strand = ele[ir_idx['Chrom']], ele[ir_idx['Strand']]
            key_coors = list(map(int, ele[ir_idx['KeyCoors']].rsplit(',')))
            # if strand.startswith('+'): 
            #     _5_site = (chrom, strand, key_coors[0])
            #     _3_site = (chrom, strand, key_coors[1])
            # else:
            #     _5_site = (chrom, strand, key_coors[1])
            #     _3_site = (chrom, strand, key_coors[0])
            # ir_5_sites[_5_site] = _5_site in anno_5_sites
            # ir_3_sites[_3_site] = _3_site in anno_3_sites
            ir[(chrom, strand, key_coors[0], key_coors[1])] = (chrom, strand, key_coors[0], key_coors[1]) in anno_intron
    # return ir_gene_dict, ir_5_sites, ir_3_sites
    return ir_gene_dict, ir


# def collect_nonir_sites(comb_out, ir_genes, anno_5_sites, anno_3_sites, ir_5_sites, ir_3_sites):
def collect_nonir_sites(comb_out, ir_genes, anno_intron, ir):
    # nonir_5_sites = dict()
    # nonir_3_sites = dict()
    nonir = dict()
    with open(comb_out) as fp:
        for line in fp:
            ele = line.rsplit()
            if ele[0] == 'isoformID':
                continue
            genenames = ele[idx['geneName']].rsplit()
            hit = -1
            for i, gene in enumerate(genenames):
                if gene in ir_genes:
                    hit = i
                    break
            if hit == -1:
                continue
            blockCount = int(ele[idx['blockCount']])
            if blockCount == 1:
                continue
            start_coor = int(ele[idx['startCoor0based']])
            sizes = list(map(int, ele[idx['blockSize']].rsplit(',')))
            starts = list(map(int, ele[idx['blockStarts']].rsplit(',')))
            chrom = ele[idx['chrom']]
            strand = ele[idx['geneStrand']].rsplit(',')[hit]
            # if strand == '+':
            #     for i in range(blockCount-1):
            #         _5_site = (chrom, strand, start_coor + sizes[i])
            #         if _5_site not in ir_5_sites:
            #             nonir_5_sites[_5_site] = _5_site in anno_5_sites
            #         _3_site = (chrom, strand, start_coor + starts[i+1]+1)
            #         if _3_site not in ir_3_sites:
            #             nonir_3_sites[_3_site] = _3_site in anno_3_sites
            # else:
            #     for i in range(blockCount-1):
            #         _3_site = (chrom, strand, start_coor + sizes[i])
            #         if _3_site not in ir_3_sites:
            #             nonir_3_sites[_3_site] = _3_site in anno_3_sites
            #         _5_site = (chrom, strand, start_coor + starts[i+1]+1)
            #         if _5_site not in ir_5_sites:
            #             nonir_5_sites[_5_site] = _5_site in anno_5_sites
            for i in range(blockCount-1):
                intron = (chrom, strand, start_coor + starts[i] + sizes[i], start_coor + starts[i+1]+1)
                if intron not in ir:
                    nonir[intron] = intron in anno_intron
    # return nonir_5_sites, nonir_3_sites
    return nonir


def get_seqs(ref_fa, chrom, strand, up_site, down_site):
    _5_name, _3_name, _5_seq, _3_seq = '', '', '', ''
    ref_seq = ref_fa[chrom]
    if strand == '+':
        _5_name = '{}:{}-{} {}'.format(chrom, up_site-2, up_site+6, strand)
        _5_seq = ref_seq[up_site-3:up_site+6].seq.upper()
        _3_name = '{}:{}-{} {}'.format(chrom, down_site-20, down_site+2, strand)
        _3_seq = ref_seq[down_site-21:down_site+2].seq.upper()
    else: # '-'
        _5_name = '{}:{}-{} {}'.format(chrom, down_site-6, down_site+2, strand)
        _5_seq = ref_seq[down_site-7:down_site+2].seq.upper()
        _5_seq = mp.revcomp(_5_seq)
        _3_name = '{}:{}-{} {}'.format(chrom, up_site-2, up_site+20, strand)
        _3_seq = ref_seq[up_site-3:up_site+20].seq.upper()
        _3_seq = mp.revcomp(_3_seq)
    if _5_seq[3:5] != 'GT' or _3_seq[-5:-3] != 'AG':
        _5_name = ''
    return _5_name, _5_seq, _3_name, _3_seq


def old_get_site_fa(ref_fa, out_pre, ir_5_sites, ir_3_sites, nonir_5_sites, nonir_3_sites):
    ir_5_fa, ir_3_fa, nonir_5_fa, nonir_3_fa = out_pre + 'ir.5.fa', out_pre + 'ir.3.fa', out_pre + 'nonir.5.fa', out_pre + 'nonir.3.fa'
    anno_ir_5_fa, anno_ir_3_fa, anno_nonir_5_fa, anno_nonir_3_fa = out_pre + 'anno.ir.5.fa', out_pre + 'anno.ir.3.fa', out_pre + 'anno.nonir.5.fa', out_pre + 'anno.nonir.3.fa'
    out_fps = []
    anno_fps = []
    out_fas = [ir_5_fa, ir_3_fa, nonir_5_fa, nonir_3_fa]
    anno_fas = [anno_ir_5_fa, anno_ir_3_fa, anno_nonir_5_fa, anno_nonir_3_fa]
    for out_fa, anno_fa in zip(out_fas, anno_fas):
        out_fps.append(open(out_fa, 'w'))
        anno_fps.append(open(anno_fa, 'w'))
    sites = [ir_5_sites, ir_3_sites, nonir_5_sites, nonir_3_sites]
    is_5s = [True, False, True, False]
    for out_fp, anno_fp, sites, is_5 in zip(out_fps, anno_fps, sites, is_5s):
        for site, is_anno in sites.items():
            name, seq = get_seqs(ref_fa, site, is_5)
            if not name: 
                continue
            if is_anno:
                anno_fp.write('>{}\n{}\n'.format(name, seq))
            out_fp.write('>{}\n{}\n'.format(name, seq))

    for _5_fa in [ir_5_fa, nonir_5_fa, anno_ir_5_fa, anno_nonir_5_fa]:
        cmd = 'perl {} {} > {}.score'.format(_5_perl, _5_fa, _5_fa)
        print(cmd)
    
    for _3_fa in [ir_3_fa, nonir_3_fa, anno_ir_3_fa, anno_nonir_3_fa]:
        cmd = 'perl {} {} > {}.score'.format(_3_perl, _3_fa, _3_fa)
        print(cmd)



def get_site_fa(ref_fa, out_pre, ir, nonir):
    ir_5_fa, ir_3_fa, nonir_5_fa, nonir_3_fa = out_pre + 'ir.5.fa', out_pre + 'ir.3.fa', out_pre + 'nonir.5.fa', out_pre + 'nonir.3.fa'
    anno_ir_5_fa, anno_ir_3_fa, anno_nonir_5_fa, anno_nonir_3_fa = out_pre + 'anno.ir.5.fa', out_pre + 'anno.ir.3.fa', out_pre + 'anno.nonir.5.fa', out_pre + 'anno.nonir.3.fa'
    n_ir, n_nonir, n_anno_ir, n_anno_nonir = 0, 0, 0, 0
    _5_dict, _3_dict = dict(), dict()
    out_fps, anno_fps = [], []
    out_fas = [ir_5_fa, ir_3_fa, nonir_5_fa, nonir_3_fa]
    anno_fas = [anno_ir_5_fa, anno_ir_3_fa, anno_nonir_5_fa, anno_nonir_3_fa]
    for out_fa, anno_fa in zip(out_fas, anno_fas):
        out_fps.append(open(out_fa, 'w'))
        anno_fps.append(open(anno_fa, 'w'))

    for intron, is_anno in ir.items():
        (chrom, strand, _up_site, _down_site) = intron
        _5_name, _5_seq, _3_name, _3_seq = get_seqs(ref_fa, chrom, strand, _up_site, _down_site)
        
        if _5_name:
            n_ir += 1
            if is_anno:
                n_anno_ir += 1
                if _5_name not in _5_dict:
                    anno_fps[0].write('>{}\n{}\n'.format(_5_name, _5_seq))
                if _3_name not in _3_dict:
                    anno_fps[1].write('>{}\n{}\n'.format(_3_name, _3_seq))
            if _5_name not in _5_dict:
                out_fps[0].write('>{}\n{}\n'.format(_5_name, _5_seq))
            if _3_name not in _3_dict:
                out_fps[1].write('>{}\n{}\n'.format(_3_name, _3_seq))
            _5_dict[_5_name] = 1
            _3_dict[_3_name] = 1

    for intron, is_anno in nonir.items():
        (chrom, strand, _up_site, _down_site) = intron
        _5_name, _5_seq, _3_name, _3_seq = get_seqs(ref_fa, chrom, strand, _up_site, _down_site)
        if _5_name:
            n_nonir += 1
            if is_anno:
                n_anno_nonir += 1
                if _5_name not in _5_dict:
                    anno_fps[2].write('>{}\n{}\n'.format(_5_name, _5_seq))
                if _3_name not in _3_dict:
                    anno_fps[3].write('>{}\n{}\n'.format(_3_name, _3_seq))
            if _5_name not in _5_dict:
                out_fps[2].write('>{}\n{}\n'.format(_5_name, _5_seq))
            if _3_name not in _3_dict:
                out_fps[3].write('>{}\n{}\n'.format(_3_name, _3_seq))
            _5_dict[_5_name] = 1
            _3_dict[_3_name] = 1

    print('# RI: {}\t#Non-RI: {}'.format(n_ir, n_nonir))
    print('# Anno-RI: {}\t#Anno-Non-RI: {}'.format(n_anno_ir, n_anno_nonir))
    
    for _5_fa in [ir_5_fa, nonir_5_fa, anno_ir_5_fa, anno_nonir_5_fa]:
        cmd = 'perl {} {} > {}.score'.format(_5_perl, _5_fa, _5_fa)
        print(cmd)
    
    for _3_fa in [ir_3_fa, nonir_3_fa, anno_ir_3_fa, anno_nonir_3_fa]:
        cmd = 'perl {} {} > {}.score'.format(_3_perl, _3_fa, _3_fa)
        print(cmd)


if __name__ == '__main__':
    if len(sys.argv) != 6:
        print('{} ref.fa IR.list comb.isocirc.out intron.bed out_pre'.format(sys.argv[0]))
        sys.exit(1)
    
    [ref_fa, ir_list, comb_out, intron_bed, out_pre] = sys.argv[1:]
    
    # anno_5_sites, anno_3_sites = get_anno_intron(intron_bed)
    anno_intron = get_anno_intron(intron_bed)
    # ir_gene_dict, ir_5_sites, ir_3_sites = get_ir_gene_and_sites(ir_list, anno_5_sites, anno_3_sites)
    ir_gene_dict, ir = get_ir_gene_and_sites(ir_list, anno_intron)
    
    # nonir_5_sites, nonir_3_sites = collect_nonir_sites(comb_out, ir_gene_dict, anno_5_sites, anno_3_sites, ir_5_sites, ir_3_sites)
    nonir = collect_nonir_sites(comb_out, ir_gene_dict, anno_intron, ir)
    # get_site_fa(Fasta(ref_fa), out_pre, ir_5_sites, ir_3_sites, nonir_5_sites, nonir_3_sites)
    get_site_fa(Fasta(ref_fa), out_pre, ir, nonir)



    