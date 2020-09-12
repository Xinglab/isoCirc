import argparse
import re,sys,os
import threading
from threading import Thread
from collections import defaultdict as dd
import pysam as ps
from pyfaidx import Fasta
from Bio.Seq import Seq

import isocirc.parse_bam as pb
import isocirc.parse_gff as pg
import isocirc.uniq_isoform as ui
import isocirc.basic_stats as bs
import isocirc.utils as ut
from isocirc.__init__ import __program__
from isocirc.__init__ import __version__
from isocirc.__init__ import whole_output_header
from isocirc.__init__ import whole_output_header_idx
from isocirc.__init__ import isoform_output_header
from isocirc.__init__ import isoform_output_header_idx

threads = 8
flank_len = 500
site_dis = 0
end_dis = 10
bedtools = 'bedtools'
bed2exonGtf = 'bed2exonGtf'  # dir_path + '/bin/bed2exonGtf'
itst_gtf_bed = 'itst_gtf_bed'  # dir_path + '/bin/itst_gtf_bed'
itst_gtf_gtf = 'itst_gtf_gtf'  # dir_path + '/bin/itst_gtf_gtf'
gtf2gene = 'gtf2gene'  # dir_path + '/bin/gtf2gene'
gtfToGenePred = 'gtfToGenePred'
genePredToBed = 'genePredToBed'
gtf2bed = 'gtf2bed'  # dir_path + '/bin/gtf2bed'


cigar_op_dict = pb.cigar_op_dict

itst_bed_name = ['CDS', 'UTR', 'lincRNA', 'antisense', 'rRNA', 'Alu', 'allRepeat']


def output_header(out_fp, herder):
    # version information
    out_fp.write('#' + __program__ + '\t' + __version__ + '\n')
    out_fp.write('\t'.join(herder) + '\n')

def output_whole_eval(out_fp, all_out, start_id, end_id):
    for id in range(start_id, end_id):
        eval_out = all_out[id]
        out_array = []
        for h in whole_output_header:
            out_array.append(str(eval_out[h]))
        out_fp.write('\t'.join(out_array) + '\n')


def circRNA_isoform_bed12(circRNA_bed, circRNA_out, iso_to_name_dict):
    with open(circRNA_bed, 'w') as bed_fp:
        bed_fp.write('#' + __program__ + '\t' + __version__ + '\n')
        for iso_id, read_ids in iso_to_name_dict.items():
            out1 = circRNA_out[read_ids[0]]
            bed_fp.write('{}\t{}\t{}\t{}{}\t0\t{}\t0\t0\t0\t{}\t{}\t{}\n'.format(out1['chrom'], out1['startCoor0based'], out1['endCoor'], __program__, iso_id, out1['canoBSJMotif'][0], out1['blockCount'], out1['blockSize'], out1['blockStarts']))


def is_fullLength(out1, all_out, read_ids):
    if int(out1['blockCount']) == 1:
            out1['isFullLength'] = 'True'
    else:
        out1['isFullLength'] = 'True'
        bsj_strand = out1['canoBSJMotif'][0]
        for sj_motif in out1['canoFSJMotif'].rsplit(','):
            if sj_motif[0] != bsj_strand:
                out1['isFullLength'] = 'False'
                return
        # all FSJ has canonical motifs and same strand with BSJ
        # if any FSJ is NOT known and NOT high-quality: False

        known_sj = []
        known_fsj_ss = out1['isKnownSS'].rsplit(',')[1:-1]
        for i in range(int(len(known_fsj_ss) / 2)):
            known_sj.append(','.join([known_fsj_ss[i*2], known_fsj_ss[i*2+1]]))
        cano_sj = out1['isCanoFSJ'].rsplit(',')
        high_sj = ['False'] * len(cano_sj)
        for read_id in read_ids:
            for i, high1 in enumerate(all_out[read_id]['isHighFSJ'].rsplit(',')):
                if high1 == 'True':
                    high_sj[i] = 'True'
        if len(known_sj) != len(high_sj):
            ut.fatal_format_time('output_isoform_eval', 'Unmatched list size.')
        for k, h in zip(known_sj, high_sj):
            if 'False' in k and 'False' in h:
                out1['isFullLength'] = 'False'
                break


# add 'isFullLength', 'isFSM', 'isNIC', 'isNNC', # FSM: full splice match, NIC: novel in catalog, NNC, novel and not in catalog
def output_isoform_eval(out_fp, all_out, itst_out_dict, all_trans, iso_to_name_dict={}):
    ut.err_format_time('output_isoform_eval', 'Writing isoform-wise evaluation result to file ...')
    for iso_id, read_ids in iso_to_name_dict.items():
        iso_name = '{}{}'.format(__program__, iso_id)
        out_array = [iso_name]
        out1 = all_out[read_ids[0]] # XXX everything except isHighSJ (high mapping quality SJ) are the same for all the reads
        # 0. write intersect info and cons_info to eval_out
        out1['geneID'], out1['geneName'], out1['geneStrand'] = itst_out_dict['geneID'][str(iso_name)], itst_out_dict['geneName'][str(iso_name)], itst_out_dict['geneStrand'][str(iso_name)]
        out1['CDS'], out1['UTR'], out1['lincRNA'], out1['antisense'] = itst_out_dict['CDS'][str(iso_name)], itst_out_dict['UTR'][str(iso_name)], itst_out_dict['lincRNA'][str(iso_name)], itst_out_dict['antisense'][str(iso_name)]
        out1['rRNA'], out1['Alu'], out1['allRepeat'], out1['upFlankAlu'], out1['downFlankAlu'] = itst_out_dict['rRNA'][str(iso_name)], itst_out_dict['Alu'][str(iso_name)], itst_out_dict['allRepeat'][str(iso_name)], itst_out_dict['upFlankAlu'][str(iso_name)], itst_out_dict['downFlankAlu'][str(iso_name)]
        # block type, anno
        out1['blockType'], out1['blockAnno'] = itst_out_dict['blockType'][str(iso_name)], itst_out_dict['blockAnno'][str(iso_name)]
        # full-length: BSJ is already high-confidenc, only check FSJ (internal SJ) XXX
        is_fullLength(out1, all_out, read_ids)
        
        # FSM/NIC/NNC for BSJ and interIso
        if 'True' in out1['isKnownBSJ']:
            bsj_cate = 'FSM'
        elif out1['isKnownSS'].startswith('True') and out1['isKnownSS'].endswith('True'):
            bsj_cate = 'NIC'
        else: bsj_cate = 'NNC'
        out1['BSJCate'] = bsj_cate
        iso_cate, start, gene_ids = None, int(out1['startCoor0based']), out1['geneID'].rsplit(',')
        if int(out1['blockCount']) == 1:
            iso_cate = 'FSM'
        else:
            start_array, size_array = out1['blockStarts'].split(','), out1['blockSize'].split(',')
            if '' in start_array: start_array.remove('')
            if '' in size_array: size_array.remove('')
            iso_coor = []
            for s, l in zip(start_array, size_array):
                iso_coor.extend([start+int(s)+1, start+int(s)+int(l)])
            int_iso_coor_str = '_'.join(map(str, iso_coor[1:-1]))
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
                iso_cate = 'NIC' if 'False' not in out1['isKnownSS'][1:-1] else 'NNC'
        out1['FSJCate'] = iso_cate

        for h in isoform_output_header[1:-2]:
            out_array.append(str(out1[h]))
        read_names = []
        for id in read_ids:
            name = all_out[id]['#readID'].rsplit('_cons')[0]
            read_names.append(name)
        out_array.extend([str(len(set(read_names))), ','.join(set(read_names))])
        out_fp.write('\t'.join(out_array) + '\n')
    ut.err_format_time('output_isoform_eval', 'Writing isoform-wise evaluation result to file done!')


# readName, readLen, consLen, copyNum, consFrac, chimInfo
def get_cons_info(cons_info_fp):
    cons_info_dict = {}
    for line in cons_info_fp:
        if line.startswith('#'):
            continue
        ele = line.rsplit()
        cons_info_dict[ele[0]] = ele[1:]
    return cons_info_dict


def bed12_to_gene_id(in_bam_bed, in_anno, site_dis, end_dis, out_five_gene_name_id, out_three_gene_name_id):
    out_dir = os.path.dirname(os.path.abspath(out_three_gene_name_id)) + '/'
    bam_five_site_bed, bam_five_site_exon_gtf = in_bam_bed + '.five.site.bed', in_bam_bed + '.five.site.exon.gtf'
    bam_three_site_bed, bam_three_site_exon_gtf = in_bam_bed + '.three.site.bed', in_bam_bed + '.three.site.exon.gtf'
    anno_five_site_bed, anno_five_site_exon_gtf = out_dir + os.path.basename(in_anno) + '.five.site.bed', out_dir + os.path.basename(in_anno) + '.five.site.exon.gtf'
    anno_three_site_bed, anno_three_site_exon_gtf = out_dir + os.path.basename(in_anno) + '.three.site.bed', out_dir + os.path.basename(in_anno) + '.three.site.exon.gtf'
    pg.bed12_to_site_bed(in_bam_bed, bam_five_site_bed, bam_three_site_bed, site_dis, end_dis)
    pg.gtf_to_site_bed(in_anno, anno_five_site_bed, anno_three_site_bed, 0, 0)
    ut.exec_cmd(sys.stderr, 'bed2exonGtf', '{} {} {}'.format(bed2exonGtf, bam_five_site_bed, bam_five_site_exon_gtf))
    ut.exec_cmd(sys.stderr, 'bed2exonGtf', '{} {} {}'.format(bed2exonGtf, bam_three_site_bed, bam_three_site_exon_gtf))
    ut.exec_cmd(sys.stderr, 'bed2exonGtf', '{} {} {}'.format(bed2exonGtf, anno_five_site_bed, anno_five_site_exon_gtf))
    ut.exec_cmd(sys.stderr, 'bed2exonGtf', '{} {} {}'.format(bed2exonGtf, anno_three_site_bed, anno_three_site_exon_gtf))
    ut.exec_cmd(sys.stderr, 'itst_gtf_gtf', '{} {} {} {}'.format(itst_gtf_gtf, bam_five_site_exon_gtf, anno_five_site_exon_gtf, out_five_gene_name_id))
    ut.exec_cmd(sys.stderr, 'itst_gtf_gtf', '{} {} {} {}'.format(itst_gtf_gtf, bam_three_site_exon_gtf, anno_three_site_exon_gtf, out_three_gene_name_id))


def get_site_gene_name_id(five_site_name_id_fn, three_site_name_id_fn, id_dict, name_dict, strand_dict):
    tot_dict = dd(lambda: dd(lambda: 0))
    with open(five_site_name_id_fn) as in_fp:
        for line in in_fp:
            if line.startswith('#'): continue
            ele = line.rsplit()  # read_id, gene_id, gene_name, gene_strand
            tot_dict[ele[0]][(ele[1], ele[2], ele[3])] += 1
    with open(three_site_name_id_fn) as in_fp:
        for line in in_fp:
            if line.startswith('#'): continue
            ele = line.rsplit()  # read_id, gene_id, gene_name, gene_strand
            tot_dict[ele[0]][(ele[1], ele[2], ele[3])] += 1

    for read_id, tup in tot_dict.items():
        max_cnt = max(tup.values())
        if max_cnt == 0: continue
        for (gene_id, gene_name, gene_strand), cnt in tup.items():
            if cnt == max_cnt:
                id_dict[read_id] = gene_id if id_dict[read_id] == 'NA' else id_dict[read_id] + ',' + gene_id 
                name_dict[read_id] = gene_name if name_dict[read_id] == 'NA' else name_dict[read_id] + ',' + gene_name 
                strand_dict[read_id] = gene_strand if strand_dict[read_id] == 'NA' else strand_dict[read_id] + ',' + gene_strand


def get_ovlp_gene_name_id(name_id_fn, id_dict, name_dict, strand_dict):
    assign_read = dd(lambda: '')
    for r in id_dict:
        assign_read[r] = ''

    with open(name_id_fn) as in_fp: # read_id, gene_id, gene_name, gene_strand
        for line in in_fp:
            if line.startswith('#'): continue
            ele = line.rsplit()
            if ele[0] not in assign_read:
                id_dict[ele[0]] = ele[1] if id_dict[ele[0]] == 'NA' else id_dict[ele[0]] + ',' + ele[1]
                name_dict[ele[0]] = ele[2] if name_dict[ele[0]] == 'NA' else name_dict[ele[0]] + ',' + ele[2]
                strand_dict[ele[0]] = ele[3] if strand_dict[ele[0]] == 'NA' else strand_dict[ele[0]] + ',' + ele[3]


def get_block_info(in_bed, ref_map_len_dict, start_coor_dict, end_coor_dict, block_count_dict, block_size_dict,
                   block_starts_dict):
    with open(in_bed) as in_fp:
        for line in in_fp:
            if line.startswith('#'): continue
            ele = line.rsplit()
            size = ele[10].split(',')
            if '' in ele[10].split(','):
                size.remove('')
            ref_map_len_dict[ele[3]] = sum(map(int, size))
            start_coor_dict[ele[3]] = int(ele[1])
            end_coor_dict[ele[3]] = int(ele[2])
            block_count_dict[ele[3]] = int(ele[9])
            block_size_dict[ele[3]] = ele[10]
            block_starts_dict[ele[3]] = ele[11]


def get_itst_len(itst_out, itst_dict):
    with open(itst_out) as in_fp:
        for line in in_fp:
            if line.startswith('#'): continue
            ele = line.rsplit()
            size = ele[10].rsplit(',')
            if '' in ele[10].rsplit(','):
                size.remove('')
            itst_dict[ele[3]] = sum(map(int, size))


def get_block_type(block_type, block_count_dict):
    block_type_dict = dd(lambda : '')
    for read_name, v in block_count_dict.items():
        type_array = ['E'] * int(v)
        for num, t in block_type[read_name].items():
            type_array[num-1] = t
        block_type_dict[read_name] = ','.join(type_array)
    return block_type_dict


# type: E/I/N exon/intron/intergenic
# based on "exon_id"
# intron_out: not exon
# intergenic_out: not exon and intron
# others: exon
def get_itst_block_type(intron_out, intergenic_out):
    block_type = dd(lambda: dd(lambda: ''))
    # if block in intron_out: block_type_dict[read_name][block_num] = 'I'
    with open(intron_out) as in_fp:
        for line in in_fp:
            if line.startswith('#'): continue
            try:
                id_idx = line.index('exon_id')
            except:
                ut.fatal_format_time('get_itst_block_type', 'No \"exon_id\" found in record.')
            id = line[id_idx+9:].split('"')[0]
            read_name, block_num = os.path.splitext(id)[0], int(os.path.splitext(id)[1][1:])
            block_type[read_name][block_num] = 'I'
    # if block in intergenic_out: block_type_dict[read_name][block_num] = 'N'
    with open(intergenic_out) as in_fp:
        for line in in_fp:
            if line.startswith('#'): continue
            try:
                id_idx = line.index('exon_id')
            except:
                ut.fatal_format_time('get_itst_block_type', 'No \"exon_id\" found in record.')
            id = line[id_idx+9:].split('"')[0]
            read_name, block_num = os.path.splitext(id)[0], int(os.path.splitext(id)[1][1:])
            block_type[read_name][block_num] = 'N'
    return block_type

# return coordinates of overlapped part on [s1,e1]
def get_block_ovlp(s1, e1, s2, e2):
    ovlp_s = 1 if s2 <= s1 else s2 - s1 + 1
    ovlp_e = e1-s1+1 if e1 <= e2 else e2 - s1 + 1#e1 - s1 + 1 - (e1 - e2)
    return ovlp_s, ovlp_e

# exon_out: 1-base coordinate
# overlap with exon => T1:E1(len)+I(len);T2:E1(len)+I(len)
# not overlap with exon => I/N, based on blockType
def get_itst_block_anno(exon_out):
    block_anno = dd(lambda: dd(lambda: dd(lambda: [])))
    with open(exon_out) as in_fp:
        for line in in_fp:
            if line.startswith('#'): continue
            ele = line.rsplit()
            block_start, block_end = int(ele[3]), int(ele[4])
            anno_start, anno_end = int(ele[19]), int(ele[20])
            ovlp_s, ovlp_e = get_block_ovlp(block_start, block_end, anno_start, anno_end)
            # print block_start, block_end, anno_start, anno_end, ovlp_s, ovlp_e
            try:
                id_idx = line.index('exon_id')
            except:
                ut.fatal_format_time('get_block_anno', 'No \"exon_id\" found in record.')
            id = line[id_idx + 9:].split('"')[0]
            read_name, block_num = os.path.splitext(id)[0], int(os.path.splitext(id)[1][1:])
            try:
                trans_id_idx = id_idx + line[id_idx:].index('transcript_id')
            except:
                ut.fatal_format_time('get_block_anno', 'No \"transcript_id\" found in record.')
            trans_id = line[trans_id_idx + 15:].split('"')[0]
            try:
                exon_num_idx = trans_id_idx + line[trans_id_idx:].index('exon_number')
            except:
                ut.fatal_format_time('get_block_anno', 'No \"exon_number\" found in record.')
            a = line[exon_num_idx + 12:].split(';')[0].split('"')
            exon_num = int(a[1]) if '' in a else int(a[0])
            block_anno[read_name][block_num][trans_id].append((ovlp_s, ovlp_e, exon_num))

    return block_anno


def get_block_anno(block_anno, block_size_dict):
    block_anno_dict = dd(lambda : '')
    for read_name, block_size in block_size_dict.items():
        size_array = list(map(int, block_size.split(',')))
        block_anno_array = ['NA'] * len(size_array)
        for block_num in block_anno[read_name]:
            block_len = size_array[block_num-1]
            anno_array = []
            for trans_id in block_anno[read_name][block_num]:
                trans_anno = block_anno[read_name][block_num][trans_id]
                trans_anno.sort()  # (ovlp_s, ovlp_e, exon_num)
                trans_anno_array = []
                # first intron
                intron_len = trans_anno[0][0] - 1
                if intron_len > 0:
                    trans_anno_array.append('I({})'.format(intron_len))
                for i, anno_block in enumerate(trans_anno):
                    # exon
                    trans_anno_array.append('E{}({})'.format(anno_block[2], anno_block[1]-anno_block[0] + 1))
                    # intron
                    intron_len = trans_anno[i + 1][0] - trans_anno[i][1] - 1 if i < len(trans_anno) - 1 else block_len - trans_anno[i][1]
                    if intron_len > 0:
                        trans_anno_array.append('I({})'.format(intron_len))
                anno_array.append('{}:{}'.format(trans_id, '+'.join(trans_anno_array)))
            block_anno_array[int(block_num)-1] = ';'.join(anno_array)
        block_anno_dict[read_name] = ','.join(block_anno_array)
    return block_anno_dict

def get_flank_bed(in_bed, flank_len, up_flank_bed, down_flank_bed):
    with open(in_bed, 'r') as in_fp, open(up_flank_bed, 'w') as up_fp, open(down_flank_bed, 'w') as down_fp:
        for line in in_fp:
            if line.startswith('#'): continue
            ele = line.rsplit()
            chr, start, end, id = ele[0], int(ele[1]), int(ele[2]), ele[3]
            up_start, up_end = start - flank_len, start
            down_start, down_end = end, end + flank_len
            if up_start >= 0:
                up_fp.write('{}\t{}\t{}\t{}\n'.format(chr, up_start, up_end, id))
            down_fp.write('{}\t{}\t{}\t{}\n'.format(chr, down_start, down_end, id))


def get_flank_Alu(flank_out, flank_dict):
    with open(flank_out, 'r') as in_fp:
        for line in in_fp:
            if line.startswith('#'): continue
            ele = line.rsplit()
            flank_dict[ele[3]] = flank_dict[ele[3]] = ele[9] + ele[7] if flank_dict[ele[3]] == 'NA' else flank_dict[ele[3]] + ',' + ele[9] + ele[7]



# TODO seperate intersecting and circRNA filtering
# input with dict{'CDS':'CDS.bed'}
def intersect_with_bed(out_dir, cons_bed12, all_anno, all_anno_bed, itst_anno_dict, flank_len, bedtools):
    itst_out_dict = dd(lambda: dd(lambda: 'NA'))
    cons_bam_exon_gtf = cons_bed12 + '.exon.gtf'
    ut.exec_cmd(sys.stderr, 'bed2exonGtf', '{} {} {}'.format(bed2exonGtf, cons_bed12, cons_bam_exon_gtf))

    # anno/CDS/UTR/antisense ...
    all_anno_exon_gtf = out_dir + os.path.basename(all_anno) + '.exon.gtf'
    all_anno_gene_bed = out_dir + os.path.basename(all_anno) + '.gene.bed'
    cds_gtf = out_dir + os.path.basename(all_anno) + '.cds.gtf'
    utr_gtf = out_dir + os.path.basename(all_anno) + '.utr.gtf'
    linc_gtf = out_dir + os.path.basename(all_anno) + '.lincRNA.gtf'
    anti_gtf = out_dir + os.path.basename(all_anno) + '.antisense.gtf'
    rRNA_gtf = out_dir + os.path.basename(all_anno) + '.rRNA.gtf'
    itst_anno_dict['CDS'], itst_anno_dict['UTR'], itst_anno_dict['lincRNA'], itst_anno_dict['antisense'], itst_anno_dict['rRNA'] = cds_gtf, utr_gtf, linc_gtf, anti_gtf, rRNA_gtf
    ut.exec_cmd(sys.stderr, 'exonGtf', 'awk -v OFS="\\t" \'($3=="exon"){print}\' ' + all_anno + ' > ' + all_anno_exon_gtf)
    ut.exec_cmd(sys.stderr, 'gtf2bed', 'awk -v OFS="\\t" \'($3=="gene"){print $1,$4-1,$5}\' ' + all_anno + ' > ' + all_anno_gene_bed)
    ut.exec_cmd(sys.stderr, 'gtf2bed', 'awk -v OFS="\\t" \'($3=="CDS"){print}\' ' + all_anno + ' > ' + cds_gtf)
    ut.exec_cmd(sys.stderr, 'gtf2bed', 'awk -v OFS="\\t" \'($3=="UTR" || $3=="five_prime_utr" || $3=="three_prime_utr"){print}\' ' + all_anno + ' > ' + utr_gtf)
    ut.exec_cmd(sys.stderr, 'gtf2bed', 'awk -v OFS="\\t" \'($3=="exon" && ($0 ~ /gene_biotype "lincRNA"/ || $0 ~ /gene_type "lincRNA"/)){print}\' ' + all_anno + ' > ' + linc_gtf)
    ut.exec_cmd(sys.stderr, 'gtf2bed', 'awk -v OFS="\\t" \'($3=="exon" && ($0 ~ /gene_biotype "antisense"/ || $0 ~ /gene_type "antisense"/)){print}\' ' + all_anno + ' > ' + anti_gtf)
    ut.exec_cmd(sys.stderr, 'gtf2bed',
                'awk -v OFS="\\t" \'($3=="exon" && ($0 ~ /gene_biotype "rRNA"/ || $0 ~ /gene_type "rRNA"/)){print}\' ' + all_anno + ' > ' + rRNA_gtf)

    # block information
    ref_map_len_dict, start_coor_dict, end_coor_dict, block_count_dict, block_size_dict, block_starts_dict = dict(), dict(), dict(), dict(), dict(), dict()
    get_block_info(cons_bed12, ref_map_len_dict, start_coor_dict, end_coor_dict, block_count_dict, block_size_dict, block_starts_dict)

    # gene bed * required
    gene_id_dict, gene_name_dict, gene_strand_dict = itst_out_dict['geneID'], itst_out_dict['geneName'], itst_out_dict['geneStrand']
    five_site_gene_name_id, three_site_gene_name_id = cons_bed12 + '.five.site.gene.out', cons_bed12 + '.three.site.gene.out'
    bed12_to_gene_id(cons_bed12, all_anno, 0, 0, five_site_gene_name_id, three_site_gene_name_id)
    get_site_gene_name_id(five_site_gene_name_id, three_site_gene_name_id, gene_id_dict, gene_name_dict, gene_strand_dict)

    ovlp_gene_name_id = cons_bed12 + '.ovlp.gene.out'
    ut.exec_cmd(sys.stderr, 'gtf2gene', '{} {} {} {}'.format(gtf2gene, cons_bam_exon_gtf, all_anno, ovlp_gene_name_id))
    get_ovlp_gene_name_id(ovlp_gene_name_id, gene_id_dict, gene_name_dict, gene_strand_dict)

    for b in itst_bed_name:
        if itst_anno_dict[b]:
            itst_out = cons_bed12 + '.{}.out'.format(b)
            ut.exec_cmd(sys.stderr, 'itst_gtf_bed', '{} {} {} {}'.format(itst_gtf_bed, cons_bam_exon_gtf, itst_anno_dict[b], itst_out))
            get_itst_len(itst_out, itst_out_dict[b])
            if b == 'Alu':
                up_flank_bed, down_flank_bed = cons_bed12 + '.up.bed', cons_bed12 + '.down.bed'
                get_flank_bed(cons_bed12, flank_len, up_flank_bed, down_flank_bed)
                up_flank_out, down_flank_out = cons_bed12 + '.up_flank_Alu.out', cons_bed12 + '.down_flank_Alu.out'
                ut.exec_cmd(sys.stderr, 'flank_Alu', '{} intersect -a {} -b {} -wb > {}'.format(bedtools, up_flank_bed, itst_anno_dict[b], up_flank_out))
                ut.exec_cmd(sys.stderr, 'flank_Alu', '{} intersect -a {} -b {} -wb > {}'.format(bedtools, down_flank_bed, itst_anno_dict[b], down_flank_out))
                up_flank_Alu_dict, down_flank_Alu_dict = itst_out_dict['upFlankAlu'], itst_out_dict['downFlankAlu']
                get_flank_Alu(up_flank_out, up_flank_Alu_dict)
                get_flank_Alu(down_flank_out, down_flank_Alu_dict)

    # block type
    if all_anno_bed and all_anno_gene_bed:
        intron_out = cons_bed12 + '.intron.out'
        ut.exec_cmd(sys.stderr, 'itst_intron', '{} intersect -v -a {} -b {} -split > {}'.format(bedtools, cons_bam_exon_gtf, all_anno_bed, intron_out))

        intergenic_out = cons_bed12 + '.intergenic.out'
        ut.exec_cmd(sys.stderr, 'itst_intergenic', '{} intersect -v -a {} -b {} > {}'.format(bedtools, cons_bam_exon_gtf, all_anno_gene_bed, intergenic_out))
        block_type = get_itst_block_type(intron_out, intergenic_out)
        block_type_dict = get_block_type(block_type, block_count_dict)
        itst_out_dict['blockType'] = block_type_dict
    # block annotation
    if all_anno_exon_gtf:
        exon_out = cons_bed12 + '.exon.out'
        ut.exec_cmd(sys.stderr, 'itst_exon', '{} intersect -a {} -b {} -wa -wb > {}'.format(bedtools, cons_bam_exon_gtf, all_anno_exon_gtf, exon_out))
        block_anno = get_itst_block_anno(exon_out)
        block_anno_dict = get_block_anno(block_anno, block_size_dict)
        itst_out_dict['blockAnno'] = block_anno_dict
    return itst_out_dict


def is_coincide_knownSS(eval_out):
    if not eval_out['isCanoBSJ']:
        return False
    dis_to_bsj = eval_out['disToCanoBSJ'].split(',')
    dis_to_known = eval_out['disToKnownFSS'].split(',')
    return dis_to_bsj[0] == dis_to_known[0] and dis_to_bsj[-1] == dis_to_known[-1]


# No xid in key 4 bases
def get_bsj_xid(cigar_str, flank_len=10, key_flank_len=2):
    if cigar_str == 'NA':
        return flank_len*2, flank_len*2
    cigar_stats = pb.cigarstring_to_cigarstats(cigar_str)
    key_cigar_tup = pb.get_spec_ref_cigar(pb.cigarstring_to_cigartuples(cigar_str), flank_len-key_flank_len, flank_len+key_flank_len)
    key_cigar_stats = pb.cigartuples_to_cigarstats(key_cigar_tup)
    return pb.get_xid_from_cigar(cigar_stats), pb.get_xid_from_cigar(key_cigar_stats)


def get_bsj_aln_score(cigar_str, flank_len=10):
    eq, xid = 1, -1
    if cigar_str == 'NA':
        return -flank_len*2, -flank_len*2
    cigar_stats = pb.cigarstring_to_cigarstats(cigar_str)
    key_cigar_tup = pb.get_spec_ref_cigar(pb.cigarstring_to_cigartuples(cigar_str), flank_len-2, flank_len+2)
    key_cigar_stats = pb.cigartuples_to_cigarstats(key_cigar_tup)
    score = cigar_stats[cigar_op_dict['=']] * eq + (cigar_stats[cigar_op_dict['X']] + cigar_stats[cigar_op_dict['D']] + cigar_stats[cigar_op_dict['I']]) * xid
    key_score = key_cigar_stats[cigar_op_dict['=']] * eq + (key_cigar_stats[cigar_op_dict['X']] + key_cigar_stats[cigar_op_dict['D']] + key_cigar_stats[cigar_op_dict['I']]) * xid
    return score, key_score


# $1=="yes"&&($2=="yes"||($3=="GT/AG"&&$5==4&&$4>=18&&($6>1||($6==1&&$7>=150))))
def circRNA_filter_core(eval_out, bsj_motifs=['GT/AG'], bsj_xid=1, bsj_key_xid=0, min_circ_dis=150):
    # filter #1
    if len(eval_out['chrom']) >= 6 or eval_out['chrom'].startswith('chrM') or eval_out['chrom'].startswith('chrUn'):
        return False
    # filter #2
    if is_coincide_knownSS(eval_out):
        return True
    else:
        # filter #3
        cano = 0
        for m in bsj_motifs:
            if eval_out['canoBSJMotif'].endswith(m):
                cano = 1
                break
        if cano == 0:
            return False
        # filter #4/5
        xid, key_xid = get_bsj_xid(eval_out['alignAroundCanoBSJ'])
        if xid > bsj_xid:
            return False
        if key_xid > bsj_key_xid:
            return False
        # filter #6/7
        if (int(eval_out['endCoor']) + int(eval_out['disToCanoBSJ'].split(',')[1])) - (int(eval_out['startCoor0based']) + int(eval_out['disToCanoBSJ'].split(',')[0])) < min_circ_dis:
            return False
    return True

def get_clu_id(d, i):
    if i not in d:
        ut.fatal_format_time('get_clu_id', 'Error: unknown i ({}).'.format(i))
    if d[i] == i:
        return i
    else:
        return get_clu_id(d, d[i])


def get_one_cons_out(mul_cons):
    if len(mul_cons) == 1:
        return mul_cons[0]
    # cluster cons based on chimInfo
    clu_id, tot_frac = dict(), dd(lambda:0)
    for cons_out in mul_cons:
        id = int(cons_out['#readID'].rsplit('cons')[-1])
        clu_id[id] = id
    for i, cons_out in enumerate(mul_cons):
        if cons_out['chimInfo'] == 'NA': continue
        # print cons_out['#readID'], cons_out['chimInfo']
        _ids = list(map(int, re.split('RC|ID', cons_out['chimInfo'])))
        ids = []
        for id in _ids:
            if id in clu_id:
                ids.append(id)
        clu_ids = [clu_id[ii] for ii in ids]
        min_clu_id = min(clu_ids)
        for iden_id in ids:
            if clu_id[iden_id] > min_clu_id:
                clu_id[iden_id] = min_clu_id
                clu_id[clu_id[iden_id]] = min_clu_id
    for i in clu_id:
        clu_id[i] = get_clu_id(clu_id, i)

    for cons_out in mul_cons:
        id = int(cons_out['#readID'].rsplit('cons')[-1])
        tot_frac[clu_id[id]] += float(cons_out['consFrac'])

    best_cons_i = 0
    best_score, best_key_score = get_bsj_aln_score(mul_cons[0]['alignAroundCanoBSJ'])
    best_cluster_frac = tot_frac[0]
    # best_cons_frac = float(mul_cons[0]['consFrac'])
    best_NM, best_AS = mul_cons[0]['mapNM'], mul_cons[0]['mapAS']
    best_chrom = mul_cons[0]['chrom']

    for i, cons_out in enumerate(mul_cons[1:]):
        id = int(cons_out['#readID'].rsplit('cons')[-1])
        cons_i = i + 1
        score, key_score = get_bsj_aln_score(cons_out['alignAroundCanoBSJ'])
        cluster_frac = tot_frac[clu_id[id]]
        # cons_frac = float(cons_out['consFrac'])
        NM, AS = cons_out['mapNM'], cons_out['mapAS']
        chrom = cons_out['chrom']
        if score > best_score:
            best_cons_i, best_score, best_key_score, best_cluster_frac, best_NM, best_AS, best_chrom = cons_i, score, key_score, cluster_frac, NM, AS, chrom
        elif score == best_score:
            if cluster_frac > best_cluster_frac:
                best_cons_i, best_score, best_key_score, best_cluster_frac, best_NM, best_AS, best_chrom = cons_i, score, key_score, cluster_frac, NM, AS, chrom
            elif cluster_frac == best_cluster_frac:
                if AS > best_AS:
                    best_cons_i, best_score, best_key_score, best_cluster_frac, best_NM, best_AS, best_chrom = cons_i, score, key_score, cluster_frac, NM, AS, chrom
                elif AS == best_AS:
                    if NM < best_NM:
                        best_cons_i, best_score, best_key_score, best_cluster_frac, best_NM, best_AS, best_chrom = cons_i, score, key_score, cluster_frac, NM, AS, chrom
                    elif NM == best_NM:
                        if key_score > best_key_score:
                            best_cons_i, best_score, best_key_score, best_cluster_frac, best_NM, best_AS, best_chrom = cons_i, score, key_score, cluster_frac, NM, AS, chrom
                        elif key_score == best_key_score:
                            if (best_chrom == 'chrM' or chrom.startswith('chrUn')) and (chrom != 'chrM' and not chrom.startswith('chrUn')):
                                best_cons_i, best_score, best_key_score, best_cluster_frac, best_NM, best_AS, best_chrom = cons_i, score, key_score, cluster_frac, NM, AS, chrom
    return mul_cons[best_cons_i]


# 1. high.out => filtered.bsj.coors
# 2. filtered.bsj.coors + all.out => filerted.out
# 3. filtered.out => circRNA.bed12
# 4. circRNA.bed12 + GTF => ovlp.bed
# 5. filtered.bsj.coors + filterted.out + ovlp.bed => circRNA.out
def circRNA_bsj_filter(all_out, cano_motifs=['GT/AG'], bsj_xid=1, key_bsj_xid=0, min_circ_dis=150):
    # filtered_coors: adjusted coors by disToCanoBSJ
    filtered_bsj_coors, all_bsj, all_bsj_stats_dict = dict(), dict(), dd(lambda:0) # ead_with_bsj_n'/'bsj_n'/'known_bsj_n'
    all_bsj_stats_dict['known_bsj_n'] = dd(lambda:0)
    last_read_name, mul_cons = '', []
    for id, eval_out in all_out.items():
        read_name = eval_out['#readID'].rsplit('_cons')[0]
        if read_name != last_read_name:
            if mul_cons:
                out1 = get_one_cons_out(mul_cons)
                if out1['isCanoBSJ']:
                    bsj = (out1['chrom'], int(out1['startCoor0based']) + int(out1['disToCanoBSJ'].split(',')[0]), int(out1['endCoor']) + int(out1['disToCanoBSJ'].split(',')[1]))

                    # all bsj
                    all_bsj_stats_dict['read_with_bsj_n'] += 1
                    if bsj not in all_bsj:
                        all_bsj[bsj] = 1
                        all_bsj_stats_dict['bsj_n'] += 1
                        for i, known_bsj in enumerate(out1['isKnownBSJ'].rsplit(',')):
                            all_bsj_stats_dict['known_bsj_n'][i] += (known_bsj == 'True')
                        all_bsj_stats_dict['known_bsj_n'][i+1] += ('False' not in out1['isKnownBSJ'])
                    if circRNA_filter_core(out1, cano_motifs, bsj_xid, key_bsj_xid, min_circ_dis):
                        # bsj coors
                        filtered_bsj_coors[bsj] = 1
            mul_cons = [eval_out]
            last_read_name = read_name
        else:
            mul_cons.append(eval_out)
    if mul_cons:
        out1 = get_one_cons_out(mul_cons)
        if out1['isCanoBSJ']:
            bsj = (out1['chrom'], int(out1['startCoor0based']) + int(out1['disToCanoBSJ'].split(',')[0]), int(out1['endCoor']) + int(out1['disToCanoBSJ'].split(',')[1]))
            # all bsj
            all_bsj_stats_dict['read_with_bsj_n'] += 1
            if bsj not in all_bsj:
                all_bsj[bsj] = 1
                all_bsj_stats_dict['bsj_n'] += 1
                for i, known_bsj in enumerate(out1['isKnownBSJ'].rsplit(',')):
                    if known_bsj == 'True':
                        all_bsj_stats_dict['known_bsj_n'][i] += 1
                if 'False' not in out1['isKnownBSJ']:
                    all_bsj_stats_dict['known_bsj_n'][i+1] += 1
            if circRNA_filter_core(out1, cano_motifs, bsj_xid, key_bsj_xid, min_circ_dis):
                # bsj coors
                filtered_bsj_coors[bsj] = 1
    return filtered_bsj_coors, all_bsj_stats_dict

def circRNA_rescue(all_out, bsj_coors, h_bam, anno_site, anno_exon):
    # filtered_coors: adjusted coors by disToCanoBSJ
    filtered_id, filtered_out, filtered_coors = 0, dict(), dict()
    last_read_name, mul_cons = '', []
    for id, eval_out in all_out.items():
        read_name = eval_out['#readID'].rsplit('_cons')[0]
        if read_name != last_read_name:
            if mul_cons:
                out1 = get_one_cons_out(mul_cons)
                if 'NA' not in out1['disToCanoBSJ']:
                    chrom, adj_start, adj_end = out1['chrom'], int(out1['startCoor0based']) + int(out1['disToCanoBSJ'].split(',')[0]), int(out1['endCoor']) + int(out1['disToCanoBSJ'].split(',')[1])
                    chrom_id = h_bam.get_tid(chrom)
                    # print(last_read_name, chrom, adj_start, adj_end)
                    if (chrom, adj_start, adj_end) in bsj_coors:
                        # coors
                        filtered_coors[filtered_id] = [chrom_id]
                        filtered_coors[filtered_id].extend(pg.get_coor_from_block(out1['startCoor0based'], out1['blockSize'], out1['blockStarts']))
                        filtered_coors[filtered_id][1] = adj_start + 1
                        filtered_coors[filtered_id][-1] = adj_end
                        # all_out
                        out1['startCoor0based'], out1['endCoor'] = adj_start, adj_end
                        start, out1['blockSize'], out1['blockStarts'] = pg.get_block_from_coor(filtered_coors[filtered_id][1:])
                        if start < 0 or min(map(int, out1['blockSize'].rsplit(','))) < 0 or min(map(int, out1['blockStarts'].rsplit(','))):
                            filtered_coors.pop(filtered_id)
                            continue
                        out1['refMapLen'] = sum(map(int, out1['blockSize'].rsplit(',')))
                        # if not is_coincide_knownSS(out1): # update
                        is_known_ss = out1['isKnownSS'].rsplit(',')
                        is_known_ss[0] = 'True' if (chrom_id, out1['canoBSJMotif'][0] == '-', True, adj_start+1) in anno_site else 'False'
                        is_known_ss[-1] = 'True' if (chrom_id, out1['canoBSJMotif'][0] == '-', False, adj_end) in anno_site else 'False'
                        out1['isKnownSS'] = ','.join(is_known_ss)

                        is_known_exon = out1['isKnownExon'].rsplit(',')
                        block_sizes = list(map(int, out1['blockSize'].rsplit(',')))
                        is_known_exon[0] = 'True' if (chrom_id, out1['canoBSJMotif'][0] == '-', adj_start+1, adj_start + block_sizes[0]) in anno_exon else 'False'
                        is_known_exon[-1] = 'True' if (chrom_id, out1['canoBSJMotif'][0] == '-', adj_end-block_sizes[-1]+1, adj_end) in anno_exon else 'False'
                        out1['isKnownExon'] = ','.join(is_known_exon)
                        filtered_out[filtered_id] = out1
                        filtered_id += 1
            mul_cons = [eval_out]
            last_read_name = read_name
        else:
            mul_cons.append(eval_out)
    if mul_cons:
        out1 = get_one_cons_out(mul_cons)
        if 'NA' not in out1['disToCanoBSJ']:
            chrom, adj_start, adj_end = out1['chrom'], int(out1['startCoor0based']) + int(out1['disToCanoBSJ'].split(',')[0]), int(out1['endCoor']) + int(out1['disToCanoBSJ'].split(',')[1])
            chrom_id = h_bam.get_tid(chrom)
            # print(last_read_name, chrom, adj_start, adj_end)
            if (chrom, adj_start, adj_end) in bsj_coors:
                # coors
                filtered_coors[filtered_id] = [chrom_id]
                filtered_coors[filtered_id].extend(pg.get_coor_from_block(out1['startCoor0based'], out1['blockSize'], out1['blockStarts']))
                filtered_coors[filtered_id][1] = adj_start + 1
                filtered_coors[filtered_id][-1] = adj_end
                # all_out
                out1['startCoor0based'], out1['endCoor'] = adj_start, adj_end
                start, out1['blockSize'], out1['blockStarts'] = pg.get_block_from_coor(filtered_coors[filtered_id][1:])
                if start < 0 or min(map(int, out1['blockSize'].rsplit(','))) < 0 or min(map(int, out1['blockStarts'].rsplit(','))):

                    filtered_coors.pop(filtered_id)
                    return filtered_out, filtered_coors
                out1['refMapLen'] = sum(map(int, out1['blockSize'].rsplit(',')))
                # if not is_coincide_knownSS(out1): # update
                is_known_ss = out1['isKnownSS'].rsplit(',')
                is_known_ss[0] = 'True' if (chrom_id, out1['canoBSJMotif'][0] == '-', True, adj_start+1) in anno_site else 'False'
                is_known_ss[-1] = 'True' if (chrom_id, out1['canoBSJMotif'][0] == '-', False, adj_end) in anno_site else 'False'
                out1['isKnownSS'] = ','.join(is_known_ss)


                is_known_exon = out1['isKnownExon'].rsplit(',')
                block_sizes = list(map(int, out1['blockSize'].rsplit(',')))
                is_known_exon[0] = 'True' if (chrom_id, out1['canoBSJMotif'][0] == '-', adj_start+1, adj_start + block_sizes[0]) in anno_exon else 'False'
                is_known_exon[-1] = 'True' if (chrom_id, out1['canoBSJMotif'][0] == '-', adj_end-block_sizes[-1]+1, adj_end) in anno_exon else 'False'
                out1['isKnownExon'] = ','.join(is_known_exon)
                filtered_out[filtered_id] = out1
                filtered_id += 1
    return filtered_out, filtered_coors


# circ_sj = [dict(), dict()] : multiple circRNA annotation files
def eval_core(id, r, cons_info_dict, ref_fa, cons_fa, all_site, all_exon, all_sj,
              circ_sj, sj_xid, key_sj_xid, site_dis, end_dis, all_out):
    r_block, coors = pb.get_block(r), []
    for b in r_block:
        coors.extend([b[2], b[3]])

    eval_out = {i: 'NA' for i in whole_output_header}  # for each BAM record
    eval_out['#readID'], eval_out['chrom'], eval_out['mapStrand'], eval_out['consMapLen'] = r.query_name, r.reference_name, '-' if r.is_reverse else '+', pb.get_aligned_read_length(r)
    eval_out['startCoor0based'], eval_out['endCoor'] = r.reference_start, r.reference_end
    eval_out['blockCount'] = len(r_block)
    _s, eval_out['blockSize'], eval_out['blockStarts'] = pg.get_block_from_coor(coors)
    eval_out['refMapLen'] = sum(map(int, eval_out['blockSize'].rsplit(',')))
    [eval_out['readLen'], eval_out['consLen'], eval_out['copyNum'], eval_out['consFrac'], eval_out['chimInfo']] = cons_info_dict[r.query_name]

    # mapping score
    eval_out['mapNM'], eval_out['mapAS'] = int(r.get_tag('NM')) if r.has_tag('NM') else 0, int(r.get_tag('AS')) if r.has_tag('AS') else 0
    ref_seq = ref_fa[r.reference_name]
    # 1. compare with whole gene annotation TODO use bedtools intersect
    is_known_site, is_known_exon, is_known_junc, is_cano_junc, dis_to_known_ss, dis_to_cano_sj, cano_splice_motif = pg.comp_with_anno(r_block, all_site, all_exon, all_sj, ref_seq, site_dis, end_dis)
    splice_strand = [m[0] if k else 'NA' for k,m in zip(is_known_junc, cano_splice_motif)]
    is_high_sj = pb.parse_sj_alignment(r, min_xid=sj_xid, key_min_xid=key_sj_xid)

    # 2. check bsj
    bsj = (r.reference_id, r.is_reverse, int(eval_out['endCoor']) + 1, int(eval_out['startCoor0based']))
    force_strand = ''
    if '+' in splice_strand and '-' not in splice_strand: force_strand = '+'
    elif '-' in splice_strand and '+' not in splice_strand: force_strand = '-'
    bsj_dis_to_known_ss = [dis_to_known_ss[0], dis_to_known_ss[-1]]
    is_known_bsj, is_cano_bsj, dis_to_cano_bsj, bsj_motif, align_bsj = pg.is_known_cano_bsj(bsj, circ_sj, ref_seq, cons_fa[r.query_name][:].seq.upper(), int(eval_out['startCoor0based']), int(eval_out['endCoor']), r.is_reverse, r.cigartuples, int(eval_out['refMapLen']), int(eval_out['consMapLen']), int(eval_out['consLen']), end_dis, force_strand, bsj_dis_to_known_ss)

    eval_out['isKnownBSJ'], eval_out['isCanoBSJ'], eval_out['disToCanoBSJ'], eval_out['canoBSJMotif'], eval_out['alignAroundCanoBSJ'] = is_known_bsj, is_cano_bsj, dis_to_cano_bsj, bsj_motif, align_bsj

    eval_out['isKnownSS'], eval_out['isKnownFSJ'], eval_out['isCanoFSJ'], eval_out['isKnownExon'], eval_out['disToKnownFSS'], eval_out['disToCanoFSJ'], eval_out['canoFSJMotif'] = ','.join(map(str, is_known_site)), ','.join(map(str, is_known_junc)),  ','.join(map(str, is_cano_junc)), ','.join(map(str, is_known_exon)), ','.join(map(str, dis_to_known_ss)), ','.join(dis_to_cano_sj), ','.join(cano_splice_motif)
    eval_out['isHighFSJ'] = ','.join(map(str, is_high_sj))

    all_out[id] = eval_out

def hcBSJ_fullIso(high_bam, low_bam, long_len, cons_info, cons,
                   ref, all_anno, circ_anno, itst_anno_dict, bedtools=bedtools,
                   flank_len=flank_len,
                   cano_motif='GT/AG', bsj_xid=1, key_bsj_xid=0, min_circ_dis=150, rescue_low = False,
                   sj_xid=1, key_sj_xid=0,
                   isoform_out_fn='iso.out', circRNA_bed='circRNA.bed', stats_out_fn='stats.out'):
    # filter criteria
    cano_motifs = ['GT/AG'] if cano_motif == 'GT/AG' else ['GT/AG', 'GC/AG', 'AT/AC']
    # read ref_fa
    ref_fa, cons_fa = Fasta(ref), Fasta(cons)
    # read whole and circRNA annotation file
    out_dir = os.path.dirname(os.path.abspath(isoform_out_fn)) + '/'

    all_anno_bed = out_dir + os.path.basename(all_anno) + '.bed'
    all_anno_gene_pred = out_dir + os.path.basename(all_anno) + '.gene_pred' 
    ut.exec_cmd(sys.stderr, 'gtfToGenePred', '{} -genePredExt -ignoreGroupsWithoutExons {} {}'.format(gtfToGenePred, all_anno, all_anno_gene_pred))
    ut.exec_cmd(sys.stderr, 'genePredToBed', '{} {} {}'.format(genePredToBed, all_anno_gene_pred, all_anno_bed))

    all_trans = pg.get_transcript_from_gene_pred(all_anno_gene_pred, high_bam)
    all_site = pg.get_splice_site_from_bed12(all_anno_bed, high_bam)
    all_sj = pg.get_splice_junction_from_bed12(all_anno_bed, False, high_bam)
    all_exon = pg.get_exon_from_bed12(all_anno_bed, high_bam)
    # multiple circ_anno input files
    circ_sj = []
    for circ_anno1 in circ_anno.rsplit(','):
        if os.path.splitext(circ_anno1)[1] == '.gtf':
            circ_anno_bed = out_dir + os.path.basename(circ_anno1) + '.bed'
            ut.exec_cmd(sys.stderr, 'gtf2bed', '{} {} {}'.format(gtf2bed, circ_anno1, circ_anno_bed))
        else:
            circ_anno_bed = circ_anno1
        circ_sj.append(pg.get_back_splice_junction_from_bed(circ_anno_bed, high_bam))

    with ps.AlignmentFile(high_bam) as h_bam, ps.AlignmentFile(low_bam) as l_bam, open(cons_info, 'r') as cons_info_fp, open(isoform_out_fn, 'w') as iso_fp:
        # write header information for eval.out
        output_header(iso_fp, isoform_output_header)
        # read cons.info
        cons_info_dict = get_cons_info(cons_info_fp)
        all_out = dict()
        processed_cnt, batch_n = 0, 100000

        # 0. read-wise evaluation
        ut.err_format_time('read_wise_eval', 'Generating read-wise evaluation result ... ')
        for r in h_bam:
            if r.is_unmapped: continue
            eval_core(processed_cnt, r, cons_info_dict, ref_fa, cons_fa, all_site, all_exon, all_sj, circ_sj, sj_xid, key_sj_xid, site_dis, end_dis, all_out)
            processed_cnt += 1
            if processed_cnt % batch_n == 0:
                ut.err_format_time('high_quality', '{} high mapping quality BAM records have been processed ... '.format(processed_cnt))
        ut.err_format_time('read_wise_eval', 'Generating read-wise evaluation result done!')

        # 1. circRNA filtering using only high-quality reads
        # all_bsj: bsj from all high-quality reads
        ut.err_format_time('filter_circRNA_read', 'Filtering back-splice-junctions ...')
        filtered_bsj_coors, all_bsj_stats_dict = circRNA_bsj_filter(all_out, cano_motifs, bsj_xid, key_bsj_xid, min_circ_dis)
        ut.err_format_time('filter_circRNA_read', 'Filtering back-splice-junctions done!')

        if rescue_low:
            high_processed_cnt = processed_cnt
            for r in l_bam:
                if r.is_unmapped: continue
                eval_core(processed_cnt, r, cons_info_dict, ref_fa, cons_fa, all_site, all_exon, all_sj, circ_sj, sj_xid, key_sj_xid, site_dis, end_dis, all_out)
                processed_cnt += 1
                if (processed_cnt - high_processed_cnt) % batch_n == 0:
                    ut.err_format_time('low_quality', '{} low mapping quality BAM records have been processed ... '.format(processed_cnt - high_processed_cnt))

        # 2. find reads with reliable bsjs
        ut.err_format_time('rescue_reads', 'Rescuing reads using reliable back-splice-junctions ...')
        circRNA_out, circRNA_coors = circRNA_rescue(all_out, filtered_bsj_coors, h_bam, all_site, all_exon)
        ut.err_format_time('rescue_reads', 'Rescuing reads using reliable back-splice-junctions done!')

        # 3. circRNA isoform
        sorted_iso_to_name_dict = ui.uniq_isoform_with_unsorted_coors(circRNA_coors, True)
        circRNA_isoform_bed12(circRNA_bed, circRNA_out, sorted_iso_to_name_dict)
        # 4. intersect with annoation
        itst_out_dict = intersect_with_bed(out_dir, circRNA_bed, all_anno, all_anno_bed, itst_anno_dict, flank_len, bedtools)
        # 5. output isoform result
        output_isoform_eval(iso_fp, circRNA_out, itst_out_dict, all_trans, sorted_iso_to_name_dict)
    bs.stats_core(long_len, cons_info, high_bam, isoform_out_fn, all_bsj_stats_dict, stats_out_fn)
    ut.exec_cmd(sys.stderr, 'remove temp', 'rm {}.* {}{}.*'.format(circRNA_bed, out_dir, os.path.basename(all_anno)))
    return


def hcBSJ_fullIso_core(args):
    itst_anno_dict = dd(lambda : '')
    itst_anno_dict['Alu'] = args.Alu
    itst_anno_dict['allRepeat'] = args.all_repeat

    hcBSJ_fullIso(args.high_bam, args.low_bam, args.long_len, args.cons_info,
        args.cons_fa, args.ref, args.all_anno, args.circRNA_anno, itst_anno_dict, args.bedtools,
        args.flank_len, args.cano_motif, 
        args.bsj_xid, args.key_bsj_xid, args.min_circ_dis, args.rescue_low,
        args.fsj_xid, args.key_fsj_xid,
        args.out, args.bed, args.stats_out)


def parser_argv():
    # parse command line arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="Evaluate circRNA long-read with gene annotation")
    parser.add_argument('long_fa', metavar='long.fa', type=str, help='Long read sequencing data generated with isoCirc.')
    parser.add_argument('long_len', metavar='long.fa.len', type=str, help='Read length file of long-read data.')
    parser.add_argument('cons_info', metavar='cons.info', type=str, help='Consensus information file.')
    parser.add_argument('cons_fa', metavar='cons.fa', type=str, help='Concatemer of two copies of consensus sequence in FASTA format.')
    parser.add_argument('high_bam', metavar='high.bam', type=str, help='Unsorted high-quality alignment file of consensus sequence.')
    parser.add_argument('low_bam', metavar='low.bam', type=str, help='Unsorted low-quality alignment file of consensus sequence.')
    parser.add_argument('ref', metavar='ref.fa', type=str, help='Reference genome sequence.')
    parser.add_argument('all_anno', metavar='all.gtf', type=str, help='Gene annotation file in GTF format.')
    parser.add_argument('circRNA_anno', metavar='circRNA.bed/gtf', type=str, help='circRNA annotation file in BED or GTF format. Use \',\' to separate multiple circRNA annotation files.')
    parser.add_argument('out', metavar='{}.out'.format(__program__), type=str, help='Isoform-wise circRNA output file.')
    parser.add_argument('bed', metavar='{}.bed'.format(__program__), type=str, help='BED12 file of \'{}.out\'.'.format(__program__))
    parser.add_argument('stats_out', metavar='{}_stats.out'.format(__program__), type=str, help='Basic stats numbers of \'{}.out\''.format(__program__))

    # parser.add_argument('--type', type=str, help='Type of sequencing data: Oxford Nanopore(ont) or Pacific Biosciences (pb).', choices=['ont', 'pb'], default='ont')
    parser.add_argument('-t', '--threads', type=int, default=threads, help='Number of thread to use.')
    parser.add_argument('--bedtools', help='Path to bedtools.', default=bedtools)


    # parser.add_argument('-s', '--site-dis', type=int, default=site_dis, help='Maximum allowed distance between circRNA internal splice-site and annoated splice-site.')
    # parser.add_argument('-S', '--end-dis', type=int, default=end_dis, help='Maximum allowed distance between circRNA back-splice-site and annoated splice-site.')
    parser.add_argument('--cano-motif', type=str, default='GT/AG', help='Canonical back-splice motif (GT/AG or all three motifs: GT/AG, GC/AG, AT/AC).', choices=['GT/AG', 'all'])
    parser.add_argument('--bsj-xid', type=int, default=1, help='Maximum allowed mis/ins/del for 20-bp exonic sequence flanking the BSJ (10-bp each side).')
    parser.add_argument('--key-bsj-xid', type=int, default=0, help='Maximum allowed mis/ins/del for 4-bp exonic sequence flanking the BSJ (2-bp each side).')
    parser.add_argument('--min-circ-dis', type=int, default=150, help='Minimum distance between the genomic coordinates of the two back-splice sites.')
    parser.add_argument('--rescue-low', default=False, action='store_true', help='Use high mapping quality reads to rescue low mapping quality reads.')

    parser.add_argument('--fsj-xid', type=int, default=1, help='Maximum allowed mis/ins/del for 20-bp exonic sequence flanking the FSJ (10-bp each side).')
    parser.add_argument('--key-fsj-xid', type=int, default=0, help='Maximum allowed mis/ins/del for 4-bp exonic sequence flanking the FSJ (2-bp each side).')

    parser.add_argument('--Alu', type=str, default='', help='Alu repetitive element annotation in BED format. ')
    parser.add_argument('--flank-len', type=int, default=flank_len, help='Length of upstream and downstream flanking sequence to search for Alu.')
    parser.add_argument('--all-repeat', type=str, default='', help='All repetitive element annotation in BED format.')
    return parser.parse_args()


def main():
    args = parser_argv()
    hcBSJ_fullIso_core(args)

if __name__ == '__main__':
    main()
