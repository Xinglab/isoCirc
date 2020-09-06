import os, sys, re
from collections import defaultdict as dd
import gffutils as gu
import pysam as ps
import itertools
from Bio.Seq import Seq

import isocirc.utils as ut
import isocirc.parse_bam as pb

def restore_gff_db(gtf_fn):
    gtf_db = None
    if gtf_fn is not None:
        gtf_db_fn = gtf_fn + '.gffdb'
        if not os.path.isfile(gtf_db_fn):
            try:
                # check if 'gene' or 'transcript' is in GTF
                disable_gene, disable_trans = False, False
                with open(gtf_fn) as fp:
                    l = 0
                    for line in fp:
                        if line[0] != '#':
                            if line.split()[2] == 'gene':
                                disable_gene = True
                            elif line.split()[2] == 'transcript':
                                disable_trans = True
                            l += 1

                        if (disable_gene and disable_trans) or l == 100: break
                ut.err_format_time('restore_gtf_db', 'Creating GTF databases for {} ...'.format(gtf_fn))
                gtf_db = gu.create_db(gtf_fn, gtf_db_fn, disable_infer_genes=disable_gene,
                                      disable_infer_transcripts=disable_trans)
                ut.err_format_time('restore_gtf_db', 'Creating GTF databases for {} done!'.format(gtf_fn))

            except:
                ut.err_format_time('restore_gtf_db',
                                   'Error in parsing {}\nCheck if annotation file format is correct'.format(gtf_fn))
                sys.exit(IOError)
        else:
            try:
                ut.err_format_time('restore_gtf_db', 'Retrieving gff database for {} ...'.format(gtf_fn))
                gtf_db = gu.FeatureDB(gtf_db_fn)
                ut.err_format_time('restore_gtf_db', 'Retrieving gff database for {} done!'.format(gtf_fn))

            except:
                ut.err_format_time('restore_gtf_db',
                                   'Error in parsing {}\nTry to remove this db file and re-run'.format(gtf_db_fn))
                sys.exit(IOError)
    return gtf_db


# return: dict{(TransID,GeneID):[(chr, is_reverse, start, end), ()...]}
def get_exon_block_from_gtf(in_gtf, is_db, bam_fn=''):
    gtf_db = in_gtf if is_db else restore_gff_db(in_gtf)
    ut.err_format_time("get_exon_block_from_gtf", "Loading exon block from {} ... ".format(in_gtf))
    exon_block = dd(lambda: [])
    bam = ps.AlignmentFile(bam_fn) if bam_fn else None
    for exon in gtf_db.features_of_type('exon', order_by='start'):
        tid = bam.get_tid(exon.chrom) if bam else exon.chrom
        ID = (exon.attributes['transcript_id'][0], exon.attributes['gene_id'][0])
        exon_block[ID].append((tid, exon.strand == '-', int(exon.start), int(exon.end)))
    ut.err_format_time("get_exon_block_from_gtf", "Loading exon block from {} done!".format(in_gtf))
    if bam: bam.close()
    return exon_block


def dict2list(d):
    array = []
    for k, v in d.items():
        a = [k] + v
        array.append(a)
    return array


# return: list[[block] [block] ... ]
# block: [(TransID,GeneID), (exon1), (exon2) ... ]
# exon: (chr, is_reverse, start, end)
# index: dict{id:(start,end), id:(start,end)}
def get_sorted_exon_block_from_gtf(in_gtf, is_db, bam_fn=''):
    if not bam_fn: ut.fatal_format_time("get_sorted_exon_block_from_gtf", 'No BAM header provided.')
    gtf_db = in_gtf if is_db else restore_gff_db(in_gtf)
    ut.err_format_time("get_sorted_exon_block_from_gtf", "Loading exon block from {} ... ".format(in_gtf))
    exon_block_dict = dd(lambda: [])

    with ps.AlignmentFile(bam_fn) as bam:
        for exon in gtf_db.features_of_type('exon', order_by='start'):
            ID = (exon.attributes['transcript_id'][0], exon.attributes['gene_id'][0])
            exon_block_dict[ID].append((bam.get_tid(exon.chrom), exon.strand == '-', int(exon.start), int(exon.end)))
    exon_block_list = dict2list(exon_block_dict)
    exon_block_list = sorted(exon_block_list, key=lambda x: (x[1][0], x[1][2], x[-1][3]))

    # index
    block_index = {}
    last_tid, tid, i_start, i_end = -1, 0, 0, 0
    for i, block in enumerate(exon_block_list):
        exon = block[1]
        tid = exon[0]
        if tid != last_tid:
            if last_tid != -1:
                i_end = i - 1
                block_index[last_tid] = (i_start, i_end)
            i_start = i
            last_tid = tid
    if last_tid != -1:
        i_end = i - 1
        block_index[last_tid] = (i_start, i_end)
    ut.err_format_time("get_sorted_exon_block_from_gtf", "Loading exon block from {} done!".format(in_gtf))
    return exon_block_list, block_index


# return: list[[block] [block] ... ]
# block: [(TransID,GeneID), (exon1), (exon2) ... ]
# exon: (chr, is_reverse, start, end)
def get_exon_block_from_bed12(in_bed, bam_fn=''):
    # chromStart: 0-base, exonStarts: 0-base
    header_ele = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand',
                  'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'exonStarts']
    bed_header = {header_ele[i]: i for i in range(len(header_ele))}
    exon_block = []
    ut.err_format_time("get_exon_block_from_bed12", "Loading exon block from {} ...".format(in_bed))
    bam = ps.AlignmentFile(bam_fn) if bam_fn else None
    with open(in_bed, 'r') as bed:
        for line in bed:
            if line.startswith('#'): continue
            ele = line[:-1].rsplit('\t')
            chrom = bam.get_tid(ele[bed_header['chrom']]) if bam else ele[bed_header['chrom']]
            strand = ele[bed_header['strand']]
            start = int(ele[bed_header['chromStart']])
            start_array, len_array = ele[bed_header['exonStarts']].split(','), ele[bed_header['blockSizes']].split(
                ',')
            if '' in start_array: start_array.remove('')
            if '' in len_array: len_array.remove('')
            exon_start = [int(i) for i in start_array]
            exon_len = [int(i) for i in len_array]
            exon_block.append([(ele[bed_header['name']], ele[bed_header['name']])])
            for s, l in zip(exon_start, exon_len):
                exon_block[-1].append((chrom, strand == '-', int(start + s + 1), int(start + s + l)))

    exon_block = sorted(exon_block, key=lambda x: (x[1][0], x[1][2], x[-1][3]))
    ut.err_format_time("get_exon_block_from_bed12", "Loading exon block from {} done!".format(in_bed))
    if bam: bam.close()
    return exon_block


# index: dict{(iD):(start,end), id:(start,end)}
def get_sorted_exon_block_from_bed12(in_bed, bam_fn=''):
    if not bam_fn: ut.fatal_format_time("get_sorted_exon_block_from_bed12", 'No BAM header provided.')
    # chromStart: 0-base, exonStarts: 0-base
    header_ele = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand',
                  'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'exonStarts']
    bed_header = {header_ele[i]: i for i in range(len(header_ele))}
    ut.err_format_time("get_sorted_exon_block_from_bed12", "Loading exon block from {} ...".format(in_bed))
    exon_block = []
    bam = ps.AlignmentFile(bam_fn) if bam_fn else None
    with open(in_bed, 'r') as bed:
        for line in bed:
            if line.startswith('#'): continue
            ele = line[:-1].rsplit('\t')
            chrom = ele[bed_header['chrom']]
            strand = ele[bed_header['strand']]
            start = int(ele[bed_header['chromStart']])
            start_array, len_array = ele[bed_header['exonStarts']].split(','), ele[bed_header['blockSizes']].split(
                ',')
            if '' in start_array: start_array.remove('')
            if '' in len_array: len_array.remove('')
            exon_start = [int(i) for i in start_array]
            exon_len = [int(i) for i in len_array]
            exon_block.append([(ele[bed_header['name']], ele[bed_header['name']])])
            for s, l in zip(exon_start, exon_len):
                exon_block[-1].append((bam.get_tid(chrom), strand == '-', int(start + s + 1), int(start + s + l)))
    if bam: bam.close()

    exon_block = sorted(exon_block, key=lambda x: (x[1][0], x[1][2], x[-1][3]))
    # index
    block_index = {}
    last_tid, tid, i_start, i_end = -1, 0, 0, 0
    for i, block in enumerate(exon_block):
        exon = block[1]
        tid = exon[0]
        if tid != last_tid:
            if last_tid != -1:
                i_end = i
                block_index[last_tid] = (i_start, i_end)
            i_start = i
            last_tid = tid
    if last_tid != -1:
        i_end = len(exon_block)
        block_index[last_tid] = (i_start, i_end)
    ut.err_format_time("get_sorted_exon_block_from_bed12", "Loading exon block from {} done!".format(in_bed))
    return exon_block, block_index


def gff2len(in_gff, out_fn):
    # print in_gff, out_fn
    with open(out_fn, 'w') as out:
        gtf = True if os.path.splitext(in_gff)[1] == '.gtf' else False
        if gtf:
            exon_block = get_exon_block_from_gtf(in_gff, False)

            for id, b in exon_block.items():
                len = 0
                for exon in b:
                    len += exon[3] - exon[2] + 1
                out.write(id + '\t' + str(len) + '\n')

        else:
            exon_block = get_exon_block_from_bed12(in_gff)
            for b in exon_block:
                id = b[0]
                len = 0
                for exon in b[1:]:
                    len += exon[3] - exon[2] + 1
                out.write(id + '\t' + str(len) + '\n')


# [(tid, is_rev, sj_start, sj_end, cov)]
def get_sorted_splice_junction_from_gtf(in_gtf, is_db, include_end, bam_fn=""):
    if not bam_fn: ut.fatal_format_time("get_sorted_splice_junction_from_gtf", 'No BAM header provided.')
    gtf_db = in_gtf if is_db else restore_gff_db(in_gtf)
    ut.err_format_time("get_sorted_splice_junction_from_gtf", "Loading splice junction from {} ... ".format(in_gtf))
    sj = dict()
    bam = ps.AlignmentFile(bam_fn) if bam_fn else None
    for trans in gtf_db.features_of_type('transcript', order_by='start'):
        sj_start, sj_end = -1, -1
        tid = bam.get_tid(trans.chrom)
        is_rev = trans.strand == '-'
        for exon in gtf_db.children(trans, featuretype='exon', order_by='start'):
            sj_end = exon.start - 1
            if sj_start > 0 and sj_end > 0:
                sj[(tid, is_rev, sj_start, sj_end)] = 1
            sj_start = exon.end + 1
        if include_end:
            sj[(tid, is_rev, trans.end + 1, trans.start - 1)] = 1
    if bam: bam.close()

    sj_list = dict2list(sj)
    sj_list = sorted(sj_list, key=lambda x: (x[0], x[2], x[3]))

    # index
    sj_index = {}
    last_tid, tid, i_start, i_end = -1, 0, 0, 0
    for i, sj1 in enumerate(sj_list):
        tid = sj1[0]
        if tid != last_tid:
            if last_tid != -1:
                i_end = i - 1
                sj_index[last_tid] = (i_start, i_end)
            i_start = i
            last_tid = tid
    if last_tid != -1:
        i_end = i - 1
        sj_index[last_tid] = (i_start, i_end)
    ut.err_format_time("get_sorted_splice_junction_from_gtf", "Loading splice junction from {} done!".format(in_gtf))
    return sj_list, sj_index


# {(tid, is_rev, sj_start, sj_end):coverage}
def get_splice_junction_from_gtf(in_gtf, is_db, include_end, bam_fn=""):
    gtf_db = in_gtf if is_db else restore_gff_db(in_gtf)
    ut.err_format_time("get_splice_junction_from_gtf", "Loading splice junction from {} ... ".format(in_gtf))
    sj = dict()
    bam = ps.AlignmentFile(bam_fn) if bam_fn else None
    for trans in gtf_db.features_of_type('transcript', order_by='start'):
        sj_start, sj_end = -1, -1
        tid = bam.get_tid(trans.chrom) if bam else trans.chrom
        is_rev = trans.strand == '-'
        for exon in gtf_db.children(trans, featuretype='exon', order_by='start'):
            sj_end = exon.start - 1
            if sj_start > 0 and sj_end > 0:
                sj[(tid, is_rev, sj_start, sj_end)] = 1
            sj_start = exon.end + 1
        if include_end:
            sj[(tid, is_rev, trans.end + 1, trans.start - 1)] = 1  # coverage count?

    ut.err_format_time("get_splice_junction_from_gtf", "Loading splice junction from {} done!".format(in_gtf))
    if bam: bam.close
    return sj

def get_back_splice_junction_from_bed(in_bed, bam_fn=""):
    ut.err_format_time("get_back_splice_junction_from_bed", "Loading splice junction from {} ... ".format(in_bed))
    # chromStart: 0-base, exonStarts: 0-base
    header_ele = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand']
    bed_header = {header_ele[i]: i for i in range(len(header_ele))}
    sj = dict()
    bam = ps.AlignmentFile(bam_fn) if bam_fn else None
    with open(in_bed, 'r') as bed:
        for line in bed:
            if line.startswith('#'): continue
            ele = line[:-1].rsplit()
            chrom = ele[bed_header['chrom']]
            start = int(ele[bed_header['chromStart']])
            end = int(ele[bed_header['chromEnd']])
            tid = bam.get_tid(chrom) if bam else chrom
            if len(ele) >= 6:
                is_rev = (ele[bed_header['strand']] == '-')
                sj[(tid, is_rev, end+1, start)] = 1
            else:
                sj[(tid, False, end+1, start)] = 1
                sj[(tid, True, end+1, start)] = 1
    ut.err_format_time("get_back_splice_junction_from_bed", "Loading splice junction from {} done!".format(in_bed))
    if bam: bam.close()
    return sj

def get_splice_junction_from_bed12(in_bed, include_end, bam_fn=""):
    ut.err_format_time("get_splice_junction_from_bed12", "Loading splice junction from {} ... ".format(in_bed))
    # chromStart: 0-base, exonStarts: 0-base
    header_ele = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand',
                  'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'exonStarts']
    bed_header = {header_ele[i]: i for i in range(len(header_ele))}
    sj = dict()
    bam = ps.AlignmentFile(bam_fn) if bam_fn else None
    with open(in_bed, 'r') as bed:
        for line in bed:
            if line.startswith('#'): continue
            ele = line[:-1].rsplit('\t')
            chrom = ele[bed_header['chrom']]
            strand = ele[bed_header['strand']]
            start = int(ele[bed_header['chromStart']])
            start_array, len_array = ele[bed_header['exonStarts']].split(','), ele[bed_header['blockSizes']].split(
                ',')
            if '' in start_array: start_array.remove('')
            if '' in len_array: len_array.remove('')
            exon_start = [int(i) for i in start_array]
            exon_len = [int(i) for i in len_array]
            tid = bam.get_tid(chrom) if bam else chrom
            is_rev = strand == '-'
            sj_start, sj_end = -1, -1
            for s, l in zip(exon_start, exon_len):
                sj_end = start + s
                if sj_start > 0 and sj_end > 0:
                    sj[(tid, is_rev, sj_start, sj_end)] = 1
                sj_start = start + s + l + 1
            if include_end:
                sj[(tid, is_rev, sj_start, start)] = 1
    ut.err_format_time("get_splice_junction_from_bed12", "Loading splice junction from {} done!".format(in_bed))
    if bam: bam.close()
    return sj

# return:
# {(tid, is_rev, is_left, site_coor) : is_start/end}
# exonic 1-base coordinate
# is_left: site is on the left of the exon
def get_splice_site_from_gtf(in_gtf, is_db, bam_fn=""):
    split_site = dict()
    gtf_db = in_gtf if is_db else restore_gff_db(in_gtf)
    ut.err_format_time("get_splice_site_from_gtf", "Loading splice site from {} ... ".format(in_gtf))
    bam = ps.AlignmentFile(bam_fn) if bam_fn else None
    for trans in gtf_db.features_of_type('transcript', order_by='start'):
        tid = bam.get_tid(trans.chrom) if bam else trans.chrom
        is_rev = trans.strand == '-'
        split_site[(tid, is_rev, True, trans.start)] = True
        split_site[(tid, is_rev, False, trans.end)] = True
        for exon in gtf_db.children(trans, featuretype='exon', order_by='start'):
            if (tid, is_rev, True, exon.start) not in split_site:
                split_site[(tid, is_rev, True, exon.start)] = False
            if (tid, is_rev, False, exon.end) not in split_site:
                split_site[(tid, is_rev, False, exon.end)] = False
    ut.err_format_time("get_splice_site_from_gtf", "Loading splice site from {} done!".format(in_gtf))
    if bam: bam.close()
    return split_site

# return
# {gene_id: ['s_e_s_e ...']}
# s/e: exonic 1-base position
# exonic 1-base coordinate
def get_transcript_from_gene_pred(in_gene_pred, bam_fn=""):
    trans = dd(lambda: [])
    ut.err_format_time('get_transcript_from_bed12', 'Loading transcript from {} ... '.format(in_gene_pred))
    header_ele = ['transID', 'chrom', 'strand', 'transStart', 'transEnd', 'cdsStart', 'cdsEnd', 'blockCount', 'exonStarts', 'exonEnds', 'score', 'geneID', 'cdsStartStatus', 'cdsEndStatus', 'exonFrame']
    pred_header = {header_ele[i]: i for i in range(len(header_ele))}
    # bam = ps.AlignmentFile(bam_fn) if bam_fn else None
    with open(in_gene_pred, 'r') as pred:
        for line in pred:
            if line.startswith('#'): continue
            ele = line[:-1].rsplit('\t')
            # chrom = ele[bed_header['chrom']]
            # strand = ele[bed_header['strand']]
            start_array, end_array = ele[pred_header['exonStarts']].split(','), ele[pred_header['exonEnds']].split(
                ',')
            gene_id = ele[pred_header['geneID']]
            if '' in start_array: start_array.remove('')
            if '' in end_array: end_array.remove('')
            exon_start = [int(i) + 1 for i in start_array]
            exon_end = [int(i) for i in end_array]
            # tid = bam.get_tid(chrom) if bam else chrom
            # is_rev = strand == '-'
            coor = []
            for s, e in zip(exon_start, exon_end):
                coor.extend([s, e])
            trans[gene_id].append('_'.join(map(str, coor)))
    ut.err_format_time('get_transcript_from_gene_pred', 'Loading transcript from {} done!'.format(in_gene_pred))
    # if bam: bam.close()
    return trans

# return:
# {(tid, is_rev, is_left, site_coor) : is_start/end}
# exonic 1-base coordinate
# is_left: site is on the left of the exon
# excluding start/end of each transcript
def get_splice_site_from_bed12(in_bed, bam_fn=""):
    splice_site = dict()
    ut.err_format_time('get_splice_site_from_bed12', 'Loading splice site from {} ... '.format(in_bed))
    # chromStart: 0-base, exonStarts: 0-base
    header_ele = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand',
                  'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'exonStarts']
    bed_header = {header_ele[i]: i for i in range(len(header_ele))}
    bam = ps.AlignmentFile(bam_fn) if bam_fn else None
    with open(in_bed, 'r') as bed:
        for line in bed:
            if line.startswith('#'): continue
            ele = line[:-1].rsplit('\t')
            chrom = ele[bed_header['chrom']]
            strand = ele[bed_header['strand']]
            start = int(ele[bed_header['chromStart']])
            end = int(ele[bed_header['chromEnd']])
            start_array, len_array = ele[bed_header['exonStarts']].split(','), ele[bed_header['blockSizes']].split(
                ',')
            if '' in start_array: start_array.remove('')
            if '' in len_array: len_array.remove('')
            exon_start = [int(i) for i in start_array]
            exon_len = [int(i) for i in len_array]
            tid = bam.get_tid(chrom) if bam else chrom
            is_rev = strand == '-'

            for i in range(len(exon_start)):
                s, l = exon_start[i], exon_len[i]
                keep_left = False if i == 0 else True
                keep_right = False if i == len(exon_start) - 1 else True
                
                if keep_left:  
                    splice_site[(tid, is_rev, True, start + 1 + s)] = True
                if keep_right:
                    splice_site[(tid, is_rev, False, start + s + l)] = True

    ut.err_format_time('get_splice_site_from_bed12', 'Loading splice site from {} done!'.format(in_bed))
    if bam: bam.close()
    return splice_site


# return:
# {(tid, is_rev, start, end) : is_start/end}
# exonic 1-base coordinate
def get_exon_from_gtf(in_gtf, is_db, bam_fn=""):
    exon_dict = dict()
    gtf_db = in_gtf if is_db else restore_gff_db(in_gtf)
    ut.err_format_time("get_exon_from_gtf", "Loading exon from {} ... ".format(in_gtf))
    bam = ps.AlignmentFile(bam_fn) if bam_fn else None
    for trans in gtf_db.features_of_type('transcript', order_by='start'):
        tid = bam.get_tid(trans.chrom) if bam else trans.chrom
        is_rev = trans.strand == '-'
        for exon in gtf_db.children(trans, featuretype='exon', order_by='start'):
            exon_dict[(tid, is_rev, exon.start, exon.end)] = True \
                if exon.start == trans.start or exon.end == trans.end else False
    ut.err_format_time("get_exon_from_gtf", "Loading exon from {} done!".format(in_gtf))
    if bam: bam.close()
    return exon_dict


def get_exon_from_bed12(in_bed, bam_fn=""):
    ut.err_format_time("get_exon_from_bed12", "Loading exon from {} ... ".format(in_bed))
    # chromStart: 0-base, exonStarts: 0-base
    header_ele = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand',
                  'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'exonStarts']
    bed_header = {header_ele[i]: i for i in range(len(header_ele))}
    exon = dict()
    bam = ps.AlignmentFile(bam_fn) if bam_fn else None
    with open(in_bed, 'r') as bed:
        for line in bed:
            if line.startswith('#'): continue
            ele = line[:-1].rsplit('\t')
            chrom = ele[bed_header['chrom']]
            strand = ele[bed_header['strand']]
            start = int(ele[bed_header['chromStart']])
            end = int(ele[bed_header['chromEnd']])
            start_array, len_array = ele[bed_header['exonStarts']].split(','), ele[bed_header['blockSizes']].split(
                ',')
            if '' in start_array: start_array.remove('')
            if '' in len_array: len_array.remove('')
            exon_start = [int(i) for i in start_array]
            exon_len = [int(i) for i in len_array]
            tid = bam.get_tid(chrom) if bam else chrom
            is_rev = strand == '-'
            for s, l in zip(exon_start, exon_len):
                exon[(tid, is_rev, start + s + 1, start + s + l)] = True \
                    if s == 0 or start + s + l == end else False
    ut.err_format_time("get_exon_from_bed12", "Loading exon from {} done!".format(in_bed))
    if bam: bam.close()
    return exon

# sort by:
# 1. abs(left_dis_to_center - right_dis_to_center)
# 2. abs(left_dis_to_center) + abs(right_dis_to_center)
def get_min_offset(left_idxs, left_center, right_idxs, right_center):
    left_dis = [i - left_center for i in left_idxs]
    right_dis = [i - right_center for i in right_idxs]
    min_delta_dis, min_sum_dis = sys.maxsize, sys.maxsize
    left_idx, right_idx = None, None
    for i, ld in enumerate(left_dis):
        for j, rd in enumerate(right_dis):
            if abs(ld-rd) < min_delta_dis or \
                    (abs(ld-rd) == min_delta_dis and abs(ld) + abs(rd) < min_sum_dis):
                left_idx, right_idx = left_idxs[i], right_idxs[j]
                min_delta_dis = abs(ld-rd)
                min_sum_dis = abs(ld) + abs(rd)
    return left_idx, right_idx

def get_cano_sj(left_seq, left_center, right_seq, right_center, strand):
    is_cano, cano_strand, left_idx, right_idx, left_motif, right_motif = False, 'NA', -1, -1, '', ''
    for_motif = [('GT', 'AG'), ('GC', 'AG'), ('AT', 'AC')]
    rev_motif = [('CT', 'AC'), ('CT', 'GC'), ('GT', 'AT')]
    cano_motif = ['GT/AG', 'GC/AG', 'AT/AC']

    if not strand:  # forward and reverse have same priority
        for fm, rm, cm in zip(for_motif, rev_motif, cano_motif):
            f_left_idx, f_right_idx, r_left_idx, r_right_idx = None, None, None, None
            if fm[0] in left_seq and fm[1] in right_seq:  # forward
                left_idxs = [i.start() for i in re.finditer(fm[0], left_seq)]  # forward
                right_idxs = [i.start() for i in re.finditer(fm[1], right_seq)]
                f_left_idx, f_right_idx = get_min_offset(left_idxs, left_center, right_idxs, right_center)
            if rm[0] in left_seq and rm[1] in right_seq:  # reverse
                left_idxs = [i.start() for i in re.finditer(rm[0], left_seq)]
                right_idxs = [i.start() for i in re.finditer(rm[1], right_seq)]
                r_left_idx, r_right_idx = get_min_offset(left_idxs, left_center, right_idxs, right_center)
            if f_left_idx is not None and r_left_idx is not None:
                if abs((r_left_idx - left_center) - (r_right_idx - right_center)) < abs((f_left_idx - left_center) - (f_right_idx - right_center)) or \
                        (abs((r_left_idx - left_center) - (r_right_idx - right_center)) == abs((f_left_idx - left_center) - (f_right_idx - right_center)) and
                         abs(r_left_idx - left_center) + abs(r_right_idx - right_center) < abs(f_left_idx - left_center) + abs(f_right_idx - right_center)):
                    return True, '-', r_left_idx, r_right_idx, cm
                else:
                    return True, '+', f_left_idx, f_right_idx, cm
            elif f_left_idx is not None:
                return True, '+', f_left_idx, f_right_idx, cm
            elif r_left_idx is not None:
                return True, '-', r_left_idx, r_right_idx, cm
    else:
        if strand == '+':
            in_motif = for_motif
        elif strand == '-':
            in_motif = rev_motif
        for motif, cm in zip(in_motif, cano_motif):
            if motif[0] in left_seq and motif[1] in right_seq:
                # find indices that are the most close to the centers
                left_idxs = [i.start() for i in re.finditer(motif[0], left_seq)]
                right_idxs = [i.start() for i in re.finditer(motif[1], right_seq)]
                left_idx, right_idx = get_min_offset(left_idxs, left_center, right_idxs, right_center)
                return True, strand, left_idx, right_idx, cm
    return is_cano, cano_strand, left_idx, right_idx, 'NA'

def get_cano_sjs(down_seq, up_seq, strand):
    canos = []
    is_cano, cano_strand, down_idx, up_idx, down_motif, up_motif = False, 'NA', -1, -1, '', ''
    for_motif = [('GT', 'AG'), ('GC', 'AG'), ('AT', 'AC')]
    rev_motif = [('CT', 'AC'), ('CT', 'GC'), ('GT', 'AT')]
    cano_motif = ['GT/AG', 'GC/AG', 'AT/AC']
    if not strand:  # forward and reverse have same priority
        for fm, rm, cm in zip(for_motif, rev_motif, cano_motif):
            f_down_idx, f_up_idx, r_down_idx, r_up_idx = None, None, None, None
            if fm[0] in down_seq and fm[1] in up_seq:  # forward
                down_idxs = [i.start() for i in re.finditer(fm[0], down_seq)]  # forward
                up_idxs = [i.start() for i in re.finditer(fm[1], up_seq)]
                for down_idx, up_idx in list(itertools.product(down_idxs, up_idxs)):
                    canos.append(('+', down_idx, up_idx, cm))
            if rm[0] in down_seq and rm[1] in up_seq:  # reverse
                down_idxs = [i.start() for i in re.finditer(rm[0], down_seq)]
                up_idxs = [i.start() for i in re.finditer(rm[1], up_seq)]
                for down_idx, up_idx in list(itertools.product(down_idxs, up_idxs)):
                    canos.append(('-', down_idx, up_idx, cm))
            if canos:
                break
    else:
        if strand == '+':
            in_motif = for_motif
        elif strand == '-':
            in_motif = rev_motif
        for motif, cm in zip(in_motif, cano_motif):
            if motif[0] in down_seq and motif[1] in up_seq:
                # find indices that are the most close to the centers
                down_idxs = [i.start() for i in re.finditer(motif[0], down_seq)]
                up_idxs = [i.start() for i in re.finditer(motif[1], up_seq)]
                for down_idx, up_idx in list(itertools.product(down_idxs, up_idxs)):
                    canos.append((strand, down_idx, up_idx, cm))
            if canos:
                break
    return canos

# return alnScore, alnCigar
def get_cano_bsj_align(up_dis, down_dis, bsj_strand, ref_seq, read_seq, start, end, end_dis, is_reverse, cigartuples, ref_map_len, cons_len):
    # ref
    up_ref_cano_pos, down_ref_cano_pos = start + up_dis, end + down_dis
    bsj_ref_seq = ref_seq[down_ref_cano_pos-end_dis:down_ref_cano_pos].seq.upper() + ref_seq[up_ref_cano_pos:up_ref_cano_pos+end_dis].seq.upper()
    # read
    # read_seq += read_seq
    if is_reverse: read_seq = str(Seq(read_seq).reverse_complement())

    # NEW CHANGE
    # corresponding query's base positions of up_ref_cano_pos+end_dis and down_ref_cano_pos-end_dis: up_read_cano_pos and down_read_cano_pos
    # bsj_read_seq = read_seq[down_read_cano_pos_start : up_read_cano_pos_end]
    if up_dis + end_dis < 0 or down_dis - end_dis > 0:
        return -sys.maxsize, 'NA'
    up_cigartuples = pb.get_spec_ref_cigar(cigartuples, 0, up_dis + end_dis)
    up_cigarstats = pb.cigartuples_to_cigarstats(up_cigartuples)
    up_len = up_cigarstats[0] + up_cigarstats[1] + up_cigarstats[7] + up_cigarstats[8]
    down_cigartuples = pb.get_spec_ref_cigar(cigartuples, ref_map_len + down_dis - end_dis, ref_map_len)
    down_cigarstats = pb.cigartuples_to_cigarstats(down_cigartuples)
    down_len = down_cigarstats[0] + down_cigarstats[1] + down_cigarstats[7] + down_cigarstats[8]

    left_soft_clip = cigartuples[0][1] if cigartuples[0][0] == 4 or cigartuples[0][0] == 5 else 0
    right_soft_clip = cigartuples[-1][1] if cigartuples[-1][0] == 4 or cigartuples[-1][0] == 5 else 0
    if left_soft_clip + up_len <= cons_len:
        up_end = left_soft_clip + up_len + cons_len
        down_start = cons_len * 2 - right_soft_clip - down_len
    else:
        up_end = left_soft_clip + up_len
        down_start = cons_len - right_soft_clip - down_len
    if down_start < up_end and down_start >= 0 and up_end <= cons_len * 2:
        bsj_read_seq = read_seq[down_start:up_end]
        return pb.pairwise_align(bsj_ref_seq, bsj_read_seq, 'g', True)
    else:
        return -sys.maxsize, 'NA'

    # left_soft_clip += cons_len if left_soft_clip < end_dis else 0
    # bsj_read_seq = read_seq[left_soft_clip - end_dis + up_dis:left_soft_clip + end_dis + up_dis]
    # return pb.pairwise_align(bsj_ref_seq, bsj_read_seq, 'g', True)


# TODO use bedtools intersect to extract known sites/junctions
# bsj[rid, is_rev, end, start]
# NEW_CHANGE
# delta_aln_len = consMapLen - consLen
# constrain: | canoBSJLen - consLen | <= dis
def is_known_cano_bsj(bsj, anno_bsj, ref_seq, read_seq, start, end, is_reverse, cigartuples, ref_map_len, cons_map_len, cons_len, end_dis, force_strand, dis_to_known_ss):
    known, cano, up_dis, down_dis, strand, cano_motif, alignBSJ = [], False, 'NA', 'NA', 'NA', 'NA', 'NA'
    for anno_bsj1 in anno_bsj:
        known.append(False)
    # TODO if cons_map_len / (cons_len+0.0) <= pb.high_repeat_ratio: # multi-copies cons

    delta_aln_cons = cons_map_len - cons_len # delta_aln_cano = s1 - s2
    delta_ref_aln_cons = int(round((cons_map_len - cons_len) * ref_map_len / cons_map_len))
    # AG | up -- circRNA -- down | GT
    up_range = [0]
    if delta_ref_aln_cons >= 0:
        # (-dis, +dis+delta), (-dis-delta, +dis)
        up_range_s = -end_dis
        up_range_e = delta_ref_aln_cons + end_dis
        down_range_s = -delta_ref_aln_cons - end_dis
        down_range_e = end_dis
        for i in range(1, end_dis + 1):
            up_range.append(i)
            # down_range.append(i)
            up_range.append(-i)
            # down_range.append(-i)
        for i in range(end_dis + 1, end_dis + delta_ref_aln_cons + 1):
            up_range.append(i)
            # down_range.append(-i)
    else:
        # TODO (-dis+delta, +dis), (-dis, +dis-delta)
        # (-dis, +dis), (-dis, +dis)
        up_range_s = -end_dis
        up_range_e = end_dis
        down_range_s = -end_dis
        down_range_e = end_dis
        for i in range(1, end_dis+1):
            up_range.append(i)
            up_range.append(-i)
            # down_range.append(i)
            # down_range.append(-i)
        # for i in range(end_dis + 1, end_dis - delta_ref_aln_cons + 1):
        # 	up_range.append(-i)
        # 	down_range.append(i)
    end_range = [0]
    for i in range(1, end_dis+1):
        end_range.append(i)
        end_range.append(-i)
    anno_start, anno_end = bsj[2], bsj[3] # down, up
    # for i in end_range:
    #     for s1 in up_range:
    #         s2 = s1 - delta_ref_aln_cons + i
    #         if s2 < down_range_s or s2 > down_range_e:
    #             continue
    #         for anno_i in range(len(anno_bsj)):
    #             if (force_strand and (bsj[0], force_strand=='-', bsj[2]+s2, bsj[3]+s1) in anno_bsj[anno_i]) or \
    #                     (not force_strand and ((bsj[0], True, bsj[2]+s2, bsj[3]+s1) in anno_bsj[anno_i] or \
    #                         (bsj[0], False, bsj[2]+s2, bsj[3]+s1) in anno_bsj[anno_i])):
    #                 known[anno_i] = True
    #                 anno_start, anno_end = bsj[2] + s2, bsj[3] + s1
    #         if True in known:
    #             break
    #     if True in known:
    #         break
    if ref_seq:
        ref_len = len(ref_seq)
        if not(anno_start + 1 > ref_len or anno_end > ref_len or anno_start - 1 < 0 or anno_end - 2 < 0):
            # S = 0 if known else S
            # if True in known:
            #     up_range_s = up_range_e = down_range_s = down_range_e = 0
            down_seq = ref_seq[max(0, anno_start - 1 + down_range_s): min(ref_len, anno_start + 1 + down_range_e)].seq.upper()
            down_start = max(0, anno_start - 1 + down_range_s) + 1
            # down_center = anno_start - dwon_start
            up_seq = ref_seq[max(0, anno_end - 2 + up_range_s): min(ref_len, anno_end + up_range_e)].seq.upper()
            up_start = max(0, anno_end - 2 + up_range_s) + 1
            # up_center = anno_end - up_start - 1

            canos = get_cano_sjs(down_seq, up_seq, force_strand)
            if canos:
                max_score = -sys.maxsize
                dis = sys.maxsize
                for strand1, down_idx, up_idx, cano_motif1 in canos:
                    up_dis1, down_dis1 = up_start + up_idx - (bsj[3]-1), down_start + down_idx - bsj[2]
                    if up_dis1 + end_dis < 0 or down_dis1-end_dis > 0:
                        continue
                    delta_cano_aln_cons = (up_dis1 - down_dis1) * cons_map_len / ref_map_len
                    if abs(delta_aln_cons - delta_cano_aln_cons) > end_dis:
                        continue
                    cano_motif1 = strand1 + cano_motif1
                    score, alignBSJ1 = get_cano_bsj_align(up_dis1, down_dis1, strand1, ref_seq, read_seq, start, end, end_dis, is_reverse, cigartuples, ref_map_len, cons_len)
                    delta_dis = abs(up_dis1 - down_dis1 - delta_ref_aln_cons)
                    if score > max_score or (score == max_score and delta_dis < dis):
                        cano = True
                        up_dis, down_dis = up_dis1, down_dis1
                        cano_motif = cano_motif1
                        alignBSJ = alignBSJ1
                        max_score = score
                        dis = delta_dis
            if cano:
                for anno_i in range(len(anno_bsj)):
                    if (force_strand and (bsj[0], force_strand=='-', bsj[2]+down_dis, bsj[3]+up_dis) in anno_bsj[anno_i]) or \
                            (not force_strand and ((bsj[0], True, bsj[2]+down_dis, bsj[3]+up_dis) in anno_bsj[anno_i] or \
                                (bsj[0], False, bsj[2]+down_dis, bsj[3]+up_dis) in anno_bsj[anno_i])):
                        known[anno_i] = True
    return ','.join(map(str, known)), cano, '{},{}'.format(up_dis, down_dis), cano_motif, alignBSJ

# sj1: 1-base intronic position
def is_known_cano_sj(sj1, anno_sj, ref_seq, site_dis, end_dis, force_strand):
    anno_start, anno_end = sj1[2], sj1[3]
    known, cano, dis_to_cano, strand, cano_motif = False, False, 'NA,NA', 'NA', 'NA'
    S = site_dis if sj1[2] < sj1[3] else end_dis  # back-splice-junction
    S_range = [0]
    for i in range(1, S + 1):
        S_range.append(i)
        S_range.append(-i)

    for s1 in S_range:
        for s2 in S_range:
            if (force_strand and (sj1[0], force_strand=='-', sj1[2]+s1, sj1[3]+s2) in anno_sj) or \
                    (not force_strand and ((sj1[0], True, sj1[2]+s1, sj1[3]+s2) in anno_sj or \
                        (sj1[0], False, sj1[2]+s1, sj1[3]+s2) in anno_sj)):
                known = True
                anno_start, anno_end = sj1[2] + s1, sj1[3] + s2
                break
    if ref_seq:
        ref_len = len(ref_seq)
        if anno_start + 1 > ref_len or anno_end > ref_len or anno_start - 1 < 0 or anno_end - 2 < 0:
            cano = False
        else:
            S = 0 if known else S
            left_seq = ref_seq[max(0, anno_start - 1 - S): min(ref_len, anno_start + S + 1)].seq.upper()
            left_start = max(0, anno_start - 1 - S) + 1
            left_center = anno_start - left_start
            right_seq = ref_seq[max(0, anno_end - 2 - S): min(ref_len, anno_end + S)].seq.upper()
            right_start = max(0, anno_end - 2 - S) + 1
            right_center = anno_end - right_start - 1

            cano, strand, left_idx, right_idx, cano_motif = get_cano_sj(left_seq, left_center, right_seq, right_center, force_strand)
            if cano:
                dis_to_cano = '{},{}'.format(left_start + left_idx - sj1[2], right_start + right_idx - (sj1[3]-1))

    return known, cano, dis_to_cano, strand, cano_motif


def is_known_exon(exon, is_left, is_right, anno_exon, site_dis, end_dis):
    S_left = end_dis if is_left else site_dis
    S_right = end_dis if is_right else site_dis
    S_left_range, S_right_range = [0], [0]
    for i in range(1, S_left + 1):
        S_left_range.append(i)
        S_left_range.append(-i)
    for i in range(1, S_right + 1):
        S_right_range.append(i)
        S_right_range.append(-i)
    for sl in S_left_range:
        for sr in S_right_range:
            if (exon[0], True, exon[2] + sl, exon[3] + sr) in anno_exon or \
                    (exon[0], False, exon[2] + sl, exon[3] + sr) in anno_exon:
                return True
    return False


# return: is_known, is_reverse, dis_to_known
def is_known_splice_site(site, is_end, anno_site, site_dis, end_dis):
    dis_lim = end_dis if is_end else site_dis
    S_range = [0]
    for i in range(1, dis_lim + 1):
        S_range.append(i)
        S_range.append(-i)

    for s in S_range:
        if (site[0], True, site[2], site[3] + s) in anno_site or (site[0], False, site[2], site[3] + s) in anno_site:
            return True, s
    return False, 'NA'


# exon_block = [(chr, is_rev, start, end), ()...]
# return: isKnown-site, isKnown-exon,  isKnown-junc, isCano-junc, disToKnown
def comp_with_anno(exon, anno_site, anno_exon, anno_sj, ref_seq, site_dis, end_dis):
    known_site, known_exon, known_junc, cano_junc, dis_to_known, dis_to_cano, splice_motif = [], [], [], [], [], [], []
    don, acc = -1, -1
    site_i, exon_i, junc_i = 0, 0, 0
    exon_n = len(exon)
    if exon_n == 1:
        known_junc = ['NA']
        cano_junc = ['NA']
        splice_motif = ['NA']
        dis_to_cano = ['NA']
    for e in exon:
        is_left = True if exon_i == 0 else False
        is_right = True if exon_i == exon_n - 1 else False
        known_exon.append(is_known_exon(e, is_left, is_right, anno_exon, site_dis, end_dis))
        exon_i += 1
        is_known_site, dis = is_known_splice_site((e[0], e[1], True, e[2]), is_left, anno_site, site_dis, end_dis)
        known_site.append(is_known_site), dis_to_known.append(dis)
        is_known_site, dis = is_known_splice_site((e[0], e[1], False, e[3]), is_right, anno_site, site_dis, end_dis)
        known_site.append(is_known_site), dis_to_known.append(dis)

        acc = e[2] - 1
        if don > 0 and acc > 0:
            # dis_to_cano is NOT used
            k_sj, c_sj, d_to_c, strand, motif = is_known_cano_sj((e[0], e[1], don, acc), anno_sj, ref_seq, site_dis, end_dis, '')
            known_junc.append(k_sj)
            cano_junc.append(c_sj)
            dis_to_cano.append(d_to_c)
            if strand != 'NA' and motif != 'NA':
                splice_motif.append('{}{}'.format(strand, motif))
            else: splice_motif.append('NA')
        don = e[3] + 1
    return known_site, known_exon, known_junc, cano_junc, dis_to_known, dis_to_cano, splice_motif


# return all 1base coordinates
def get_coor_from_block(start0base, block_size, block_starts):
    size = block_size.split(',')
    if '' in size:
        size.remove('')
    starts = block_starts.split(',')
    if '' in starts:
        starts.remove('')
    isize = list(map(int, size))
    istarts = list(map(int, starts))
    coor = []
    for l, s in zip(isize, istarts):
        coor.append(int(start0base) + s + 1)
        coor.append(int(start0base) + s + l)
    return coor


# input 1base coordinates
# return start0base, block_size, block_starts
def get_block_from_coor(coor=[]):
    if len(coor) % 2:
        ut.err_format_time('Error in coor.')
        exit(1)
    start0base = coor[0] - 1
    block_size, block_starts = [], []
    for i, c in enumerate(coor):
        if i % 2:
            block_size.append(c - coor[i - 1] + 1)
        else:
            block_starts.append(c - start0base - 1)
    return start0base, ','.join(map(str, block_size)), ','.join(map(str, block_starts))


# TODO 5'/3' site
# site: exonic base
def bed12_to_site_bed(in_bed_fn, out_five_site_fn, out_three_site_fn, site_dis, end_dis):
    with open(in_bed_fn) as in_bed, open(out_five_site_fn, 'w') as five_out_bed, open(out_three_site_fn, 'w') as three_out_bed:
        for line in in_bed:
            if line.startswith('#'): continue
            ele = line.rsplit()
            start0base = int(ele[1])
            block_size = ele[10]
            block_starts = ele[11]
            coor = get_coor_from_block(start0base, block_size, block_starts)

            five_site_coor, three_site_coor = [], []
            for i, c in enumerate(coor):
                if i % 2:
                    if i == len(coor) - 1:
                        five_site_coor.extend([max(1, c - end_dis), c + end_dis])
                    else:
                        five_site_coor.extend([max(1, c - site_dis), c + site_dis])
                else:
                    if i == 0:
                        three_site_coor.extend([max(1, c - end_dis), c + end_dis])
                    else:
                        three_site_coor.extend([max(1, c - site_dis), c + site_dis])
            # five site
            site_start0base, site_block_size, site_block_starts = get_block_from_coor(five_site_coor)
            ele[1], ele[2] = str(site_start0base), str(five_site_coor[-1])
            ele[6], ele[7] = ele[1], ele[2]
            ele[9] = str(len(site_block_size.split(',')))
            ele[10], ele[11] = site_block_size, site_block_starts
            five_out_bed.write('\t'.join(ele) + '\n')
            # three site
            site_start0base, site_block_size, site_block_starts = get_block_from_coor(three_site_coor)
            ele[1], ele[2] = str(site_start0base), str(three_site_coor[-1])
            ele[6], ele[7] = ele[1], ele[2]
            ele[9] = str(len(site_block_size.split(',')))
            ele[10], ele[11] = site_block_size, site_block_starts
            three_out_bed.write('\t'.join(ele) + '\n')


# site: exonic base
def gtf_to_site_bed(in_gtf_fn, out_five_site_bed, out_three_site_bed, site_dis, end_dis):
    with open(in_gtf_fn) as in_fp, open(out_five_site_bed, 'w') as out_five_fp, open(out_three_site_bed, 'w') as out_three_fp:
        five_site_dict = dd(lambda: set())
        three_site_dict = dd(lambda: set())
        for line in in_fp:
            ele = line.rsplit()
            if len(ele) > 6 and ele[2] == 'exon':  # exon line
                chr, start, end, strand = ele[0], int(ele[3]), int(ele[4]), ele[6]
                id_mat = re.search('gene_id \"(\S+)\"', line)
                name_mat = re.search('gene_name \"(\S+)\"', line)
                gene_id, gene_name = id_mat.group(1) if id_mat else 'NA', name_mat.group(1) if name_mat else 'NA'
                three_site_dict[('{},{},{}'.format(gene_id, gene_name, strand), chr, strand)].add(start)
                five_site_dict[('{},{},{}'.format(gene_id, gene_name, strand), chr, strand)].add(end)
        # write site to bed
        for (id_name, chr, strand), site_set in five_site_dict.items():
            site_block_list = []
            site_list = sorted(site_set)
            for s in site_list[:-1]:
                site_block_list.extend([max(1, s - site_dis), s + site_dis])
            site_block_list.extend([max(1, site_list[-1]) - end_dis, site_list[-1] + end_dis])

            start0base, block_size, block_starts = get_block_from_coor(site_block_list)
            end = site_block_list[-1]
            out_five_fp.write(
                '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    chr, start0base, end, id_name, 0, strand,
                    start0base, end, 0, len(block_size.split(',')),
                    block_size, block_starts))
        for (id_name, chr, strand), site_set in three_site_dict.items():
            site_block_list = []
            site_list = sorted(site_set)
            site_block_list.extend([max(1, site_list[0] - end_dis), site_list[0] + end_dis])
            for s in site_list[1:]:
                site_block_list.extend([max(1, s - site_dis), s + site_dis])

            start0base, block_size, block_starts = get_block_from_coor(site_block_list)
            end = site_block_list[-1]
            out_three_fp.write(
                '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    chr, start0base, end, id_name, 0, strand,
                    start0base, end, 0, len(block_size.split(',')),
                    block_size, block_starts))


if __name__ == '__main__':
    # bed12_to_site_bed(sys.argv[1], sys.argv[2], 0, 0)
    gtf_to_site_bed(sys.argv[1], sys.argv[2], sys.argv[3], 0, 0)
