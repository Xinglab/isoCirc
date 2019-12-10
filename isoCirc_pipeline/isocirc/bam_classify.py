import argparse
import pysam as ps
import sys

import isocirc.parse_bam as pb
import isocirc.utils as ut

high_max_ratio = 1.1
high_min_ratio = 0.9
high_iden_ratio = 0.75 #TODO old: 0.75
high_repeat_ratio = 0.6
low_repeat_ratio = 1.9

max_ins_len = 50 # max ins/del length within a single block
max_del_len = 50


# map / cons
# iden / align
# best-sec / best
# map / cons ~= 2 => increase copy number to get full alignment
# chimeric
# partial mapped, map / cons < 0.9

def get_iden_ratio(r):
    new_block = 1
    ins_len, del_len = 0, 0
    for tuples in r.cigartuples:
        if tuples[0] == pb.BAM_CINS:
            ins_len += tuples[1]
            if ins_len > max_ins_len: return -1.0
        elif tuples[0] == pb.BAM_CDEL:
            del_len += tuples[1]
            if del_len > max_del_len: return -1.0
        elif tuples[0] == pb.BAM_CREF_SKIP:
            ins_len, del_len = 0, 0
    map_len = pb.get_aligned_read_length(r)
    if not r.has_tag('NM'): ut.fatal_format_time('bam_classify', 'No NM tag found.\n')
    NM = int(r.get_tag('NM'))
    del_len = int(pb.get_cigar_len(r, pb.BAM_CDEL))
    iden_len = map_len - NM + del_len
    return iden_len / (map_len + 0.0)


#if any returnable record is aligned to non-primary chromosome, discard all the records
#TODO secondary / best < ratio ???
def high_qual_record(r_array, high_max_ratio=high_max_ratio, high_min_ratio=high_min_ratio, high_iden_ratio=high_iden_ratio):
    if not r_array: return None
    primary_r = r_array[0]
    # primary_r = ps.AlignedSegment()
    if primary_r.is_secondary or primary_r.is_supplementary or primary_r.is_unmapped:
        ut.fatal_format_time('high_qual_record', 'Error: input SAM file is sorted or modified.')
    primary_start = primary_r.reference_start + 1
    primary_end = primary_start + pb.get_ref_op_length(primary_r) - 1
    rlen = 0.0 + pb.get_read_op_length(primary_r)
    cons_len = rlen / 2

    best_i = -1
    best_r = None
    best_AS = -1
    best_iden_ratio = -1.0
    primary_is_high = False
    for i, r in iter(enumerate(r_array)):
        map_len = pb.get_aligned_read_length(r)
        mc = map_len / cons_len
        iden_ratio = get_iden_ratio(r)
        AS = int(r.get_tag('AS'))
        if high_min_ratio <= mc <= high_max_ratio and iden_ratio >= high_iden_ratio:
            if len(r.reference_name) >= 6 or r.reference_name.startswith('chrM') or r.reference_name.startswith('chrUn'):
                return None
            if AS > best_AS:
                if i == 0:
                    primary_is_high = True
                    best_r, best_i, best_AS = r, i, AS
                # if r is not primary record, r has to NOT overlap with primary r
                elif r.reference_name != primary_r.reference_name or (r.reference_start + 1 > primary_end or r.reference_start + pb.get_ref_op_length(r) < primary_start):
                    best_r, best_i, best_AS = r, i, AS

    if best_i == -1:
        return None
    else:
        return primary_r if primary_is_high else best_r

# cons is 2/3... copies of the real cons
# 2 or more parts are mapped to one same mapped positions and same mapped length
def cons_repeat_record(r_array, high_iden_ratio=high_iden_ratio, high_repeat_ratio=high_repeat_ratio):
    map_pos, map_len = [], []
    rlen = 0.0 + pb.get_read_op_length(r_array[0])
    cons_len = rlen / 2
    for r in r_array:
        map_pos.append(r.reference_start + 1)
        map_len.append(pb.get_aligned_read_length(r))

    for i in range(len(map_pos) - 1):
        for j in range(i + 1, len(map_pos)):
            if abs(map_pos[i] - map_pos[j]) <= 10 and \
                    abs(map_len[i] - map_len[j]) <= 10 and map_len[i] / cons_len <= high_repeat_ratio:
                if get_iden_ratio(r_array[i]) >= high_iden_ratio and get_iden_ratio(r_array[j]) >= high_iden_ratio:
                    return r_array[i]
    return None


def cons_self_tandem_record(r_array, low_repeat_ratio=low_repeat_ratio):
    rlen = 0.0 + pb.get_read_op_length(r_array[0])
    cons_len = rlen / 2
    for r in r_array:
        map_len = pb.get_aligned_read_length(r)
        if map_len / cons_len >= low_repeat_ratio:
            return r
    return None


# same part of read is mapped to two or more genomic positions
def is_cons_partial(r_array):
    query_pos, map_len = [], []
    for r in r_array:
        qstart = (0 if r.cigartuples[-1][0] != 4 else r.cigartuples[-1][1] + 1) if r.is_reverse else r.query_alignment_start + 1
        query_pos.append(qstart)
        map_len.append(pb.get_aligned_read_length(r))
    for i in range(1, len(query_pos)):
        if abs(query_pos[i] - query_pos[i - 1]) > 10:
            return False
        elif abs(map_len[i] - map_len[i - 1]) > 10:
            return False
    return True


# r_array: records of one read
def cla_record(r_array, high_bam, low_bam, high_max_ratio=high_max_ratio, high_min_ratio=high_min_ratio, high_iden_ratio=high_iden_ratio, high_repeat_ratio=high_repeat_ratio, low_repeat_ratio=low_repeat_ratio):
    if not r_array: return

    if len(r_array) == 1:  # single record
        rlen = 0.0 + pb.get_read_op_length(r_array[0])
        cons_len = rlen / 2

        map_len = pb.get_aligned_read_length(r_array[0])
        mc = map_len / cons_len
        if high_min_ratio <= mc <= high_max_ratio and get_iden_ratio(r_array[0]) >= high_iden_ratio:
            high_bam.write(r_array[0])
        elif mc >= low_repeat_ratio:
            low_bam.write(r_array[0])
        else:
            low_bam.write(r_array[0])
    else:  # multiple records
        high_r = high_qual_record(r_array, high_max_ratio, high_min_ratio, high_iden_ratio)
        if high_r:
            high_bam.write(high_r)
        else:
            rep_r = cons_repeat_record(r_array, high_iden_ratio, high_repeat_ratio)
            if rep_r:
                # high_bam.write(rep_r) # TODO
                low_bam.write(rep_r)
                return
            self_r = cons_self_tandem_record(r_array, low_repeat_ratio)
            if self_r:
                low_bam.write(self_r)
            else:
                for r in r_array:
                    if not (r.is_secondary or r.is_supplementary):
                        low_bam.write(r)
                        break


def bam_classify(in_bam_fn, high_bam_fn, low_bam_fn, 
                high_max_ratio=high_max_ratio, high_min_ratio=high_min_ratio,
                high_iden_ratio=high_iden_ratio, high_repeat_ratio=high_repeat_ratio, low_repeat_ratio=low_repeat_ratio):
    ut.err_format_time('classify_bam_core', 'Processing {} ... '.format(in_bam_fn))
    with ps.AlignmentFile(in_bam_fn) as in_bam, ps.AlignmentFile(high_bam_fn, 'wb', template=in_bam) as high_bam, \
            ps.AlignmentFile(low_bam_fn, 'wb', template=in_bam) as low_bam:
        cnt = 0
        r_array = []
        r_name = ''
        for r in in_bam:
            if pb.is_unmapped(r) : continue
            if r.query_name != r_name:
                cla_record(r_array, high_bam, low_bam, high_max_ratio, high_min_ratio, high_iden_ratio, high_repeat_ratio, low_repeat_ratio)
                r_name = r.query_name
                r_array = [r]
            else:
                r_array.append(r)
            cnt += 1
            if cnt % 100000 == 0:
                ut.err_format_time('classify_bam_core', '{} BAM records done ... '.format(cnt))
        if r_array:
            cla_record(r_array, high_bam, low_bam, high_max_ratio, high_min_ratio, high_iden_ratio, high_repeat_ratio, low_repeat_ratio)
    ut.err_format_time('classify_bam_core', 'Processing {} done.'.format(in_bam_fn))
    return high_bam_fn, low_bam_fn


def bam_classify_core(args):
    bam_classify(args.cons_all_sam, args.high_bam, args.low_bam,
        args.high_max_ratio, args.high_min_ratio, args.high_iden_ratio, args.high_repeat_ratio, args.low_repeat_ratio)


def parser_argv():
    # parse command line arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="Classify alignment file of consensus sequence into high quality and low quality")
    parser.add_argument("cons_all_sam", metavar='cons.sam', type=str, help='Unsorted SAM alignment file of consensus sequence.')
    parser.add_argument('high_bam', metavar='high.bam', type=str, help='High-quality BAM alignment.')
    parser.add_argument('low_bam', metavar='low.bam', type=str, help='Low-quality BAM alignment.')

    parser.add_argument('-F', '--high-max-ratio', type=float, help='Maximum mappedLen / consLen ratio for high-quality alignment.', default=high_max_ratio)
    parser.add_argument('-f', '--high-min-ratio', type=float, help='Minimum mappedLen /consLen ratio for high-quality alignment.', default=high_min_ratio)
    parser.add_argument('-d', '--high-iden-ratio', type=float, help='Minimum identicalBases/ consLen ratio for high-quality alignment.', default=high_iden_ratio)
    parser.add_argument('-s', '--high-repeat-ratio', type=float, help='Maximum mappedLen / consLen ratio for high-quality self-tandem consensus.', default=high_repeat_ratio)
    parser.add_argument('-t', '--low-repeat-ratio', type=float, help='Minimum mappedLen / consLen ratio for low-quality self-tandem alignment.', default=low_repeat_ratio)

    return parser.parse_args()

def main():
    args = parser_argv()
    bam_classify_core(args)

if __name__ == '__main__':
    main()

