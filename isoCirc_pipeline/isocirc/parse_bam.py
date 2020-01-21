import sys, re
import pysam as ps
from collections import defaultdict as dd
import copy
from Bio import pairwise2

import isocirc.utils as ut

#### cigar operation:
BAM_CMATCH = 0  # M
BAM_CINS = 1  # I
BAM_CDEL = 2  # D
BAM_CREF_SKIP = 3  # N
BAM_CSOFT_CLIP = 4  # S
BAM_CHARD_CLIP = 5  # H
BAM_CPAD = 6  # P
BAM_CEQUAL = 7  # =
BAM_CDIFF = 8  # X
BAM_CBACK = 9  # B

cigar_op_dict = {'M': 0, 'I': 1, 'D': 2, 'N': 3, 'S': 4, 'H': 5, 'P': 6, '=': 7, 'X': 8, 'B': 9,
                 0: 'M', 1: 'I', 2: 'D', 3: 'N', 4: 'S', 5: 'H', 6: 'P', 7: '=', 8: 'X', 9: 'B'}

def is_cigar_M(c):
    if c == BAM_CMATCH or c == BAM_CEQUAL or c == BAM_CDIFF:
        return True
    else:
        return False
# cigar stats:
# M	BAM_CMATCH	0
# I	BAM_CINS	1
# D	BAM_CDEL	2
# N	BAM_CREF_SKIP	3
# S	BAM_CSOFT_CLIP	4
# H	BAM_CHARD_CLIP	5
# P	BAM_CPAD	6
# =	BAM_CEQUAL	7
# X	BAM_CDIFF	8
# B	BAM_CBACK	9
# NM	NM tag	10

# for pairwise alignment
# The match parameters are:
#
# CODE  DESCRIPTION
# x     No parameters. Identical characters have score of 1, otherwise 0.
# m     A match score is the score of identical chars, otherwise mismatch
#       score.
# d     A dictionary returns the score of any pair of characters.
# c     A callback function returns scores.
# The gap penalty parameters are:
#
# CODE  DESCRIPTION
# x     No gap penalties.
# s     Same open and extend gap penalties for both sequences.
# d     The sequences have different open and extend gap penalties.
# c     A callback function returns the gap penalties.


# align1: ref, align2: query
# return: ref_pos,cigar
def get_cigar_from_pairwise_res(r):
    (align1, align2, score, start, end) = r
    ref_pos = start - align1[:start].count('-')  # 0-base
    read_left_clip = start - align2[:start].count('-')
    read_right_clip = len(align2) - end - align2[end:].count('-')
    cigartuples = [(cigar_op_dict['S'], read_left_clip)] if read_left_clip else []
    align_str1 = align1[start:end]
    align_str2 = align2[start:end]
    for (s1, s2) in zip(align_str1, align_str2):
        if s1 == '-' and s2 != '-':  # insertion
            op = 'I'
        elif s1 != '-' and s2 == '-':  # deletion
            op = 'D'
        elif s1 == s2:  # match
            op = '='
        else:  # mismatch
            op = 'X'
        if cigartuples and cigartuples[-1][0] == cigar_op_dict[op]:
            cigartuples[-1] = (cigar_op_dict[op], cigartuples[-1][1] + 1)
        else:
            cigartuples.append((cigar_op_dict[op], 1))
    if read_right_clip: cigartuples.append((cigar_op_dict['S'], read_right_clip))
    cigarstring = cigartuples_to_cigarstring(cigartuples)
    return cigarstring

# return: (align1, align2, score, start, end)
# return: ref_pos,cigar
def pairwise_align(seq1='', seq2='', align_mode='g', cigar=False):  # match:1, mismatch:-1, no gap penalties
    # global, local
    match, mismatch, gap_open, gap_ext = 2, -4, -6, -2
    res = pairwise2.align.globalms(seq1, seq2, match, mismatch, gap_open, gap_ext) if align_mode == 'g' else pairwise2.align.localms(seq1, seq2, match, mismatch, gap_open, gap_ext)
    if len(res) == 0: return 0, 'NA'
    r = res[0]
    if cigar:
        # print r[2], get_cigar_from_pairwise_res(r)
        return r[2], get_cigar_from_pairwise_res(r)
    else:
        return r[2], r


def cigartuples_to_cigarstring(cigartuples=[]):
    cigarstring = ''
    for t in cigartuples:
        cigarstring += '{}{}'.format(t[1], cigar_op_dict[t[0]])
    return cigarstring


def cigarstring_to_cigartuples(cigarstring=''):
    S = re.sub(r'(?<=[MIDSH=X])(?=[^\s])', r' ', cigarstring)
    cigartuples = []
    for r in S.split():
        cigartuples.append((cigar_op_dict[r[-1]], int(r[:-1])))
    return cigartuples

def cigartuples_to_cigarstats(cigartuples=[]):
    cigarstats = dd(lambda: 0)
    for t in cigartuples:
        cigarstats[t[0]] += t[1]
    return cigarstats

def cigarstring_to_cigarstats(cigar=''):
    cigarstats = dd(lambda: 0)
    for c in 'MIDNHS=X':
        if c in cigar:
            count = sum(map(int, (re.findall(r'(\d+)' + c, cigar))))
            cigarstats[cigar_op_dict[c]] = count
    return cigarstats


def get_ref_op_length(r=ps.AlignedSegment):
    cigar_stats = r.get_cigar_stats()[0]
    # get op length for MDNP=X
    op_len = cigar_stats[0] + cigar_stats[2] + cigar_stats[3] + cigar_stats[6] + cigar_stats[7] + cigar_stats[8]
    return op_len


def get_read_op_length(r=ps.AlignedSegment):
    cigar_stats = r.get_cigar_stats()[0]
    # get op length for MISH=X
    op_len = cigar_stats[0] + cigar_stats[1] + cigar_stats[4] + cigar_stats[5] + cigar_stats[7] + cigar_stats[8]
    return op_len


def get_aligned_read_length(r=ps.AlignedSegment):
    cigar_stats = r.get_cigar_stats()[0]
    # get op length for MI=X
    op_len = cigar_stats[0] + cigar_stats[1] + cigar_stats[7] + cigar_stats[8]
    return op_len


def manipulate_cigar(r=ps.AlignedSegment, old='', new=''):
    r.cigarstring = re.sub(r'%s' % old, new, r.cigarstring)


def get_cigar_count(cigartuples, c):
    cnt = 0
    for tuples in cigartuples:
        if tuples[0] == c:
            cnt += 1
    return cnt


def get_cigar_len(r, c):
    cnt = 0
    cigartuples = r.cigartuples
    for tuples in cigartuples:
        if tuples[0] == c:
            cnt += tuples[1]
    return cnt


def is_unmapped(r=ps.AlignedSegment):
    return r.is_unmapped


# '20^T3AC20T', 10, 30 => ' 10^T3AC5
def get_spec_MD(mdstr='', start=0, end=0):
    mSub = re.sub(r'([ACGTNacgtn^])', ' \\1 ', mdstr)
    mSplit = re.split('[ ]+', mSub)
    if '' in mSplit:
        mSplit.remove('')
    start_remain_len = start
    end_remain_len = end - start
    ret_md = []
    mi = 0
    is_del = False
    # print mSplit, start_remain_len, end_remain_len
    while start_remain_len > 0:
        if mSplit[mi].isdigit():
            is_del = False
            if int(mSplit[mi]) > start_remain_len:
                mSplit[mi] = str(int(mSplit[mi]) - start_remain_len)
                break
            else:
                start_remain_len -= int(mSplit[mi])
        elif mSplit[mi] == '^':
            is_del = True
        else:  # isalpha()
            if not is_del:
                start_remain_len -= 1
        mi += 1
    while end_remain_len > 0:
        if mSplit[mi].isdigit():
            is_del = False
            if int(mSplit[mi]) >= end_remain_len:
                ret_md.append(str(end_remain_len))
                break
            else:
                end_remain_len -= int(mSplit[mi])
                ret_md.append(mSplit[mi])
        elif mSplit[mi] == '^':
            is_del = True
            ret_md.append(mSplit[mi])
        else:  # isalpha()
            if not is_del:
                end_remain_len -= 1
            ret_md.append(mSplit[mi])
        mi += 1
    # print ret_md
    return ret_md
    # return ''.join(ret_md)

# '15M1D5M2I10M', 10, 25 => '5M1D5M2I4M'
def get_spec_ref_cigar(cigartuples=[], start=0, end=0):  # Aligned bases
    tuples = copy.copy(cigartuples)
    start_remain_len = start
    end_remain_len = end - start
    ret_cigar = []
    mi = 0
    # print start_remain_len, end_remain_len
    # print start, end, tuples
    while start_remain_len > 0:
        if is_cigar_M(tuples[mi][0]) or tuples[mi][0] == BAM_CDEL:
            if tuples[mi][1] > start_remain_len:
                tuples[mi] = (tuples[mi][0], tuples[mi][1] - start_remain_len)
                break
            else:
                start_remain_len -= tuples[mi][1]

        mi += 1
    if end_remain_len > 0 and tuples[mi][0] == BAM_CINS:
        mi += 1
    # print tuples
    while end_remain_len > 0:
        # print end_remain_len
        if is_cigar_M(tuples[mi][0]) or tuples[mi][0] == BAM_CDEL:
            if tuples[mi][1] > end_remain_len:
                ret_cigar.append((tuples[mi][0], end_remain_len))
                break
            else:
                end_remain_len -= tuples[mi][1]
                ret_cigar.append(tuples[mi])
        else:
            ret_cigar.append(tuples[mi])
        mi += 1
    # print ret_cigar
    return ret_cigar

# '15M2D5M1I10M', 10, 25 => '5M2D5M1I4M'
def get_spec_read_cigar(cigartuples=[], start=0, end=0):  # Aligned bases
    tuples = copy.copy(cigartuples)
    start_remain_len = start
    end_remain_len = end - start
    ret_cigar = []
    mi = 0
    # print start_remain_len, end_remain_len
    # print tuples, start, end
    while start_remain_len > 0:
        # print mi
        if is_cigar_M(tuples[mi][0]) or tuples[mi][0] == BAM_CINS:
            if tuples[mi][1] > start_remain_len:
                tuples[mi] = (tuples[mi][0], tuples[mi][1] - start_remain_len)
                break
            else:
                start_remain_len -= tuples[mi][1]

        mi += 1
    # print tuples
    while end_remain_len > 0:
        # print end_remain_len
        if is_cigar_M(tuples[mi][0]) or tuples[mi][0] == BAM_CINS:
            if tuples[mi][1] > end_remain_len:
                ret_cigar.append((tuples[mi][0], end_remain_len))
                break
            else:
                end_remain_len -= tuples[mi][1]
                ret_cigar.append(tuples[mi])
        else:
            ret_cigar.append(tuples[mi])
        mi += 1
    # print ret_cigar
    return ret_cigar


# MISMATCH: read_pos(first), ref_pos(first), len, read_base, ref_base
# INSERTION: ins_read_pos(first), ins_ref_pos(left),  len, ins_base
# DELETION: del_read_pos(left), del_ref_pos(first), len, del_base
def get_error_from_MD(cigartuples=[], mdstr='', full_query_seq='', ref_start=0):
    mis, ins, dele = [], [], []
    last_error = ''
    md_i, m_pos = 0, 0
    mdSub = re.sub(r'([\\^][ACGTNacgtn]+)[0]*', ' \\1 ', mdstr)
    mdSplit = mdSub.rsplit()
    ref_pos, query_pos = ref_start, 0

    for tuples in cigartuples:
        if tuples[0] == BAM_CMATCH:
            m = mdSplit[md_i]

            if m.startswith('^'):
                ut.format_time(sys.stderr, 'get_error_from_MD', 'Unexpected MD string: {}\n'.format(mdstr))
                sys.exit(1)
            mSub = re.sub(r'([ACGTNacgtn])', ' \\1 ', m)
            m_len = sum(map(int, (re.sub(r'([ACGTNacgtn])', '1', mSub)).rsplit()))

            # from m_pos to m_pos + tuples[1]
            sub_ms = get_spec_MD(m, m_pos, m_pos + tuples[1])

            for ms in sub_ms:
                if ms.isalpha():  # MISMATCH
                    if full_query_seq[query_pos] != ms:
                        if last_error != 'MIS' or mis[-1][0] + mis[-1][2] != query_pos:
                            mis_error = [query_pos, ref_pos, 1, full_query_seq[query_pos], ms]
                            mis.append(mis_error)
                        else:  # last_error == 'MIS' and  mis[-1][2] == ap[0] - 1:
                            mis[-1][-3] += 1
                            mis[-1][-2] += full_query_seq[query_pos]
                            mis[-1][-1] += ms
                        last_error = 'MIS'
                    else:
                        ut.fatal_format_time('get_error_from_MD', 'MIS error: {} v.s {}.'.format(full_query_seq[query_pos], ms))
                    query_pos += 1
                    ref_pos += 1
                elif ms.isdigit():  # MATCH
                    query_pos += int(ms)
                    ref_pos += int(ms)

            if m_pos + tuples[1] == m_len:
                md_i += 1
                m_pos = 0
            elif m_pos + tuples[1] < m_len:
                m_pos += tuples[1]
            else:  #
                ut.format_time(sys.stderr, 'get_error_from_MD', 'Unexpected MD string: {}\n'.format(mdstr))
                sys.exit(1)
        elif tuples[0] == BAM_CDEL:
            m = mdSplit[md_i]
            if not m.startswith('^'):
                ut.format_time(sys.stderr, 'get_error_from_MD', 'Unexpected MD string: {}\n'.format(mdstr))
                sys.exit(1)
            del_error = [query_pos - 1, ref_pos, tuples[1], m[1:]]
            dele.append(del_error)
            ref_pos += tuples[1]
            last_error = 'DEL'
            md_i += 1
        elif tuples[0] == BAM_CINS:
            ins_error = [query_pos, ref_pos - 1, tuples[1], full_query_seq[query_pos:query_pos + tuples[1]]]
            ins.append(ins_error)
            query_pos += tuples[1]
            last_error = 'INS'
        elif tuples[0] == BAM_CSOFT_CLIP or tuples[0] == BAM_CHARD_CLIP:
            query_pos += tuples[1]
        elif tuples[0] == BAM_CREF_SKIP:
            ref_pos += tuples[1]
        else:
            ut.format_time(sys.stderr, 'get_error_from_MD', 'Unexpected cigar: {}\n'.format(cigartuples))
            sys.exit(1)

    return mis, ins, dele

def old_get_error_from_MD(cigartuples=[], mdstr='', full_query_seq='', ref_start=0):
    mis, ins, dele = [], [], []
    last_error = ''
    md_i, m_pos = 0, 0
    mdSub = re.sub(r'([\\^][ACGTNacgtn]+)[0]*', ' \\1 ', mdstr)
    mdSplit = mdSub.rsplit()
    ref_pos, query_pos = ref_start, 0

    for tuples in cigartuples:
        if tuples[0] == BAM_CMATCH:
            m = mdSplit[md_i]

            if m.startswith('^'):
                ut.format_time(sys.stderr, 'get_error_from_MD', 'Unexpected MD string: {}'.format(mdstr))
                sys.exit(1)
            mSub = re.sub(r'([ACGTNacgtn])', ' \\1 ', m)
            m_len = sum(map(int, (re.sub(r'([ACGTNacgtn])', '1', mSub)).rsplit()))

            # from m_pos to m_pos + tuples[1]
            sub_ms = get_spec_MD(m, m_pos, m_pos + tuples[1])

            for ms in sub_ms:
                if ms.isalpha():  # MISMATCH
                    if last_error != 'MIS' or mis[-1][0] != query_pos - 1:
                        mis_error = [query_pos, ref_pos, 1, full_query_seq[query_pos], ms]
                        mis.append(mis_error)
                    else:  # last_error == 'MIS' and  mis[-1][1] == ap[0] - 1:
                        mis[-1][-3] += 1
                        mis[-1][-2] += full_query_seq[query_pos]
                        mis[-1][-1] += ms
                    query_pos += 1
                    ref_pos += 1
                    last_error = 'MIS'
                elif ms.isdigit():  # MATCH
                    query_pos += int(ms)
                    ref_pos += int(ms)

            if m_pos + tuples[1] == m_len:
                md_i += 1
                m_pos = 0
            elif m_pos + tuples[1] < m_len:
                m_pos += tuples[1]
            else:  #
                ut.format_time(sys.stderr, 'get_error_from_MD', 'Unexpected MD string: {}'.format(mdstr))
                sys.exit(1)
        elif tuples[0] == BAM_CDEL:
            m = mdSplit[md_i]
            if not m.startswith('^'):
                ut.format_time(sys.stderr, 'get_error_from_MD', 'Unexpected MD string: {}'.format(mdstr))
                sys.exit(1)
            del_error = [query_pos - 1, ref_pos, tuples[1], m[1:]]
            dele.append(del_error)
            ref_pos += tuples[1]
            last_error = 'DEL'
            md_i += 1
        elif tuples[0] == BAM_CINS:
            ins_error = [query_pos, ref_pos - 1, tuples[1], full_query_seq[query_pos:query_pos + tuples[1]]]
            ins.append(ins_error)
            query_pos += tuples[1]
            last_error = 'INS'
        elif tuples[0] == BAM_CSOFT_CLIP or tuples[0] == BAM_CHARD_CLIP:
            query_pos += tuples[1]
        elif tuples[0] == BAM_CREF_SKIP:
            ref_pos += tuples[1]
        else:
            ut.format_time(sys.stderr, 'get_error_from_MD', 'Unexpected cigar: {}'.format(cigartuples))
            sys.exit(1)

    return mis, ins, dele


def get_error_from_Cigar(cigartuples=[], full_query_seq='', align_ref_seq='', ref_start=0):
    mis, ins, dele = [], [], []
    last_error = ''
    ref_pos, query_pos = ref_start, 0
    for tuples in cigartuples:
        if is_cigar_M(tuples[0]):
            for q, r in zip(full_query_seq[query_pos:query_pos + tuples[1]],
                            align_ref_seq[ref_pos - ref_start:ref_pos - ref_start + tuples[1]]):
                if q != r:  # MISMATCH
                    if last_error != 'MIS' or mis[-1][0] != query_pos - 1:
                        mis_error = [query_pos, ref_pos, 1, q, r]
                        mis.append(mis_error)
                    else:  # last_error == 'MIS' and  mis[-1][1] == ap[0] - 1:
                        mis[-1][-3] += 1
                        mis[-1][-2] += q
                        mis[-1][-1] += r
                    last_error = 'MIS'
                ref_pos += 1
                query_pos += 1
        elif tuples[0] == BAM_CDEL:
            del_error = [query_pos - 1, ref_pos, tuples[1],
                         align_ref_seq[ref_pos - ref_start:ref_pos - ref_start + tuples[1]]]
            dele.append(del_error)
            last_error = 'DEL'
            ref_pos += tuples[1]
        elif tuples[0] == BAM_CINS:
            ins_error = [query_pos, ref_pos - 1, tuples[1], full_query_seq[query_pos:query_pos + tuples[1]]]
            ins.append(ins_error)
            last_error = 'INS'
            query_pos += tuples[1]
        elif tuples[0] == BAM_CHARD_CLIP or tuples[0] == BAM_CSOFT_CLIP:
            query_pos += tuples[1]
        elif tuples[0] == BAM_CREF_SKIP:
            ref_pos += tuples[1]
        else:
            ut.format_time(sys.stderr, 'get_error_from_Cigar', 'Unexpected cigar: {}'.format(cigartuples))
            sys.exit(1)
    return mis, ins, dele


def get_XD_from_MD(mdstr='', start=0, end=0):
    # print mdstr, start, end
    mSub = re.sub(r'([ACGTNacgtn^])', ' \\1 ', mdstr)
    mSplit = re.split('[ ]+', mSub)
    if '' in mSplit:
        mSplit.remove('')
    start_remain_len = start
    end_remain_len = end - start
    xi, mi = [0, 0], 0

    is_del = False
    # print mSplit, start_remain_len, end_remain_len
    while start_remain_len > 0:
        if mSplit[mi].isdigit():
            is_del = False
            if int(mSplit[mi]) > start_remain_len:
                mSplit[mi] = str(int(mSplit[mi]) - start_remain_len)
                break
            else:
                start_remain_len -= int(mSplit[mi])
        elif mSplit[mi] == '^':
            is_del = True
        else:  # isalpha()
            if not is_del:
                start_remain_len -= 1
        mi += 1
    while end_remain_len > 0:
        if mSplit[mi].isdigit():
            is_del = False
            if int(mSplit[mi]) >= end_remain_len:
                break
            else:
                end_remain_len -= int(mSplit[mi])
        elif mSplit[mi] == '^':
            is_del = True
        else:  # isalpha()
            if not is_del:
                end_remain_len -= 1
                xi[0] += 1
            else:
                xi[1] += 1
        mi += 1
    return xi


def get_I_from_cigar(cigartuples=[], start=0, end=0):  # Aligned bases
    tuples = copy.copy(cigartuples)
    start_remain_len = start
    end_remain_len = end - start
    ins, mi = 0, 0

    # print start_remain_len, end_remain_len
    while start_remain_len > 0:
        if is_cigar_M(tuples[mi][0]) or tuples[mi][0] == BAM_CINS:
            if tuples[mi][1] > start_remain_len:
                tuples[mi] = (BAM_CMATCH, tuples[mi][1] - start_remain_len)
                break
            else:
                start_remain_len -= tuples[mi][1]

        mi += 1
    # print tuples
    while end_remain_len > 0:
        # print end_remain_len
        if is_cigar_M(tuples[mi][0]):
            if tuples[mi][1] > end_remain_len:
                break
            else:
                end_remain_len -= tuples[mi][1]
        elif tuples[mi][0] == BAM_CINS:
            if tuples[mi][1] > end_remain_len:
                ins += end_remain_len
                break
            else:
                end_remain_len -= tuples[mi][1]
                ins += tuples[mi][1]
        mi += 1
    # print ret_cigar
    return ins


# MISMATCH: read_pos(first), ref_pos(first), len, read_base, ref_base
# INSERTION: ins_read_pos(first), ins_ref_pos(left),  len, ins_base
# DELETION: del_read_pos(left), del_ref_pos(first), len, del_base

def get_block(r=ps.AlignedSegment):
    tid = r.reference_id
    is_rev = r.is_reverse
    cigartuples = r.cigartuples
    ref_start = r.reference_start + 1

    block = []
    start = ref_start
    end = ref_start - 1
    for tuples in cigartuples:
        if is_cigar_M(tuples[0]) or tuples[0] == BAM_CDEL:
            end += tuples[1]
        elif tuples[0] == BAM_CREF_SKIP:
            block.append((tid, is_rev, start, end))
            start = end + tuples[1] + 1
            end = start - 1
    block.append((tid, is_rev, start, end))
    return block


def get_exon_block(tid, is_rev, cigartuples=[], ref_start=0):
    exon_block = []
    start = ref_start
    end = ref_start - 1
    for tuples in cigartuples:
        if is_cigar_M(tuples[0]) or tuples[0] == BAM_CDEL:
            end += tuples[1]
        elif tuples[0] == BAM_CREF_SKIP:
            exon_block.append((tid, is_rev, start, end))
            start = end + tuples[1] + 1
            end = start - 1
    exon_block.append((tid, is_rev, start, end))
    return exon_block


# sj = [tid, is_rev, 1-base_intron_start, 1-base_intron_end]
def get_splice_junction(r=ps.AlignedSegment, include_bsj=False):
    sj = []
    tid, is_rev, _start = r.reference_id, r.is_reverse, r.reference_start + 1
    start = _start
    for tuples in r.cigartuples:
        if tuples[0] == BAM_CREF_SKIP:
            end = start + tuples[1] - 1
            sj.append((tid, is_rev, start, end))
            start = end + 1
        elif is_cigar_M(tuples[0]) or tuples[0] == BAM_CDEL:
            start += tuples[1]
    if include_bsj:
        end = _start - 1
        sj.append((tid, is_rev, start, end))
    return sj


def get_error_rate(in_sam_fn=''):
    tot_mapped_n, tot_unmapped_n, tot_mapped_base, tot_ins, tot_del, tot_mis, tot_match = 0, 0, 0, 0, 0, 0, 0
    with ps.AlignmentFile(in_sam_fn) as in_sam:
        for r in in_sam:
            # r = ps.AlignedSegment
            if r.is_unmapped:
                tot_unmapped_n += 1
                continue
            elif r.is_secondary or r.is_supplementary:
                continue
            tot_mapped_n += 1
            tot_mapped_base += r.query_alignment_length
            cigar_stats = cigarstring_to_cigarstats(r.cigarstring)
            tot_match += cigar_stats[cigar_op_dict['=']]
            tot_ins += cigar_stats[cigar_op_dict['I']]
            tot_del += cigar_stats[cigar_op_dict['D']]
            tot_mis += cigar_stats[cigar_op_dict['X']]
    return tot_mapped_n, tot_mapped_base, '{0:.1f}%'.format((tot_ins+tot_del+tot_mis) / (tot_ins+tot_mis+tot_match+0.0) * 100)

def get_xid_from_cigar(cigar_stats):
    return cigar_stats[cigar_op_dict['X']] + cigar_stats[cigar_op_dict['I']] + cigar_stats[cigar_op_dict['D']]

def parse_sj_alignment(r=ps.AlignedSegment, flank_len=10, min_xid=1, key_flank_len=2, key_min_xid=0):
    is_high_sj = []
    last_read_end, read_end, aln_read_len = 0, 0, get_aligned_read_length(r)
    for tuples in r.cigartuples:
        if tuples[0] == BAM_CREF_SKIP:
            # left
            if read_end - last_read_end < flank_len:
                is_high_sj.append(False)
                last_read_end = read_end
                continue
            start, key_start, end = read_end - flank_len, read_end - key_flank_len, read_end
            left_cigar_stats = cigartuples_to_cigarstats(get_spec_read_cigar(r.cigartuples, start, end))
            key_left_cigar_stats = cigartuples_to_cigarstats(get_spec_read_cigar(r.cigartuples, key_start, end))
            # right
            start, end, key_end = read_end, read_end + flank_len, read_end + key_flank_len
            if end > aln_read_len:
                is_high_sj.append(False)
                last_read_end = read_end
                continue
            right_cigar_stats = cigartuples_to_cigarstats(get_spec_read_cigar(r.cigartuples, start, end))
            key_right_cigar_stats = cigartuples_to_cigarstats(get_spec_read_cigar(r.cigartuples, start, key_end))

            if get_xid_from_cigar(left_cigar_stats) + get_xid_from_cigar(right_cigar_stats) > min_xid or\
                    get_xid_from_cigar(key_left_cigar_stats) + get_xid_from_cigar(key_right_cigar_stats) > key_min_xid:
                is_high_sj.append(False)
            else:
                is_high_sj.append(True)
            last_read_end = read_end
        elif is_cigar_M(tuples[0]) or tuples[0] == BAM_CINS:
            read_end += tuples[1]
    if not is_high_sj:
        is_high_sj.append('NA')
    return is_high_sj
