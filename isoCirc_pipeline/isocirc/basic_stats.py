import argparse
import sys
import os
import parse_bam as pb
import pysam
import utils as ut
import isocirc_stats as ps
import isocirc
from collections import defaultdict as dd
from __init__ import __program__
from __init__ import __version__

isoform_output_header = isocirc.isoform_output_header
isoform_output_header_idx = isocirc.isoform_output_header_idx

def get_error_rate(in_sam_fn=''):
    map_read_name = dict()
    tot_mapped_read_n, tot_mapped_cons_n, tot_unmapped_n, tot_mapped_base, tot_ins, tot_del, tot_mis, tot_match = 0, 0, 0, 0, 0, 0, 0, 0
    with pysam.AlignmentFile(in_sam_fn) as in_sam:
        for r in in_sam:
            if r.is_unmapped:
                tot_unmapped_n += 1
                continue
            elif r.is_secondary or r.is_supplementary:
                continue
            read_name = r.query_name.rsplit('_cons')[0]
            if read_name not in map_read_name:
                map_read_name[read_name] = 1
                tot_mapped_read_n += 1
            tot_mapped_cons_n += 1
            tot_mapped_base += r.query_alignment_length
            cigar_stats = pb.cigarstring_to_cigarstats(r.cigarstring)
            tot_match += cigar_stats[pb.cigar_op_dict['=']]
            tot_ins += cigar_stats[pb.cigar_op_dict['I']]
            tot_del += cigar_stats[pb.cigar_op_dict['D']]
            tot_mis += cigar_stats[pb.cigar_op_dict['X']]
    return tot_mapped_read_n, tot_mapped_cons_n, tot_mapped_base, '{0:.1f}%'.format((tot_ins+tot_del+tot_mis) / (tot_ins+tot_mis+tot_match+0.0) * 100)
# Read number
# Base-pair number:
# Cons. read number:
# Cons. sequence number:
# Mappable cons. number:
# Error rate:

# TODO only for pacbio
def pb_stats_core(long_read_len, cons_info, cons_bam, isoform_out, all_bsj_stats_dict, stats_out):
    ut.err_format_time('basic_stats_core', 'Writing basic stats to file ... ')
    # basic stats of read/cons
    tot_subread_n, tot_polyread_n, tot_base, tot_cons_read_n, tot_cons_n, tot_cons_base, tot_map_cons_n, tot_map_cons_base, error_rate = 0, 0, 0, 0, 0, 0, 0, 0, '0.0%'
    poly_names = dict()
    with open(long_read_len, 'r') as in_fp:
        for line in in_fp:
            [name, read_len] = line.rsplit()
            read_len = int(read_len)
            tot_subread_n += 1
            tot_base += read_len
            if len(name.rsplit('/')) > 1:
                poly_name = name.rsplit('/')[0]+ '/' + name.rsplit('/')[1]
                poly_names[poly_name] = 1
    tot_polyread_n = len(poly_names)
    tot_map_cons_n, tot_map_cons_base, error_rate = pb.get_error_rate(cons_bam)
    with open(cons_info) as in_fp:
        cons_names = dict()
        for line in in_fp:
            ele = line.rsplit()
            cons_name = ele[0].rsplit('_cons')[0]
            cons_names[cons_name] = 1
            tot_cons_n += 1
            tot_cons_base += int(ele[2])
    tot_cons_read_n = len(cons_names)

    # detailed stats of circRNA
    # tot_known_circRNA
    tot_circRNA_read_n, tot_isoform, tot_known_bsj, tot_known_bsj_subread_n, tot_known_bsj_polyread_n, tot_cano_bsj, tot_cano_bsj_subread_n, tot_cano_bsj_polyread_n = 0, 0, 0, 0, 0, 0, 0, 0
    with open(isoform_out) as in_fp:
        for line in in_fp:
            if line.startswith('#'): continue
            tot_isoform += 1
            ele = line.rsplit()
            tot_circRNA_read_n += int(ele[isoform_output_header_idx['readCount']])
            if ele[isoform_output_header_idx['isKnownBSJ']] == 'True':
                tot_known_bsj += 1
                subreads_cnt = ps.get_subreads_cnt(ele[isoform_output_header_idx['readIDs']])
                polyreads_cnt = ps.get_polyreads_cnt(ele[isoform_output_header_idx['readIDs']])
                tot_known_bsj_subread_n += subreads_cnt
                tot_known_bsj_polyread_n += polyreads_cnt
            elif ele[isoform_output_header_idx['isCanoBSJ']] == 'True':
                tot_cano_bsj += 1
                subreads_cnt = ps.get_subreads_cnt(ele[isoform_output_header_idx['readIDs']])
                polyreads_cnt = ps.get_polyreads_cnt(ele[isoform_output_header_idx['readIDs']])
                tot_cano_bsj_subread_n += subreads_cnt
                tot_cano_bsj_polyread_n += polyreads_cnt

    with open(stats_out, 'w') as out:
        out.write('Total_subread_number\t{}\n'.format(tot_subread_n))
        out.write('Total_polyread_number\t{}\n'.format(tot_polyread_n))
        out.write('Total_base\t{}\n'.format(tot_base))
        out.write('Total_cons_read_number\t{}\n'.format(tot_cons_read_n))
        out.write('Total_cons_seq_number\t{}\n'.format(tot_cons_n))
        out.write('Total_cons_seq_base\t{}\n'.format(tot_cons_base))
        out.write('Total_mappable_cons_number\t{}\n'.format(tot_map_cons_n))
        out.write('Total_mappable_cons_base\t{}\n'.format(tot_map_cons_base))
        out.write('Error_rate\t{}\n'.format(error_rate))
        out.write('Total_circRNA_read_number\t{}\n'.format(tot_circRNA_read_n))
        out.write('Total_isoform\t{}\n'.format(tot_isoform))
        out.write('Total_isoform_with_known_BSJ\t{}\n'.format(tot_known_bsj))
        out.write('Total_subread_with_known_BSJ\t{}\n'.format(tot_known_bsj_subread_n))
        out.write('Total_polyread_with_known_BSJ\t{}\n'.format(tot_known_bsj_polyread_n))
        out.write('Total_isoform_with_unknown_cano_BSJ\t{}\n'.format(tot_cano_bsj))
        out.write('Total_subread_with_unknown_cano_BSJ\t{}\n'.format(tot_cano_bsj_subread_n))
        out.write('Total_polyread_cano_BSJ\t{}\n'.format(tot_cano_bsj_polyread_n))
    ut.err_format_time('basic_stats_core', 'Writing basic stats to file done!')

def stats_core(long_read_len, cons_info, cons_bam, isoform_out, all_bsj_stats_dict, stats_out):
    # basic stats of read/cons
    tot_read_n, tot_base, tot_read_cons_n, tot_cons_n, tot_cons_base, tot_map_cons_n, tot_map_cons_base, error_rate = 0, 0, 0, 0, 0, 0, 0, '0.0%'
    with open(long_read_len, 'r') as in_fp:
        for line in in_fp:
            [name, read_len] = line.rsplit()
            read_len = int(read_len)
            tot_read_n += 1
            tot_base += read_len
    tot_map_read_n, tot_map_cons_n, tot_map_cons_base, error_rate = get_error_rate(cons_bam)
    with open(cons_info) as in_fp:
        cons_names = dict()
        for line in in_fp:
            ele = line.rsplit()
            cons_name = ele[0].rsplit('_cons')[0]
            cons_names[cons_name] = 1
            tot_cons_n += 1
            tot_cons_base += int(ele[2])
    tot_read_cons_n = len(cons_names)

    # detailed stats of circRNA
    # tot_known_circRNA
    tot_isoform, tot_bsj, tot_known_bsj, tot_circRNA_read_n, tot_iso_with_known_bsj, tot_known_bsj_read_n, tot_iso_with_cano_bsj, tot_cano_bsj_read_n = 0, 0, dd(lambda:0), 0, 0, 0, 0, 0
    tot_iso_with_cano_sj, tot_read_with_cano_sj = 0, 0
    tot_iso_with_high_sj, tot_iso_with_known_ss, tot_iso_with_high_sj_known_ss = 0, 0, 0
    tot_read_with_high_sj, tot_read_with_known_ss, tot_read_with_high_sj_known_ss = 0, 0, 0
    tot_full_iso, tot_full_read = 0, 0
    tot_full_iso_bsj_fsm_iso, tot_full_iso_bsj_fsm_read, tot_full_iso_bsj_nic_iso, tot_full_iso_bsj_nic_read, tot_full_iso_bsj_nnc_iso, tot_full_iso_bsj_nnc_read = 0,0,0,0,0,0
    tot_full_iso_int_fsm_iso, tot_full_iso_int_fsm_read, tot_full_iso_int_nic_iso, tot_full_iso_int_nic_read, tot_full_iso_int_nnc_iso, tot_full_iso_int_nnc_read = 0,0,0,0,0,0
    bsj_dict = dict()
    with open(isoform_out) as in_fp:
        for line in in_fp:
            if line.startswith('#'): continue
            tot_isoform += 1
            ele = line.rsplit()
            read_cnt = int(ele[isoform_output_header_idx['readCount']])
            bsj = (ele[isoform_output_header_idx['chrom']], ele[isoform_output_header_idx['startCoor0base']], ele[isoform_output_header_idx['endCoor']])
            if bsj not in bsj_dict:
                bsj_dict[bsj] = 1
                tot_bsj += 1
                for i, known_bsj in enumerate(ele[isoform_output_header_idx['isKnownBSJ']].rsplit(',')):
                    tot_known_bsj[i] += (known_bsj == 'True')
                tot_known_bsj[i+1] += ('False' not in ele[isoform_output_header_idx['isKnownBSJ']])
            tot_circRNA_read_n += read_cnt
            if 'False' not in ele[isoform_output_header_idx['isKnownSS']]:
                tot_read_with_known_ss += read_cnt
                tot_iso_with_known_ss += 1
            full_len = False
            if 'False' not in ele[isoform_output_header_idx['isCanoSJ']]:
                tot_iso_with_cano_sj += 1
                tot_read_with_cano_sj += read_cnt
                if 'False' not in ele[isoform_output_header_idx['isHighSJ']]:
                    tot_iso_with_high_sj += 1
                    tot_read_with_high_sj += read_cnt
                    full_len = True
                    if 'False' not in ele[isoform_output_header_idx['isKnownSS']]:
                        tot_iso_with_high_sj_known_ss += 1
                        tot_read_with_high_sj_known_ss += read_cnt
            if 'False' not in ele[isoform_output_header_idx['isKnownSS']]:
                full_len = True
            if full_len:
                tot_full_iso += 1
                tot_full_read += read_cnt
                if ele[isoform_output_header_idx['BSJCate']] == 'FSM':
                    tot_full_iso_bsj_fsm_iso += 1
                    tot_full_iso_bsj_fsm_read += read_cnt
                elif ele[isoform_output_header_idx['BSJCate']] == 'NIC':
                    tot_full_iso_bsj_nic_iso += 1
                    tot_full_iso_bsj_nic_read += read_cnt
                else:
                    tot_full_iso_bsj_nnc_iso += 1
                    tot_full_iso_bsj_nnc_read += read_cnt

                if ele[isoform_output_header_idx['interIsoCate']] == 'FSM':
                    tot_full_iso_int_fsm_iso += 1
                    tot_full_iso_int_fsm_read += read_cnt
                elif ele[isoform_output_header_idx['interIsoCate']] == 'NIC':
                    tot_full_iso_int_nic_iso += 1
                    tot_full_iso_int_nic_read += read_cnt
                else:
                    tot_full_iso_int_nnc_iso += 1
                    tot_full_iso_int_nnc_read += read_cnt

    ut.err_format_time('basic_stats_core', 'Writing basic stats to file ... ')
    with open(stats_out, 'w') as out:
        out.write('#' + __program__ + '\t' + __version__ + '\n')
        out.write('# cons=consensus sequence, high=high-confidence, cano=canonical, BSJ=back-splice junction, SJ=splice junction, SS=splice site\n')
        out.write('1_Total_reads\t{:,}\n'.format(tot_read_n))
        # out.write('Total_base\t{:,}\n'.format(tot_base))
        out.write('2_Total_reads_with_cons\t{:,}\n'.format(tot_read_cons_n))
        out.write('3_Total_mappable_reads_with_cons\t{:,}\n'.format(tot_map_read_n))
        out.write('4_Total_reads_with_candidate_BSJs\t{:,}\n'.format(all_bsj_stats_dict['read_with_bsj_n']))
        out.write('5_Total_candidate_BSJs\t{:,}\n'.format(all_bsj_stats_dict['bsj_n']))
        if len(all_bsj_stats_dict['known_bsj_n']) > 2:
            for i in all_bsj_stats_dict['known_bsj_n']:
                idx = 'all' if i == len(all_bsj_stats_dict['known_bsj_n']) - 1 else i
                out.write('6_Total_candidate_BSJs_known_in_{}\t{:,}\n'.format(idx, all_bsj_stats_dict['known_bsj_n'][i]))
        else:
            out.write('6_Total_known_candidate_BSJs\t{:,}\n'.format(all_bsj_stats_dict['known_bsj_n'][0]))
        out.write('7_Total_reads_with_high_confidence_BSJs\t{:,}\n'.format(tot_circRNA_read_n))
        out.write('8_Total_high_confidence_BSJs\t{:,}\n'.format(tot_bsj))
        if len(tot_known_bsj) > 2:
            for i in tot_known_bsj:
                idx = 'all' if i == len(tot_known_bsj) - 1 else i
                out.write('9_Total_high_confidence_BSJs_known_in_{}\t{:,}\n'.format(idx, tot_known_bsj[i]))
        else:
            out.write('9_Total_known_high_confidence_BSJs\t{:,}\n'.format(tot_known_bsj[0]))
        out.write('10_Total_isoforms_with_high_BSJs\t{:,}\n'.format(tot_isoform))
        out.write('11_Total_isoforms_with_high_BSJs_cano_SJs\t{:,}\n'.format(tot_iso_with_cano_sj))
        out.write('12_Total_isoforms_with_high_BSJs_high_SJs\t{:,}\n'.format(tot_iso_with_high_sj))
        out.write('13_Total_isoforms_with_high_BSJ_known_SSs\t{:,}\n'.format(tot_iso_with_known_ss))
        out.write('14_Total_isoforms_with_high_BSJs_high_SJs_known_SSs\t{:,}\n'.format(tot_iso_with_high_sj_known_ss))
        out.write('15_Total_full_length_isoforms\t{:,}\n'.format(tot_full_iso))
        out.write('16_Total_reads_for_full_length_isoforms\t{:,}\n'.format(tot_full_read))
        # FSM/NIC/NNC
        out.write('17_Total_full_length_isoforms_with_FSM_BSJ\t{:,}\n'.format(tot_full_iso_bsj_fsm_iso))
        out.write('18_Total_reads_for_full_length_isoforms_with_FSM_BSJ\t{:,}\n'.format(tot_full_iso_bsj_fsm_read))
        out.write('19_Total_full_length_isoforms_with_NIC_BSJ\t{:,}\n'.format(tot_full_iso_bsj_nic_iso))
        out.write('20_Total_reads_for_full_length_isoforms_with_NIC_BSJ\t{:,}\n'.format(tot_full_iso_bsj_nic_read))
        out.write('21_Total_full_length_isoforms_with_NNC_BSJ\t{:,}\n'.format(tot_full_iso_bsj_nnc_iso))
        out.write('22_Total_reads_for_full_length_isoforms_with_NNC_BSJ\t{:,}\n'.format(tot_full_iso_bsj_nnc_read))

        out.write('23_Total_full_length_and_internal_FSM_isoforms\t{:,}\n'.format(tot_full_iso_int_fsm_iso))
        out.write('24_Total_reads_full_length_and_internal_FSM_isoforms\t{:,}\n'.format(tot_full_iso_int_fsm_read))
        out.write('25_Total_full_length_and_internal_NIC_isoforms\t{:,}\n'.format(tot_full_iso_int_nic_iso))
        out.write('26_Total_reads_for_full_length_and_internal_NIC_isoforms\t{:,}\n'.format(tot_full_iso_int_nic_read))
        out.write('27_Total_full_length_and_internal_NNC_isoforms\t{:,}\n'.format(tot_full_iso_int_nnc_iso))
        out.write('28_Total_reads_for_full_length_and_internal_NNC_isoforms\t{:,}\n'.format(tot_full_iso_int_nnc_read))
    ut.err_format_time('basic_stats_core', 'Writing basic stats to file done!')

def basic_stats(args):
    all_bsj_stats_dict = dd(lambda:0)
    all_bsj_stats_dict['known_bsj_n'] = [0]
    # if args.type == 'ont':
    stats_core(args.long_read_len, args.cons_info, args.cons_bam, args.out, all_bsj_stats_dict, args.stats_out)
    # elif args.type == 'pb':
    #     pb_stats_core(args.long_read_len, args.cons_info, args.cons_bam, args.out, all_bsj_stats_dict, args.stats_out)

def parser_argv():
    # parse command line arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Generate stats result of {} result".format(__program__))
    parser.add_argument('long_read_len', metavar='long_read.len', type=str, help='Read length file of long read.')
    parser.add_argument('cons_info', metavar='cons.info', type=str, help='Information file of extracted consensus sequence.')
    parser.add_argument('cons_bam', metavar='cons.bam', type=str, help='Alignment file of consensus sequence.')
    parser.add_argument('out', metavar='{}.out'.format(__program__), type=str, help='Isoform-wise circRNA output file of {} result.'.format(__program__))
    parser.add_argument('stats_out', metavar='{}_stats.out'.format(__program__), type=str, help='Stats output file of {} result.'.format(__program__))
    # parser.add_argument('--type', type=str, help='Type of sequencing data: Oxford Nanopore(ont) or Pacific Biosciences (pb).', choices=['ont', 'pb'], default='ont')

    return parser.parse_args()

if __name__ == '__main__':
    args = parser_argv()
    basic_stats(args)