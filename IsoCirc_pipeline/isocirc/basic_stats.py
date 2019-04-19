import argparse
import sys
import os
import parse_bam as pb
import utils as ut
import isocirc_stats as ps
import isocirc
from __init__ import __program__

isoform_output_header = isocirc.isoform_output_header
isoform_output_header_idx = isocirc.isoform_output_header_idx

# Read number
# Base-pair number:
# Cons. read number:
# Cons. sequence number:
# Mappable cons. number:
# Error rate:

# TODO only for pacbio
def pb_stats_core(long_read_len, cons_info, cons_bam, isoform_out, stats_out):
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
			tot_circRNA_read_n += int(ele[isoform_output_header_idx['readCount']])
			ele = line.rsplit()
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
		out.write('Total_isoform_with_cano_BSJ\t{}\n'.format(tot_cano_bsj))
		out.write('Total_subread_with_cano_BSJ\t{}\n'.format(tot_cano_bsj_subread_n))
		out.write('Total_polyread_cano_BSJ\t{}\n'.format(tot_cano_bsj_polyread_n))
	ut.err_format_time('basic_stats_core', 'Writing basic stats to file done!')	

def stats_core(long_read_len, cons_info, cons_bam, isoform_out, stats_out):
	ut.err_format_time('basic_stats_core', 'Writing basic stats to file ... ')
	# basic stats of read/cons
	tot_read_n, tot_base, tot_cons_read_n, tot_cons_n, tot_cons_base, tot_map_cons_n, tot_map_cons_base, error_rate = 0, 0, 0, 0, 0, 0, 0, '0.0%'
	read_names = dict()
	with open(long_read_len, 'r') as in_fp:
		for line in in_fp:
			[name, read_len] = line.rsplit()
			read_len = int(read_len)
			tot_read_n += 1
			tot_base += read_len
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
	tot_isoform, tot_known_bsj, tot_known_bsj_read_n, tot_cano_bsj, tot_cano_bsj_read_n = 0, 0, 0, 0, 0
	with open(isoform_out) as in_fp:
		for line in in_fp:
			if line.startswith('#'): continue
			tot_isoform += 1
			ele = line.rsplit()
			if ele[isoform_output_header_idx['isKnownBSJ']] == 'True':
				tot_known_bsj += 1
				tot_known_bsj_read_n += int(ele[isoform_output_header_idx['readCount']])
			elif ele[isoform_output_header_idx['isCanoBSJ']] == 'True':
				tot_cano_bsj += 1
				tot_cano_bsj_read_n += int(ele[isoform_output_header_idx['readCount']])

	with open(stats_out, 'w') as out:
		out.write('Total_read_number\t{}\n'.format(tot_read_n))
		out.write('Total_base\t{}\n'.format(tot_base))
		out.write('Total_cons_read_number\t{}\n'.format(tot_cons_read_n))
		out.write('Total_cons_seq_number\t{}\n'.format(tot_cons_n))
		out.write('Total_cons_seq_base\t{}\n'.format(tot_cons_base))
		out.write('Total_mappable_cons_number\t{}\n'.format(tot_map_cons_n))
		out.write('Total_mappable_cons_base\t{}\n'.format(tot_map_cons_base))
		out.write('Error_rate\t{}\n'.format(error_rate))
		out.write('Total_isoform\t{}\n'.format(tot_isoform))
		out.write('Total_isoform_with_known_BSJ\t{}\n'.format(tot_known_bsj))
		out.write('Total_read_with_known_BSJ\t{}\n'.format(tot_known_bsj_read_n))
		out.write('Total_isoform_with_cano_BSJ\t{}\n'.format(tot_cano_bsj))
		out.write('Total_read_with_cano_BSJ\t{}\n'.format(tot_cano_bsj_read_n))
	ut.err_format_time('basic_stats_core', 'Writing basic stats to file done!')	

def basic_stats(args):
	if args.type == 'ont':
		stats_core(args.long_read_len, args.cons_info, args.cons_bam, args.out, args.stats_out)
	elif args.type == 'pb':
		pb_stats_core(args.long_read_len, args.cons_info, args.cons_bam, args.out, args.stats_out)

def parser_argv():
	# parse command line arguments
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Stats result of {} result".format(__program__))
	parser.add_argument('long_read_len', metavar='long_read.len', type=str, help='Read length file of long read.')
	parser.add_argument('cons_info', metavar='cons.info', type=str, help='Information file of extracted consensus sequence.')
	parser.add_argument('cons_bam', metavar='cons.bam', type=str, help='Alignment file of consensus sequence.')
	parser.add_argument('out', metavar='{}.out'.format(__program__), type=str, help='Isoform-wise circRNA output file of {} result.'.format(__program__))
	parser.add_argument('stats_out', metavar='{}_stats.out'.format(__program__), type=str, help='Stats output file of {} result.'.format(__program__))
	parser.add_argument('--type', type=str, help='Type of sequencing data: Oxford Nanopore(ont) or Pacific Biosciences (pb).', choices=['ont', 'pb'], default='ont')

	return parser.parse_args()

if __name__ == '__main__':
	args = parser_argv()
	basic_stats(args)
