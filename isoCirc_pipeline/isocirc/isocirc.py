#!/usr/bin/env python
threads = 8
whole_output_header = ['#readID', 'chrom', 'startCoor0base', 'endCoor', 'mapStrand',
					   'geneStrand', 'geneID', 'geneName',  # 'transID', 'transName',
					   'blockCount', 'blockSize', 'blockStarts', 'refMapLen',  # mapping information
					   'blockType', 'blockAnno',  # evaluation with whole annotation for each block
					   'readLen', 'consLen', 'consMapLen', 'copyNum', 'consFrac', 'chimInfo',
					   # original read and consensus sequence information
					   'isKnownBSJ', 'disToKnownBSJ', 'isCanoBSJ', 'disToCanoBSJ', 'canoBSJMotif', 'alignAroundCanoBSJ',
					   # back-splice junction
					   'isKnownSS', 'isKnownSJ', 'isKnownExon', #
					   'isCanoSJ', 'canoSJMotif', 'isHighSJ', # isHighSJ: high-confidence SJ based on alignment around SJ
					   'CDS', 'UTR', 'lincRNA', 'antisense',  # gene_type/biotype
					   'rRNA', 'Alu', 'allRepeat', 'upFlankAlu', 'downFlankAlu']  # repeat element
whole_output_header_idx = {h: i for i, h in enumerate(whole_output_header)}

isoform_output_header = ['#isoformID', 'chrom', 'startCoor0base', 'endCoor', # 1-4
						 'geneStrand', 'geneID', 'geneName', # 5-7
						 'blockCount', 'blockSize', 'blockStarts', 'refMapLen', # 8-11
						 'blockType', 'blockAnno', # 12-13
						 'isKnownSS', 'isKnownSJ', 'isCanoSJ', 'isHighSJ', 'isKnownExon', # 14-18
						 'isKnownBSJ', 'isCanoBSJ', 'canoBSJMotif', # 19-21
					     'isFullLength', 'BSJCate', 'interIsoCate', # FSM: full splice match, NIC: novel in catelog, NNC, novel and not in catelog
						 'CDS', 'UTR', 'lincRNA', 'antisense',  # 22-25 gene_type/biotype 
						 'rRNA', 'Alu', 'allRepeat', 'upFlankAlu', 'downFlankAlu', # 26-30 repeat element
						 'readCount', 'readIDs'] # 31-32
isoform_output_header_idx = {h: i for i, h in enumerate(isoform_output_header)}

import sys
import error_corr as ec
import tandem_repeat as tr
import cons_align as ca
import bam_classify as bc
import eval_with_anno as ea
import utils as ut
import subprocess as sp
#import isocirc_plot as pp
import basic_stats as bs
from __init__ import __version__
from __init__ import __program__

import os
from collections import defaultdict as dd
import argparse

bedtools_v = '2.25.0'
minimap2_v = '2.11'

def satisfied_version(res, tool, version):
	v = []
	if tool == 'bedtools':
		v = list(map(int ,res.rsplit()[1].rsplit('v')[1].rsplit('.')))
	elif tool == 'minimap2':
		v = list(map(int, res.rsplit('-')[0].rsplit('.')))
	min_v = list(map(int, version.rsplit('.')))
	if len(v) != len(min_v):
		return False
	for i, j in zip(v, min_v):
		if i > j:
			return True
		elif i < j:
			return False
	return True

def check_depen(bedtools='bedtools', minimap2='minimap2'):
	names = ['bedtools', 'minimap2']
	urls = ['https://bedtools.readthedocs.io/', 'https://github.com/lh3/minimap2']
	ut.err_format_time('check_dependencies', 'Checking dependencies ...')
	for t, n, v, u in zip([bedtools, minimap2], names, [bedtools_v, minimap2_v], urls):
		try:
			r = sp.check_output([t, '--version'])
			# print t, r
		except:
			ut.fatal_format_time('check_dependencies', '{} ({}, >= v{}) is not installed.'.format(n, u, v))
		if not satisfied_version(r, n, v):
			ut.fatal_format_time('check_dependencies', 'Version of {} ({}) is too low, >= v{} is needed.'.format(n, u, v))
	ut.err_format_time('check_dependencies', 'Checking dependencies done!')


def isocirc_core(args):
	if not os.path.exists(args.out_dir):
		os.makedirs(args.out_dir)
	#if not os.path.exists(args.out_dir + '/plot'):
	#	os.makedirs(args.out_dir + '/plot')
	#plot_out_dir = args.out_dir + '/plot'

	# 0. error correction with short-read
	if args.short_read:
		ut.err_format_time('Error-correction', 'Hybrid error correction using {} ...'.format(args.short_read))
		corr_long = args.out_dir + '/long_corrected.fa'
		ec.error_corr(args.in_long_read, args.short_read, corr_long, args.lordec, args.kmer, args.solid, args.threads)
		ut.err_format_time('Error-correction', 'Hybrid error correction done!')
	else:
		corr_long = args.in_long_read

	cons_fa, cons_info = args.out_dir + '/cons.fa', args.out_dir + '/cons.info'
	if args.use_tidehunter:
		# 1. generate cons with TideHunter
		ut.err_format_time('TideHunter', 'Finding tandem repeats with TideHunter ...')
		th_out = args.out_dir + '/th.out'
		tr.run_tidehunter(corr_long, cons_fa, cons_info, th_out, tidehunter=args.tidehunter, cons_min_len=args.min_len, cons_min_copy=args.min_copy, cons_min_frac=args.min_frac, threads=args.threads)
		ut.err_format_time('TideHunter', 'Finding tandem repeats with TideHunter done!')
	else:
		# 1. generate cons with TRF
		ut.err_format_time('Tandem-Repeats-Finder', 'Finding tandem repeats with TRF ...')
		trf_out = args.out_dir + '/trf.out'
		if args.threads > 1:
			tr.run_trf_parall(corr_long, cons_fa, cons_info, trf_out, args.trf, args.match, args.mismatch, args.indel, args.match_frac, args.indel_frac, args.min_score, args.max_period, args.min_len, args.min_copy, args.min_frac, args.threads)
		else:
			tr.run_trf(corr_long, cons_fa, cons_info, trf_out, args.trf, args.match, args.mismatch, args.indel, args.match_frac, args.indel_frac, args.min_score, args.max_period, args.min_len, args.min_copy, args.min_frac)
		ut.err_format_time('Tandem-Repeats-Finder', 'Finding tandem repeats with TRF done!')

	# 2. extract consensus and align it to genome
	ut.err_format_time('Mapping', 'Mapping consensus sequence to genome ...')
	cons_all_sam = cons_fa + '.sam'
	ca.align_cons(cons_fa, args.ref, cons_all_sam, args.minimap2, args.threads)
	ut.err_format_time('Mapping', 'Mapping consensus sequence to genome done!')

	# 3. classify alignment
	ut.err_format_time('Classifying', 'Classifying consensus alignment ...')
	high_bam, low_bam = args.out_dir + '/high.bam', args.out_dir + '/low.bam'
	bc.bam_classify(cons_all_sam, high_bam, low_bam, args.high_max_ratio, args.high_min_ratio, args.high_iden_ratio, args.high_repeat_ratio, args.low_repeat_ratio)
	ut.err_format_time('Classifying', 'Classifying consensus alignment done!')
	# 4. evaluate alignment with annotations
	itst_bed_dict = dd(lambda: '')
	itst_bed_dict['Alu'], itst_bed_dict['allRepeat'] = args.Alu, args.all_repeat

	if os.path.exists(corr_long + '.len'):
		long_len_fn = corr_long + '.len'
	else:
		long_len_fn = args.out_dir + '/' + os.path.basename(corr_long) + '.len'
		if not os.path.exists(long_len_fn):
			ut.exec_cmd(sys.stderr, 'fxtools', '{} lp {} > {} 2> /dev/null'.format(tr.fxtools, corr_long,long_len_fn))

	# ut.err_format_time('Filtering', 'Filtering circRNA reads ...')
	isoform_out, bed_out, stats_out = args.out_dir + '/{}.out'.format(__program__), args.out_dir + '/{}.bed'.format(__program__), args.out_dir + '/{}_stats.out'.format(__program__)
	ea.eval_with_anno(high_bam, low_bam, long_len_fn, cons_info, cons_fa,
		args.ref, args.gene_anno, args.circRNA_anno, itst_bed_dict, args.bedtools, args.flank_len, args.site_dis, args.end_dis,
		args.cano_motif, args.bsj_xid, args.key_bsj_xid, args.min_circ_dis, args.rescue_low,
		args.sj_xid, args.key_sj_xid,
		isoform_out, bed_out, stats_out)
	# ut.err_format_time('Filtering', 'Filtering circRNA reads done!')
	# 5. generate stats plot
	# ut.err_format_time('Plotting', 'Generating stats figures ... ')
	# block_max, iso_per_gene_max = 20, 10
	# pp.stats_plot_core(isoform_out, max(args.end_dis, args.site_dis), block_max, iso_per_gene_max, plot_out_dir + '/')
	# ut.err_format_time('Plotting', 'Generating stats figures done!')


def parser_argv():
	# parse command line arguments
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="{}: circular RNA profiling and analysis using long-read sequencing".format(__program__))
	parser.add_argument("in_long_read", metavar='long.fa', type=str, help='Long read data generated from long-read circRNA sequencing technique.')
	parser.add_argument("ref", metavar='ref.fa', type=str, help='Reference genome sequence file.')
	parser.add_argument('gene_anno', metavar='anno.gtf', type=str, help='Whole gene annotation file in GTF format.')
	parser.add_argument('circRNA_anno', metavar='circRNA.bed/gtf', type=str, help='circRNA annotation file in BED12 or GTF format. Use \',\' to separate multiple circRNA annotation files.')
	parser.add_argument('out_dir', metavar='out_dir', type=str, help='Output directory for final result and temporary files.')

	parser.add_argument('-v', '--version', action='version', version=__program__ + ' ' + __version__)

	general_par = parser.add_argument_group('General options')
	# general_par.add_argument('--type', type=str, help='Type of sequencing data: Oxford Nanopore(ont) or Pacific Biosciences (pb).', choices=['ont', 'pb'], default='ont')
	general_par.add_argument('-t', '--threads', type=int, default=threads, help='Number of thread to use.')
	general_par.add_argument('--bedtools', help='Path to bedtools.', default=ea.bedtools)
	general_par.add_argument('--minimap2', help='Path to minimap2.', default=ca.minimap2)

	corr_par = parser.add_argument_group('Hybrid error-correction with short-read data (LoRDEC)')
	corr_par.add_argument("--short-read", metavar='short.fa', type=str, default='', help="Short-read data for error correction. Use \',\' to connect multiple or paired-end short read data.")
	corr_par.add_argument('--lordec', type=str, help='Path to lordec-correct.', default=ec.lordec)
	corr_par.add_argument('--kmer', type=int, help='k-mer size.', default=ec.kmer)
	corr_par.add_argument('--solid', type=int, help='Solid k-mer abundance threshold.', default=ec.solid)

	tr_par = parser.add_argument_group('Detecting tandem-repeats with TRF(Tandem Repeats Finder)')
	tr_par.add_argument('-T', '--use-tidehunter', default=False, action='store_true', help='Use TideHunter as tandem repeats detection tool.')
	tr_par.add_argument('--trf', type=str, help='Path to trf program.', default=tr.trf)
	tr_par.add_argument('--tidehunter', type=str, help='Path to TideHunter.', default=tr.tidehunter)

	tr_par.add_argument('--match', type=int, help='Match score.', default=tr.match)
	tr_par.add_argument('--mismatch', type=int, help='Mismatch penalty.', default=tr.mismatch)
	tr_par.add_argument('--indel', type=int, help='Indel penalty.', default=tr.indel)
	tr_par.add_argument('--match-frac', type=int, help='Match probability.', default=tr.match_frac)
	tr_par.add_argument('--indel-frac', type=int, help='Indel probability.', default=tr.indel_frac)
	tr_par.add_argument('--min-score', type=int, help='Minimum alignment score to report.', default=tr.min_score)
	tr_par.add_argument('--max-period', type=int, help='Maximum period size to report.', default=tr.max_period)

	aln_par = parser.add_argument_group('Extracting and aligning consensus sequence to genome (minimap2)')
	aln_par.add_argument('--min-len', type=int, help='Minimum consensus length to keep.', default=tr.cons_min_len)
	aln_par.add_argument('--min-copy', type=float, help='Minimum copy number of consensus to keep.', default=tr.cons_min_copy)
	aln_par.add_argument('--min-frac', type=float, help='Minimum fraction of original long read to keep.', default=tr.cons_min_frac)

	# aln_par.add_argument('-f', '--do-classify', default=False, action='store_true', help="Classify circRNA alignment into high-quality and low-quality.")
	aln_par.add_argument('--high-max-ratio', type=float, help='Maximum mappedLen / consLen ratio for high-quality alignment.', default=bc.high_max_ratio)
	aln_par.add_argument('--high-min-ratio', type=float, help='Minimum mappedLen /consLen ratio for high-quality alignment.', default=bc.high_min_ratio)
	aln_par.add_argument('--high-iden-ratio', type=float, help='Minimum identicalBases/ consLen ratio for high-quality alignment.', default=bc.high_iden_ratio)
	aln_par.add_argument('--high-repeat-ratio', type=float, help='Maximum mappedLen / consLen ratio for high-quality self-tandem consensus.', default=bc.high_repeat_ratio)
	aln_par.add_argument('--low-repeat-ratio', type=float, help='Minimum mappedLen / consLen ratio for low-quality self-tandem alignment.', default=bc.low_repeat_ratio)

	eval_par = parser.add_argument_group('Evaluating circRNA with annotation')
	eval_par.add_argument('--Alu', type=str, default='', help='Alu repetitive element annotation in BED format. ')
	eval_par.add_argument('--flank-len', type=int, default=ea.flank_len, help='Length of upstream and downstream flanking sequence to search for Alu.')
	eval_par.add_argument('--all-repeat', type=str, default='', help='All repetitive element annotation in BED format.')

	eval_par.add_argument('-s', '--site-dis', type=int, default=ea.site_dis, help='Maximum allowed distance between circRNA internal-splice-site and annoated splice-site.')
	eval_par.add_argument('-S', '--end-dis', type=int, default=ea.end_dis, help='Maximum allowed distance between circRNA back-splice-site and annoated splice-site.')

	circ_par = parser.add_argument_group('circRNA filtering criteria')
	circ_par.add_argument('--cano-motif', type=str, default='GT/AG', help='Canonical back-splice motif (GT/AG or all three motifs: GT/AG, GC/AG, AT/AC).', choices=['GT/AG', 'all'])

	circ_par.add_argument('--bsj-xid', type=int, default=1, help='Maximum allowed mis/ins/del for 20 bp sequence alignment around the back-splice junction.')
	circ_par.add_argument('--key-bsj-xid', type=int, default=0, help='Maximum allowed mis/ins/del for 4 bp sequence alignment around the back-splice junction.')
	circ_par.add_argument('--min-circ-dis', type=int, default=150, help='Minimum distance of the start and end coordinates of circRNA.')
	circ_par.add_argument('--rescue-low', default=False, action='store_true', help='Use high mapping quality reads to rescue low mapping quality reads.')

	circ_par.add_argument('--sj-xid', type=int, default=1, help='Maximum allowed mis/ins/del for 20 bp sequence alignment around the internal splice junction.')
	circ_par.add_argument('--key-sj-xid', type=int, default=0, help='Maximum allowed mis/ins/del for 4 bp sequence alignment around the internal splice junction.')
	return parser.parse_args()


def main():
	args = parser_argv()
	check_depen(args.bedtools, args.minimap2)
	isocirc_core(args)

if __name__ == '__main__':
	main()
