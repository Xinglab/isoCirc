import sys
import os
from collections import defaultdict as dd
import pysam as ps

header = ['#isoformID', 'chrom', 'startCoor0base', 'endCoor', # 1-4
          'geneStrand', 'geneID', 'geneName', # 5-7
          'blockCount', 'blockSize', 'blockStarts', 'refMapLen', # 8-11
          'blockType', 'blockAnno', # 12-13
          'isKnownSS', 'isKnownSJ', 'isCanoSJ', 'isHighSJ', 'isKnownExon', # 14-18
          'isKnownBSJ', 'isCanoBSJ', 'canoBSJMotif', # 19-21
          'isFullLength', 'BSJCate', 'FSJCate', # 22-24 FSM: full splice match, NIC: novel in catelog, NNC, novel and not in catelog
          'CDS', 'UTR', 'lincRNA', 'antisense',  # 25-28 gene_type/biotype
          'rRNA', 'Alu', 'allRepeat', 'upFlankAlu', 'downFlankAlu', # 29-33 repeat element
          'readCount', 'readIDs'] # 34-35
isocirc_idx = {h: i for i, h in enumerate(header)}

start_n=0
tot_n=18
cons_fns=[
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_43/cons.info',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_44/cons.info',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_45/cons.info',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_58/cons.info',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_59/cons.info',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_60/cons.info',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_adipose/all.cons.info',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_adrenal/all.cons.info',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_blood_peripheral_leukocytes/all.cons.info',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_brain/all.cons.info',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_heart/all.cons.info',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_kidney/all.cons.info',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_liver/all.cons.info',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_lung/all.cons.info',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_prostate/all.cons.info',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_skeletal_muscle/all.cons.info',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_smooth_muscle/all.cons.info',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_testis/all.cons.info'
]

high_bams=[
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_43/high.bam',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_44/high.bam',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_45/high.bam',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_58/high.bam',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_59/high.bam',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_60/high.bam',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_adipose/high.bam',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_adrenal/high.bam',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_blood_peripheral_leukocytes/high.bam',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_brain/high.bam',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_heart/high.bam',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_kidney/high.bam',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_liver/high.bam',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_lung/high.bam',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_prostate/high.bam',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_skeletal_muscle/high.bam',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_smooth_muscle/high.bam',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_testis/high.bam'
]

low_bams=[
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_43/low.bam',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_44/low.bam',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_45/low.bam',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_58/low.bam',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_59/low.bam',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_60/low.bam',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_adipose/low.bam',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_adrenal/low.bam',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_blood_peripheral_leukocytes/low.bam',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_brain/low.bam',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_heart/low.bam',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_kidney/low.bam',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_liver/low.bam',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_lung/low.bam',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_prostate/low.bam',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_skeletal_muscle/low.bam',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_smooth_muscle/low.bam',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_testis/low.bam'
]

isocirc_fns=[
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/bio_reps/nano_43/isocirc.out',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/bio_reps/nano_44/isocirc.out',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/bio_reps/nano_45/isocirc.out',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/bio_reps/nano_58/isocirc.out',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/bio_reps/nano_59/isocirc.out',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/bio_reps/nano_60/isocirc.out',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_adipose/isocirc.out',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_adrenal/isocirc.out',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_blood_peripheral_leukocytes/isocirc.out',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_brain/isocirc.out',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_heart/isocirc.out',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_kidney/isocirc.out',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_liver/isocirc.out',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_lung/isocirc.out',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_prostate/isocirc.out',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_skeletal_muscle/isocirc.out',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_smooth_muscle/isocirc.out',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/tissue/nano_testis/isocirc.out'
]

names=['Nano-43', 'Nano-44', 'Nano-45', 'Nano-58', 'Nano-59', 'Nano-60',
'Adipose', 'Adrenal', 'Blood', 'Brain', 'Heart', 'Kidney', 'Liver', 'Lung', 'Prostate', 'SkeletalMuscle', 'SmoothMuscle', 'Testis']

cons_info_header = ['name', 'readLen', 'consLen', 'copyNum', 'frac', 'label']
idx = {j:i for i,j in enumerate(cons_info_header)}

cate_group_dict = {'all': ['all', 'low', 'high', 'bsj', 'iso'], 
				   'map': ['high', 'low', 'iso', 'bsj'],
				   'bsj': ['iso', 'bsj'],
				   'iso': ['iso']}
# cate:
#	all: all read
#   low: cons with low mapping quality 
# 	high: cons with high mapping quality
#   bsj: cons with BSJ
# 	iso: cons with full-length isoform
def get_cate_dict(high_bam, low_bam, isocirc_fn):
	cate_dict = dd(lambda:'All') # read_name : isFullLength
	with ps.AlignmentFile(high_bam) as in_high, ps.AlignmentFile(low_bam) as in_low:
		for r in in_low:
			name = r.query_name.rsplit('_cons')[0]
			cate_dict[name] = 'low'
		for r in in_high:
			name = r.query_name.rsplit('_cons')[0]
			cate_dict[name] = 'high'
	with open(isocirc_fn) as in_fp:
		for line in in_fp:
			if line.startswith('#'): continue
			ele = line.rsplit()
			for name in ele[isocirc_idx['readIDs']].rsplit(','):
				cate_dict[name] = 'iso' if ele[isocirc_idx['isFullLength']] == 'True' else 'bsj'
	return cate_dict


def gen_len(cate, read_fig_fn, cons_fig_fn):
	tmp_dir = os.path.dirname(os.path.abspath(read_fig_fn))
	read_out = tmp_dir + '/{}_read.len'.format(cate)
	cons_out = tmp_dir + '/{}_cons.len'.format(cate)
	with open(read_out, 'w') as read_fp, open(cons_out, 'w') as cons_fp:
		read_fp.write('Name\tType\tLen\n')
		cons_fp.write('Name\tType\tLen\n')
		for (cons_fn, high_bam, low_bam, isocirc_fn, samp_name) in zip(cons_fns, high_bams, low_bams, isocirc_fns, names):
			if names.index(samp_name) < start_n:
				continue
			elif names.index(samp_name) >= tot_n:
				break
			cate_dict = get_cate_dict(high_bam, low_bam, isocirc_fn)
			last_name = ''
			largest_frac = 0.0
			last_readLen, last_consLen = 0, 0
			with open(cons_fn) as fp:
				for line in fp:
					ele = line.rsplit()
					cons_name, readLen, consLen, frac = ele[idx['name']], int(ele[idx['readLen']]), int(ele[idx['consLen']]), float(ele[idx['frac']])
					readName = cons_name.rsplit('_')[0]
					if readName != last_name:
						if last_name != '':
							if cate_dict[last_name] in cate_group_dict[cate]:
								read_fp.write('{}\treadLen\t{}\n'.format(samp_name, last_readLen))
								cons_fp.write('{}\tconsLen\t{}\n'.format(samp_name, last_consLen))
						largest_frac = frac
						last_name = readName
						last_readLen = readLen
						last_consLen = consLen
					else:
						if frac > largest_frac:
							largest_frac = frac
							last_readLen = readLen
							last_consLen = consLen
				if cate_dict[last_name] in cate_group_dict[cate]:
					read_fp.write('{}\treadLen\t{}\n'.format(samp_name, last_readLen))
					cons_fp.write('{}\tconsLen\t{}\n'.format(samp_name, last_consLen))
	cmd= 'Rscript /home/gaoy1/program/circ_plot/all_read_len.R {} {}'.format(read_out, read_fig_fn)
	print(cmd)
	os.system(cmd)
	cmd= 'Rscript /home/gaoy1/program/circ_plot/all_cons_len.R {} {}'.format(cons_out, cons_fig_fn)
	print(cmd)
	os.system(cmd)

if __name__ == '__main__':
	if len(sys.argv) != 4:
		print('{} cate read.fig cons.fig'.format(sys.argv[0]))
		sys.exit(1)
	cate = sys.argv[1]
	if cate.upper() == 'ALL':
		cate = 'all'
	elif cate.upper() == 'MAP':
		cate = 'map'
	elif cate.upper() == 'ISO':
		cate = 'iso'
	elif cate.upper() == 'BSJ':
		cate ='bsj'
	else:
		print('\'cate\' should be one of followings:')
		print('\tall')
		print('\tmap: mappable consensus')
		print('\thigh: high-quality mappable consensus')
		print('\tlow: low-quality mappable consensus')
		print('\tbsj: high-confidence BSJ')
		print('\tiso: full-length isoform')

	gen_len(cate, sys.argv[2], sys.argv[3])
