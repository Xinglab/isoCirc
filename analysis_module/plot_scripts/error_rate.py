import os, sys
import pysam as ps
from collections import defaultdict as dd
import mappy as mp
from pyfaidx import Fasta

ep_header = ['#READ_NAME',      'READ_LEN', 'UNMAP',   'INS',     'DEL',     'MIS',     'MATCH',   'CLIP',    'SKIP',    'ERR_RATE']
ep_idx = {j: i for i,j in enumerate(ep_header)}

info_header = ['CONS_NAME', 'READ_LEN', 'CONS_LEN', 'COPY_NUM', 'ID_RATE', 'FLAG']
info_idx = {j: i for i,j in enumerate(info_header)}

read_fas=[
'/home/gaoy1/circRNA/raw_long_reads/nano_43.fa',
'/home/gaoy1/circRNA/raw_long_reads/nano_44.fa',
'/home/gaoy1/circRNA/raw_long_reads/nano_45.fa',
'/home/gaoy1/circRNA/raw_long_reads/nano_58.fa',
'/home/gaoy1/circRNA/raw_long_reads/nano_59.fa',
'/home/gaoy1/circRNA/raw_long_reads/nano_60.fa'
]

cons_info_fns = [
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_43/cons.info',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_44/cons.info',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_45/cons.info',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_58/cons.info',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_59/cons.info',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_60/cons.info'
]

cons_ep_fns = [
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_43/cons.ep',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_44/cons.ep',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_45/cons.ep',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_58/cons.ep',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_59/cons.ep',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_60/cons.ep'

#'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_43/high.bam.ep',
#'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_44/high.bam.ep',
#'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_45/high.bam.ep',
#'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_58/high.bam.ep',
#'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_59/high.bam.ep',
#'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_60/high.bam.ep'
]

ref_fns = [
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_43/cons.ref.fa',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_44/cons.ref.fa',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_45/cons.ref.fa',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_58/cons.ref.fa',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_59/cons.ref.fa',
'/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_60/cons.ref.fa'
# '/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_43/high.ref.fa',
# '/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_44/high.ref.fa',
# '/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_45/high.ref.fa',
# '/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_58/high.ref.fa',
# '/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_59/high.ref.fa',
# '/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/nano_60/high.ref.fa'
]

samples = ['Rep1_1', 'Rep1_2', 'Rep1_3', 'Rep2_1', 'Rep2_2', 'Rep2_3']

def get_mp_error_rate(ref_seq, read_seq):
    a = mp.Aligner(seq=ref_seq)
    nms = []
    mis = []
    for h in a.map(read_seq):
        mi = 0
        for c in h.cigar:
            if c[1] == 0 or c[1] == 1:
                mi += c[0]
        # if not h.is_primary: continue
        nms.append(h.NM)
        mis.append(mi)
        # print(h.cigar_str, h.NM, h.mlen, h.blen, len(ref_seq), len(read_seq))
    # if not errors: return 0.5
    # else: return sum(errors)/len(errors)
    if not nms: return 0.5
    else: return sum(nms) / (sum(mis)+0.0)

def get_copy_bin(copy_num):
    if copy_num < 6:
        return copy_num
    elif copy_num <= 10:
        return 6
    else:
        return 11

def raw_error_rate(n, fig_fn):
    tmp_out = os.path.dirname(os.path.abspath(fig_fn)) + '/raw_cons_error_{}.out'.format(n)

    with open(tmp_out, 'w') as out_fp:
        if n == 0:
            out_fp.write('Sample\tCopyNum\tType\tError\n') # Type: raw/2/3/4/5/6-10/>10
        sample, read_fn, ref_fn, info_fn, cons_ep_fn = samples[n], read_fas[n], ref_fns[n], cons_info_fns[n], cons_ep_fns[n]
        read_fa = Fasta(read_fn)
        ref_fa = Fasta(ref_fn)
        with open(cons_ep_fn) as cons_ep_fp, open(info_fn) as info_fp:
            last_name = ''
            min_cons_error = 100.0
            for cons_name in ref_fa.keys():
                read_name = cons_name.rsplit('_')[0]
                if read_name != last_name and last_name != '':
                    out_fp.write('{}\t{}\t{}\t{}\n'.format(sample, copy_num, 'Raw', raw_error))
                    out_fp.write('{}\t{}\t{}\t{}\n'.format(sample, copy_num, get_copy_bin(copy_num), min_cons_error))
                    min_cons_error = 100.0
                last_name = read_name
                cons_error = 0, 0, 0
                for eline in cons_ep_fp:
                    if eline.startswith('#'): continue
                    ele = eline.rsplit()
                    name, error = ele[ep_idx['#READ_NAME']], float(ele[ep_idx['ERR_RATE']][:-1])/100.0
                    # print('Cons', name, error, cons_name)
                    if name == cons_name:
                        cons_error = error
                        break
                    else:
                        continue
                if cons_error < min_cons_error:
                    min_cons_error = cons_error

                    ref_seq = ref_fa[cons_name][:].seq.upper()
                    read_seq = read_fa[read_name][:].seq.upper()
                    raw_error = get_mp_error_rate(ref_seq, read_seq)
                    # print('Raw', cons_name, raw_error)
                    # if raw_error < 0:
                        # last_name = ''
                        # continue

                    for sline in info_fp:
                        ele = sline.rsplit()
                        name, num = ele[info_idx['CONS_NAME']], float(ele[info_idx['COPY_NUM']])
                        # print(name, num)
                        if name == cons_name:
                            copy_num = int(num)
                            break
                        else:
                            continue
                
            if last_name != '':
                out_fp.write('{}\t{}\t{}\t{}\n'.format(sample, copy_num, 'Raw', raw_error))
                out_fp.write('{}\t{}\t{}\t{}\n'.format(sample, copy_num, get_copy_bin(copy_num), min_cons_error))
                    

    cmd = 'Rscript /home/gaoy1/program/circ_plot/error_rate.R {} {}'.format(tmp_out, fig_fn)
    print(cmd)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('{} n out.fig')
        print('\tn= 0 ~ 5')
        sys.exit(1)
    raw_error_rate(int(sys.argv[1]), sys.argv[2])