import argparse
from collections import defaultdict as dd
from collections import OrderedDict as od
import sys
import math
import copy

from isocirc.__init__ import isoform_output_header_idx
from isocirc.__init__ import isoform_output_header
from isocirc.plot import percAreaPlot as pap
from isocirc.plot import disToKnownPlot as dkp
from isocirc.plot import stackedBarPlot as sbp
from isocirc.plot import histogramPlot as hp
import isocirc.utils as ut

# genomic region category
hier_order = ['rRNA', 'Alu', 'OtherRepeat', 'Exon', 'Exon&Intron', 'Intron', 'Intergenic']
plot_name_order = ['Exon', 'Exon&Intron', 'Intron', 'Intergenic', 'OtherRepeat', 'Alu', 'rRNA']

subread_t = "Subread"
polyread_t = "PolymeraseRead"

# remove '_cons'
def get_subreads_cnt(ids):
    ele = ids.rsplit(',')
    read_ids = []
    for cons_id in ele:
        read_id = cons_id.rsplit('_cons')[0]
        if read_id not in read_ids:
            read_ids.append(read_id)
    return len(read_ids)

# remove '/start_end'
def get_polyreads_cnt(ids):
    ele = ids.rsplit(',')
    read_ids = []
    for cons_id in ele:
        if len(name.rsplit('/')) > 1:
            read_id = cons_id.rsplit('/')[0]+ '/' + cons_id.rsplit('/')[1]
            if read_id not in read_ids:
                read_ids.append(read_id)
    return len(read_ids)

def get_genomic_cate1(ele, order, cov_rat=0.75):
    block_type = ele[isoform_output_header_idx['blockType']].split(',')
    is_exon, is_exon_intron, is_intron, is_intergenic = False, False, False, True
    if 'E' in block_type and 'I' in block_type:
        is_exon_intron = True
        is_intergenic = False
    elif 'E' in block_type:
        is_exon = True
        is_intergenic = False
    elif 'I' in block_type:
        is_intron = True
        is_intergenic = False
    # non_primary = True if len(ele[isoform_output_header_idx['chrom']]) > 5 else False
    alu_len = 0 if ele[isoform_output_header_idx['Alu']] == 'NA' else int(ele[isoform_output_header_idx['Alu']])
    rRNA_len = 0 if ele[isoform_output_header_idx['rRNA']] == 'NA' else int(ele[isoform_output_header_idx['rRNA']])
    all_TE_len = 0 if ele[isoform_output_header_idx['allRepeat']] == 'NA' else int(ele[isoform_output_header_idx['allRepeat']])
    thres_len = cov_rat * int(ele[isoform_output_header_idx['refMapLen']]) # 0
    for cate in order:
        if (cate == 'Exon' and is_exon) or (cate == 'Exon&Intron' and is_exon_intron) or \
            (cate == 'Intron' and is_intron) or (cate == 'Intergenic' and is_intergenic) or \
            (cate == 'Alu' and alu_len > thres_len) or  (cate == 'rRNA' and rRNA_len > thres_len) or (cate == 'OtherRepeat' and all_TE_len > thres_len):
            return cate
    ut.fatal_format_time('get_genomic_cate1', 'No category is assigned.')

def pb_get_all_stats(samp, isoform_res, stats_out):
    with open(isoform_res, 'r') as in_fp, open(stats_out, 'w') as out_fp:
        samp = isoform_res if not samp else samp
        out_fp.write('Sample\tReadType\tCount\tCategory\tConsLen\tMapLen\tKnownSS\n')
        for line in in_fp:
            if line.startswith('#'): continue
            ele = line.rsplit()
            subreads_cnt = get_subreads_cnt(ele[isoform_output_header_idx['readIDs']])
            polyreads_cnt = get_polyreads_cnt(ele[isoform_output_header_idx['readIDs']])
            cons_len = int(ele[isoform_output_header_idx['consLen']])
            map_len = int(ele[isoform_output_header_idx['consMapLen']])
            # circ_len = int(ele[isoform_output_header_idx['circLen']])
            cate = get_genomic_cate1(ele, hier_order, 0.75)
            known_ss = "False" if 'False' in ele[isoform_output_header_idx['isKnownSS']] else "True"
            out_fp.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(samp, subread_t, subreads_cnt, cate, cons_len, map_len, known_ss))
            out_fp.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(samp, polyread_t, polyreads_cnt, cate, cons_len, map_len, known_ss))

def nano_get_all_stats(samp, isoform_res, stats_out):
    with open(isoform_res, 'r') as in_fp, open(stats_out, 'w') as out_fp:
        samp = isoform_res if not samp else samp
        out_fp.write('Sample\tCount\tCategory\tMapLen\tKnownSS\n')
        for line in in_fp:
            if line.startswith('#'): continue
            ele = line.rsplit()
            read_cnt = int(ele[isoform_output_header_idx['readCount']])
            map_len = int(ele[isoform_output_header_idx['refMapLen']])
            cate = get_genomic_cate1(ele, hier_order, 0.75)
            known_ss = 'False' if 'False' in ele[isoform_output_header_idx['isKnownSS']] else 'True'
            out_fp.write('{}\t{}\t{}\t{}\t{}\n'.format(samp, read_cnt, cate, map_len, known_ss))

def classify_read_by_mapped_region(samp, isoform_res, out_fn):
    with open(isoform_res, 'r') as in_fp, open(out_fn, 'w') as out_fp:
        out_fp.write('Sample\tReadType\tCount\tCategory\n')
        for line in in_fp:
            if line.startswith('#'): continue
            ele = line.rsplit()
            subreads_cnt = get_subreads_cnt(ele[isoform_output_header_idx['readIDs']])
            polyreads_cnt = get_polyreads_cnt(ele[isoform_output_header_idx['readIDs']])
            cate = get_genomic_cate1(ele, hier_order, 0.75)
            out_fp.write('{}\t{}\t{}\t{}}\n'.format(samp, subread_t, subreads_cnt, cate))
            out_fp.write('{}\t{}\t{}\t{}}\n'.format(samp, polyread_t, polyreads_cnt, cate))

# start of dis to known site plot
def get_dis_to_known(ele, int_dict, back_dict):
    dis = ele[isoform_output_header_idx['disToKnownSS']].split(',')
    read_cnt = int(ele[isoform_output_header_idx['readCount']])
    # internal site
    for d in dis[1:-1]:
        if d != 'NA':
            int_dict[int(d)] += read_cnt
        # else: int_dict[max_dis+1] += read_cnt
    bd =[dis[0], dis[-1]]
    for d in bd:
        if d != 'NA':
            back_dict[int(d)] += read_cnt
        # else: back_dict[max_dis+1] += read_cnt

def dis_to_known_site_plot(in_isoform_res, max_dis=30, out_plot='dis_to_known.png'):
    ut.err_format_time('dis_to_known_site_plot', 'Plotting dis_known_site_plot for {} ... '.format(in_isoform_eval_out))
    internal_dis_dict = dd(lambda: 0)
    back_dis_dict = dd(lambda: 0)
    dis_list = [i for i in range(-max_dis, max_dis+1)]
    xlabel = 'Distance'
    ylabel = 'log10(Read count + 1)'

    with open(in_isoform_res) as in_fp:
        for line in in_fp:
            if line.startswith('#'): continue
            ele = line.rsplit()
            get_dis_to_known(ele, internal_dis_dict, back_dis_dict)
    for i in dis_list:
        internal_dis_dict[i] += 0
        back_dis_dict[i] += 0

    inter_od_dict = od([(i, math.log10(internal_dis_dict[i]+1)) for i in dis_list])
    back_od_dict = od([(i, math.log10(back_dis_dict[i]+1)) for i in dis_list])
    dis_dict = {'Internal splice-site': list(inter_od_dict.values()), 'Back-splice-site': list(back_od_dict.values())}
    # print dis_dict
    dkp.dis_to_known_plot(out_fig=out_plot, in_dict=dis_dict, subgroup=[1, 1], title='Distance to known splice-site',
                          xticks=dis_list, xlabel=xlabel, ylabel=ylabel)
    ut.err_format_time('dis_to_known_site_plot', 'Plotting dis_known_site_plot for {} done!'.format(in_isoform_eval_out))
    return
# end of dis to known site plot


# for each circRNA isoform
# return dict{
# known: [1_num, 2_num, ...]
# known_bsj: [1_num, 2_num, ...]
# novel_bsj: [1_num, 2_num, ...]
# }
def exon_block_number_plot(in_isoform_res, cnt_max=20, out_plot='exon_block_num.png'):
    ut.err_format_time('exon_block_number_plot', 'Plotting exon_block_number_plot for {} ... '.format(in_isoform_eval_out))
    inc_cnt = 1 # for isoform
    l = []
    for i in range(cnt_max):
        l.append(0)
    cnt_dict = dd(lambda : copy.copy(l))

    cnt_dict['Known circRNA'][0] = 0
    cnt_dict['Known Back-splice-junction'][0] = 0
    cnt_dict['Novel Back-splice-junction'][0] = 0

    with open(in_isoform_res) as in_fp:
        for line in in_fp:
            if line.startswith('#'): continue
            ele = line.rsplit()
            # block count
            cnt = int(ele[isoform_output_header_idx['blockCount']])
            cnt_idx = cnt_max-1 if cnt > cnt_max else cnt-1
            # back-splice-site
            # is_known_circRNA:0, is_known_bsj:1, novel_bsj:2
            if ele[isoform_output_header_idx['novelFlag']] == '0':  # TODO known circRNA
                cnt_dict['Known circRNA'][cnt_idx] += inc_cnt
            elif ele[isoform_output_header_idx['isKnownBSJ']] == 'True':
                cnt_dict['Known Back-splice-junction'][cnt_idx] += 0#inc_cnt
            else:
                cnt_dict['Novel Back-splice-junction'][cnt_idx] += 0#inc_cnt

    subgroup=[1,1,1]
    title='circRNA block number'
    xticks = list(map(str, range(1, cnt_max)))
    xticks.append('>= {}'.format(cnt_max))
    xlabel = 'Block number'
    ylabel = 'circRNA isoform number'
    sbp.stackedBarPlot_core(out_fig=out_plot, in_dict=cnt_dict, subgroup=subgroup, title=title, xticks=xticks, xlabel=xlabel, ylabel=ylabel)

    ut.err_format_time('exon_block_number_plot', 'Plotting exon_block_number_plot for {} done!'.format(in_isoform_eval_out))


def circ_len_plot(in_isoform_eval_out, out_plot):
    ut.err_format_time('circ_len_plot', 'Plotting circ_len_plot for {} ... '.format(in_isoform_eval_out))
    circ_len = dd(lambda :[])
    
    circ_len['Known circRNA'] = []
    # circ_len['Known Back-splice-junction'] = []
    # circ_len['Novel Back-splice-junction'] = []

    with open(in_isoform_eval_out) as in_fp:
        for line in in_fp:
            if line.startswith('#'): continue
            else:
                ele = line.rsplit()
                if ele[isoform_output_header_idx['novelFlag']] == '0':  # TODO known circRNA
                    circ_len['Known circRNA'].append(int(ele[isoform_output_header_idx['circLen']]))
                # elif ele[isoform_output_header_idx['isKnownBSJ']] == 'True':
                #     circ_len['Known Back-splice-junction'].append(int(ele[isoform_output_header_idx['consMapLen']]))
                # else:
                #     circ_len['Novel Back-splice-junction'].append(int(ele[isoform_output_header_idx['consMapLen']]))

    bins = [i*100 for i in range(0,21)]
    xticks = bins
    title = 'circRNA length'
    xlabel = 'circRNA length'
    ylabel = 'circRNA isoform number'
    hp.histogramPlot_core(out_fig=out_plot, bins=bins, in_dict=circ_len,  title=title, xticks=xticks, xlabel=xlabel, ylabel=ylabel)
    ut.err_format_time('circ_len_plot', 'Plotting circ_len_plot for {} done!'.format(in_isoform_eval_out))

def isoform_num_per_gene_plot(in_isoform_res, iso_max=10, out_plot='iso_per_gene.png'):
    ut.err_format_time('isoform_num_per_gene_plot', 'Plotting isoform_num_per_gene_plot for {} ... '.format(in_isoform_eval_out))
    l = []
    for i in range(iso_max):
        l.append(0)
    iso_num_dict = dd(lambda: copy.copy(l))
    gene_iso_dict = dd(lambda: 0)
    with open(in_isoform_res) as in_fp:
        for line in in_fp:
            if line.startswith('#'): continue
            else:
                ele = line.rsplit()
                if ele[isoform_output_header_idx['novelFlag']] == '0':  # TODO known circRNA
                    for g in ele[isoform_output_header_idx['geneID']].split(','):
                        gene_iso_dict[g] += 1

    iso_num_dict['Known circRNA'][0] = 0
    for g, cnt in gene_iso_dict.items():
        cnt_idx = iso_max-1 if cnt > iso_max else cnt-1
        iso_num_dict['Known circRNA'][cnt_idx] += 1
    subgroup = [1]
    title = 'circRNA isoform number per gene'
    xticks = list(map(str, range(1, iso_max)))
    xticks.append('>= {}'.format(iso_max))
    xlabel = 'circRNA isoform number'
    ylabel = 'Gene number'
    # print iso_num_dict
    sbp.stackedBarPlot_core(out_fig=out_plot, in_dict=iso_num_dict, subgroup=subgroup, title=title, xticks=xticks,
                            xlabel=xlabel, ylabel=ylabel)

    ut.err_format_time('isoform_num_per_gene_plot', 'Plotting isoform_num_per_gene_plot for {} done!'.format(in_isoform_eval_out))



def inc_dict(d, start, end, inc):
    for i in range(start, end+1): d[i] += inc
# ovlp_isoform:
# 1: ovlp(A,B)
# 2: A
# blank
# 3: B
def get_shared_read_with_known_BSJ(samp, ovlp_isoform_res, out_fn):
    subreads_with_known_bsj = dd(lambda: (dd(lambda : 0), dd(lambda : 0))) # samp: {num of min reads : (shared, diff)}
    polyreads_with_known_bsj = dd(lambda: (dd(lambda : 0), dd(lambda : 0))) # num of min reads : (shared, diff)
    samp = ovlp_isoform_res if not samp else samp
    header_n = len(isoform_output_header)
    A_shared_ids, B_shared_ids = [], []
    with open(ovlp_isoform_res, 'r') as in_fp, open(out_fn, 'w') as out_fp:
        # out_fp.write('Sample\tReadType\tMinCount\tShared\n')
        is_A = True
        for line in in_fp:
            if line.startswith('#'): continue
            ele = line.rsplit()
            if len(ele) <= 1: 
                is_A = False
                continue
            if len(ele) == header_n * 2:
                if ele[isoform_output_header_idx['readIDs']] not in A_shared_ids:
                    A_shared_ids.append(ele[isoform_output_header_idx['readIDs']])
                    subreads_cnt1 = get_subreads_cnt(ele[isoform_output_header_idx['readIDs']]) 
                    polyreads_cnt1 = get_polyreads_cnt(ele[isoform_output_header_idx['readIDs']]) 
                    inc_dict(subreads_with_known_bsj['SampA'][0], 1, subreads_cnt1, 1)
                    inc_dict(polyreads_with_known_bsj['SampA'][0], 1, polyreads_cnt1, 1)
                if ele[isoform_output_header_idx['readIDs'] + header_n] not in B_shared_ids:
                    B_shared_ids.append(ele[header_n + isoform_output_header_idx['readIDs']])
                    subreads_cnt2 = get_subreads_cnt(ele[header_n+isoform_output_header_idx['readIDs']]) 
                    polyreads_cnt2 = get_polyreads_cnt(ele[header_n+isoform_output_header_idx['readIDs']]) 
                    inc_dict(subreads_with_known_bsj['SampB'][0], 1, subreads_cnt2, 1)
                    inc_dict(polyreads_with_known_bsj['SampB'][0], 1, polyreads_cnt2, 1)
                # out_fp.write('{}\t{}\t{}\t{}\n'.format(samp, subread_t, min(subreads_cnt1,subreads_cnt2), "Shared"))
                # out_fp.write('{}\t{}\t{}\t{}\n'.format(samp, polyread_t, min(polyreads_cnt1, polyreads_cnt2), "Shared"))
            elif is_A: # A
                if ele[isoform_output_header_idx['readIDs']] not in A_shared_ids:
                    subreads_cnt1 = get_subreads_cnt(ele[isoform_output_header_idx['readIDs']]) 
                    polyreads_cnt1 = get_polyreads_cnt(ele[isoform_output_header_idx['readIDs']]) 
                    inc_dict(subreads_with_known_bsj['SampA'][1], 1, subreads_cnt1, 1)
                    inc_dict(polyreads_with_known_bsj['SampA'][1], 1, polyreads_cnt1, 1)
                    # out_fp.write('{}\t{}\t{}\t{}\n'.format(samp, subread_t, subreads_cnt1, "Different"))
                    # out_fp.write('{}\t{}\t{}\t{}\n'.format(samp, polyread_t, polyreads_cnt1, "Different"))
            else:
                if ele[isoform_output_header_idx['readIDs']] not in B_shared_ids:
                    subreads_cnt1 = get_subreads_cnt(ele[isoform_output_header_idx['readIDs']]) 
                    polyreads_cnt1 = get_polyreads_cnt(ele[isoform_output_header_idx['readIDs']]) 
                    inc_dict(subreads_with_known_bsj['SampB'][1], 1, subreads_cnt1, 1)
                    inc_dict(polyreads_with_known_bsj['SampB'][1], 1, polyreads_cnt1, 1)

                    # out_fp.write('{}\t{}\t{}\t{}\n'.format(samp, subread_t, subreads_cnt1, "Different"))
                    # out_fp.write('{}\t{}\t{}\t{}\n'.format(samp, polyread_t, polyreads_cnt1, "Different"))

        all_dict = [subreads_with_known_bsj, polyreads_with_known_bsj]
        all_type = [subread_t, polyread_t]
        max_cnt = 10
        out_fp.write('Sample\tReadType\tCountCutOff\tShared\tRatio\tTotal\n')
        for i in range(1, max_cnt+1):
            for s in ('SampA', 'SampB'):
                for ds, t in zip(all_dict, all_type):
                    d = ds[s]
                    out_fp.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(s, t, i, d[0][i], d[0][i]/(d[0][i]+d[1][i]+0.0), d[0][i]+d[1][i]))


def stats_core(args):
    if args.type == 'all':
        if args.read_type == 'nano':
            nano_get_all_stats(args.samp_name, args.isoform_res, args.stats_out)
        else:
            pb_get_all_stats(args.samp_name, args.isoform_res, args.stats_out)
    elif args.type == 'regionArea':
        classify_read_by_mapped_region(args.samp_name, args.isoform_res, args.stats_out)
    elif args.type == 'sharedKnownBSJ':
        get_shared_read_with_known_BSJ(args.samp_name, args.isoform_res, args.stats_out)
    else:
        ut.fatal_format_time('stats_core', 'Unknown stats type: {}.'.format(type))

def parser_argv():
    # parse command line arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="Generate stats result for {} reslut.".format(__program__))
    parser.add_argument('isoform_res', metavar='isoform.res', type=str, help='Isoform-wise reslut file of {}.'.format(__program__))
    parser.add_argument('stats_out', metavar='stats.out', type=str, help='Output file.')
    parser.add_argument('--read-type', '-r', type=str, default='nano', help='Read type.', choices=['nano', 'pb'])
    parser.add_argument('--samp-name', '-n', type=str, default='', help='Sample name.')
    parser.add_argument('--type', '-t', type=str, default='all', help='Stats type.', choices=['all', 'sharedKnownBSJ', 'regionArea'])
    return parser.parse_args()

if __name__ == '__main__':
    args = parser_argv()
    stats_core(args)
