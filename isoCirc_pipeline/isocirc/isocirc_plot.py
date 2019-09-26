import argparse
from plot import percAreaPlot as pap
from plot import disToKnownPlot as dkp
from plot import stackedBarPlot as sbp
from plot import histogramPlot as hp
from collections import defaultdict as dd
from collections import OrderedDict as od
import sys
import math
import utils as ut
import copy
import isocirc
from __init__ import __program__

eval_header_idx = isocirc.isoform_output_header_idx


def increase_dict(d, start, end, inc_cnt):
    for i in range(start, end + 1, 1):
        d[i] += inc_cnt

# start of genomic distribution plot
# based on blockType(E:exon, I:intron, N:intergenic)
def get_genomic_len(ele, order, all_dict, max_cnt_thd=10, cov_rat=0.75):
    block_type = ele[eval_header_idx['blockType']].split(',')
    is_exon, is_intron, is_intergenic = False, False, True
    if 'E' in block_type:
        is_exon = True
    elif 'I' in block_type:
        is_intron = True
    non_primary = True if len(ele[eval_header_idx['chrom']]) > 5 else False
    alu_len = 0 if ele[eval_header_idx['Alu']] == 'NA' else int(ele[eval_header_idx['Alu']])
    rRNA_len = 0 if ele[eval_header_idx['rRNA']] == 'NA' else int(ele[eval_header_idx['rRNA']])
    all_TE_len = 0 if ele[eval_header_idx['allRepeat']] == 'NA' else int(ele[eval_header_idx['allRepeat']])
    # other_TE_len = 0 if ele[eval_header_idx['allRepeat']] == 'NA' else int(ele[eval_header_idx['allRepeat']]) - alu_len - rRNA_len
    read_cnt = int(ele[eval_header_idx['readCount']])
    cov_cnt = min(read_cnt, max_cnt_thd)
    inc_cnt = cov_cnt

    for cate in order:
        increase_dict(all_dict[cate], 1, cov_cnt, 0)
    for cate in order:
        if cate == 'Exon':  #
            if is_exon:
                increase_dict(all_dict[cate], 1, cov_cnt, inc_cnt)
                break
        elif cate == 'Intron':
            if is_intron:
                increase_dict(all_dict[cate], 1, cov_cnt, inc_cnt)
                break
        elif cate == 'Intergenic':
            if is_intergenic:
                increase_dict(all_dict[cate], 1, cov_cnt, inc_cnt)
                break
        elif cate == 'Alu':
            if alu_len > 0: # / map_len >= cov_rat:
                increase_dict(all_dict[cate], 1, cov_cnt, inc_cnt)
                break
        elif cate == 'rRNA':
            if rRNA_len > 0: # / map_len >= cov_rat:
                increase_dict(all_dict[cate], 1, cov_cnt, inc_cnt)
                break
        elif cate == 'OtherRepeat':
            if all_TE_len > 0: # / map_len >= cov_rat:
                increase_dict(all_dict[cate], 1, cov_cnt, inc_cnt)
                break
        elif cate == 'Non-primaryAssembly':
            if non_primary:
                increase_dict(all_dict[cate], 1, cov_cnt, inc_cnt)
                break
        else:
            increase_dict(all_dict[cate], 1, cov_cnt, inc_cnt)
            break


def genomic_distribution_plot(in_isoform_eval_out, out_plot='genomic_distribution.png'):
    hier_order = ['Exon', 'Alu', 'rRNA', 'OtherRepeat', 'Intron', 'Intergenic', 'Non-primaryAssembly', 'Others']
    plot_name_order = ['Exon', 'Intron', 'Intergenic', 'Non-primaryAssembly', 'Others', 'OtherRepeat', 'Alu', 'rRNA']
    xticks = [str(i) for i in range(1, 11)]
    xlabel = 'Read count >='
    ylabel = 'Percentage (%)'

    all_dict = dd(lambda: dd(lambda: 0))

    with open(in_isoform_eval_out) as in_fp:
        for line in in_fp:
            if line.startswith('#'): continue
            ele = line.rsplit()
            get_genomic_len(ele, hier_order, all_dict, 10, 0.75)
    for name in plot_name_order:
        all_dict[name] = list(all_dict[name].values())
    for name in all_dict:
        if not all_dict[name]:
            ut.err_format_time('genomic_distribution_plot', 'Not enough annotation provided for plotting genomic_distribution_plot for {}.'.format(in_isoform_eval_out))
            return
    ut.err_format_time('genomic_distribution_plot', 'Plotting genomic_distribution_plot for {} ... '.format(in_isoform_eval_out))
    pap.per_area_plot(out_fig=out_plot, in_dict=all_dict, group_name_list=plot_name_order, subgroup=[4,1,3], title='Genomic distribution',
                      xticks=xticks, xlabel=xlabel, ylabel=ylabel)
    ut.err_format_time('genomic_distribution_plot', 'Plotting genomic_distribution_plot for {} done!'.format(in_isoform_eval_out))
# end of genomic distribution plot

# start of dis to known site plot
def get_dis_to_known(ele, int_dict, back_dict):
    dis = ele[eval_header_idx['disToKnownSS']].split(',')
    read_cnt = int(ele[eval_header_idx['readCount']])
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



def dis_to_known_site_plot(in_isoform_eval_out, max_dis=30, out_plot='dis_to_known.png'):
    ut.err_format_time('dis_to_known_site_plot', 'Plotting dis_known_site_plot for {} ... '.format(in_isoform_eval_out))
    internal_dis_dict = dd(lambda: 0)
    back_dis_dict = dd(lambda: 0)
    dis_list = [i for i in range(-max_dis, max_dis+1)]
    xlabel = 'Distance'
    ylabel = 'log10(Read count + 1)'

    with open(in_isoform_eval_out) as in_fp:
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
def exon_block_number_plot(in_isoform_eval_out, cnt_max=20, out_plot='exon_block_num.png'):
    ut.err_format_time('exon_block_number_plot', 'Plotting exon_block_number_plot for {} ... '.format(in_isoform_eval_out))
    inc_cnt = 1 # for isoform
    l = []
    for i in range(cnt_max):
        l.append(0)
    cnt_dict = dd(lambda : copy.copy(l))

    cnt_dict['Known circRNA'][0] = 0
    cnt_dict['Known Back-splice-junction'][0] = 0
    cnt_dict['Novel Back-splice-junction'][0] = 0

    with open(in_isoform_eval_out) as in_fp:
        for line in in_fp:
            if line.startswith('#'): continue
            ele = line.rsplit()
            # block count
            cnt = int(ele[eval_header_idx['blockCount']])
            cnt_idx = cnt_max-1 if cnt > cnt_max else cnt-1
            # back-splice-site
            # is_known_circRNA:0, is_known_bsj:1, novel_bsj:2
            if ele[eval_header_idx['novelFlag']] == '0':  # TODO known circRNA
                cnt_dict['Known circRNA'][cnt_idx] += inc_cnt
            elif ele[eval_header_idx['isKnownBSJ']] == 'True':
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
                if ele[eval_header_idx['novelFlag']] == '0':  # TODO known circRNA
                    circ_len['Known circRNA'].append(int(ele[eval_header_idx['circLen']]))
                # elif ele[eval_header_idx['isKnownBSJ']] == 'True':
                #     circ_len['Known Back-splice-junction'].append(int(ele[eval_header_idx['consMapLen']]))
                # else:
                #     circ_len['Novel Back-splice-junction'].append(int(ele[eval_header_idx['consMapLen']]))

    bins = [i*100 for i in range(0,21)]
    xticks = bins
    title = 'circRNA length'
    xlabel = 'circRNA length'
    ylabel = 'circRNA isoform number'
    hp.histogramPlot_core(out_fig=out_plot, bins=bins, in_dict=circ_len,  title=title, xticks=xticks, xlabel=xlabel, ylabel=ylabel)
    ut.err_format_time('circ_len_plot', 'Plotting circ_len_plot for {} done!'.format(in_isoform_eval_out))


def isoform_num_per_gene_plot(in_isoform_eval_out, iso_max=10, out_plot='iso_per_gene.png'):
    ut.err_format_time('isoform_num_per_gene_plot', 'Plotting isoform_num_per_gene_plot for {} ... '.format(in_isoform_eval_out))
    l = []
    for i in range(iso_max):
        l.append(0)
    iso_num_dict = dd(lambda: copy.copy(l))
    gene_iso_dict = dd(lambda: 0)
    with open(in_isoform_eval_out) as in_fp:
        for line in in_fp:
            if line.startswith('#'): continue
            else:
                ele = line.rsplit()
                if ele[eval_header_idx['novelFlag']] == '0':  # TODO known circRNA
                    for g in ele[eval_header_idx['geneID']].split(','):
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


def stats_plot_core(in_isoform_eval_out, max_dis, block_max, iso_max, plot_prefix):
    dis_to_known_site_plot(in_isoform_eval_out, max_dis, plot_prefix + '.dis_to_known.png')
    circ_len_plot(in_isoform_eval_out, plot_prefix + '.circRNA_len.png')
    genomic_distribution_plot(in_isoform_eval_out, plot_prefix + '.genom_distribu.png')
    exon_block_number_plot(in_isoform_eval_out, block_max, plot_prefix + '.exon_block_num.png')
    isoform_num_per_gene_plot(in_isoform_eval_out, iso_max, plot_prefix + '.isoform_num_per_gene.png')

def plot_core(args):
    whole_out = args.whole_out
    iso = args.iso_out
    pre = args.plot_prefix
    plot_type = args.plot_type
    if plot_type == 'knownBSJ':
        known_BSJ_plot(iso, pre + '.known_BSJ.png')
    elif plot_type == 'regionArea':
        genomic_distribution_plot(iso, pre + '.genom_distribu.png')
    else:
        ut.fatal_format_time(__func__, 'Unknown plot type: {}.'.format(plot_type))


def parser_argv():
    # parse command line arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="Generate plots for {} reslut.".format(__program__))
    parser.add_argument('whole_out', metavar='whole.out', type=str, help='All reads\' output file of PARRIS result.')
    parser.add_argument('isoform_out', metavar='isoform.out', type=str, help='All isofroms\' output file of PARRIS result.')
    parser.add_argument('plot_prefix', metavar='prefix', type=str, help='Prefix of generated plots.')
    parser.add_argument('plot_type', metavar='type', type=str, help='Plot type', choices=['knownBSJ', 'regionArea'])
    return parser.parse_args()

if __name__ == '__main__':
    # stats_plot_core(sys.argv[1], 30, 20, 10, sys.argv[1])
    # circ_len_plot(sys.argv[1], sys.argv[1] + '.circRNA_len.png')
    args = parser_argv()
    plot_core(args)
