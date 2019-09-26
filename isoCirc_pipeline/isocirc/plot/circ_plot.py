import sys, os, math, re
import pysam as ps
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import matplotlib.lines as mlines
from matplotlib.patches import Wedge
from matplotlib.patches import Arc
from matplotlib.patches import FancyArrow
from matplotlib.patches import FancyArrowPatch
from matplotlib.patches import Rectangle
import utils as ut
import parse_bam as pb
import parse_gff as pg


font_family = 'monospace'  # , 'serif', 'sans-serif', 'cursive', 'fantasy', 'monospace'
font_weight = 'bold'  # 'normal', 'bold', 'large', 'black'
font_color = 'black'


# return:
# {circRNA_ID:[(chr, strand, start, end), ()...]}
def get_circRNA_from_gtf(in_gtf):
    gtf_db = pg.restore_gff_db(in_gtf)

    circRNA = {}
    for exon in gtf_db.features_of_type('exon', order_by='start'):
        if exon.attributes['transcript_id'][0] not in circRNA:
            circRNA[exon.attributes['transcript_id'][0]] = [(exon.chrom, exon.strand, exon.start, exon.end)]
        else:
            circRNA[exon.attributes['transcript_id'][0]].append((exon.chrom, exon.strand, exon.start, exon.end))
    return circRNA


# return:
# {circRNA_ID:[(chr, strand, start, end), ()...]}
def get_circRNA_from_bed12(in_bed):
    # chr1    9991948 9994918 hsa_circ_0000014        1000    -       9991948 9994918 0,0,255 3       82,96,99        0,912,2871
    # chromStart: 0-base
    # exonStarts: 0-base
    header_ele = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand',
                  'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'exonStarts']
    bed_header = {header_ele[i]: i for i in range(len(header_ele))}
    circRNA = {}
    with open(in_bed, 'r') as bed:
        for line in bed:
            if line.startswith('#'): continue
            ele = line.rsplit('\t')
            chrom = ele[bed_header['chrom']]
            strand = ele[bed_header['strand']]
            start = int(ele[bed_header['chromStart']])
            exon_start = [int(i) for i in ele[bed_header['exonStarts']].split(',')]
            exon_len = [int(i) for i in ele[bed_header['blockSizes']].split(',')]
            circRNA[ele[bed_header['name']]] = []
            for s, l in zip(exon_start, exon_len):
                circRNA[ele[bed_header['name']]].append((chrom, strand, start + s + 1, start + s + l))
    return circRNA


# MISMATCH: read_pos(first), ref_pos(first), len, read_base, ref_base
# INSERTION: ins_read_pos(first), ins_ref_pos(left),  len, ins_base
# DELETION: del_read_pos(left), del_ref_pos(first), len, del_base
def get_align_detial(r, ref_fa):
    if r.cigartuples[0][0] == pb.BAM_CSOFT_CLIP or r.cigartuples[0][0] == pb.BAM_CHARD_CLIP:
        left_clip = r.cigartuples[0][1]
    else:
        left_clip = 0
    if r.cigartuples[-1][0] == pb.BAM_CSOFT_CLIP or r.cigartuples[-1][0] == pb.BAM_CHARD_CLIP:
        right_clip = r.cigartuples[-1][1]
    else:
        right_clip = 0

    if r.has_tag('MD'):
        mdstr = r.get_tag('MD')
        mis_err, ins_err, dele_err = pb.get_error_from_MD(r.cigartuples, mdstr, r.query_sequence, r.reference_start)
    else:
        ref = ps.FastaFile(ref_fa)
        align_ref_seq = ref.fetch(r.reference_name, r.reference_start, r.reference_start + r.reference_length)
        mis_err, ins_err, dele_err = pb.get_error_from_Cigar(r.cigartuples, r.query_sequence, align_ref_seq,
                                                             r.reference_start)
    return [r.is_reverse, r.infer_read_length(), r.reference_start, r.reference_length, left_clip, right_clip, mis_err,
            ins_err, dele_err]


def generate_wedge(init_degree, start, end, init_x, init_y, init_r, width, dis, color):
    patches = []
    if init_degree >= 0 and init_degree < 180:  # [0, 180)
        round = int(start) / 180
    else:  # [180, 360)
        round = (int(start) - 180) / 180

    remain_degree = end - start
    last_degree = start % 360.0

    while remain_degree > 0:
        start = last_degree
        if last_degree >= 0 and last_degree < 180:  # [0, 180)
            if last_degree + remain_degree >= 180:
                remain_degree -= (180 - last_degree)
                end = 180
                last_degree = 180
            else:
                end = last_degree + remain_degree
                remain_degree = 0
        else:  # [180, 360)
            if last_degree + remain_degree >= 360:
                remain_degree -= (360 - last_degree)
                end = 360
                last_degree = 0
            else:
                end = last_degree + remain_degree
                remain_degree = 0
        if init_degree >= 0 and init_degree < 180:  # [0, 180)
            cur_x, cur_y, cur_r = init_x + round % 2 * (width + dis) / 2, init_y, init_r + round * (width + dis) / 2
        else:
            cur_x, cur_y, cur_r = init_x - round % 2 * (width + dis) / 2, init_y, init_r + round * (width + dis) / 2

        round += 1

        patches.append(Wedge((cur_x, cur_y), cur_r, start, end, width=width, color=color))
    return patches


def generate_donut_ring(ax, strand, start_degree, length, x, y, radius, width, edge_width, face_color, edge_color):
    patches = []
    PI = math.pi
    len_frac = [0.0]
    N = len(length)
    tot_len = sum(length) + 0.0
    for l in length:
        len_frac.append(len_frac[-1] + l / tot_len)

    # 1. donut ring
    font_size = 12 + 20 * width
    for i in range(len(len_frac) - 1):
        start = start_degree + len_frac[i] * 360.0
        end = start_degree + len_frac[i + 1] * 360.0
        patches.append(Wedge((x, y), radius, start, end, width=width, facecolor=face_color,
                             edgecolor=edge_color, linewidth=edge_width))
        # 2. exon number
        text, offset = 'Exon' + (str(i + 1) if strand == '+' else str(N - i)), width / 3
        ax.text(x - offset + (radius - width / 2) * math.cos((end + (start - end) / 2) * PI / 180),
                y + (radius - width / 2) * math.sin((end + (start - end) / 2) * PI / 180), text, fontsize=font_size,
                family=font_family, weight=font_weight, color=edge_color)
    if N == 1:  # single-exon circRNA, add backsplicing line
        ax.add_line(mlines.Line2D([x + radius - width, x + radius], [y, y], linewidth=edge_width, color=edge_color))
    ax.text(x + radius - 0.95 * width, y + width / 10, 'Backsplicing', fontsize=font_size - 5, color=edge_color,
            family=font_family, weight=font_weight)
    # 3. arc and arrow
    arc_radius = radius / 3
    arc_width = 6.0
    arrow_length = 0.3 * width
    arrow_angle_rad = 30 * PI / 180

    if strand == '-':
        arc_start_angle = 30
        arc_end_angle = 360 - 15
        arc_start_angle_rad = arc_start_angle * PI / 180  # degrees to radians

        arrow = FancyArrow(x + arc_radius * math.cos(arc_start_angle_rad),
                           y + arc_radius * math.sin(arc_start_angle_rad),
                           arrow_length * math.sin(arrow_angle_rad), -arrow_length * math.cos(arrow_angle_rad),
                           head_width=arrow_length, head_length=arrow_length, length_includes_head=True,
                           color=face_color)
    else:
        arc_start_angle = 15
        arc_end_angle = 360 - arc_start_angle * 2
        arc_end_angle_rad = arc_end_angle * PI / 180  # degrees to radians

        arrow = FancyArrow(x + arc_radius * math.cos(arc_end_angle_rad), y + arc_radius * math.sin(arc_end_angle_rad),
                           arrow_length * math.sin(arrow_angle_rad), arrow_length * math.cos(arrow_angle_rad),
                           head_width=arrow_length, head_length=arrow_length, length_includes_head=True,
                           color=face_color)
    arc = Arc((x, y), arc_radius * 2, arc_radius * 2,  # ellipse width and height
              theta1=arc_start_angle, theta2=arc_end_angle, linewidth=arc_width, color=face_color)
    patches.append(arc)
    patches.append(arrow)

    return patches


def generate_alignment_wedge(start_degree, init_degree, init_pos, read_tot_degree, align_detail, exon_tot_len,
                             mis_color, ins_color, del_color, clip_color, read_x, read_y, read_r, read_w, read_dis):
    [read_is_reverse, read_tot_len, ref_align_start, ref_align_len, left_clip, right_clip, mis_err, ins_err,
     del_err] = align_detail
    patches = []
    # 3.1 left_clipping: [ref_start-left_clip, ref_start]
    left_clip_degree = 360.0 * left_clip / exon_tot_len
    patches.extend(
        generate_wedge(init_degree, init_degree, init_degree + left_clip_degree, read_x, read_y, read_r, read_w,
                       read_dis, clip_color))
    # 3.2 mismatch
    mis_degree = 360.0 * 1 / exon_tot_len
    for mis in mis_err:
        # mis = [read_pos(first), ref_pos(first), len, read_base, ref_base]
        for i in range(mis[2]):
            mis_start_degree = start_degree + init_degree + 360.0 * (mis[1] + i - init_pos) / exon_tot_len
            patches.extend(
                generate_wedge(init_degree, mis_start_degree, mis_start_degree + mis_degree, read_x, read_y, read_r,
                               read_w,
                               read_dis, mis_color[mis[3][i]]))
    # 3.3 insertion
    ins_degree = 360.0 * 1 / exon_tot_len
    for ins in ins_err:
        # ins: [ins_read_pos(first), ins_ref_pos(left),  len, ins_base]
        ins_start_degree = start_degree + init_degree + 360.0 * (ins[1] - init_pos) / exon_tot_len
        patches.extend(
            generate_wedge(init_degree, ins_start_degree, ins_start_degree + ins_degree, read_x, read_y, read_r, read_w,
                           read_dis, ins_color))
    # 3.4 deletion
    for dele in del_err:
        # dele: [del_read_pos(left), del_ref_pos(first), len, del_base]
        del_degree = 360.0 * dele[2] / exon_tot_len
        del_start_degree = start_degree + init_degree + 360.0 * (dele[1] - init_pos) / exon_tot_len
        patches.extend(
            generate_wedge(init_degree, del_start_degree, del_start_degree + del_degree, read_x, read_y, read_r, read_w,
                           read_dis, del_color))
    # 3.5 right_clipping
    right_clip_degree = 360.0 * right_clip / exon_tot_len
    # print right_clip, exon_tot_len
    patches.extend(
        generate_wedge(init_degree, init_degree + read_tot_degree - right_clip_degree, init_degree + read_tot_degree,
                       read_x, read_y, read_r, read_w, read_dis, clip_color))
    return patches


# TODO: add match/misatch/in-del/unmapped statistics
def add_legend(ax, xlim, ylim, circRNA, circRNA_name, read_tot_len, exon_face_color, exon_r, exon_w, read_color, read_w,
               mis_color, ins_color, del_color, clip_color):
    patches = []
    start_x, start_y, dis = 0.6 * xlim, 0.85 * ylim, 0.05

    # 1. color legend
    exon_rect_width = exon_w / 4
    exon_rect_len = 3 * exon_rect_width
    read_rect_width = read_w
    base_rect_len = 0.7 * read_rect_width
    font_size = 14 + 20 * exon_w

    exon_rect = Rectangle((start_x, start_y), exon_rect_len, exon_rect_width, color=exon_face_color)
    patches.append(exon_rect)
    ax.text(start_x + exon_rect_len + dis / 4, start_y, 'CircRNA exon', fontsize=font_size, color=font_color,
            family=font_family, weight=font_weight)

    read_rect = Rectangle((start_x, start_y - dis), exon_rect_len, read_rect_width, color=read_color)
    patches.append(read_rect)
    ax.text(start_x + exon_rect_len + dis / 4, start_y - dis, 'Match', fontsize=font_size, color=font_color,
            family=font_family, weight=font_weight)

    unmaped_read_rect = Rectangle((start_x, start_y - 2 * dis), exon_rect_len, read_rect_width, color=clip_color)
    patches.append(unmaped_read_rect)
    ax.text(start_x + exon_rect_len + dis / 4, start_y - 2 * dis, 'Unmapped', fontsize=font_size, color=font_color,
            family=font_family, weight=font_weight)

    for i, b in enumerate('ACGT'):
        rect = Rectangle((start_x + i * dis, start_y - 3 * dis), base_rect_len, read_rect_width, color=mis_color[b])
        patches.append(rect)
        ax.text(start_x + i * dis + base_rect_len + dis / 4, start_y - 3 * dis, b, fontsize=font_size, color=font_color,
                family=font_family, weight=font_weight)

    ins_rect = Rectangle((start_x, start_y - 4 * dis), base_rect_len, read_rect_width, color=ins_color)
    patches.append(ins_rect)
    ax.text(start_x + base_rect_len + dis / 4, start_y - 4 * dis, 'Insertion', fontsize=font_size, color=font_color,
            family=font_family, weight=font_weight)

    del_rect = Rectangle((start_x + 2 * dis, start_y - 4 * dis), base_rect_len, read_rect_width, facecolor=del_color,
                         edgecolor=read_color)
    patches.append(del_rect)
    ax.text(start_x + 2 * dis + base_rect_len + dis / 4, start_y - 4 * dis, 'Deletion', fontsize=font_size,
            color=font_color, family=font_family, weight=font_weight)

    # 2. add linear circRNA exon structure
    strand = circRNA[0][1]
    exon_len = [l[3] - l[2] + 1 for l in circRNA]
    exon_N = len(exon_len)

    exon_tot_len = sum(exon_len) + 0.0
    intron_linear_len = 2 * exon_rect_width
    exon_circumference = math.pi * 2 * (exon_r - exon_w) - (exon_N - 1) * intron_linear_len
    exon_linear_len = [0.0]
    for l in exon_len:
        exon_linear_len.append(l / exon_tot_len * exon_circumference)
    for i, l in enumerate(exon_linear_len[:-1]):
        exon_linear_rect = Rectangle(
            (start_x + sum(exon_linear_len[:i + 1]) + i * intron_linear_len, start_y - 5 * dis),
            exon_linear_len[i + 1], exon_rect_width, color=exon_face_color)
        patches.append(exon_linear_rect)
        ax.text(start_x + sum(exon_linear_len[:i + 1]) + i * intron_linear_len + exon_linear_len[i + 1] / 2,
                start_y - 5 * dis + exon_rect_width / 8, str(i + 1) if strand == '+' else str(exon_N - i),
                fontsize=font_size, color='white',
                family=font_family, weight=font_weight)
    line_width = 3
    ax.add_line(
        mlines.Line2D([start_x - intron_linear_len, start_x + sum(exon_linear_len) + exon_N * intron_linear_len],
                      [start_y - 5 * dis + exon_rect_width / 2, start_y - 5 * dis + exon_rect_width / 2],
                      linewidth=line_width, color=exon_face_color))
    arrow_dis = exon_rect_width
    style = '->' if strand == '+' else '<-'
    arrow_y = start_y - 5 * dis + exon_rect_width / 2
    for i in range(1 + exon_N):
        arrow_start_x = start_x - intron_linear_len + sum(exon_linear_len[:i+1]) + i * intron_linear_len
        for j in range(int(intron_linear_len / arrow_dis)):
            arrow = FancyArrowPatch((arrow_start_x + j * arrow_dis, arrow_y),
                                    (arrow_start_x + (j + 1) * arrow_dis, arrow_y),
                                    arrowstyle=style, linewidth=line_width - 1.0, mutation_scale=30,
                                    color=exon_face_color)
            patches.append(arrow)

    # 3. read and circRNA information
    font_size = 16 + 20 * exon_w
    ax.text(start_x, start_y - 7 * dis, u'\u2022 circRNA long read: {:,}bp'.format(read_tot_len), fontsize=font_size,
            color=font_color,
            family=font_family, weight=font_weight)
    # {circRNA_ID:[(chr, strand, start, end), ()...]}
    ax.text(start_x, start_y - 8 * dis,
            u'\u2022 {}: {:,}bp\n\u2022 {} {} : {:,} - {:,}'.format(circRNA_name, exon_tot_len, circRNA[0][1],
                                                                circRNA[0][0], circRNA[0][2], circRNA[-1][3]),
            fontsize=font_size, color=font_color,
            family=font_family, weight=font_weight)

    return patches


def circ_plot(align_detail, read_name, circRNA, circRNA_name, out_dir):
    read_name = re.sub('/', '_', read_name)
    fig_name = '{}/{}_{}.pdf'.format(out_dir, circRNA_name, read_name)
    ut.err_format_time('ggplot2', 'Ploting ' + fig_name)
    mpl.style.use('ggplot')
    patches = []
    start_degree = 0.0
    # rectangle
    fig = plt.figure(figsize=(20, 10))
    x_width, xlim, y_width, ylim = 1.0, 1.0, 2.0, 0.5
    ax = fig.add_axes([0.0, 0.0, x_width, y_width], frameon=True, aspect=1)
    default_x, default_y, default_r, default_w, default_edge_width = 0.3, 0.25, 0.1, 0.05, 1

    # square
    # fig = plt.figure(figsize=(10, 10))
    # ax = fig.add_axes([0.0, 0.0, 1.0, 1.0], frameon=True, aspect=1)
    # x, y, r, w, edge_width = 0.5, 0.5, 0.2, 0.1, 1.5

    # 0. init plot
    exon_len = [l[3] - l[2] + 1 for l in circRNA]
    exon_tot_len = sum(exon_len) + 0.0
    [read_is_reverse, read_tot_len, ref_align_start, ref_align_len, left_clip, right_clip, mis_err, ins_err,
     del_err] = align_detail
    read_cov_tot_len = ref_align_len + left_clip + right_clip + 0.0
    cov = int(read_cov_tot_len / exon_tot_len) + 1

    x, y, r, w, edge_width = default_x, default_y, default_r, default_w, default_edge_width
    if r * (1 + 0.25 + cov * 13 / 120.0) > 0.25:
        r = 0.25 / (1 + 0.25 + cov * 13 / 120.0)
        w = r / 2

    # 1. donut ring of circRNA exons
    exon_x, exon_y, exon_r, exon_w, exon_edge_width = x, y, r, w, edge_width
    exon_edge_color, exon_face_color = 'white', 'black'
    strand = circRNA[0][1]
    patches.extend(
        generate_donut_ring(ax, strand, start_degree, exon_len, exon_x, exon_y, exon_r, exon_w, exon_edge_width,
                            exon_face_color, exon_edge_color))

    # 2. wedge of read
    # [ref_align_start, ref_align_len, left_clip, right_clip, mis_err, ins_err, del_err] = align_detail

    read_x, read_y, read_r, read_w, read_dis, read_color = x, y, r + exon_w / 2, exon_w / 6, exon_w / 20, 'C1' #'red' # '#1874CD' #'dodgerblue3', #'#009ACD' #'deepskyblue3'

    read_is_reverse = align_detail[0]
    init_pos = ref_align_start - left_clip
    ref_start = (ref_align_start - left_clip) % exon_tot_len
    read_cov_tot_len = ref_align_len + left_clip + right_clip

    read_tot_degree = 360.0 * read_cov_tot_len / exon_tot_len
    init_degree = (start_degree + 360.0 * ref_start / exon_tot_len) % 360.0

    patches.extend(
        generate_wedge(init_degree, init_degree, init_degree + read_tot_degree, read_x, read_y, read_r, read_w,
                       read_dis, read_color))
    # 3. wedge of alignment errors
    ins_color, del_color, clip_color = 'black', 'white', 'C3'
    mis_color = {'A': 'C2', 'a': 'C2', 'C': 'C0', 'c': 'C0', 'G': 'C5', 'g': 'C5', 'T': 'C4', 't': 'C4'}
    patches.extend(
        generate_alignment_wedge(start_degree, init_degree, init_pos, read_tot_degree, align_detail, exon_tot_len,
                                 mis_color,
                                 ins_color, del_color, clip_color, read_x, read_y, read_r,
                                 read_w, read_dis))
    # 4. legend
    exon_r, exon_w, read_w = default_r, default_w, default_w / 6
    patches.extend(
        add_legend(ax, xlim, ylim, circRNA, circRNA_name, read_tot_len, exon_face_color, exon_r, exon_w, read_color,
                   read_w, mis_color, ins_color, del_color, clip_color))
    for i in patches:
        ax.add_patch(i)
    # ax.grid(color='r', linestyle='-', linewidth=2)
    plt.axis('off')
    # plt.show()
    # plt.savefig('{}/{}_{}.png'.format(out_dir, circRNA_name, read_name), format='png', dpi=500)
    plt.savefig(fig_name, format='pdf', dpi=500)
    return


def plot_main(in_circRNA, gtf, in_bam, circRNA_fa, out_dir):
    # 1. extract circRNA information
    all_circRNA = get_circRNA_from_gtf(in_circRNA) if gtf else get_circRNA_from_bed12(in_circRNA)

    with ps.AlignmentFile(in_bam) as bam:
        for r in bam:
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue

            # 2. extract alignment details
            align_detail = get_align_detial(r, circRNA_fa)

            # 3. plot donut and wedge for circRNA and read
            # circRNA = all_circRNA[r.reference_name[:r.reference_name.find('_copy')]]
            circRNA_name = r.reference_name[:r.reference_name.find('::')]
            circRNA = all_circRNA[circRNA_name]
            circ_plot(align_detail, r.query_name, circRNA, circRNA_name, out_dir)
    return


def align_to_circRNA_fa(ref_fa, read_fa, in_circRNA, gtf, cons_list, out_dir):
    circRNA_fa = out_dir + '/circRNA_multi_copies.fa'
    circRNA_read_fa = out_dir + '/circRNA_read.fa'
    align_sam = out_dir + '/circRNA_read.fa.mem.sam'
    tmp_circRNA, tmp_fa, tmp_dup_fa = out_dir + '/tmp.circRNA', out_dir + '/tmp.fa', out_dir + '/tmp_dup.fa'

    if os.path.exists(circRNA_fa): ut.exec_cmd(sys.stderr, 'Generate circRNA read sequences',
                                               'rm {}'.format(circRNA_fa))
    if os.path.exists(circRNA_read_fa): ut.exec_cmd(sys.stderr, 'Generate circRNA read sequences',
                                                    'rm {}'.format(circRNA_read_fa))

    # 1. generate circRNA.fa
    with open(cons_list, 'r') as cons_fp:
        for line in cons_fp:
            if line.startswith('#'): continue
            ele = line.split()
            read_name, circRNA_name, copy_n = ele[0], ele[1], float(ele[2]) + 1
            circRNA_name = re.sub('/', '_', circRNA_name)
            cmd = 'grep {} {}> {}'.format(circRNA_name, in_circRNA, tmp_circRNA)
            ut.exec_cmd(sys.stderr, 'Generate circRNA exon sequences', cmd)

            if gtf:
                cmd = 'gffread {} -g {} -w {}; fxtools df {} {} > {}; cat {} >> {}'.format(tmp_circRNA, ref_fa, tmp_fa,tmp_fa, copy_n, tmp_dup_fa, tmp_dup_fa, circRNA_fa)
                ut.exec_cmd(sys.stderr, 'Generate circRNA exon sequences', cmd)
            else:
                cmd = 'bedtools getfasta -split -name -fi {} -bed {} -fo {}; fxtools df {} {} > {};  cat {} >> {}'.format(
                    ref_fa, tmp_circRNA, tmp_fa, tmp_fa, copy_n, tmp_dup_fa, tmp_dup_fa, circRNA_fa)
                ut.exec_cmd(sys.stderr, 'Generate circRNA exon sequences', cmd)
            cmd = 'fxtools sd {} {} 1 -1 >> {}'.format(read_fa, read_name, circRNA_read_fa)
            ut.exec_cmd(sys.stderr, 'Generate circRNA read sequences', cmd)
    ut.exec_cmd(sys.stderr, 'Generate circRNA read sequences', 'rm {} {} {}'.format(tmp_circRNA, tmp_fa, tmp_dup_fa))

    # 2. align read to circRNA.fa
    cmd = 'bwa index {}'.format(circRNA_fa)
    ut.exec_cmd(sys.stderr, 'Align to circRNA sequences', cmd)
    cmd = 'bwa mem -t 8 {} {} > {}'.format(circRNA_fa, circRNA_read_fa, align_sam)
    ut.exec_cmd(sys.stderr, 'Align to circRNA sequences', cmd)


"""
#circRNA.bam: alignment of original circRNA long read and multi-copies of corresponding genome sequence
ref.fa: genome sequence in fasta format
circRNA.bed/gtf: detailed information of circRNA in bed12 or gtf format
cons.list: read name | corresponding circRNA name | copy number
output: plot of circRNA and long read sequence alignment 
"""
# TODO without cons.list: align all read to all circRNA exon sequences
if __name__ == '__main__':
    if len(sys.argv) != 6:
        print('Usage:')
        print('{} ref.fa read.fa circRNA.gtf/bed cons.lit output_dir'.format(sys.argv[0]))
        # print '{} circRNA.bam/sam circRNA.bed/gtf circRNA.fa output_dir'.format(sys.argv[0])
        sys.exit(1)

    ref_fa = sys.argv[1]
    read_fa = sys.argv[2]
    in_circRNA = sys.argv[3]
    cons_list = sys.argv[4]
    out_dir = sys.argv[5]

    gtf = True if os.path.splitext(in_circRNA)[1] == '.gtf' else False

    circRNA_fa = out_dir + '/circRNA_multi_copies.fa'
    align_sam = out_dir + '/circRNA_read.fa.mem.sam'
    align_to_circRNA_fa(ref_fa, read_fa, in_circRNA, gtf, cons_list, out_dir)
    #
    plot_main(in_circRNA, gtf, align_sam, circRNA_fa, out_dir)
