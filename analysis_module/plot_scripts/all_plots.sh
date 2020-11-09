#!/bin/bash
root_dir=/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out # raw
out_dir=$root_dir/all_figures
# out_dir=/scr1/users/gaoy1/isocirc_figs

# HEK293 6 reps result
nano_43_out=$root_dir/bio_reps/nano_43/isocirc.out
nano_44_out=$root_dir/bio_reps/nano_44/isocirc.out
nano_45_out=$root_dir/bio_reps/nano_45/isocirc.out
nano_58_out=$root_dir/bio_reps/nano_58/isocirc.out
nano_59_out=$root_dir/bio_reps/nano_59/isocirc.out
nano_60_out=$root_dir/bio_reps/nano_60/isocirc.out
rep1_out=$root_dir/bio_reps/rep1/isocirc.out
rep2_out=$root_dir/bio_reps/rep2/isocirc.out
comb_6rep_out=$root_dir/bio_reps/all.isocirc.out
comb_12tissues_out=/mnt/isilon/xing_lab/aspera/isocirc_datasets/raw_long_reads/raw_isoCirc_out/20200805/read_count_2/12Tissues.isocirc.out.withID.tsv
comb_12tissues_out_cnt1=/mnt/isilon/xing_lab/aspera/isocirc_datasets/raw_long_reads/raw_isoCirc_out/20200805/read_count_0/12Tissues.isocirc.out.withID.tsv
all_comb_out=/mnt/isilon/xing_lab/aspera/isocirc_datasets/raw_long_reads/raw_isoCirc_out/20200805/read_count_2/combined.isocirc.out.withID.tsv

ciri_short_out=/home/gaoy1/circRNA/CIRI_results/sample_BSJ_reads_PE_all.txt

src_dir=/home/gaoy1/program/circ_plot
tissue_dir=$root_dir/tissue
non_fig=$out_dir/tmp.pdf

# 0 combined 12 tissues
# src=$src_dir/isoCircCombine.py
# out=$root_dir/tissue/12tissues_cnt2_isocirc.out
# echo "python3 $src 2 $tissue_dir/nano_*/isocirc.out $out"
# python3 $src 2 $tissue_dir/nano_*/isocirc.out $out
# 0 combined 6 reps
# out=$root_dir/bio_reps/HEK293_6reps_cnt1_isocirc.out
# echo "python3 $src 1 $root_dir/bio_reps/nano_*/isocirc.out $out"
# python3 $src 1 $root_dir/bio_reps/nano_*/isocirc.out $out

#1. 6 replicates nano_43/44/45/58/59/60
#Fig A/B/S1/S2. pair-wise comparison of BSJ/isoform
src=$src_dir/circRepPairHeatMap.py
fig1=$out_dir/fig_2a_s3_6rep_bsj_pairwise.pdf
fig2=$out_dir/fig_2b_s4_6rep_iso_pairwise.pdf
echo "python $src 3 $nano_43_out $nano_44_out $nano_45_out $nano_58_out $nano_59_out $nano_60_out $fig1 $fig2"
# python $src 3 $nano_43_out $nano_44_out $nano_45_out $nano_58_out $nano_59_out $nano_60_out $fig1 $fig2

#2. 2 bio-reps nano_43/44/45 vs nano_58/59/60
#Fig 2c. anno stack plot of 2 bio-reps
src=$src_dir/circHighBSJRep.py
r_src=$src_dir/circHighBSJRepAveCnt.R
min_cnt=1
fig1=$out_dir/fig_2c_2bioreps_bsj_anno_cnt1.pdf
fig2=$out_dir/fig_2c_2bioreps_bsj_cnt_cnt1.pdf
fig3=$out_dir/fig_2c_2bioreps_bsj_ave_cnt_cnt1.pdf
echo "python $src $min_cnt $rep1_out $rep2_out $fig1 $fig2 $non_fig"
echo "Rscript $r_src $out_dir/cnt1_bsj_read_cnt.out $fig3"

min_cnt=2
fig1=$out_dir/fig_2c_2bioreps_bsj_anno_cnt2.pdf
fig2=$out_dir/fig_2c_2bioreps_bsj_cnt_cnt2.pdf
fig3=$out_dir/fig_s11a_2bioreps_bsj_ave_cnt_cnt2.pdf
echo "python $src $min_cnt $rep1_out $rep2_out $fig1 $fig2 $non_fig"
echo "Rscript $r_src $out_dir/cnt2_bsj_read_cnt.out $fig3"
# python $src 2 $rep1_out $rep2_out $out_dir/1 $fig 3 $non_fig
min_cnt=3
fig1=$out_dir/fig_s10_2bioreps_bsj_anno_cnt3.pdf
fig2=$out_dir/fig_2c_2bioreps_bsj_cnt_cnt3.pdf
fig3=$out_dir/fig_s11b_2bioreps_bsj_ave_cnt_cnt3.pdf
echo "python $src $min_cnt $rep1_out $rep2_out $fig1 $fig2 $non_fig"
echo "Rscript $r_src $out_dir/cnt3_bsj_read_cnt.out $fig3"

#Fig 2d/S13. BSJ-internal catalog of full-length isoform
src=$src_dir/circCateTable.py
fig1=$out_dir/fig_2d_6rep_cate_table_cnt2.pdf
fig2=$out_dir/fig_s13_6rep_cate_table_cnt3.pdf
fig3=$out_dir/fig_2d_6rep_cate_table_cnt1.pdf
echo "python $src 3 $comb_6rep_out $fig1"
src=$src_dir/circCateTable_hek293.R
echo "Rscript $src $out_dir/bsj_iso_cate.out $fig1"
# python $src 3 $comb_6rep_out $fig1
src=$src_dir/circCateTableCnt3_hek293.R
echo "Rscript $src $out_dir/bsj_iso_cate.out $fig2"
src=$src_dir/circCateTableCnt1_hek293.R
echo "Rscript $src $out_dir/bsj_iso_cate.out $fig3"
# Rscript $src $out_dir/bsj_iso_cate.out $fig2

#Fig S10A/S10B. isoform count per gene and per (actual/log2 scale)
# src=$src_dir/circIsoNum.py
# fig1=$out_dir/fig_s10a_iso_num_per_gene.pdf
# fig2=$out_dir/fig_s10b_iso_num_per_bsj.pdf
# echo "python $src $comb_6rep_out 2 $fig1"
# # python $src $comb_6rep_out 2 $fig1
# src=$src_dir/circIsoNumPerBSJCDF.R
# echo "Rscript $src $out_dir/iso_num_per_gene_bsj.dat $fig2"
# Rscript $src $out_dir/iso_num_per_gene_bsj.dat $fig2

src=$src_dir/tissue_iso_number.py
fig=$out_dir/fig_3a_tissue_iso_num_bar.pdf
cnt=3 # 2
echo "python $src $cnt $comb_12tissues_out_cnt1 $fig"

# fig 3b
src=/home/yan/program/circ_plot/shortLongTissueComp.R
dat=/home/yan/Dropbox/Submission/IsoCirc/Re-submission/tissue_short_long.text
fig=/home/yan/data/raw_isocirc_figs/fig_3b_tissue_short_long.pdf
echo "/usr/bin/Rscript $src $dat $fig"

#Fig 3c/d CDF plot of exon count/length of full-length isoform for each catalog
src=$src_dir/tissue_exon_len.py
fig1=$out_dir/fig_3c_exon_num_each_cate.pdf
fig2=$out_dir/fig_3d_length_each_cate.pdf
echo "python $src $comb_12tissues_out $fig1 $fig2"
# python $src 2 $comb_6rep_out $out_dir/1 2 3 $non_fig $fig2 $fig1
# fig1=$out_dir/fig_s_FSM_block_each_cate.pdf
# fig2=$out_dir/fig_s_FSM_length_each_cate.pdf
# src1=$src_dir/circFullIsoRepBlockCntFSM.R
# src2=$src_dir/circFullIsoRepLenFSM.R
# echo "Rscript $src1 $out_dir/cnt1_iso_len.out $fig1"
# # Rscript $src1 $out_dir/cnt1_iso_len.out $fig1
# echo "Rscript $src2 $out_dir/cnt1_iso_len.out $fig2"
# Rscript $src2 $out_dir/cnt1_iso_len.out $fig2

# 12 tissue result
#4. 12 tissues
	# Fig 3A. isoform number of each tissue
# src=$src_dir/isonumBarPlot.py

# python $src $tissue_dir/nano_*/isocirc.out $fig
	# Fig 3C. heatmap table of tissues-specific gene 
	# Fig 3D. heatmap table tissues-specific isoform
	# Fig S11. ratio pari-wise comparison of tissue-specific isoforms
src=$src_dir/circSpecificIso.py
fig1=$out_dir/fig_3e_tissue_spe_gene.pdf
fig2=$out_dir/fig_3f_tissue_spe_iso.pdf
fig3=$out_dir/fig_s18_tissue_spe_pairwise.pdf
echo "python3 $src 2 2 $comb_12tissues_out $fig2 $fig3"
# python3 $src 2 2 $tissue_dir/nano_*/isocirc.out $fig2 $fig3
src=$src_dir/circSpecificGene.R
echo "Rscript $src $out_dir/spe.out $fig1"
# Rscript $src $out_dir/spe.out $fig1

# read/cons error rate
# TODO
# bedtools bamtobed -split -bed12 -i cons.fa.sam > cons.bed
# bedtools getfasta -fi hg19.fa -bed cons.bed -fo cons.ref.fa -name -split
src=$src_dir/error_rate.py
fig=$out_dir/fig_s1_read_cons_err_rate.pdf
echo "python $src $fig"
src=$src_dir/error_rate.R
echo "Rscript $src $out_dir/raw_cons_error.out $fig"

### rebuttal figures
# read/cons length
src=$src_dir/all_read_cons_len.py
fig1=$out_dir/all_read_len.pdf
fig2=$out_dir/all_cons_len.pdf
echo "python $src all $fig1 $fig2"
fig1=$out_dir/fig_s2a_bsj_read_len.pdf
fig2=$out_dir/bsj_s2b_cons_len.pdf
echo "python $src bsj $fig1 $fig2"
fig1=$out_dir/iso_read_len.pdf
fig2=$out_dir/iso_cons_len.pdf
echo "python $src iso $fig1 $fig2"

#Fig S5 pair-wise comparison of BSJ read count
src=$src_dir/circHighBSJRep.py
fig=$out_dir/fig_s5_6reps_bsj_readCnt_pairwise.pdf
echo "python $src 1 $nano_43_out $nano_44_out $nano_45_out $nano_58_out $nano_59_out $nano_60_out $out_dir/1 2 $fig"
# python $src 1 $nano_43_out $nano_44_out $nano_45_out $nano_58_out $nano_59_out $nano_60_out $out_dir/1 $non_fig 3 $fig
#Fig S6 pair-wise comparison of full-length isoform read count
src=$src_dir/circFullIsoRep.py
fig=$out_dir/fig_s6_6reps_iso_readCnt_pairwise.pdf
echo "python $src 1 $nano_43_out $nano_44_out $nano_45_out $nano_58_out $nano_59_out $nano_60_out $fig"

fig=/home/yan/data/raw_isocirc_figs/fig_s7_isocirc_illumina_ratio.pdf
src=/home/yan/program/circ_plot/circShortLongRatio.R
echo "Rscript $src /home/yan/Dropbox/Submission/IsoCirc/circRNA_read_ratio_6_isoCirc/circRNA_read_ratio_6_isoCirc_replicate_all.txt $fig"

# s8 heat map of illumina short-reads
src=$src_dir/circRepPairHeatMap.py
fig1=$out_dir/fig_s8_3short_bsj_pairwise.pdf
fig2=$non_fig
echo "python $src 3 $ciri_short_out $fig1 $fig2"
r_src=$src_dir/pairHeatMap_illumina.R
echo "Rscript $r_src $out_dir/bsj_cnt.out $fig1"

#3. combined 6 reps
#Fig S9. BSJ read count comparison of short-read and 6 reps 
src=$src_dir/circShortLong.py 
fig=$out_dir/fig_s9_short_long_bsj.pdf
echo "python $src $comb_6rep_out $ciri_short_out 2 $fig"
# python $src $comb_6rep_out $ciri_short_out 2 $fig



src=$src_dir/bsjAluPair.py
alu_bed=/home/gaoy1/data/genome/hg19_RepeatMasker_Alu.sorted.bed
gene_pred=/home/gaoy1/data/genome/Homo_sapiens.GRCh37.87.chr_patch_hapl_scaff_UCSC.gtf.gene_pred
circbase=/home/gaoy1/data/genome/hsa_hg19_circRNA.bed
mionco=/home/gaoy1/data/genome/mioncocirc_hg19_sorted_uniq.bed
fig=$out_dir/fig_s12_alu_pair.pdf
min_cnt=2
for min_cnt in 2 3; do for flank_len in 1000 2000; do echo "python $src $min_cnt $flank_len $rep1_out $rep2_out $alu_bed $gene_pred $circbase $mionco $fig"; done done
cat $out_dir/*bsj_alu_pair.out > $out_dir/bsj_alu_pair_all.out
r_src=$src_dir/bsjAluPair.R
echo "Rscript $r_src $out_dir/bsj_alu_pair_all.out $fig"

src=$src_dir/tissueCircCateTable.py
fig1=$out_dir/fig_s16_12tissue_cate_table_cnt2.pdf
fig2=$out_dir/fig_s16_12tissue_cate_table_cnt3.pdf
echo "python $src 3 $comb_12tissues_out $fig1 $fig2"
# python $src 3 $tissue_comb_out $fig1 $fig2

src=$src_dir/tissueIsoPerGeneBSJ.py
fig1=$out_dir/fig_s17a_12tissue_iso_num_per_gene.pdf
fig2=$out_dir/fig_s17b_12tissue_iso_num_per_bsj.pdf
echo "python $src $comb_12tissues_out $fig1 $fig2"
# python $src $tissue_comb_out $fig1 $fig2

# alternative splicing figs
# TODO
src=/home/gaoy1/program/circ_plot/alt_splice.py
comb_dir=/home/gaoy1/circRNA/raw_long_reads/raw_isoCirc_out/20200805/read_count_2/
all_sep_out=$comb_dir/combined.isocirc.out.withSepID.sorted.tsv 
gene_pred=/home/gaoy1/data/genome/Homo_sapiens.GRCh37.87.chr_patch_hapl_scaff_UCSC.gtf.gene_pred
alt_out=$comb_dir/alt_splice.list
bsj_out=$comb_dir/bsj_dis.list
echo "python $src $all_sep_out $gene_pred $alt_out $bsj_out"

src=$src_dir/as_cate_pie.R
fig=$out_dir/as_cate.pdf
echo "Rscript $src $fig"

# all
for i in ES A5S A3S IR; do echo $i; grep "^$i" $alt_out > $comb_dir/${i}_as.list; awk '{print $4,$5}' $comb_dir/${i}_as.list | sort | uniq | wc -l; done
# FSM/NIC
for i in ES A5S A3S IR; do echo $i; grep "^$i" $alt_out | awk '($10 == "FSM" && $11 == "FSM")' > $comb_dir/FSM_${i}_as.list; awk '{print $4,$5}' $comb_dir/FSM_${i}_as.list | sort | uniq | wc -l; done
# FSM
for i in ES A5S A3S IR; do echo $i; grep "^$i" $alt_out | awk '(($10 == "FSM" || $10 == "NIC") && ($11 == "FSM"|| $11 == "NIC"))' >  $comb_dir/FSM_NIC_${i}_as.list; awk '{print $4,$5}' $comb_dir/FSM_NIC_${i}_as.list | sort | uniq | wc -l; done

# minor_read_cnt
for i in ES A5S A3S IR; do echo $i; awk '($9>=2){print $4,$5}' $comb_dir/${i}_as.list | sort | uniq | wc -l; done
for i in ES A5S A3S IR; do echo $i; awk '($9>=5){print $4,$5}' $comb_dir/${i}_as.list | sort | uniq | wc -l; done
for i in ES A5S A3S IR; do echo $i; awk '($9>=10){print $4,$5}' $comb_dir/${i}_as.list | sort | uniq | wc -l; done

# maxEnt score of known/novel BSJs
src=$src_dir/known_novel_maxEntScore.py
ref_fa=/home/gaoy1/data/genome/hg19.fa
score_out_pre=$out_dir/known_novel_maxent_
echo "python $src $ref_fa $rep1_out $rep2_out $score_out_pre"
# combine scores
src=$src_dir/combineSiteScores.py
echo "python $src $score_out_pre*score $out_dir/known_novel.comb.score known"
# Boxplot
fig1=$out_dir/known_novel_maxent_box.pdf
fig2=$out_dir/known_novel_maxent_cdf.pdf
src=$src_dir/known_novel_maxEntBox.R
echo "Rscript $src $out_dir/known_novel.comb.score $fig1"
src=$src_dir/known_novel_maxEntCDF.R
echo "Rscript $src $out_dir/known_novel.comb.score $fig2"


# maxEnt score of RIs and non-RIs
# all sites of introns that are from circRNAs
ir_list=$comb_dir/IR.as.list
intron_bed=/home/gaoy1/data/genome/hg19_intron.bed
score_out_pre=$comb_dir/ri_
src=$src_dir/ir_nonir_anno_collect_fa.py
echo "python $src $ref_fa $ir_list $all_sep_out $intron_bed $out_pre"

ir_5_fa=$out_pre/ir.5.fa
ir_3_fa=$out_pre/ir.3.fa
nonir_5_fa=$out_pre/nonir.5.fa
nonir_3_fa=$out_pre/nonir.3.fa
anno_ir_5_fa=$out_pre/anno.ir.5.fa
anno_ir_3_fa=$out_pre/anno.ir.3.fa
anno_nonir_5_fa=$out_pre/anno.nonir.5.fa
anno_nonir_3_fa=$out_pre/anno.nonir.3.fa

src=$src_dir/combineSiteScores.py
echo "python $src $score_out_pre*score $out_dir/ri.comb.score ir"
src=$src_dir/RI_maxEntBox.R
fig1=$out_dir/ri_box.pdf
fig2=$out_dir/ri_cdf.pdf
echo "Rscript $src $out_dir/ri.comb.score $fig1"
src=$src_dir/RI_maxEntCDF.R
echo "Rscript $src $out_dir/ri.comb.score $fig2"
fig1=$out_dir/anno_ri_box.pdf
fig2=$out_dir/anno_ri_cdf.pdf
src=$src_dir/RI_anno_maxEntBox.R
echo "Rscript $src $out_dir/ri.comb.score fig1"
src=$src_dir/RI_anno_maxEntCDF.R
echo "Rscript $src $out_dir/ri.comb.score fig2"
