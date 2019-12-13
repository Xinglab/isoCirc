1) all_samples_isocirc0320.pl 
This script summarizes outputs of the isoCirc pipeline for the R script (overlap_datasets_example.R) in the next step. 
Its usage is: perl all_samples_isocirc0320.pl /path/to/isoCirc_output_folders/

For example: 'perl all_samples_isocirc0320.pl /mnt/isilon/xing_lab/aspera/isocirc_datasets/raw_long_reads/raw_isoCirc_out/ >
all_samples_isocirc4overlap_nanopore_raw_0417.txt' 

will generate a summary of isoCirc outputs for raw data.

2) R script overlap_datasets_example.R
The input file in this script is 'all_samples_isocirc4overlap_nanopore_raw_0417.txt' generated in the above step. 
Please alter this input file and labels/output files to generate a new Venn graph.

3) miRNA_binding_site_isoforms_isocirc.pl 
This script outputs miRNA seeds in all detected circRNA isoforms.
The three input files are /path/to/separate_chromosomes/, /path/to/isoCirc_output, and /path/to/miRNA_list.
The output file first records the chromosome length and number of circRNAs in each.
Next, the file lists each circRNA ID, sequence, and miRNA seed found in the sequence. 
Only miRNA seeds found >3 times in a circRNA are listed.
An example output of this script is 'miRNA_binding_site_isoforms_isocirc_pb3.txt'. 

4) isoform_analysis.tar.gz is a zipped TAR file containing various scripts for analyzing circRNA isoforms detected by isoCirc. 
Functions of these scripts include:
* assignment of IDs to circRNA isoforms found in the isoCirc output file 
* visualization of various statistics on circRNA isoforms (e.g. exon number per isoform, transcript length per isoform, 
read count per isoform, etc.)
* identification of tissue-specific and tissue-stable circRNA isoforms using isoCirc output for sequencing of 12 human tissues
* visualization of proportions of circRNA isoforms across 12 human tissues for a user-specific gene of interest
Details regarding each script can be found in README files contained within the TAR
Decompress the file by running: tar -zxvf isoform_analysis.tar.gz
