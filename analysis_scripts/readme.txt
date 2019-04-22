1) all_samples_isocirc0320.pl is a script used to summarize outputs of isoCirc pipeline for the R script (overlap_datasets_example.R) in next step. Its usage is:
perl all_samples_isocirc0320.pl /path/to/isoCirc_output_folders/
For example: 'perl all_samples_isocirc0320.pl /mnt/isilon/xing_lab/aspera/isocirc_datasets/raw_long_reads/raw_isoCirc_out/ > all_samples_isocirc4overlap_nanopore_raw_0417.txt' will generate a summary of isoCirc outputs for raw data.

2) R script overlap_datasets_example.R
The input file in this script is 'all_samples_isocirc4overlap_nanopore_raw_0417.txt' generated in the above step. Please alter this input file and labels/output files to generate new Venn graph.

3) miRNA_binding_site_isoforms_isocirc.pl is a script to output miRNA seeds in all detected circRNA isoforms.
The three input files are /path/to/seperate_chromosomes/, /path/to/isoCirc_output, and /path/to/miRNA_list.
The output file first recorded length of chromosomes and number of circRNAs in each.
Next it list each circRNA ID, sequence and miRNA seed found in its sequence. Only miRNA seeds found more than 3 times in a circRNA are listed.
An example output of this script is 'miRNA_binding_site_isoforms_isocirc_pb3.txt'. I also told Yan about the details of the output format.
