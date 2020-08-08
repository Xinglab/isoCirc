This is a module for conducting downstream analysis of circular RNA isoforms detected
using isoCirc.

Please contact [Robert Wang](mailto:robwang@pennmedicine.upenn.edu) for feedback, questions,
or bugs. 

## TABLE OF CONTENTS
* [Dependencies](#dependencies)
* [Module Overview](#module-overview)
* [Usage](#usage)


## DEPENDENCIES

#### R
The analysis module requires R (version 4.0.2), which is available at their 
[website](https://www.r-project.org/). It also requires argparse (version 2.0.1),
ggplot2 (version 3.3.2), dplyr (version 1.0.0), gridExtra (version 2.3), 
tidyr (version 1.1.1), reshape2 (version 1.4.4), data.table (version 1.13.0),
RColorBrewer (version 1.1.2), and stringr (version 1.4.0).


#### Python
The analysis module requires Python3 (version 3.8.3). It also requires
pandas (version 1.1.0).


## MODULE OVERVIEW
### Workflow
The workflow of the analysis module:

![](https://raw.githubusercontent.com/Xinglab/isoCirc/master/analysis_module/Workflow.png?token=AFTHNXRJZWEH5AKO7VOHYCC7HBGZ2 "Workflow")

## USAGE
As shown in the [Module Overview](#module-overview), the analysis module can be used
to (i) visualize various statistics related to detected circular RNA isoforms (e.g.
exon number per isoform, transcript length per isoform, read count per isoform per
tissue/cell line, relative proportions of isoforms for a gene across all tissues/
cell lines) and (ii) detect tissue-specific and tissue-stable circular RNA isoforms
relative to all tissues/cell lines. Inputs for the R scripts that carry out these
applications are generated from a set of file pre-processing scripts (written in 
Python3). In brief, the pre-processing steps involve merging full-length circular
RNA isoforms detected across all tested tissues/cell lines and assigning IDs to
each isoform in the merged file based on read count abundances. 


#### Merging individual isoCirc output files
Raw isoCirc output files for tested tissues/cell lines will be merged together into 
a single file that contains the union of full-length circular RNA isoforms detected
across all tissues/cell lines. Specifically, the file will only include isoforms
for which there exists at least one tissue/cell line in which the sum of read counts
across all respective biological replicates is greater than or equal to a user-defined
threshold. To generate this file:
```bash
	python merge_isocirc.py -c /path/to/configuration/file \
		-r [read count thresold] \
		-o /path/to/output/file
```
Note that the configuration file is a user-supplied comma-separated file (CSV) with the
following field information: (i) tissue/cell line name, (ii) biological replicate number,
(iii) path name to corresponding isoCirc output file.

#### Assigning isoform IDs
IDs are assigned to each full-length circular RNA isoform in the merged file based on
the highest median read count abundance across all tissues/cell lines. To generate
this file:
```bash
	python assign_id.py -i /path/to/input/merged/file \
		 -o /path/to/output/file
```

#### Visualizing statistics related to detected circular RNA isoforms
Once the merged file of full-length circular RNA isoforms with IDs has been generated, a
set of R scripts have been written to allow the user to visualize various statistics 
related to their merged file of isoforms. 

To visualize the distribution of exon number per isoform, run:
```bash
	Rscript exon_number_CDF.R -i /path/to/input/merged/ID/file \
		-o /path/to/output/PDF
```

To visualize the distribution of transcript length per isoform, run:
```bash
	Rscript transcript_length_CDF.R -i /path/to/input/merged/ID/file \
		-o /path/to/output/PDF
```

To visualize the distribution of read count per isoform per tissue/cell line, 
run:
```bash
	Rscript read_count_CDF.R -i /path/to/input/merged/ID/file \
		-r [read count threshold] \
		-o /path/to/output/PDF
```
The option '-r' allows the user to only keep isoforms with read counts 
greater than or equal to a specified threshold for a given tissue/cell line.

If you are interested in a particular gene (e.g. KDM1A), you can also visualize
the relative proportions of circular RNA isoforms for the gene across all 
tissues/cell lines as a stacked barplot by running:
```bash
	Rscript stacked_barplot.R -i /path/to/input/merged/ID/file \
		-g [gene name] \
		-o /path/to/output/PDF
```

#### Detecting tissue-specific and tissue-stable circular RNA isoforms
To detect tissue-specific and tissue-stable circular RNA isoforms among the
merged set of full-length isoforms relative to all tested tissues/cell lines, you
can run the following:
```bash
	Rscript tissue_specificity.R -i /path/to/input/merged/ID/file \
		-p [number of processors] \
		-o /path/to/output/prefix
```
The script will output two files: (i) a file containing detected tissue-stable
circular RNA isoforms, (ii) a file containing detected tissue-specific circular
RNA isoforms. Methods used for detecting tissue-specific and tissue-stable circular
RNA isoforms can be found in our paper "isoCirc catalogs full-length circular RNA 
isoforms in human transcriptomes" (insert citation)
