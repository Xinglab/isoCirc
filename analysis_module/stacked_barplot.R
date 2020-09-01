#!/usr/bin/env Rscript

# isoCirc (1.0.0)
# Author: Robert Wang, PhD Student (Xing Lab)
# Date: 2020.08.06
# E-mail: robwang@pennmedicine.upenn.edu

# This is a script that generates stacked barplots to display relative proportions
# of circular RNA isoforms of a gene across tissues/cell lines. Only isoforms where
# the back-splice junction and forward-splice junction are not classified as 'NNC'
# will be plotted. For definitions of the different categories describing back-splice
# junctions and forward-splice junctions, refer to our GitHub page at
# https://github.com/Xinglab/isoCirc

# Usage:
#	Rscript stacked_barplot.R -i /path/to/input/file \
#		-g [gene name] \
#		-o /path/to/output/PDF

# Dependencies:
#	* argparse
#	* ggplot2
#	* dplyr
#	* tidyr
#	* reshape2
#	* data.table
#	* RColorBrewer
#	* stringr

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(stringr))

format_columns <- function(geneDF) {
	# Establish a data frame with total read counts (sum over all isoforms)
	# for each tissue/cell line
	countDF <- as.data.frame(setDT(as.data.frame(geneDF[,-1] %>% colSums(.)), 
		keep.rownames = TRUE)[])
	colnames(countDF) <- c('group', 'readCount')
	
	# Add total read counts to the column names of geneDF
	colnames(geneDF) <- c('isoformID', paste(paste(paste(countDF$group, ' (', sep=''), 
		as.character(countDF$readCount), sep=''), ')', sep=''))
	
	# Return geneDF
	return(geneDF)
}

get_matrix <- function(inDF, gene) {
	# Create a list of lists for genes assigned to each circRNA isoform
	geneLL <- strsplit(as.character(inDF$geneName), ',')
	
	# Retrieve row indices for isoforms assigned to the given gene
	rowIdx <- which(unlist(lapply(geneLL, function(x) gene %in% x)))
	
	# Extract rows from the input data frame that match rowIdx and 
	# select for isoformID and sample columns
	geneDF <- inDF[rowIdx,-c(2:4)]
	
	# Remove isoform IDs containing gene names that are not the current gene
	geneDF <- as.data.frame(geneDF %>% 
		mutate(isoformID = strsplit(as.character(isoformID), ',')) %>%
		unnest(cols=c(isoformID)) %>%
		subset(., grepl(paste(gene, 'circRNA', sep='.'), isoformID))
	)
	
	# Format column names to include sum of read counts over all isoforms
	geneDF <- format_columns(geneDF)
	
	# Return the geneDF object
	return(geneDF)
}

read_df <- function(infile) {
	# Read in input file as a data frame
	inDF <- read.table(infile, sep='\t', header=TRUE)
	
	# Select the columns: isoformID, geneName, BSJCate, FSJCate, tissue/cell-line sample columns
	inDF <- inDF[,c(1,7,23,24,34:ncol(inDF))]
	
	# Filter out isoforms with 'NA' gene assignment
	inDF <- inDF[!(is.na(inDF$geneName)),]
	
	# Filter out isoforms where either BSJCate or FSJCate is 'NNC'
	inDF <- inDF[!(inDF$BSJCate == 'NNC' | inDF$FSJCate == 'NNC'),]
	
	# Collapse read counts by tissue/cell-line sample columns
	samples <- colnames(inDF)[-c(1:4)]
	groups <- unlist(lapply(samples, function(x) unlist(strsplit(x, '_'))[1]))
	groupSet <- unique(groups)
	
	# Establish a data frame mapping group names to sample names
	groupDF <- data.frame(groups=groups, samples=samples)
	
	# Prepare an output data frame
	outDF <- inDF[,c(1:4)]
	
	# Iterate through groups and appending to the out data frame
	for(grp in groupSet) {
		# Retrieve sample names corresponding to each group
		sampleNames <- groupDF[groupDF$groups == grp,]$samples
		
		# Get matrix of read counts for each group
		grpRC <- inDF[,sampleNames]
		
		# Collapse matrix by rowSums function if there is more than one sample
		if(length(sampleNames) > 1){
			grpRC <- rowSums(grpRC)
		}
		
		# Append grpRC to outDF as a new column specified by group name
		outDF <- cbind(outDF, data.frame(tmp=grpRC))
		colnames(outDF)[ncol(outDF)] <- grp
	}
	
	# Return resulting data frame
	return(outDF)
}

get_plot_df <- function(geneMatrix) {
	# Compute the relative abundance of isoforms in each group
	for(grp in colnames(geneMatrix)[-1]) {
		total <- sum(geneMatrix[[grp]])
		if(total == 0) {
			geneMatrix[[grp]] <- 0
		}
		else {
			geneMatrix[[grp]] <- geneMatrix[[grp]]/total
		}
	}
	
	# Melt the matrix into a data frame for plotting
	# Suppress the warning message here (reshape2 is deprecated)
	plotDF <- suppressWarnings(melt(geneMatrix, 
		id.vars=c('isoformID'), 
		variable.name='group', 
		value.name='proportion'))
	
	# Re-order factors in plotDF
	plotDF$isoformID <- factor(plotDF$isoformID, levels=str_sort(unique(plotDF$isoformID), numeric=TRUE))
	
	# Return the plotDF object
	return(plotDF)
}

theme_plot <- function() {
	# Establish plot theme
	theme_gray(base_size=12, base_family = '') %+replace% theme(
		panel.background = element_rect(fill=NULL),
		panel.grid.minor.y = element_line(size=0),
		panel.grid.minor.x = element_line(size=0),
		panel.grid.major = element_line(colour = 'grey80',size=0.2)
	)
}

make_plot <- function(plotDF, outfile, gene, numIsoforms) {
	# Create a color palette for stacked barplot
	getPalette <- colorRampPalette(brewer.pal(9, 'Set1'))

	# Set plot theme
	theme_set(theme_plot())
	
	# Establish axis limits:
	ymin <- 0
	ymax <- 1
	
	# Establish axis/legend titles:
	ylabel <- 'circRNA isoform proportion'
	legend_title <- paste(gene, 'circRNA isoforms', sep=" ")
	
	# Establish line widths and font features
	line_size <- 1
	font_size <- 15
	font_family <- 'sans'
	font_color <- 'black'
	
	# Establish color palette
	set.seed(123)
	colorChoices <- sample(getPalette(numIsoforms))
	
	# Create stacked barplot
	p <- ggplot(plotDF, aes(x=group, y=proportion, fill=isoformID)) +
		geom_bar(stat='identity', color='black') +
		ylab(ylabel) +
		scale_y_continuous(limit=c(ymin,ymax), breaks=seq(ymin,ymax, by=0.25)) +
		theme(
			text=element_text(size=font_size, family=font_family),
			legend.title = element_text(face = 'bold'),
			legend.position = 'right',
			axis.title.x=element_blank(),
			axis.text.x=element_text(angle=45, hjust=1, face='bold', color='black')) +
		scale_fill_manual(values = colorChoices) +
		guides(fill=guide_legend(legend_title, ncol=1)) 
	
	# Save ggplot object to outfile
	ggsave(outfile, plot=p, dpi=300, width=9.5, height=6.5)
}

main <- function() {
	moduleSummary <- 'Creates a stacked barplot to display relative proportions of circRNA isoforms of a gene across tissues/cell lines.'
	parser <- ArgumentParser(description = moduleSummary)
	
	# Add options to the argument parser
	parser$add_argument('-i', metavar='/path/to/input/file', required=TRUE, help='path to input file')
	parser$add_argument('-g', metavar='<string>', required=TRUE, help='gene name')
	parser$add_argument('-o', metavar='/path/to/output/PDF', required=TRUE, help='path to output PDF')
	
	# Parse command line arguments
	pargs <- parser$parse_args()
	infile <- pargs$i
	gene <- pargs$g
	outfile <- pargs$o
	
	# Read input file as data frame, filtering out 'NA' isoforms and isoforms where
	# the BSJCate or FSJCate is classified as 'NNC'
	inDF <- read_df(infile)
	
	# Isolate gene matrix associated with input gene
	geneMatrix <- get_matrix(inDF, gene)
	
	# Count number of isoforms associated with gene
	numIsoforms <- nrow(geneMatrix)
	
	# Get data frame for creating stacked barplot
	plotDF <- get_plot_df(geneMatrix)
	
	# Generate stacked barplot and save to outfile
	make_plot(plotDF, outfile, gene, numIsoforms)
}

main()