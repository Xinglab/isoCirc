#!/usr/bin/env Rscript

# isoCirc (1.0.0)
# Author: Robert Wang, PhD Student (Xing Lab)
# Date: 2020.08.06
# E-mail: robwang@pennmedicine.upenn.edu

# This is a script that creates a cumulative distribution plot for the number of exons
# per circular RNA isoform. Isoforms are grouped based on categorizations of their 
# back-splice junctions and forward-splice junctions respectively as follows: 
# (i) FSM/NIC - FSM/NIC, (ii) FSM/NIC - NNC, (iii) NNC - FSM/NIC, (iv) NNC - NNC, 
# (v) combined. For definitions of the different categories describing back-splice
# junctions and forward-splice junctions, refer to our GitHub page at
# https://github.com/Xinglab/isoCirc

# Usage:
#	Rscript exon_number_CDF.R -i /path/to/input/file \
#		-o /path/to/output/PDF

# Dependencies:
#	* argparse
#	* ggplot2
#	* dplyr

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))

format_label <- function(groupString, counts){
	# Retrieve number of isoforms associated with groupString
	n <- counts[counts$group == groupString,]$count
		
	# Format number for output
	outString <- format(n, big.mark = ',', scientific = FALSE)
	outString <- paste(outString, ')', sep = '')
	outString <- paste('(', outString, sep = '')
	outString <- paste(groupString, outString, sep = ' ')
		
	# Return the reformatted groupString
	return(outString)
}

get_counts <- function(inDF) {
	# Summarize counts for all isoform groups
	counts <- data.frame(table(inDF$group))
	colnames(counts) <- c('group', 'count')
	
	# Re-format group labels to include isoform counts
	inDF$group <- unlist(lapply(inDF$group, format_label, counts=counts))
	
	# Cast re-formatted group labels as factors
	inDF$group <- as.factor(inDF$group)
	inDF$group <- factor(inDF$group, levels = c(levels(inDF$group)[-1], levels(inDF$group)[1]))
	
	# Return resulting data frame
	return(inDF)
}

filter_df <- function(inDF, groupString) {
	# Filter out an input data frame based on an input group string
	if(groupString == 'FSM/NIC-FSM/NIC') {
		outDF <- inDF[inDF$BSJCate != 'NNC' & inDF$FSJCate != 'NNC',]
	}
	else if(groupString == 'FSM/NIC-NNC') {
		outDF <- inDF[inDF$BSJCate != 'NNC' & inDF$FSJCate == 'NNC',]
	}
	else if(groupString == 'NNC-FSM/NIC') {
		outDF <- inDF[inDF$BSJCate == 'NNC' & inDF$FSJCate != 'NNC',]
	}
	else if(groupString == 'NNC-NNC') {
		outDF <- inDF[inDF$BSJCate == 'NNC' & inDF$FSJCate == 'NNC',]
	}
	else {
		# All-All circRNA isoforms
		outDF <- inDF
	}
	
	# Append a column describing BSJ-FSJ category (specified by groupString)
	outDF$group <- rep(groupString, nrow(outDF))
	
	# Return the data frame
	return(outDF)
}

create_df <- function(infile) {
	# Read in input file as data frame and select BSJCate, FSJCate, and blockCount columns
	inDF <- read.table(infile, sep='\t', header=TRUE)
	inDF <- inDF %>% select(BSJCate, FSJCate, blockCount)
	
	# Create separate data frame for plotting CDF
	plotDF <- data.frame(BSJCate=c(), FSJCate=c(), blockCount=c(), group=c())
	
	# Establish groupings
	groups <- c('All-All', 'FSM/NIC-FSM/NIC', 'FSM/NIC-NNC', 'NNC-FSM/NIC', 'NNC-NNC')
	
	# Augment plotDF with sub data frames specific to each groupString
	for(groupString in groups) {
		plotDF <- rbind(plotDF, filter_df(inDF, groupString))
	}
	
	# Annotate each group label with formatted counts
	plotDF <- get_counts(plotDF)
	
	# Return resulting data frame
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

make_plot <- function(plotDF, outfile) {
	# Set plot theme
	theme_set(theme_plot())
	
	# Establish axis limits:
	xmin <- 0
	xmax <- (floor(max(plotDF$blockCount)/5)+1)*5
	ymin <- 0
	ymax <- 1
	
	# Establish axis/plot/legend titles:
	xlabel <- 'Exon number'
	ylabel <- 'Cumulative fraction'
	legend_title <- 'BSJ-FSJ category'
	
	# Establish line widths and font features
	line_size <- 1
	font_size <- 15
	font_family <- 'sans'
	font_color <- 'black'
	
	# Create plot
	p <- ggplot(plotDF, aes(blockCount, color=group)) + 
		stat_ecdf(size=line_size) +
		xlab(xlabel) +
		ylab(ylabel) +
		scale_x_continuous(limit=c(xmin, xmax), breaks=seq(xmin, xmax, by=5)) +
		scale_y_continuous(limit=c(ymin,ymax), breaks=seq(ymin,ymax, by=0.25)) +
		theme(
			text=element_text(size=font_size, family=font_family),
			legend.position= 'right', 
			axis.text=element_text(color=font_color)
		) +
		scale_color_manual(
			values=c('#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#000000'),
			name=legend_title
		)
	
	# Save plot to output file
	ggsave(outfile, plot=p, dpi=300, width=12.2, height=6)
}

main <- function() {
	moduleSummary <- 'Creates CDF plot for exon number per circRNA isoform.'
	parser <- ArgumentParser(description = moduleSummary)
	
	# Add options to the argument parser
	parser$add_argument('-i', metavar='/path/to/input/file', required=TRUE, help='path to input file')
	parser$add_argument('-o', metavar='/path/to/output/PDF', required=TRUE, help='path to output PDF')
	
	# Parse command line arguments
	pargs <- parser$parse_args()
	infile <- pargs$i
	outfile <- pargs$o
	
	# Create data frame for CDF plots
	plotDF <- create_df(infile)
	
	# Generate CDF plot and save to outfile
	make_plot(plotDF, outfile)
}

main()