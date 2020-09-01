#!/usr/bin/env Rscript

# isoCirc (1.0.0)
# Author: Robert Wang, PhD Student (Xing Lab)
# Date: 2020.08.06
# E-mail: robwang@pennmedicine.upenn.edu

# This is a script that creates a cumulative distribution plot for the read counts
# per circular RNA isoform per tissue/cell-line. Isoforms are grouped based on 
# categorizations of their back-splice junctions and forward-splice junctions 
# respectively as follows: (i) FSM/NIC - FSM/NIC, (ii) FSM/NIC - NNC, (iii) NNC - FSM/NIC, 
# (iv) NNC - NNC, (v) combined. For definitions of the different categories describing 
# back-splice junctions and forward-splice junctions, refer to our GitHub page at
# https://github.com/Xinglab/isoCirc

# Usage:
#	Rscript read_count_CDF.R -i /path/to/input/file \
#		-r [read count threshold] \
#		-o /path/to/output/PDF

# Dependencies:
#	* argparse
#	* ggplot2
#	* gridExtra

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))

# Set penalty for using scientific notation to false
options(scipen = 0)

get_maxRC <- function(infile) {
	# Read in input file as a data frame
	inDF <- read.table(infile, sep='\t', header=TRUE)
	
	# Compute sum of read counts across all samples
	sumRC <- rowSums(inDF[,-(1:33)])
	
	# Return maximum value
	return(max(sumRC))
}

get_groups <- function(infile) {
	# Read in header of the input file as data frame
	header <- read.table(infile, sep='\t', header=FALSE, nrows=1)
	colnames(header) <- NULL
	header <- unlist(as.list(header))[-(1:33)]
	
	# Extract groups from the header
	groups <- unlist(lapply(header, function(x) unlist(strsplit(x, '_'))[1]))
	
	# Return data frame summarizing group to column name mappings
	return(data.frame(groups=groups, samples=header))
}

collapse_df <- function(inDF) {
	# Construct a subset of the input data frame with just the
	# categories for backsplice junction and forward splice junction
	tmpDF <- inDF[,c('BSJCate', 'FSJCate')]
	
	# Create a vector to store read count sums across remaining columns
	grpRC <- inDF[, !(names(inDF) %in% c('BSJCate', 'FSJCate'))]
	
	# Check if grpRC contains more than 1 column
	if(ncol(inDF) > 3){
		grpRC <- rowSums(grpRC)
	}
	
	# Return collapsed data frame
	return(data.frame(BSJCate=tmpDF$BSJCate, FSJCate=tmpDF$FSJCate, grpRC=grpRC))
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

theme_plot <- function() {
	# Establish plot theme
	theme_gray(base_size=12, base_family = '') %+replace% theme(
		panel.background = element_rect(fill=NULL),
		panel.grid.minor.y = element_line(size=0),
		panel.grid.minor.x = element_line(size=0),
		panel.grid.major = element_line(colour = 'grey80',size=0.2)
	)
}

create_plot <- function(infile, sampleNames, maxRC, threshold) {
	# Determine group name from sampleNames
	grpName <- unlist(strsplit(sampleNames[1], '_'))[1]
	
	# Read in input file as a data frame
	inDF <- read.table(infile, sep='\t', header=TRUE)
	
	# Selection of BSJCate, FSJCate, and columns specified by sampleNames
	inDF <- inDF[,c('BSJCate', 'FSJCate', sampleNames)]
	
	# Collapse the read counts across the columns dictated by sampleNames
	inDF <- collapse_df(inDF)
	
	# Filter isoforms by a user-defined read count threshold
	inDF <- inDF[inDF$grpRC >= threshold,]
	
	# Create separate data frame for plotting CDF
	plotDF <- data.frame(BSJCate=c(), FSJCate=c(), grpRC=c(), group=c())
	
	# Establish groupings
	groups <- c('All-All', 'FSM/NIC-FSM/NIC', 'FSM/NIC-NNC', 'NNC-FSM/NIC', 'NNC-NNC')
	
	# Augment plotDF with sub data frames specific to each groupString
	for(groupString in groups) {
		plotDF <- rbind(plotDF, filter_df(inDF, groupString))
	}
	
	# Annotate each group label with formatted counts
	plotDF <- get_counts(plotDF)
	
	# Set plot theme
	theme_set(theme_plot())
	
	# Establish axis limits:
	xmin <- 10^floor(log(threshold, 10))
	xmax <- 10^ceiling(log(maxRC, 10)) 
	ymin <- 0
	ymax <- 1
	
	# Establish plot/legend titles:
	plot_title <- grpName
	legend_title <- 'BSJ-FSJ category'
	
	# Establish line widths and font features
	line_size <- 1
	font_size <- 15
	font_family <- 'sans'
	font_color <- 'black'
	
	# Create plot
	p <- ggplot(plotDF, aes(grpRC, color=group)) + 
		stat_ecdf(size=line_size) +
		ggtitle(plot_title) +
		scale_x_continuous(limits=c(xmin, 0.2*xmax),
			breaks=c(1, 10^seq(log10(xmin), log10(xmax)-1), 0.2*xmax),
			trans='log10',
			labels=function(x) format(x, scientific = TRUE)) +
		scale_y_continuous(limit=c(ymin,ymax), breaks=seq(ymin,ymax, by=0.25)) +
		theme(
			text=element_text(size=font_size, family=font_family),
			plot.title = element_text(hjust = 0.5),
			legend.position= 'right', 
			axis.text=element_text(color=font_color),
			axis.text.x=element_text(angle=45, hjust=1),
			axis.title.x=element_blank(),
			axis.title.y=element_blank()
		) +
		scale_color_manual(
			values=c('#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#000000'),
			name=legend_title
		)

	# Return ggplot object
	return(p)
}

main <- function() {
	moduleSummary <- 'Creates CDF plot for read count per circRNA isoform per tissue.'
	parser <- ArgumentParser(description = moduleSummary)
	
	# Add options to the argument parser
	parser$add_argument('-i', metavar='/path/to/input/file', required=TRUE, help='path to input file')
	parser$add_argument('-r', metavar='###', default=2, help='read count threshold')
	parser$add_argument('-o', metavar='/path/to/output/PDF', required=TRUE, help='path to output PDF')
	
	# Parse command line arguments
	pargs <- parser$parse_args()
	infile <- pargs$i
	threshold <- as.numeric(pargs$r)
	outfile <- pargs$o
	
	# Create a data frame mapping group name to column names
	groupDF <- get_groups(infile)
	
	# Determine maximum read count sum across all samples to use as an upper limit
	maxRC <- get_maxRC(infile)
	
	# Retrieve list of group names
	groupNames <- unique(groupDF$groups)
	
	# Initialize a list object to store group-specific CDF plots for read counts
	cdfList <- list()
	
	# Iterate through group names
	for(i in 1:length(groupNames)) {
		grp <- groupNames[i]
		
		# Retrieve column names associated with the group as a vector
		sampleNames <- groupDF[groupDF$groups == grp,]$samples
		
		# Create ggplot object for CDF of read count sum over samples
		p <- create_plot(infile, sampleNames, maxRC, threshold)
		
		# Append ggplot object to cdfList
		cdfList[[i]] <- p
	}
	
	# Use grid.arrange to assemble multiple plots
	g <- do.call('grid.arrange', c(cdfList, ncol=2))
	
	# Save plot to outfile
	ggsave(outfile, plot=g, dpi=300, width=17, height=22)
}

main()