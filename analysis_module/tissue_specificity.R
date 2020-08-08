#!/usr/bin/env Rscript

# isoCirc (1.0.0)
# Author: Robert Wang, PhD Student (Xing Lab)
# Date: 2020.08.06
# E-mail: robwang@pennmedicine.upenn.edu

# This is a script to identify tissue-specific and tissue-stable circular
# RNA isoforms among full-length circular RNA isoforms detected across 
# various tissues/cell-lines using isoCirc.

# Usage:
#	Rscript tissue_specificity.R -i /path/to/input/file \
#		-p [number of processors] \
#		-o /path/to/output/prefix

# Dependencies:
#	* argparse
#	* stringr
#	* tidyr
#	* dplyr
#	* parallel

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(parallel))

get_BSJ_FSJ <- function(id, isoformLL, inDF) {
	# Get row index matching the input isoform ID
	rowIdx <- which(unlist(lapply(isoformLL, function(x) id %in% x)))
	
	# Retrieve categorizations for BSJ and FSJ
	BSJ <- as.character(inDF[rowIdx, 3])
	FSJ <- as.character(inDF[rowIdx, 4])
	
	return(c(BSJ, FSJ))
}

get_gene_matrix <- function(gene, geneLL, inDF) {
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
	
	# Create a matrix object from the data frame
	geneMatrix <- data.matrix(geneDF[,-1])
	rownames(geneMatrix) <- as.character(geneDF$isoformID)
	
	# Return the resulting matrix object
	return(geneMatrix)
}

normalize <- function(geneMatrix) {
	# Get group names
	groups <- colnames(geneMatrix)
	
	# Iterate through group names
	for(grp in groups) {
		# Get sum of read counts across all isoforms in a given column
		total <- sum(geneMatrix[,grp])
		if(total == 0) {
			geneMatrix[,grp] <- rep(0, nrow(geneMatrix))
		}
		else {
			geneMatrix[,grp] <- geneMatrix[,grp]/total
		}
	}
	
	# Return the resulting matrix of normalized counts
	return(geneMatrix)
}

get_stable_isoforms <- function(gene, geneLL, isoformLL, inDF) {
	# Retrieve a matrix of read counts for all isoforms associated with gene
	geneMatrix <- get_gene_matrix(gene, geneLL, inDF)
	
	# Normalize read count matrix by sum of reads in each group
	normMatrix <- normalize(geneMatrix)
	
	# Check if there exists an isoform where the read count proportion is > 50% across all tissues
	isStable <- apply(normMatrix, 1, min) > 0.5
	
	if(any(isStable)) {
		# Pull out corresponding isoform ID of stable isoform along w/ other statistics
		id <- rownames(normMatrix)[which(isStable)]
		numIso <- nrow(normMatrix)
		minGrpPct <- min(normMatrix[id,])
		minGrpCnt <- min(geneMatrix[id,])
		sjCate <- get_BSJ_FSJ(id, isoformLL, inDF)
		
		# Impose final requirement that the group read count of the candidate isoform is 
		# greater than or equal to 2 reads across all groups
		if(minGrpCnt >= 2) {
			info <- c(gene, numIso, id, minGrpPct, minGrpCnt, sjCate)
			return(paste(info, collapse=','))
		}
		else {
			return(NA)
		}
	}
	else {
		return(NA)
	}
}

filter_genes <- function(geneLL) {
	# Remove repeated gene names for each isoform
	geneLL <- lapply(geneLL, unique)
	
	# Generate a data frame summarizing the number of isoforms per gene
	isoformCounts <- as.data.frame(table(unlist(geneLL)))
	colnames(isoformCounts) <- c('gene', 'isoformCnt')
	
	# Keep genes with at least two isoforms
	geneList <- as.character(isoformCounts[isoformCounts$isoformCnt >= 2,]$gene)
	
	# Return the filtered gene list
	return(geneList)
}

filter_df <- function(infile) {
	# Read in input file as a data frame
	inDF <- read.table(infile, sep='\t', header=TRUE)
	
	# Select the columns: isoformID, geneName, BSJCate, FSJCate, tissue/cell-line sample columns
	inDF <- inDF[,c(1,7,23,24,34:ncol(inDF))]
	
	# Filter out isoforms with 'NA' gene assignment
	inDF <- inDF[!(is.na(inDF$geneName)),]
	
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

write_stable <- function(results, outfile) {
	# Filter out results that are NA
	results <- results[!is.na(results)]
	
	# Write results to the outfile if there are any
	if(length(results) > 0) {
		outDF <- data.frame(do.call('rbind', strsplit(unlist(results), ',')))
		
		# The outfile should have the following fields:
		#	1. gene name
		#	2. number of detected isoforms
		#	3. isoform ID of tissue-stable isoform
		#	4. minimum read count proportion of the tissue-stable isoform across all samples
		#	5. minimum read count of the tissue-stable isoform across all samples
		#	6. BSJ category for the tissue-stable isoform
		#	7. FSJ category for the tissue-stable isoform
		colnames(outDF) <- c('gene', 'numIsoforms', 'isoformID', 'minGroupPct', 
			'minGrpReadCounts', 'BSJCate', 'FSJCate')
		outDF[,'gene'] <- as.character(outDF[,'gene'])
		outDF[,'numIsoforms'] <- as.numeric(as.character(outDF[,'numIsoforms']))
		outDF[,'isoformID'] <- as.character(outDF[,'isoformID'])
		outDF[,'minGroupPct'] <- as.numeric(as.character(outDF[,'minGroupPct']))
		outDF[,'minGrpReadCounts'] <- as.numeric(as.character(outDF[,'minGrpReadCounts']))
		outDF[,'BSJCate'] <- as.character(outDF[,'BSJCate'])
		outDF[,'FSJCate'] <- as.character(outDF[,'FSJCate'])
		
		# Write data frame to the output file
		write.table(outDF, file=outfile, quote=FALSE, sep='\t', row.names=FALSE)
	}
}

get_specific_genes <- function(gene, geneLL, isoformLL, inDF) {
	# Retrieve a matrix of read counts for all isoforms associated with gene
	geneMatrix <- get_gene_matrix(gene, geneLL, inDF)
	
	# Retrieve column indices where the column read count sum is greater than 0
	viableIdx <- which(colSums(geneMatrix) > 0)
	
	# Filter out columns where the read count sum across all isoforms is 0
	geneMatrix <- geneMatrix[, viableIdx]
	
	# Identify genes with tissue-enriched isoforms
	if(length(viableIdx) > 1) {
		# Perform Chi-Square Test of Homogeneity
		pVal <- chisq.test(geneMatrix)$p.value
		return(paste(c(gene, pVal), collapse=','))
	}
	# Identify genes with tissue-exclusive isoforms
	else {
		return(paste(c(gene, 0), collapse=','))
	}
}

filter_specific_genes <- function(specificGenesDF) {
	# Remove genes where the resulting p-value was NA
	specificGenesDF <- specificGenesDF[!is.na(specificGenesDF)]
	specificGenesDF <- data.frame(do.call('rbind', strsplit(unlist(specificGenesDF), ',')))
	colnames(specificGenesDF) <- c('gene', 'p.value')
	specificGenesDF[,'gene'] <- as.character(specificGenesDF[,'gene'])
	specificGenesDF[,'p.value'] <- as.numeric(as.character(specificGenesDF[,'p.value']))
	
	# Selection of genes harboring potentially tissue-specific isoforms using FDR < 5%
	specificGenesDF[,'p.value'] <- p.adjust(specificGenesDF$p.value, method='hochberg')
	filtered_genes <- specificGenesDF[specificGenesDF$p.value < 0.05,]$gene

	# Return filtered list of genes
	return(filtered_genes)
}

binomial_test <- function(geneMatrix) {
	# Compute expected proportions for each isoform assuming no tissue-specificity (null)
	nullProp <- rowSums(geneMatrix)/sum(geneMatrix)
		
	# Compute column sums for group read counts
	grpRC <- colSums(geneMatrix)
		
	# Initialize an empty matrix to store binomial test p-values
	pValMatrix <- matrix(, nrow = nrow(geneMatrix), ncol = ncol(geneMatrix))
		
	# Compute binomial test p-values
	for(i in 1:nrow(geneMatrix)) {
		for(j in 1:ncol(geneMatrix)) {
			pValMatrix[i,j] <- binom.test(geneMatrix[i,j], grpRC[j], p=nullProp[i], alternative='greater')$p.value
		}
	}
	
	# Correct p-values for false-discovery rate
	pValMatrix <- matrix(p.adjust(pValMatrix, method='hochberg'), nrow=nrow(geneMatrix), ncol=ncol(geneMatrix))
	
	# Annotate p-value matrix
	rownames(pValMatrix) <- rownames(geneMatrix)
	colnames(pValMatrix) <- colnames(geneMatrix)
	
	# Return the resulting p-value matrix
	return(pValMatrix)
}

get_specific_isoforms <- function(gene, geneLL, isoformLL, inDF) {
	# Retrieve a matrix of read counts for all isoforms associated with gene
	geneMatrix <- get_gene_matrix(gene, geneLL, inDF)
	
	# Normalize read count matrix by sum of reads in each group
	normMatrix <- normalize(geneMatrix)
	
	# Retrieve column indices where the column read count sum is greater than 0
	viableIdx <- which(colSums(geneMatrix) > 0)
	
	# Identification of tissue-enriched isoforms
	if(length(viableIdx) > 1) {
		# Filter geneMatrix and normMatrix by viableIdx
		geneMatrix <- geneMatrix[,viableIdx]
		normMatrix <- normMatrix[,viableIdx]
		
		# Perform a binomial test for each cell of the matrix and return FDR-corrected p-values
		pValMatrix <- binomial_test(geneMatrix)
		
		# Identify cells indices with FDR < 5%
		idx <- which(pValMatrix < 0.05, arr.ind = TRUE)
		
		# Output tissue-enriched isoforms
		returnString <- c()
		if(nrow(idx) > 0) {
			for(i in 1:nrow(idx)) {
				# Extract information specific to each tissue-enriched isoform
				id <- rownames(pValMatrix)[idx[i,1]]
				tissue <- colnames(pValMatrix)[idx[i,2]]
				isoPct <- normMatrix[id, tissue]
				isoRC <- geneMatrix[id, tissue]
				sjCate <- get_BSJ_FSJ(id, isoformLL, inDF)
				
				# Require that the tissue-enriched isoform has at least 2 reads in given sample
				if(isoRC >= 2) {
					# Require the difference in isoform proprtion for the tissue-enriched isoform against 
					# all other samples be >= 5%
					if(all(isoPct - normMatrix[id,][-idx[i,2]] >= 0.05)) {
						# Report adjusted p-value from Benjamini Hochberg correction
						pValue <- pValMatrix[id, tissue]
						info <- c(gene, nrow(geneMatrix), id, tissue, isoPct, isoRC, pValue, sjCate, 'ENRICHED')
						returnString <- c(returnString, paste(info, collapse=','))
					}
				}
			}
		}
		if(length(returnString) == 0) {
			return(NA)
		}
		else {
			return(paste(returnString, collapse=';'))
		}	
	}
	# Identification of tissue-exclusive isoforms
	else {
		tissue <- colnames(geneMatrix)[viableIdx]
		
		# Filter columns of geneMatrix and propMatrix by viableIdx
		geneMatrix <- geneMatrix[,viableIdx]
		normMatrix <- normMatrix[,viableIdx]
		
		# Output tissue-exclusive isoforms
		returnString <- c()
		for(id in rownames(geneMatrix)) {
			isoPct <- normMatrix[id]
			isoRC <- geneMatrix[id]
			sjCate <- get_BSJ_FSJ(id, isoformLL, inDF)
			info <- c(gene, nrow(geneMatrix), id, tissue, isoPct, isoRC, 0, sjCate, 'EXCLUSIVE')
			returnString <- c(returnString, paste(info, collapse=','))
		}
		return(paste(returnString, collapse=';'))
	}
}

write_specific <- function(results, outfile) {
	# Filter out results that are NA
	results <- results[!is.na(results)]
	
	# Write results to the outfile if there are any
	if(length(results) > 0) {
		results <- unlist(strsplit(unlist(results), ';'))
		outDF <- data.frame(do.call('rbind', strsplit(results, ',')))
		
		# The output file contains the following columns describing tissue-specific isoforms:
		#	1. gene name
		#	2. number of detected circRNA isoforms
		#	3. isoform ID
		#	4. tissue/cell-line group associated with the tissue-specific isoform
		#	5. the read count proportion of the tissue-specific isoform
		#	6. the read count of the tissue-specific isoform
		#	7. adjusted p-value
		#	8. BSJ category of the tissue-specific isoform
		#	9. FSJ category of the tissue-specific isoform
		#	10. classification of tissue-specific isoform as either: (i) tissue-exclusive or 
		#		(ii) tissue-enriched (see definitions in Methods)
		
		colnames(outDF) <- c('gene', 'numIsoforms', 'isoformID', 'group', 'isoformPct', 
			'isoformReadCounts', 'pAdj', 'BSJCate', 'FSJCate', 'status')
		outDF[,'gene'] <- as.character(outDF[,'gene'])
		outDF[,'numIsoforms'] <- as.numeric(as.character(outDF[,'numIsoforms']))
		outDF[,'isoformID'] <- as.character(outDF[,'isoformID'])
		outDF[,'group'] <- as.character(outDF[,'group'])
		outDF[,'isoformPct'] <- as.numeric(as.character(outDF[,'isoformPct']))
		outDF[,'isoformReadCounts'] <- as.numeric(as.character(outDF[,'isoformReadCounts']))
		outDF[,'pAdj'] <- as.numeric(as.character(outDF[,'pAdj']))
		outDF[,'BSJCate'] <- as.character(outDF[,'BSJCate'])
		outDF[,'FSJCate'] <- as.character(outDF[,'FSJCate'])
		outDF[,'status'] <- as.character(outDF[,'status'])
		
		# Write data frame to the output file
		write.table(outDF, file=outfile, quote=FALSE, sep='\t', row.names=FALSE)
	}
}

main <- function() {
	moduleSummary <- 'Identifies tissue-specific and tissue-stable circRNA isoforms.'
	parser <- ArgumentParser(description = moduleSummary)
	
	# Add options to the argument parser
	parser$add_argument('-i', metavar='/path/to/input/file', required=TRUE, help='path to input file')
	parser$add_argument('-p', metavar='###', default=1, help='number of processors')
	parser$add_argument('-o', metavar='/path/to/output/prefix', required=TRUE, help='path to output prefix')
	
	# Parse command line arguments
	pargs <- parser$parse_args()
	infile <- pargs$i
	nproc <- as.numeric(pargs$p)
	outprefix <- pargs$o
	
	# Create names for output files using the outprefix
	stableOut <- paste(outprefix,'tissue_stable.tsv',sep='.')
	specificOut <- paste(outprefix,'tissue_specific.tsv',sep='.')
	
	# Read input file as a data frame, filter for (i) isoformID, (ii) geneName, (iii) BSJCate,
	# (iv) FSJCate, (v) tissue/cell-line columns with group read count sums. Isoforms with no
	# gene assignment are filtered out.
	inDF <- filter_df(infile)
	
	# Create a list of lists for genes and isoform IDs assigned to each circRNA isoform
	geneLL <- strsplit(as.character(inDF$geneName), ',')
	isoformLL <- strsplit(as.character(inDF$isoformID), ',')
	
	# Based on geneLL, establish a list of genes that have at least 2 detected circRNA isoforms
	geneList <- filter_genes(geneLL)
	
	# Identify tissue-stable isoforms
	stableResults <- mclapply(geneList, 
		get_stable_isoforms,
		geneLL=geneLL,
		isoformLL=isoformLL,
		inDF=inDF, 
		mc.cores=nproc
	)
	
	# Write results to corresponding output file
	write_stable(stableResults, stableOut)
	
	# Identify genes harboring candidate tissue-specific isoforms
	specificGenesDF <- mclapply(geneList, 
		get_specific_genes,
		geneLL=geneLL,
		isoformLL=isoformLL,
		inDF=inDF, 
		mc.cores=nproc
	)
	
	# Filter current gene list for genes that are likely to have
	# tissue-specific isoforms
	filteredGeneList <- filter_specific_genes(specificGenesDF)
	
	# Identify tissue-specific isoforms
	specificResults <- mclapply(filteredGeneList, 
		get_specific_isoforms,
		geneLL=geneLL,
		isoformLL=isoformLL,
		inDF=inDF, 
		mc.cores=nproc
	)
	
	# Write results to corresponding output file
	write_specific(specificResults, specificOut)
	
}

main()