#!/usr/bin/env python3

# isoCirc (1.0.0)
# Author: Robert Wang, PhD Student (Xing Lab)
# Date: 2020.08.04
# E-mail: robwang@pennmedicine.upenn.edu

# This is a script for merging raw isocirc.out files into a single file containing only
# full-length circular RNA isoforms. This script only keeps full-length isoforms for which
# there exists at least one tissue/cell-line in which the sum of read counts across all
# respective biological replicates is greater than or equal to a specified read count
# threshold.

# Usage:
#   python merge_isocirc.py -c /path/to/configuration/file \
#       -r [read count thresold] \
#       -o /path/to/output/file

# The configuration file is a user-written comma-separated file with the following structure:
#   Field 1: tissue OR cell-line name
#   Field 2: biological replicate number  
#   Field 3: path name to corresponding isocirc.out file

# Example:
#   Hek293,1,./nano_43/isocirc.out
#   Hek293,2,./nano_44/isocirc.out
#   Hek293,3,./nano_45/isocirc.out
#   Brain,1,./nano_brain/isocirc.out
#   Blood,1,./nano_blood_peripheral_leukocytes/isocirc.out
#   SkeletalMuscle,1,./nano_skeletal_muscle/isocirc.out

# Dependencies:
#   * pandas
#   * argparse

import pandas as pd
import argparse

def write_output(mergedict, outfile):
    f = open(outfile, 'a')
    
    # Iterate through all keys in the merge dictionary
    for id in list(mergedict.keys()):
        f.write('\t'.join(map(str,list(id) + mergedict[id]))+'\n')
    
    f.close()

def meet_threshold(rcArr, groupdict, threshold):
    grpList = []
    for grp in list(groupdict.keys()):
        ctr = 0
        for idx in groupdict[grp]:
            ctr += rcArr[idx]
        grpList.append(ctr)
    if max(grpList) >= int(threshold):
        return True
    else:
        return False
        

def update_dict(fldict, rcdict, line, n, k):
    # Remove fields "isoformID" and "readIDs"
    info = line.split('\t')[1:-1]

    # Set isHighSJ to "NA" (alignment information is not important here)
    info[15] = 'NA'
    
    # Create a temporary identifier for the isoform
    id = tuple(info[0:20] + info[21:32])
    
    # Parse isFullLength field from input line
    if info[20] == 'True':
        isFullLength = True
    else:
        isFullLength = False
    
    # Parse readCount field from input line
    readCount = int(info[32])
                    
    # Update full-length status dictionary
    if id in fldict:
        # Update by taking union 
        fldict[id] = fldict[id] or isFullLength
    else:
        # Create new key
        fldict[id] = isFullLength
    
    # Update read count dictionary
    if id in rcdict:
        # Update read counts
        rcdict[id][k] += readCount
    else:
        # Create new key
        rcdict[id] = [0] * n
        rcdict[id][k] += readCount
    
    # Return the updated dictionaries
    return fldict, rcdict
    
def merge_file(config, threshold):
    # Parse the configuration file
    df = pd.read_csv(config, header=None)
    nrow = df.shape[0]
    files = list(df.loc[0:nrow,2])
    
    # Determine tissue/cell-line groups based on configuration file
    samples = list(df.loc[0:nrow,0])
    groups = list(set(samples))
    
    # Create a dictionary to store a list of indices for each tissue/cell-line group
    groupdict = {}
    for grp in groups:
        groupdict[grp] = [i for i, x in enumerate(samples) if x == grp]
    
    # Establish dictionaries to map isoforms to full length status and read counts
    fldict = {}
    rcdict = {}
    
    # Create a counter to keep track of file number
    ctr = 0
    
    # Iterate through each file
    for tsv in files:
        # Open each file and read each line
        with open(tsv) as infile:
            for line in infile:
                line = line.strip()
                # Skip header lines in the file
                if not line.startswith('#'):
                    # Update dictionaries with each line
                    fldict, rcdict = update_dict(fldict, rcdict, line, nrow, ctr)
        
        # Update counter
        ctr += 1
    
    # Create dictionary to store output
    outdict = {}
    
    # Filter for isoforms for which there exists at least one tissue/cell-line where
    # the sum of read counts across all biological replicates is greater than or equal
    # to a threshold value
    for id in list(fldict.keys()):
        # If isoform is full-length:
        if fldict[id]:
            # If isoform meets read-length threshold
            if meet_threshold(rcdict[id], groupdict, threshold):
                # Create a new ID from existing ID that includes isFullLength field
                finalID = tuple(list(id[0:20]) + ['True'] + list(id[20:31]))
                outdict[finalID] = rcdict[id]
    
    return outdict
                
def write_header(config, outfile):
    header = []
    
    # Add existing column names for isocirc.out files, removing "isoformID", "readCount", "readIDs"
    header += ['chrom', 'startCoor0based', 'endCoor', 'geneStrand', 'geneID', 'geneName', 'blockCount', 'blockSize',
        'blockStarts', 'refMapLen', 'blockType', 'blockAnno', 'isKnownSS', 'isKnownSJ', 'canoSJMotif', 'isHighSJ', 
        'isKnownExon', 'isKnownBSJ', 'isCanoBSJ', 'canoBSJMotif', 'isFullLength', 'BSJCate', 'FSJCate', 'CDS', 'UTR',
        'lincRNA', 'antisense', 'rRNA', 'Alu', 'allRepeat', 'upFlankAlu', 'downFlankAlu']
    
    # Parse the configuration file and extract sample and replicate columns
    df = pd.read_csv(config, header=None)
    nrow = df.shape[0]
    names = list(df.loc[0:nrow,0])
    reps = list(df.loc[0:nrow,1])
    
    # Add new column names for storing read counts associated with each biological sample
    header += ['{}_{}'.format(names_, reps_) for names_, reps_ in zip(names, reps)]
    
    # Print header to output file using tab-separated fields
    f = open(outfile, 'w')
    f.write('\t'.join(map(str,header))+'\n')
    f.close()

def main():
    moduleSummary = 'Merges isocirc.out files into single file with full-length circRNA isoforms.'
    parser = argparse.ArgumentParser(description=moduleSummary)
    
    # Add arguments
    parser.add_argument('-c', metavar='/path/to/configuration/file', required=True, help='path to configuration file')
    parser.add_argument('-r', metavar='###', default=2, help='read count threshold')
    parser.add_argument('-o', metavar='/path/to/output/file', required=True, help='path to output file')
    
    # Parse command line arguments
    args = parser.parse_args()
    config, threshold, outfile = args.c, args.r, args.o
    
    # Write a header for the output file based on sample information in the configuration file
    write_header(config, outfile)
    
    # Merge isocirc.out files into a single file with full-length circRNA isoforms
    mergedict = merge_file(config, threshold)
    
    # Print mergedict to outfile
    write_output(mergedict, outfile)
    
if __name__ == '__main__':
	main()