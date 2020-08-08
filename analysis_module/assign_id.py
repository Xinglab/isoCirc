#!/usr/bin/env python3

# isoCirc (1.0.0)
# Author: Robert Wang, PhD Student (Xing Lab)
# Date: 2020.08.05
# E-mail: robwang@pennmedicine.upenn.edu

# This is a script for assigning IDs to each circular RNA isoform of a given gene using
# the following format: [gene name].circRNA.[rank]. The "rank" is determined by the 
# following tier list of rules: (i) highest median isoform read count sum (taken over
# all sequencing libraries of a given tissue/cell-line) across all tissues/cell-lines,
# (ii) highest mean isoform read count sum (taken over all sequencing libraries of a 
# given tissue/cell-line) across all tissues/cell-lines. For intergenic circular RNA
# isoforms, the script assigns IDs using the following format: NA.[chromosome #].circRNA.
# [sort index]. The "sort index" is determined by sorting the genomic coordinates of the
# back-splice junction. 

# Usage:
#   python assign_id.py -i /path/to/input/file \
#       -o /path/to/output/file

# Dependencies:
#   * statistics
#   * argparse

import statistics
import argparse

def get_stats(a):
    # Input is a list of integers
    # Returns a tuple (median, mean)
    return tuple([statistics.median(a), statistics.mean(a)]) 

def sum_rc(arr, idx):
    # Gathers elements in array at specified indices
    # and returns the sum
    return sum([arr[i] for i in idx])

def combine_rc(arr, groups):
    # Returns a list of read counts collapsed by the
    # group indices in the groups dictionary
    return [sum_rc(arr, groups[grp]) for grp in list(groups.keys())]

def get_bsj(arr):
    # Returns the start and end coordinates of the 
    # back-splice junction as a tuple
    start = int(arr[1])
    end = int(arr[2])
    return tuple([start, end])

def assign_id_na(outdict, nadict):
    # Iterate through all chromosomes in the NA dictionary
    for chr in list(nadict.keys()):
        # Retrieve list of NA isoforms on the given chromosome, where each isoform
        # is characterized by a list
        isoforms = nadict[chr]
        
        # Isolate coordinates of the back-splice junction as a tuple
        bsjList = [get_bsj(info) for info in isoforms]
        sortidx = sorted(range(len(bsjList)), key=lambda k: bsjList[k])
        
        # Update outdict
        for i in range(len(sortidx)):
            idx = str(i+1)
            
            # Construct isoform ID
            id = 'NA.' + chr + '.circRNA.' + idx
            info = tuple(isoforms[sortidx[i]])
            
            # Add isoform to outdict
            outdict[info] = id
    
    # Return updated dictionary
    return outdict
            
def assign_id_gene(outdict, genedict, groups):
    # Iterate through all genes in the gene dictionary
    for gene in list(genedict.keys()):
        # Retrieve list of isoforms, where each isoform is characterized by a list
        isoforms = genedict[gene]
        
        # Isolate read count matrix across samples
        rcMatrix = [list(map(int,info[32:])) for info in isoforms]
        
        # Group together read counts belonging to the same group
        grpMatrix = [combine_rc(a, groups) for a in rcMatrix]
        
        # Create a list of tuples where each tuple is (median, mean)
        stats = [get_stats(a) for a in grpMatrix]
        sortidx = list(reversed(sorted(range(len(stats)), key=lambda k: stats[k])))
        
        # Update outdict
        # Check if the isoform annotation already exists in dictionary (likely due to
        # isoforms with multiple gene assignments)
        for i in range(len(sortidx)):
            idx = str(i+1)
            
            # Construct isoform ID
            id = gene + '.circRNA.' + idx
            info = tuple(isoforms[sortidx[i]])
            
            # Check if info has already been added to dictionary
            if info in outdict:
                outdict[info] = outdict[info] + ',' + id
            else:
                outdict[info] = id
    
    # Return updated dictionary    
    return outdict
            
def update_gene(genedict, line):
    # Isoforms may be assigned to multiple genes, so split on them
    info = line.split('\t')
    genearr = info[5].split(',')
    
    # Iterate through list of genes
    for gene in genearr:
        # Check if gene already exists in dictionary
        if gene in genedict:
            # Check if info is already stored in genedict[gene]
            # to avoid adding the same annotation
            if info not in genedict[gene]:
                genedict[gene].append(info)
        else:
            # Create new instance and update with line information
            genedict[gene] = list()
            genedict[gene].append(info)

    # Return the updated dictionary
    return genedict

def update_na(nadict, line):
    info = line.split('\t')
    chr = info[0]
    
    # Check if chromosome already exists in dictionary
    if chr in nadict:
        # Append info to existing list
        nadict[chr].append(info)
    else:
        # Create new instance and update with line info
        nadict[chr] = list()
        nadict[chr].append(info)

    # Return the updated dictionary
    return nadict
    
def get_groups(header):
    groupdict = {}
    
    # Isolate columns referring to tissues/cell-lines
    samples = [s.split('_')[0] for s in header[32:]]
    groups = list(set(samples))
    for grp in groups:
        # Augment indices to represent column number in input file
        groupdict[grp] = [i for i, x in enumerate(samples) if x == grp]
    
    return groupdict

def check_intergenic(line):
    info = line.split('\t')
    gene = info[5]
    
    # Check if the circRNA isoform has a gene assignment
    if gene == 'NA':
        return True
    else:
        return False

def assign_id(infile):
    genedict = {}
    nadict = {}
    outdict = {}
    
    # Read contents of the input file
    with open(infile) as f:
        # Save header line to a header list
        hline = f.readline()
        hline = hline.strip()
        header = hline.split('\t')
        groups = get_groups(header)
        
        # Process remaining lines of file
        for line in f:
            line = line.strip()
            # Check if the isoform has a gene assignment or not
            isNA = check_intergenic(line)
            if isNA:
                # Isoform does not have a gene assignment
                # Update NA dictionary
                nadict = update_na(nadict, line)
            else:
                # Isoform has a gene assignment
                # Update gene dictionary
                genedict = update_gene(genedict, line)

    # Assign IDs to isoforms with gene assignments
    outdict = assign_id_gene(outdict, genedict, groups)
    
    # Assign IDs to intergenic isoforms
    outdict = assign_id_na(outdict, nadict)
    
    # Return the updated dictionary
    return outdict
                
def write_header(infile, outfile):
    header = ['isoformID']
    
    # Read header line from the infile
    with open(infile, 'r') as f:
        fline = f.readline()
    
    # Process contents from header line of infile
    fline = fline.strip()
    header += fline.split('\t')
    
    # Print contents of updated header to outfile
    of = open(outfile, 'w')
    of.write('\t'.join(map(str,header))+'\n')
    of.close()

def write_output(outdict, outfile):
    f = open(outfile, 'a')
    
    # Iterate through the keys in outdict
    for info in list(outdict.keys()): 
        f.write('\t'.join(map(str,[outdict[info]] + list(info)))+'\n')

    f.close()

def main():
    moduleSummary = 'Assigns IDs to full-length circRNA isoforms of each gene based on read count abundance.'
    parser = argparse.ArgumentParser(description=moduleSummary)
    
    # Add arguments
    parser.add_argument('-i', metavar='/path/to/input/file', required=True, help='path to input file')
    parser.add_argument('-o', metavar='/path/to/output/file', required=True, help='path to output file')
    
    # Parse command line arguments
    args = parser.parse_args()
    infile, outfile = args.i, args.o
    
    # Write a header for the output file based on column names in the input file
    write_header(infile, outfile)
    
    # Create dictionary mapping individual isoforms to isoform IDs
    outdict = assign_id(infile)
    
    # Print contents of outdict to the outfile
    write_output(outdict, outfile)
    
if __name__ == '__main__':
	main()
