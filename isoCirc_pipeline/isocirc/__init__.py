__program__ = "isocirc"
__version__ = "1.0.2"

whole_output_header = ['#readID', 'chrom', 'startCoor0based', 'endCoor', 'mapStrand',
                       'geneStrand', 'geneID', 'geneName',  # 'transID', 'transName',
                       'blockCount', 'blockSize', 'blockStarts', 'refMapLen',  # mapping information
                       'blockType', 'blockAnno',  # evaluation with whole annotation for each block
                       'readLen', 'consLen', 'consMapLen', 'copyNum', 'consFrac', 'chimInfo',
                       # original read and consensus sequence information
                       'isKnownBSJ', 'disToKnownBSJ', 'isCanoBSJ', 'disToCanoBSJ', 'canoBSJMotif', 'alignAroundCanoBSJ',
                       # back-splice junction
                       'isKnownSS', 'isKnownFSJ', 'isKnownExon', #
                       'isCanoFSJ', 'canoFSJMotif', 'isHighFSJ', # isHighSJ: high-confidence SJ based on alignment around SJ
                       'CDS', 'UTR', 'lincRNA', 'antisense',  # gene_type/biotype
                       'rRNA', 'Alu', 'allRepeat', 'upFlankAlu', 'downFlankAlu']  # repeat element
whole_output_header_idx = {h: i for i, h in enumerate(whole_output_header)}

isoform_output_header = ['#isoformID', 'chrom', 'startCoor0based', 'endCoor', # 1-4
                         'geneStrand', 'geneID', 'geneName', # 5-7
                         'blockCount', 'blockSize', 'blockStarts', 'refMapLen', # 8-11
                         'blockType', 'blockAnno', # 12-13
                         'isKnownSS', 'isKnownFSJ', 'canoFSJMotif', 'isHighFSJ', 'isKnownExon', # 14-18
                         'isKnownBSJ', 'isCanoBSJ', 'canoBSJMotif', # 19-21
                         'isFullLength', 'BSJCate', 'FSJCate', # FSM: full splice match, NIC: novel in catalog, NNC, novel and not in catalog
                         'CDS', 'UTR', 'lincRNA', 'antisense',  # 22-25 gene_type/biotype 
                         'rRNA', 'Alu', 'allRepeat', 'upFlankAlu', 'downFlankAlu', # 26-30 repeat element
                         'readCount', 'readIDs'] # 31-32
isoform_output_header_idx = {h: i for i, h in enumerate(isoform_output_header)}
