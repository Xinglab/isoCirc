# IsoCirc

## <a name="IsoCirc"></a>What is IsoCirc ?
xxx
xxx

## Table of Contents

- [What is IsoCirc?](#IsoCirc)
- [Installation](#install)
  - [Install IsoCirc with `pip`](#pip)
  - [Dependencies](#depen)
- [Getting started](#start)
- [Input and output](#input_output)
  - [Input files](#input)
  - [Output files](#output)
  - [Column explanation of detailed information](#detailed)
- [FAQ](#FAQ)
- [Contact](#contact)
- [Changelog](#change)

## <a name="install"></a>Installation
### <a name="pip"></a>Install IsoCirc with `pip`
**IsoCirc** is written with `python`, please use `pip` to install **IsoCirc**:
```
pip install isocirc           # first time installation
pip install isocirc --upgrade # update to the latest version
```
### <a name="depen"></a>Dependencies 
IsoCirc is dependent on two open-source software: [`bedtools`(>= v2.25.0)](https://bedtools.readthedocs.io/) and minimap2 [`minimap2`(>= 2.11)](https://github.com/lh3/minimap2).
Please make sure they are installed before running IsoCirc.
 
## <a name="start"></a>Getting started
Example command #1:
```
isocirc -t 8 long_circRNA.fa reference.fa gene_anno.gtf circRNA.bed output_folder
```
Example command #2:
```
isocirc -t 8 long_circRNA.fa reference.fa gene_anno.gtf circRNA.bed output_folder \
  --short-read short_read.fa \
  --Alu ./anno/hg19/alu.bed   \
  --all-repeat ./anno/hg19/all_repeat.bed
```

Detailed arguments:
```
usage: isocirc [-h] [-v] [-t THREADS] [--bedtools BEDTOOLS]
               [--minimap2 MINIMAP2] [--short-read short.fa]
               [--lordec LORDEC] [--kmer KMER] [--solid SOLID] [--trf TRF]
               [--match MATCH] [--mismatch MISMATCH] [--indel INDEL]
               [--match-frac MATCH_FRAC] [--indel-frac INDEL_FRAC]
               [--min-score MIN_SCORE] [--max-period MAX_PERIOD]
               [--min-len MIN_LEN] [--min-copy MIN_COPY]
               [--min-frac MIN_FRAC] [--high-max-ratio HIGH_MAX_RATIO]
               [--high-min-ratio HIGH_MIN_RATIO]
               [--high-iden-ratio HIGH_IDEN_RATIO]
               [--high-repeat-ratio HIGH_REPEAT_RATIO]
               [--low-repeat-ratio LOW_REPEAT_RATIO] [--Alu ALU]
               [--flank-len FLANK_LEN] [--all-repeat ALL_REPEAT]
               [-s SITE_DIS] [-S END_DIS] [--cano-motif {GT/AG,all}]
               [--bsj-score BSJ_SCORE] [--key-bsj-score KEY_BSJ_SCORE]
               [--min-circ-dis MIN_CIRC_DIS] [--rescue-low]
               long.fa ref.fa anno.gtf circRNA.bed/gtf out_dir

IsoCirc: Profiling and Annotating ciRcular RNA with Iso-Seq

positional arguments:
  long.fa               Long read data generated from long-read circRNA
                        sequencing technique.
  ref.fa                Reference genome sequence file.
  anno.gtf              Whole gene annotation file in GTF format.
  circRNA.bed/gtf       circRNA annotation file in BED12 or GTF format.
  out_dir               Output directory for final result and temporary files.

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

General options:
  -t THREADS, --threads THREADS
                        Number of thread to use. (default: 8)
  --bedtools BEDTOOLS   Path to bedtools. (default: bedtools)
  --minimap2 MINIMAP2   Path to minimap2. (default: minimap2)

Hybrid error-correction with short-read data (LoRDEC):
  --short-read short.fa
                        Short-read data for error correction. Use ',' to
                        connect multiple or paired-end short read data.
                        (default: )
  --lordec LORDEC       Path to lordec-correct. (default: lordec-correct)
  --kmer KMER           k-mer size. (default: 21)
  --solid SOLID         Solid k-mer abundance threshold. (default: 3)

Detecting tandem-repeat with TRF(Tandem Repeat Finder):
  --trf TRF             Path to trf program. (default: trf409.legacylinux64)
  --match MATCH         Match score. (default: 2)
  --mismatch MISMATCH   Mismatch penalty. (default: 7)
  --indel INDEL         Indel penalty. (default: 7)
  --match-frac MATCH_FRAC
                        Match probability. (default: 80)
  --indel-frac INDEL_FRAC
                        Indel probability. (default: 10)
  --min-score MIN_SCORE
                        Minimum alignment score to report. (default: 100)
  --max-period MAX_PERIOD
                        Maximum period size to report. (default: 2000)

Extracting and aligning consensus sequence to genome (minimap2):
  --min-len MIN_LEN     Minimum consensus length to keep. (default: 30)
  --min-copy MIN_COPY   Minimum copy number of consensus to keep. (default:
                        2.0)
  --min-frac MIN_FRAC   Minimum fraction of original long read to keep.
                        (default: 0.0)
  --high-max-ratio HIGH_MAX_RATIO
                        Maximum mappedLen / consLen ratio for high-quality
                        alignment. (default: 1.1)
  --high-min-ratio HIGH_MIN_RATIO
                        Minimum mappedLen /consLen ratio for high-quality
                        alignment. (default: 0.9)
  --high-iden-ratio HIGH_IDEN_RATIO
                        Minimum identicalBases/ consLen ratio for high-quality
                        alignment. (default: 0.75)
  --high-repeat-ratio HIGH_REPEAT_RATIO
                        Maximum mappedLen / consLen ratio for high-quality
                        self-tandem consensus. (default: 0.6)
  --low-repeat-ratio LOW_REPEAT_RATIO
                        Minimum mappedLen / consLen ratio for low-quality
                        self-tandem alignment. (default: 1.9)

Evaluating circRNA with annotation:
  --Alu ALU             Alu repetitive element annotation in BED format.
                        (default: )
  --flank-len FLANK_LEN
                        Length of upstream and downstream flanking sequence to
                        search for Alu. (default: 500)
  --all-repeat ALL_REPEAT
                        All repetitive element annotation in BED format.
                        (default: )
  -s SITE_DIS, --site-dis SITE_DIS
                        Maximum allowed distance between circRNA internal-
                        splice-site and annoated splice-site. (default: 0)
  -S END_DIS, --end-dis END_DIS
                        Maximum allowed distance between circRNA back-splice-
                        site and annoated splice-site. (default: 10)

circRNA filtering criteria:
  --cano-motif {GT/AG,all}
                        Canonical back-splicing motif (GT/AG or all three
                        motifs: GT/AG, GC/AG, AT/AC). (default: GT/AG)
  --bsj-score BSJ_SCORE
                        Minimum alignment score of sequence around the back-
                        splicing junction (<=20). (default: 18)
  --key-bsj-score KEY_BSJ_SCORE
                        Minimum alignment score of the back-splicing motif
                        (<=4). (default: 4)
  --min-circ-dis MIN_CIRC_DIS
                        Minimum distance of the start and end coordinates of
                        circRNA. (default: 150)
  --rescue-low          Use high mapping quality reads to rescue low mapping
                        quality reads. (default: False)
```
## <a name="input_output"></a>Input and output
### <a name="input_file"></a>Input files
IsoCirc takes long-read that contains muliple copies of circRNA sequence as input.

It also require reference genome sequence and gene annotation to be provided.

### <a name="output_file"></a>Output files
IsoCirc outputs three result file in a user specified directory:

| No. | File name         |  Explanation | 
|:---:|   :---            | ---        |
|  1  | IsoCirc.bed       | coordinates of circRNA with reliabe back-splice-junctions in bed12 format |
|  2  | IsoCirc.out       | detailed information for each circRNA isoform in tabular format |
|  3  | IsoCirc_stats.out | stats numbers |

### <a name="detailed"></a>Explanation of detailed information
For detailed information output file, 30 columns are generated for each circRNA isoform:

| No. | Column name         |  Explanation | 
|:---:|   :---         | ---        |
|  1  | isoformID      | assigned isoform ID |
|  2  | chrom          | chromosome ID |
|  3  | startCoor0base | start coordinate of circRNA, 0-base |
|  4  | endCoor        | end coordinate of circRNA |
|  5  | geneStrand     | gene strand (+/-) |
|  6  | geneID         | gene ID  |
|  7  | geneName       | gene name  |
|  8  | blockCount     | number of block |
|  9  | blockSize      | size of each block, separated by ',' |
|  10 | blockStarts    | relative start coordinates of each block, separated by ‘,’, refer to ‘bed12’ format for further details |
|  11 | refMapLen      | total length of circRNA |
|  12 | blockType    | category for each block, E:exon, I:intron, N:intergenic |
|  13 | blockAnno      | detailed annotation for each block, format: "TransID:E1(100)+I(50)+E2(30)", TransID is the id of corresponding transcript, E1 and E2 are the 1st and 2nd exon of the transcript.Multiple blocks are seperated by ',', multiple transcripts of one block are seperated by ';' |
|  14 | isKnownSS      | `True` if splice-site is known in whole gene annotation, `False` if not, separated by ‘,’ |
|  15 | isKnownSJ      | `True` if splice-junction is known in whole gene annotation, `False` if not, separated by ‘,’ | 
|  16 | isKnownExon    | `True` if block is known exon in whole gene annotation, `False` if not, separated by ‘,’ | 
|  17 | isKnownBSJ     | `True` if back-splice-junction is known in circRNA annotation, `False` if not | 
|  18 | isCanoBSJ      | `True` if back-splice-junction has canonical motif (GT/AG), `False` if not | 
|  19 | canoBSJMotif   | strand and motif of back-splice-junction: `-GT/AG`, `NA` if back-splice-junction is not canonical | 
|  20 | CDS            | number of bases that are mapped to CDS region |
|  21 | UTR            | number of bases that are mapped to UTR region |
|  22 | lincRNA        | number of bases that are mapped to lincRNA region |
|  23 | antisense    | number of bases that are mapped to antisense region |
|  24 | rRNA         | number of bases that are mapped to rRNA region |
|  25 | Alu            | number of bases that are mapped to Alu element, `NA` if Alu annotation is not provided |
|  26 | allRepeat      | number of bases that are mapped to all repeat element, `NA` if repeat annotation is not provided |
|  27 | upFlankAlu     | flanking alu element in upstream, `NA` if none or Alu annotation is not provided |
|  28 | downFlankAlu   | flanking alu element in downstream, ‘NA’ if none or Alu annotation is not provided |
|  29 | readCount      | number of reads that come from this circRNA isoform |
|  30 | readIDs        | ID of reads that come from this circRNA isoform, separated by ‘,’  |

### <a name="stats"></a>Stats numbers
| No. | Name         |  Explanation | 
|:---:|   :---       | ---          |
|  1  | Total_read_number            | Total number of reads |
|  2  | Total_base                   | Total number of bases |
|  3  | Total_cons_read_number       | Total number of reads that have a consensus sequence |
|  4  | Total_cons_seq_base          | Total number of bases of consensus sequence |
|  5  | Total_mappable_cons_number   | Total number of mappable consensus sequence |
|  6  | Error_rate                   | Error rate of mappable consensus sequence |
|  7  | Total_isoform                | Total number of circRNA isoforms with reliable back-splice-junction |
|  8  | Total_isoform_with_known_BSJ | Total number of circRNA isoforms with known back-splice-junction |
|  9  | Total_read_with_known_BSJ    | Total number of reads with known back-splice-junction |
|  10 | Total_isoform_with_cano_BSJ  | Total number of circRNA isoforms with canonical back-splice-junction |
|  11 | Total_read_with_cano_BSJ     | Total number of reads with canonical back-splice-junction |

## <a name="FAQ"></a>FAQ
## <a name="contact"></a>Contact

## <a name="change"></a>Changelog (v1.5.12)
1. Only output reads/isoforms with highly reliable back-splice-junctions