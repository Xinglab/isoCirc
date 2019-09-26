# IsoCirc: circRNA profiling and analysis using long-read sequencing 

## <a name="isocirc"></a>What is IsoCirc ?
xxx
xxx

## Table of Contents

- [What is IsoCirc?](#isocirc)
- [Installation](#install)
  - [Dependencies](#depen)
  - [Install IsoCirc with `pip`](#pip)
  - [Install IsoCirc from source](#src)
- [Getting started](#start)
- [Input and output](#input_output)
  - [Input files](#input)
  - [Output files](#output)
  - [Column explanation of detailed information](#detailed)
- [FAQ](#FAQ)
- [Contact](#contact)
- [Changelog](#change)

## <a name="install"></a>Installation
### <a name="depen"></a>Dependencies 
IsoCirc is dependent on two open-source software: [`bedtools`(>= v2.25.0)](https://bedtools.readthedocs.io/) and minimap2 [`minimap2`(>= 2.11)](https://github.com/lh3/minimap2).
Please make sure they are installed before running IsoCirc.

### <a name="pip"></a>Install IsoCirc with `pip` (Note: this does not work for now, will do after being published)
**IsoCirc** is written with `python`, please use `pip` to install **IsoCirc**:
```
pip install isocirc            # first time installation
pip install isocirc --upgrade  # update to the latest version
```
### <a name="src"></a>Install IsoCirc from source
Alternatively, you can install **IsoCirc** from source:
```
git clone https://github.com/Xinglab/IsoCirc.git
cd IsoCirc/IsoCirc_pipeline
pip install -r requirements.txt  # install dependencies
export PATH=$PATH:$PWD/bin       # To permanently modify your PATH, you need to add it to your ~/.profile or ~/.bashrc file. 
python setup.py install          # install main package
```
 
## <a name="start"></a>Getting started with toy example in `test_data`
```
cd IsoCirc/test_data
isocirc -t 1 read_toy.fa chr16_toy.fa chr16_toy.gtf chr16_circRNA_toy.bed output
```

Command example:
```
isocirc -t 8 long_circRNA.fa reference.fa gene_anno.gtf circRNA.bed output_folder \
  --short-read short_read.fa --Alu alu.bed --all-repeat all_repeat.bed
```

Detailed arguments:
```
usage: isocirc.py [-h] [-v] [--type {ont,pb}] [-t THREADS]
                  [--bedtools BEDTOOLS] [--minimap2 MINIMAP2]
                  [--short-read short.fa] [--lordec LORDEC] [--kmer KMER]
                  [--solid SOLID] [-T] [--trf TRF] [--tidehunter TIDEHUNTER]
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
                  [--bsj-xid BSJ_XID] [--key-bsj-xid KEY_BSJ_XID]
                  [--min-circ-dis MIN_CIRC_DIS] [--rescue-low]
                  [--sj-xid SJ_XID] [--key-sj-xid KEY_SJ_XID]
                  long.fa ref.fa anno.gtf circRNA.bed/gtf out_dir

isocirc: circular RNA profiling and analysis using long-read sequencing

positional arguments:
  long.fa               Long read data generated from long-read circRNA
                        sequencing technique.
  ref.fa                Reference genome sequence file.
  anno.gtf              Whole gene annotation file in GTF format.
  circRNA.bed/gtf       circRNA annotation file in BED12 or GTF format. Use
                        ',' to separate multiple circRNA annotation files.
  out_dir               Output directory for final result and temporary files.

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

General options:
  --type {ont,pb}       Type of sequencing data: Oxford Nanopore(ont) or
                        Pacific Biosciences (pb). (default: ont)
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

Detecting tandem-repeats with TRF(Tandem Repeats Finder):
  -T, --use-tidehunter  Use TideHunter as tandem repeats detection tool.
                        (default: False)
  --trf TRF             Path to trf program. (default: trf409.legacylinux64)
  --tidehunter TIDEHUNTER
                        Path to TideHunter. (default: TideHunter)
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
                        Canonical back-splice motif (GT/AG or all three
                        motifs: GT/AG, GC/AG, AT/AC). (default: GT/AG)
  --bsj-xid BSJ_XID     Maximum allowed mis/ins/del for 20 bp sequence
                        alignment around the back-splice junction. (default:
                        1)
  --key-bsj-xid KEY_BSJ_XID
                        Maximum allowed mis/ins/del for 4 bp sequence
                        alignment around the back-splice junction. (default:
                        0)
  --min-circ-dis MIN_CIRC_DIS
                        Minimum distance of the start and end coordinates of
                        circRNA. (default: 150)
  --rescue-low          Use high mapping quality reads to rescue low mapping
                        quality reads. (default: False)
  --sj-xid SJ_XID       Maximum allowed mis/ins/del for 20 bp sequence
                        alignment around the internal splice junction.
                        (default: 1)
  --key-sj-xid KEY_SJ_XID
                        Maximum allowed mis/ins/del for 4 bp sequence
                        alignment around the internal splice junction.
                        (default: 0)
```

## <a name="input_output"></a>Input and output
### <a name="input_file"></a>Input files
IsoCirc takes long-read that contains muliple copies of circRNA sequence as input.

It also require reference genome sequence and gene annotation to be provided.

### <a name="output_file"></a>Output files
IsoCirc output three result file in a user specified directory:

| No. | File name         |  Explanation | 
|:---:|   :---            | ---        |
|  1  | isocirc.bed       | coordinates of circRNA with reliabe back-splice-junctions in bed12 format |
|  2  | isocirc.out       | detailed information for each circRNA isoform in tabular format |
|  3  | isocirc_stats.out | summary stats numbers for each sample |

### <a name="detailed"></a>Explanation of detailed information
For detailed information output file, 30 columns are generated for each circRNA isoform:

| No. | Column name    |  Explanation | 
|:---:|   :---         | ---        |
|  1  | isoformID      | assigned isoform ID |
|  2  | chrom          | chromosome ID |
|  3  | startCoor0base | start coordinate of circRNA, 0-base |
|  4  | endCoor        | end coordinate of circRNA |
|  5  | geneStrand     | gene strand (+/-) |
|  6  | geneID         | gene ID  |
|  7  | geneName       | gene name  |
|  8  | blockCount     | number of block |
|  9  | blockSize      | size of each block, separated by `,` |
|  10 | blockStarts    | relative start coordinates of each block, separated by `,`, refer to `bed12` format for further details |
|  11 | refMapLen      | total length of circRNA |
|  12 | blockType      | category for each block, E:exon, I:intron, N:intergenic |
|  13 | blockAnno      | detailed annotation for each block, format: "TransID:E1(100)+I(50)+E2(30)", TransID is the id of corresponding transcript, E1 and E2 are the 1st and 2nd exon of the transcript.Multiple blocks are seperated by `,`, multiple transcripts of one block are seperated by `;` |
|  14 | isKnownSS      | `True` if splice-site is known in whole gene annotation, `False` if not, separated by `,` |
|  15 | isKnownSJ      | `True` if splice-junction is known in whole gene annotation, `False` if not, separated by `,` | 
|  16 | isCanoSJ       | `True` if splice-junction is canonical (GT-AG/GC-AG/AT-AC), `False` if not, separated by `,`|
|  17 | isHighSJ       | `True` if alignment around the splice-junction is high-quality, `False` if not, separated by `,`  |
|  18 | isKnownExon    | `True` if block is known exon in whole gene annotation, `False` if not, separated by ‘,’ | 
|  19 | isKnownBSJ     | `True` if back-splice-junction is known in circRNA annotation, `False` if not, separated by `,` if multiple circRNA annotations are provided | 
|  20 | isCanoBSJ      | `True` if back-splice-junction has canonical motif (GT/AG), `False` if not | 
|  21 | canoBSJMotif   | strand and motif of back-splice-junction, e.g., `-GT/AG`, `NA` if back-splice-junction is not canonical | 
|  22 | isFullLength   | `True` if the isoform is considered as `full-length isoform`, `False` if not |
|  23 | BSJCate        | catelogs of the back-splice-junction: `FSM`/`NIC`/`NNC`, see explanation below. |
|  24 | interIsoCate   | catelogs of the internal isoform: `FSM`/`NIC`/`NNC` |
|  25 | CDS            | number of bases that are mapped to CDS region |
|  26 | UTR            | number of bases that are mapped to UTR region |
|  27 | lincRNA        | number of bases that are mapped to lincRNA region |
|  28 | antisense      | number of bases that are mapped to antisense region |
|  29 | rRNA           | number of bases that are mapped to rRNA region |
|  30 | Alu            | number of bases that are mapped to Alu element, `NA` if Alu annotation is not provided |
|  31 | allRepeat      | number of bases that are mapped to all repeat element, `NA` if repeat annotation is not provided |
|  32 | upFlankAlu     | flanking alu element in upstream, `NA` if none or Alu annotation is not provided |
|  33 | downFlankAlu   | flanking alu element in downstream, `NA` if none or Alu annotation is not provided |
|  34 | readCount      | number of reads that come from this circRNA isoform |
|  35 | readIDs        | ID of reads that come from this circRNA isoform, separated by `,`  |

### <a name="stats"></a>Stats numbers
| No. | Name           |  Explanation | 
|:---:|   :---         | ---          |
|  1  | Total reads                                      | Number of raw reads in the sample |
|  2  | Total reads with cons                            | Number of reads with consensus sequence called |
|  3  | Total mappable reads with cons                   | Number of reads with consensus sequence called, mappable to the genome |
|  4  | Total reads with candidate BSJ                   | Number of reads with consensus sequence called, mappable to the genome, and with BSJs ("candidate BSJs") |
|  5  | Total candidate BSJs                             | Number of candidate BSJs (circRNA species) |
|  6  | Total known candidate BSJs                       | Number of candidate BSJs that are known (reported in existing circRNA BSJ database - circBase / MiOncoCirc) |
|  7  | Total reads with high BSJs                       | Number of reads with consensus sequence called, mappable to the genome, and with high-confidence BSJs (based on additional inspection of alignment around BSJs) |
|  8  | Total high BSJs                                  | Number of high-confidence BSJs |
|  9  | Total known high BSJs                            | Number of high-confidence BSJs that are known |
|  10 | Total isoforms with high BSJs                    | Number of unique circRNA isoforms with high-confidence BSJs |
|  11 | Total isoforms with high BSJs cano SJs           | Number of unique circRNA isoforms with high-confidence BSJs, and all internal splice sites are canonical |
|  12 | Total isoforms with high BSJs high SJs           | Number of unique circRNA isoforms with high-confidence BSJs, and all internal splice sites are high-confidence (canonical, high-quality alignment around internal splice sites) |
|  13 | Total isoforms with high BSJ known SSs           | Number of unique circRNA isoforms with with high-confidence BSJs, and all internal splice sites are known (known is based in existing transcript GTF annotations for splice sites in linear RNA) |
|  14 | Total isoforms with high BSJs high SJs known SSs | Number of unique circRNA isoforms with high-confidence BSJs, and all internal splice sites are high-confidence and known |
|  15 | Total full-length isoforms                       | Number of unique circRNA isoforms with high-confidence BSJs, and all internal splice sites are high-confidence or known |
|  16 | Total reads for full-length isoforms             | Number of reads for unique circRNA isoforms with high-confidence BSJs, and all internal splice sites are high-confidence or known |

 * BSJ:  back-splice junction
 * SJ:   splice junction
 * SS:   splice site
 * cons: consensus sequence
 * cano: canonical
 * high: high-confidence (canonical and high-quality alignment around SJ/BSJ)
 * FSM: full splice match
 * NIC: novel in catelog
 * NNC: novel and not in catelog

## <a name="FAQ"></a>FAQ
## <a name="contact"></a>Contact

Yan Gao gaoy1@email.chop.edu

Yi Xing yi.xing@pennmedicine.upenn.edu

[github issues](https://github.com/Xinglab/IsoCirc/issues)
