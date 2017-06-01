# CMSeq #
 
* Provides interface for .bam files
* reference free consensus
* Breadth and Depth of coverage
* Coverage plots

**Requires:**

* samtools (> 1.x)
* numpy
* pysam
* matplotlib
* seaborn


**Note: CMSeq can be used [as python module](README_class.md) as well**

## Usage as Python Program ##


```
usage: cmseq.py [-h] {bd,consensus,coverageplot} ...
```

### Subcommand `bd` (Breadth-Depth) ###

Provides breadth and depth of coverage for the contigs in a BAM alignment file, in tabular format. The file must be indexed and sorted (alternatively, --sortindex can be used)

```
usage: cmseq.py bd [-h] [-c CONTIG] [-f] [--sortindex] BAMFILE

calculate the Breadth and Depth of coverage of BAMFILE. Focuses only on
covered regions (i.e. depth >= 1)

positional arguments:
  BAMFILE               The file on which to operate

optional arguments:
  -h, --help            show this help message and exit
  -c REFERENCE ID, --contig REFERENCE ID
                        Gets the breadth and depth of a specific reference
                        within a BAM Can be a string or a list of strings
                        separated by comma.
  -f                    If set unmapped (FUNMAP), secondary (FSECONDARY), qc-
                        fail (FQCFAIL) and duplicate (FDUP) are excluded. If
                        unset ALL reads are considered (bedtools genomecov
                        style). Default: unset
  --sortindex           Sort and index the file
```

Examples:

```
# extract breadth and depth of coverage from a sorted and indexed bam file
cmseq.py bd mybam.sorted.bam

# extract breadth and depth of coverage from an unsorted bam file
cmseq.py bd --sortindex mybam.sorted.bam 

# extract breadth and depth of coverage from an unsorted bam file, only for reads aligning against genome_1 or genome_2
cmseq.py bd --sortindex -c genome_1,genome_2 mybam.sorted.bam
```

### Subcommand `consensus` ###

Provides the Reference Free consensus for the contigs in a BAM alignment file, in FASTA format to standard output. The file must be indexed and sorted (alternatively, --sortindex can be used)

```
usage: cmseq.py consensus [-h] [-c CONTIG] [-f] [--mincov MINCOV] [--sortindex] BAMFILE

outputs the consensus in FASTA format

positional arguments:
  BAMFILE               The file on which to operate

optional arguments:
  -h, --help            show this help message and exit
  -c CONTIG, --contig CONTIG
                        Gets the breadth and depth of a specific reference
                        within a BAM Can be a string or a list of strings
                        separated by comma.
  -f                    If set unmapped (FUNMAP), secondary (FSECONDARY), qc-
                        fail (FQCFAIL) and duplicate (FDUP) are excluded. If
                        unset ALL reads are considered (bedtools genomecov
                        style). Default: unset
  --sortindex           Sort and index the file
  --mincov MINCOV       Minimum read coverage (on single position) to call the
                        consensus
```

Examples:

```
# extract the consensus from all the references from a sorted and indexed BAM file, in FASTA format
cmseq.py consensus mybam.sorted.bam

# extract the consensus from all the references from an unsorted BAM file, in FASTA format
cmseq.py consensus --sortindex mybam.sorted.bam 

# extract the consensus of genome_1 and genome_2 from a BAM file. Positions with coverage lower than 5 are ignored (N is reported instead of base-call).

cmseq.py consensus --mincov 5 -c genome_1,genome_2 mybam.sorted.bam
```

### Subcommand `coverageplot` ###

Generates a linear or polar coverage plot for the contigs in a BAM alignment file.  The file must be indexed and sorted (alternatively, --sortindex can be used)

```
usage: cmseq.py coverageplot [-h] [-c CONTIG] [-f] [--smooth SMOOTH]
                             [--l_avoid] [--l_color L_COLOR]
                             [--s_color S_COLOR] [--flavour FLAVOUR]
                             [--s_avoid] [--sortindex]
                             BAMFILE

Plot the coverage of the contigs in a BAM alignment file. The file must be indexed and sorted (alternatively, --sortindex can be used).

positional arguments:
  BAMFILE               The file on which to operate

optional arguments:
  -h, --help            show this help message and exit
  

  -c CONTIG, --contig CONTIG
                        Get the breadth of a specific contig
  -f                    If set unmapped (FUNMAP), secondary (FSECONDARY), qc-
                        fail (FQCFAIL) and duplicate (FDUP) are excluded. If
                        unset ALL reads are considered (bedtools genomecov
                        style). Default: unset
  --flavour {polar,linear}
                        choose from linear or polar plot
  --sortindex           Sort and index the file
  --format {svg,pdf,png}
                        Output format
  --smooth SMOOTH       Smooth factor. Default: 0 (no smooth)
  --l_avoid             Suppresses line-plot
  --l_color L_COLOR     Line color for matplotlib, in HTML format (#XXXXXX)
  --s_color S_COLOR     Scatter color in HTML format (#XXXXXX)
  --s_avoid             Suppresses scatter-plot
```