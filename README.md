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
#!python
usage: cmseq.py [-h] {bd,consensus,coverageplot} ...
```

### Subcommand bd (Breadth-Depth) ###

Provides breadth and depth of coverage for the contigs in a BAM alignment file, in tabular format. The file must be indexed and sorted (alternatively, --sortindex can be used)

```
#!python
usage: cmseq.py bd [-h] [-c CONTIG] [-f] [--sortindex] BAMFILE

calculate the Breadth and Depth of coverage of BAMFILE

positional arguments:
  BAMFILE               The file on which to operate

optional arguments:
  -c CONTIG, --contig CONTIG
                        Get the breadth and depth of a specific contig
  -f                    If set unmapped (FUNMAP), secondary (FSECONDARY), qc-
                        fail (FQCFAIL) and duplicate (FDUP) are excluded. If
                        unset ALL reads are considered (bedtools genomecov
                        style). Default: unset
  --sortindex           Sort and index the file
```

### Subcommand consensus ###

Provides the Reference Free consensus for the contigs in a BAM alignment file, in FASTA format to standard output. The file must be indexed and sorted (alternatively, --sortindex can be used)

```
usage: cmseq.py consensus [-h] [-c CONTIG] [-f] [--sortindex] BAMFILE

outputs the consensus in FASTA format

positional arguments:
  BAMFILE               The file on which to operate

optional arguments:
  -h, --help            show this help message and exit
  -c CONTIG, --contig CONTIG
                        Get the consensus of a specific contig
  -f                    If set unmapped (FUNMAP), secondary (FSECONDARY), qc-
                        fail (FQCFAIL) and duplicate (FDUP) are excluded. If
                        unset ALL reads are considered (bedtools genomecov
                        style). Default: unset
  --sortindex           Sort and index the file
```

### Subcommand coverageplot ###

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