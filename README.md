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
usage: cmseq.py bd [-h] [-c REFERENCE ID] [-f] [--sortindex] [--minlen MINLEN]
                   [--minqual MINQUAL] [--mincov MINCOV]
                   BAMFILE

calculate the Breadth and Depth of coverage of BAMFILE.

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
  --minlen MINLEN       Minimum Reference Length for a reference to be
                        considered
  --minqual MINQUAL     Minimum base quality. Bases with quality score lower
                        than this will be discarded. This is performed BEFORE
                        --mincov. Default: 0
  --mincov MINCOV       Minimum position coverage to perform the polymorphism
                        calculation. Position with a lower depth of coverage
                        will be discarded (i.e. considered as zero-coverage
                        positions). This is calculated AFTER --minqual.
                        Default: 1
```

Examples:

```
# extract breadth and depth of coverage from a sorted and indexed bam file
cmseq.py bd mybam.sorted.bam

# extract breadth and depth of coverage from an unsorted bam file
cmseq.py bd --sortindex mybam.sorted.bam 

# extract breadth and depth of coverage from an unsorted bam file, counting only bases with minimum quality of 30 and a minimum position-coverage of 10
cmseq.py bd --sortindex --mincoverage 10 --minqual 30 mybam.sorted.bam


# extract breadth and depth of coverage from an unsorted bam file, only for reads aligning against genome_1 or genome_2
cmseq.py bd --sortindex -c genome_1,genome_2 mybam.sorted.bam
```

### Subcommand `poly` ###

Provides the Polymorphic-rate of each reference in a sorted and indexed BAMFILE. The polymorphic rate is defined as: number_of_polymorhpic_sites / number_of_total_nucleotides. Beware that number_of_total_nucleotides depends on --minqual and --mincov, as if a position is not covered (e.g. coverage = 0) will not be counted in the denominator.


```
usage: cmseq.py poly [-h] [-c REFERENCE ID] [-f] [--sortindex]
                     [--minlen MINLEN] [--minqual MINQUAL] [--mincov MINCOV]
                     BAMFILE

Reports the polymorpgic rate of each reference (polymorphic bases / total
bases). Focuses only on covered regions (i.e. depth >= 1)

positional arguments:
  BAMFILE               The file on which to operate

optional arguments:
  -h, --help            show this help message and exit
  -c REFERENCE ID, --contig REFERENCE ID
                        Gets the polymorhic rate a specific reference
                        within a BAM Can be a string or a list of strings
                        separated by comma.
  -f                    If set unmapped (FUNMAP), secondary (FSECONDARY), qc-
                        fail (FQCFAIL) and duplicate (FDUP) are excluded. If
                        unset ALL reads are considered (bedtools genomecov
                        style). Default: unset
  --sortindex           Sort and index the file
  --minlen MINLEN       Minimum Reference Length for a reference to be
                        considered
  --minqual MINQUAL     Minimum base quality. Bases with quality score lower
                        than this will be discarded. This is performed BEFORE
                        --mincov. Default: 0
  --mincov MINCOV       Minimum position coverage to perform the polymorphism
                        calculation. Position with a lower depth of coverage
                        will be discarded (i.e. considered as zero-coverage
                        positions). This is calculated AFTER --minqual.
                        Default: 1
```

Examples:

```
# extract polymorphic rate from a sorted and indexed bam file
cmseq.py poly mybam.sorted.bam

# extract polymorphic rate from an unsorted bam file
cmseq.py poly --sortindex mybam.sorted.bam 

# extract polymorphic rate from an unsorted bam file, counting only bases with minimum quality of 30 and minimum position-coverage of 10
cmseq.py poly --sortindex --mincoverage 10 --minqual 30 mybam.sorted.bam


# extract polymorphic rate from an unsorted bam file, only for reads aligning against genome_1 or genome_2
cmseq.py poly --sortindex -c genome_1,genome_2 mybam.sorted.bam
```


### Subcommand `consensus` ###

Provides the Reference Free consensus for the contigs in a BAM alignment file, in FASTA format to standard output. The file must be indexed and sorted (alternatively, --sortindex can be used)

```
usage: cmseq.py consensus [-h] [-c CONTIG] [-f] [--sortindex]
                          [--minqual MINQUAL] [--mincov MINCOV]
                          [--minlen MINLEN]
                          BAMFILE

outputs the consensus in FASTA format. Non covered positions (or quality-
trimmed positions) are reported as a dashes: -

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
  --minqual MINQUAL     Minimum base quality. Bases with quality score lower
                        than this will be discarded. This is performed BEFORE
                        --mincov. Default: 0
  --mincov MINCOV       Minimum position coverage to perform the polymorphism
                        calculation. Position with a lower depth of coverage
                        will be discarded (i.e. considered as zero-coverage
                        positions). This is calculated AFTER --minqual.
                        Default: 1
  --minlen MINLEN       Minimum Reference Length for a reference to be
                        considered

```

Examples:

```
# extract the consensus from all the references from a sorted and indexed BAM file, in FASTA format
cmseq.py consensus mybam.sorted.bam

# extract the consensus from all the references from an unsorted BAM file, in FASTA format
cmseq.py consensus --sortindex mybam.sorted.bam 

# extract the consensus of genome_1 and genome_2 from a BAM file. Positions with coverage lower than 5 are ignored (- is reported instead of base-call).

cmseq.py consensus --mincov 5 -c genome_1,genome_2 mybam.sorted.bam

# extract the consensus of genome_1 and genome_2 from a BAM file. Positions with coverage lower than 5 "high quality" bases are ignored (- is reported instead of base-call).

cmseq.py consensus --mincov 5 --minqual 30 -c genome_1,genome_2 mybam.sorted.bam
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