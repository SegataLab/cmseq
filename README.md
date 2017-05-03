# CMSeq #

 
* Provides interface for .bam files
* reference free consensus
* Breadth and Depth of coverage
* Coverage plots

Requires samtools (> 1.x), numpy, pysam, matplotlib and seaborn

## Usage as Python Program ##
```
#!python
usage: cmseq.py [-h] {bd,consensus,coverageplot} ...
```

### Subcommand bd (Breadth-Depth) ###
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

### Subcommand: consensus ###
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

### Subcommand: coverageplot ###
```
usage: cmseq.py coverageplot [-h] [-c CONTIG] [-f] [--smooth SMOOTH]
                             [--l_avoid] [--l_color L_COLOR]
                             [--s_color S_COLOR] [--flavour FLAVOUR]
                             [--s_avoid] [--sortindex]
                             BAMFILE

Plot the coverage of a contig

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
  --smooth SMOOTH       Smooth factor. Default: 0 (no smooth)
  --l_avoid             Suppresses line-plot
  --l_color L_COLOR     Line color for matplotlib, in HTML format (#XXXXXX)
  --s_color S_COLOR     Scatter color in HTML format (#XXXXXX)
  --s_avoid             Suppresses scatter-plot

```

## Usage as Python Module ##

### class BamFile ###

Represents a collection of contig/reference of a bam file

To create a new BamContig from an unsorted BAM file:
```
#!python
collection = cmseq.BamFile(BAM_FILE_PATH,sort=True,index=True)
```

To start from a pre-sorted and indexed bam file:
```
#!python
collection = cmseq.BamFile(BAM_FILE_PATH)
```

To set the pysam stepper to a custom value (e.g. 'all', that avoids secondary alignments or 'nofilter', that includes secondary alignments):
```
#!python
#Impose a custom stepper for all the contigs of the BAMFILE
collection = cmseq.BamFile(BAM_FILE_PATH,stepper='all')
```

### class BamContig ###

Represents a contig to which some reads map against

To create a new BamContig:
*Note*: this is NOT needed if a BamFile instance has been created before, as this is done automatically for each contig within the bamfile

```
#!python
contig = cmseq.BamContig(bamHandle,contigName,contigLength)
```

* bamHandle: a pysam AlignmentFile instance, pointing to the original bam file (sorted and indexed)
* contigName: the name of the contig/reference in the bam file
* contigLength: the length of that contig/reference

**Refernece Free Consensus**

reference_free_consensus(): returns a string, long as the reference, with the consensus.

Optionally, a **custom consensus function** can be passed:

* The function takes as input a python dictionary and is applied to each column of the samtools pileup. 
* The dictionary has this structre: {'A':0,'T':0,'C':0,'G':0,'N':0} and stores the counts (coverages) for each position in each nucleotide ("N" = anything else) 
* The function must returns a char
* The default function is: *lambda array: max(array, key=array.get)* (pure majority rule).
* The function is applied only to positions with at least one covering read: other positions are reported as an "N"

Example: use a custom consensus rule: return X for each position
```
#!python
print a.get_contig_by_label("CONTIG_NAME").reference_free_consensus(consensus_rule=lambda array: 'X')
```

**Depth of Coverage**

BamContig.depth_of_coverage(): returns a tuple, with the (mean_coverage,median_coverage) values, calculated over the positions that have a coverage of at least 1 (at least one mapping read on that position)

**Breadth of Coverage**

BamContig.breadth_of_coverage: returns a float, with the percentage of the total reference length covered by at least one read

**plot_coverage**
BamContig.plot_coverage() produces a PDF file with a coverage plot. 

```
#!python
BamContig.plot_coverage(flavour='polar',path='./out.pdf',smooth=0,l_avoid=False,s_avoid=False,l_color='#000000',s_color='#000000')
```

**flavour {'polar','linear'}**: changes from polar to linear coverage plot
**path**: the path of the output image file
**smooth**: convolution window size for smoothing (default 0 = no-smoothing)
**l_avoid**: do not plot line
**s_avoid**: do not plot scatter points
**l_color**: line_color (in HTML format, string)
**s_color**: scatter color (in HTML format, string)

**set_stepper**

BamContig.set_stepper({'all','nofilter'}): resets the pysam stepper for the contig. By default the stepper is set to 'nofilter' (bedtools style).

### Examples ###

Create a new instance of a BamFile. An unsorted, unindexed bam file can be provided and will be sorted and indexed within the module:

```
#!python
import cmseq
collection = cmseq.BamFile("CONTIG_NAME",sort=True,index=True)
```

Iterate over each contig represented in the BAM/SAM file:

```
#!python
for i in collection.get_contigs():
 	print i,collection.get_contig_by_label(i).reference_free_consensus()
 	print collection.get_contig_by_label(i).depth_of_coverage()  #(mean,median)
 	print collection.get_contig_by_label(i).breadth_of_coverage()
```

Select a custom contig and get a custom consensus sequence, with "+" where coverage is higher or equal 2, - otherwise:

```
#!python
print collection.get_contig_by_label("CONTIG_NAME").reference_free_consensus(consensus_rule=lambda array: '+' if sum(array.values()) >= 2 else '-')
```

Do the same as before, without using the BamFile class, but with pysam only. The bam file needs to be sorted and indexed!

```
#!python
import pysam,cmseq
bamHandle = pysam.AlignmentFile(BAM_PATH, "rb")
lengths = dict((r,l) for r,l in zip(bamHandle.references,bamHandle.lengths))
contig = cmseq.BamContig(bamHandle,TARGET_CONTIG,lengths[TARGET_CONTIG])

print contig.reference_free_consensus(consensus_rule=lambda array: '+' if sum(array.values()) >= 2 else '-')

```