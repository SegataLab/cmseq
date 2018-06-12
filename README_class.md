# CMSeq #

 
* Provides interface for .bam files
* reference free consensus
* Breadth and Depth of coverage 

Requires samtools (> 1.x), numpy, pysam, matplotlib and seaborn
 
## Usage as Python Module ##

### class BamFile ###

Represents a collection of contig/reference of a bam file

To create a new BamContig from an unsorted BAM file:
```
#!python
collection = cmseq.BamFile(BAM_FILE_PATH,sort=True,index=True,minlen=0)
```

an optional argument ``filterInputList`` can be passed to BamFile, to filter only some references. ``filterInputList`` can be:
* a string of comma-separated IDs
* the path to a FASTA file with the to-be-filtered IDs as FASTA IDs

To start from a pre-sorted and indexed bam file:
```
#!python
collection = cmseq.BamFile(BAM_FILE_PATH)
```

To set the pysam stepper to a custom value (e.g. `all`, that avoids secondary alignments or `nofilter`, that includes secondary alignments):
```
#!python
#Chose a custom stepper for all the contigs of the BAMFILE
collection = cmseq.BamFile(BAM_FILE_PATH,stepper='all')
```

To take into accounts only references (/contigs) longer than N, use `minlen`:
```
#!python
#Build the collection only on contigs / references longer than 5000
collection = cmseq.BamFile(BAM_FILE_PATH,minlen=5000)
```

### class BamContig ###

Represents a reference to which some reads map against

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

The function can use the optional parameters:

* `minqual`: the consensus will be based only on those nucleotides with a mapping-quality higher than minqual. **Default: 0**, meaning everything is used
* `mincov`: the consensus will be based only on those positions with at least MINCOV coverage (after the quality filtering of `minqual`). **Default: 1**, meaning everything is used.

* `consensus_rule`: a custom consensus function that: 
takes as input a python dictionary. The function is applied to each column of the samtools pileup.
The dictionary has this structre: {'A':0,'T':0,'C':0,'G':0,'N':0} and stores the counts (coverages) for each position in each nucleotide ("N" = anything else). The function must return a char
The default function is: `lambda array: max(array, key=array.get)` (pure majority rule).
The function is applied only to positions that meet the requirements of `minqual` and `mincov`. Other positions are reported as "-"

Examples
```
#!python
# Get the simplest majority rule (default) consensus of REFERENCE_NAME:
print a.get_contig_by_label("REFERENCE_NAME").reference_free_consensus()

# Get the simplest majority rule (default) consensus of REFERENCE_NAME considering positions covered by at least 5 reads with qualities higher than 33:
print a.get_contig_by_label("REFERENCE_NAME").reference_free_consensus(mincov=5,minqual=33)

# Use a custom consensus rule: return X for each position
print a.get_contig_by_label("REFERENCE_NAME").reference_free_consensus(consensus_rule=lambda array: 'X')
```

**Depth of Coverage**

BamContig.**depth_of_coverage()**: returns a tuple, with the (mean_coverage,median_coverage) values, calculated over the positions that have a coverage of at least 1 (at least one mapping read on that position). Optionally, can take:

* `minqual`: the nucleotides considered are only those that have a quality score higher than MINQUAL. **Default: 0**, meaning everything is used
* `mincov`: the depth is based only on those positions with at least MINCOV coverage (after the quality filtering of `minqual`). **Default: 1**, meaning everything is used.

**Breadth of Coverage**

BamContig.**breadth_of_coverage**: returns a float, with the percentage of the total reference length covered by reads. It takes as optional parameters `mincov` and `minqual` as *depth_of_coverage*

**Polymorphic Rate**

BamContig.**polymorphism_rate**: returns a DataFrame, with the statistics of polymorphic positions, over the total number of reconstructable positions. It takes as optional parameters `mincov` and `minqual` as *depth_of_coverage*. 
 
**Set the Pysam stepper**

BamContig.**set_stepper(VALUE)**: resets the pysam stepper for the reference. VALUE can be `all` or `nofilter`, as of the pysam specifications. By default the stepper is set to 'nofilter' (bedtools style).

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
Select a custom contig and get its consensus sequence by majoriy rule:
```
#!python
print collection.get_contig_by_label("REFERENCE_NAME").reference_free_consensus()
```

Select a custom contig and plot its coverage
```
#!python
collection.get_contig_by_label("REFERENCE_NAME").plot_coverage('out.pdf')
```

Select a custom contig and get its consensus sequence by majoriy rule, only for positions covered by at least 10 high quality reads:

```
#!python
print collection.get_contig_by_label("REFERENCE_NAME").reference_free_consensus(mincov=10,minqual=33)
```

Select a custom contig and get a custom consensus sequence, with "+" where coverage is higher or equal 2, - otherwise:

```
#!python
print collection.get_contig_by_label("REFERENCE_NAME").reference_free_consensus(consensus_rule=lambda array: '+' if sum(array.values()) >= 2 else '-')
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