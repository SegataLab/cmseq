# CMSeq #


## Consensus ##

* Provides interface for .bam files
* reference free consensus

Requires samtools (> 1.x), numpy, pysam


### class BamFile ###

Represents a collection of contig/reference of a bam file

To create a new BamContig from an unsorted BAM file:
```
#!python
collection = consensus.BamFile(BAM_FILE_PATH,sort=True,index=True)
```


### class BamContig ###

Represents a contig to which some reads map against

To create a new BamContig:
*Note*: this is NOT needed if a BamFile instance has been created before, as this is done automatically for each contig within the bamfile

```
#!python
contig = consensus.BamContig(bamHandle,contigName,contigLength)
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

contig.depth_of_coverage(): returns a tuple, with the (mean_coverage,median_coverage) values, calculated over the positions that have a coverage of at least 1 (at least one mapping read on that position)

**Breadth of Coverage**

contg.breadth_of_coverage: returns a float, with the percentage of the total reference length covered by at least one read

### Examples ###

Create a new instance of a BamFile. An unsorted, unindexed bam file can be provided and will be sorted and indexed within the module:

```
#!python
import consensus
collection = consensus.BamFile("CONTIG_NAME",sort=True,index=True)
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
import pysam,consensus
bamHandle = pysam.AlignmentFile(BAM_PATH, "rb")
lengths = dict((r,l) for r,l in zip(bamHandle.references,bamHandle.lengths))
contig = consensus.BamContig(bamHandle,TARGET_CONTIG,lengths[TARGET_CONTIG])

print contig.reference_free_consensus(consensus_rule=lambda array: '+' if sum(array.values()) >= 2 else '-')

```