# CMSeq #
 
* Provides interface for .bam files
* Reference free consensus
* Breadth and Depth of coverage

**Requires:**

* samtools (> 1.x)
* numpy
* pysam
* pandas

**Note: CMSeq can be used [as python module](README_class.md) as well**

## breadth_depth.py - Breadth and Depth of coverage

Provides breadth and depth of coverage the references of BAM alignment file, in tabular format. The file must be indexed and sorted (alternatively, --sortindex can be used).


```
usage: breadth_depth.py [-h] [-c REFERENCE ID] [-f] [--sortindex]
                        [--minlen MINLEN] [--minqual MINQUAL]
                        [--mincov MINCOV] [--truncate TRUNCATE]
                        BAMFILE
```


```
positional arguments:
  BAMFILE               The file on which to operate

optional arguments:
  -h, --help            show this help message and exit
  -c REFERENCE ID, --contig REFERENCE ID
                        Focus on a subset of references in the BAM file. Can
                        be a list of references separated by commas or a FASTA
                        file (the IDs are used to subset)
  -f                    If set unmapped (FUNMAP), secondary (FSECONDARY), qc-
                        fail (FQCFAIL) and duplicate (FDUP) are excluded. If
                        unset ALL reads are considered (bedtools genomecov
                        style). Default: unset
  --sortindex           Sort and index the file
  --minlen MINLEN       Minimum Reference Length for a reference to be
                        considered
  --minqual MINQUAL     Minimum base quality. Bases with quality score lower
                        than this will be discarded. This is performed BEFORE
                        --mincov. Default: 30
  --mincov MINCOV       Minimum position coverage to perform the polymorphism
                        calculation. Position with a lower depth of coverage
                        will be discarded (i.e. considered as zero-coverage
                        positions). This is calculated AFTER --minqual.
                        Default: 1
  --truncate TRUNCATE   Number of nucleotides that are truncated at either
                        contigs end before calculating coverage values.

```
 

Breadh and Depth of coverage outputs a table with the breadth of coverage, average and median depth-of-coverage of each reference. Values are calculated only on the covered portion of the reference:

|contig|Breadth|Depth (avg)|Depth (median)|
|------|-------|-----------|--------------|
|EF401177.1.1491|0.101274312542|1.0|1.0|
|EF405039.1.1494|0.101070950469|2.69536423841|3.0|
|all_contigs|-|1.84768211921|1.0|

The last line is a summary line calculated as if all the reads were coming from the same (big) contig.

### Examples: ###

Extract breadth and depth of coverage for all the references within a sorted and indexed `BAM` file


```
breadth_depth.py mybam.sorted.bam
```

Extract breadth and depth of coverage for all the references within an unsorted `BAM` file 


```
breadth_depth.py --sortindex mybam.sorted.bam 
```

Extract breadth and depth of coverage for all the references within a sorted `BAM` file, count only the reads with minimum quality of 25 and positions with a minimum coverage of 10 


```
breadth_depth.py --mincov 10 --minqual 20 mybam.bam
```

Extract breadth and depth of coverage for the references: genome_1 and genome_2 within a sorted `BAM` file


```
breadth_depth.py -c genome_1,genome_2 mybam.bam
```

Extract breadth and depth of coverage for the references present in MYFASTA.fasta, within a sorted `BAM` file


```
breadth_depth.py -c MYFASTA.fasta mybam.sorted.bam
```

## poly.py - Polymorphic Rate

Provides the Polymorphic-rate of each reference in a sorted and indexed BAMFILE. The polymorphic rate is defined as: number_of_polymorhpic_sites / number_of_total_nucleotides. Beware that *number_of_total_nucleotides* depends on --minqual and --mincov, as if a position is not covered (e.g. coverage = 0) will not be counted in the denominator.


```
usage: poly.py [-h] [-c REFERENCE ID] [-f] [--sortindex] [--minlen MINLEN]
               [--minqual MINQUAL] [--mincov MINCOV] [--pvalue PVALUE]
               [--seq_err SEQ_ERR] [--dominant_frq_thrsh DOMINANT_FRQ_THRSH]
               BAMFILE

Reports the polymorpgic rate of each reference (polymorphic bases / total
bases). Focuses only on covered regions (i.e. depth >= 1)

positional arguments:
  BAMFILE               The file on which to operate

optional arguments:
  -h, --help            show this help message and exit
  -c REFERENCE ID, --contig REFERENCE ID
                        Focus on a subset of references in the BAM file. Can
                        be a list of references separated by commas or a FASTA
                        file (the IDs are used to subset)
  -f                    If set unmapped (FUNMAP), secondary (FSECONDARY), qc-
                        fail (FQCFAIL) and duplicate (FDUP) are excluded. If
                        unset ALL reads are considered (bedtools genomecov
                        style). Default: unset
  --sortindex           Sort and index the file
  --minlen MINLEN       Minimum Reference Length for a reference to be
                        considered. Default: 0
  --minqual MINQUAL     Minimum base quality. Bases with quality score lower
                        than this will be discarded. This is performed BEFORE
                        --mincov. Default: 30
  --mincov MINCOV       Minimum position coverage to perform the polymorphism
                        calculation. Position with a lower depth of coverage
                        will be discarded (i.e. considered as zero-coverage
                        positions). This is calculated AFTER --minqual.
                        Default:1
  --pvalue PVALUE       Binomial p-value threshold for the binomal-polymorphic
                        test. Default: 0.01
  --seq_err SEQ_ERR     Sequencing error rate. Default: 0.001
  --dominant_frq_thrsh DOMINANT_FRQ_THRSH
                        Cutoff for degree of `allele dominance` for a position
                        to be considered polymorphic. Default: 0.8

```

The output is strucutred as follows:


|referenceID|dominant_allele_distr_mean|dominant_allele_distr_perc_10|…|dominant_allele_distr_sd|tot_covered_bases|tot_polymorphic_bases|polymorphic_rate|
|----|----|----|----|----|----|----|----|----|
|EF401177.1.1491|-|-…|-|151.00|0.00|0.00|
|EF405039.1.1494|0.65|0.67|…|0.04|151.00|13.00|0.09|
|-GENOME-|0.65|0.67|…|0.04|302.00|13.00|0.04|

As for ``breadh_depth.py``, also the polymorphic rate analyisis is subjected to ``mincov``, ``minqual``, and ``minlen``. Additionally, two parameters can be set to decide when a site is polymorphic:


* ``dominant_frq_thrsh`` is a percentage: if the majoritary allele frequency at position x is greater than the threshold, x is considered non-polymorphic. Otherwise, a binomial test is performed to assure that x is polymorpfic (polymorphic if p < ``pvalue``)


### Examples ###

Extract polymorphic rate from a sorted and indexed bam file 


```
poly.py mybam.sorted.bam
```


Extract polymorphic rate from an unsorted bam file 


```
poly.py --sortindex mybam.sorted.bam 
```


Extract polymorphic rate from an unsorted bam file, counting only bases with minimum quality of 30 and minimum position-coverage of 10


```
poly.py --sortindex --mincov 10 --minqual 30 mybam.unsorted.bam
```


Extract polymorphic rate from an unsorted bam file, only for reads aligning against genome_1 or genome_2. Consider polymorphic only sites with majoritary-allele-freq < 70%


```
poly.py --sortindex -c genome_1,genome_2 --dominant_frq_thrsh 0.7 mybam.unsorted.bam
```

 
## consensus.py - Reference Free (but guided) consensis

Provides the Reference Free consensus for the references in a BAM alignment file, reconstructing the sequence from the raw reads, in FASTA format to standard output. The file must be indexed and sorted (alternatively, --sortindex can be used). Note that the length of the reconstructed sequence is bound to the original length of the reference. On that length, not all the positions may be covered. This can happen because:

* there are no reads mapping to the position
* there are too few reads (*i.e < ``mincov``*) mapping to the position
* the reads that map to the position have a low quality (*i.e. < ``minqual``*)
* the distribution of nucleotides at that position is potentially problematic (*i.e. dominant_allele_frequency < ``dominant_frq_thrsh``*): in this case, the position is excluded to reduce noise.


```
usage: consensus.py [-h] [-c REFERENCE ID] [-f] [--sortindex]
                    [--minqual MINQUAL] [--mincov MINCOV]
                    [--dominant_frq_thrsh DOMINANT_FRQ_THRSH]
                    [--minlen MINLEN]
                    BAMFILE

outputs the consensus in FASTA format. Non covered positions (or quality-
trimmed positions) are reported as a dashes: -

positional arguments:
  BAMFILE               The file on which to operate

optional arguments:
  -h, --help            show this help message and exit
  -c REFERENCE ID, --contig REFERENCE ID
                        Focus on a subset of references in the BAM file. Can
                        be a list of references separated by commas or a FASTA
                        file (the IDs are used to subset)
  -f                    If set unmapped (FUNMAP), secondary (FSECONDARY), qc-
                        fail (FQCFAIL) and duplicate (FDUP) are excluded. If
                        unset ALL reads are considered (bedtools genomecov
                        style). Default: unset
  --sortindex           Sort and index the file
  --minqual MINQUAL     Minimum base quality. Bases with quality score lower
                        than this will be discarded. This is performed BEFORE
                        --mincov. Default: 30
  --mincov MINCOV       Minimum position coverage to perform the polymorphism
                        calculation. Position with a lower depth of coverage
                        will be discarded (i.e. considered as zero-coverage
                        positions). This is calculated AFTER --minqual.
                        Default: 0
  --dominant_frq_thrsh DOMINANT_FRQ_THRSH
                        Cutoff for degree of `allele dominance` for a position
                        to be considered polymorphic. Default: 0.8
  --minlen MINLEN       Minimum Reference Length for a reference to be
                        considered. Default: 0


```



Note that positions with a majoritary allele frequency lower than dominant_frq_thrsh will be considered "problematic" and substituted with a "-", even with sufficient coverage and quality.



```
consensus.py ~/tmp.bam.sorted -c EF401177.1.1491,EF405039.1.1494 --mincov 1 --dominant_frq_thrsh 0.5
>EF401177.1.1491_consensus
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
-----------------------------------TACGTAGGGGGCAAGCGTTATCCGG
ATTTACTGGGTGTAAAGGGAGCGTAGACGGCGAGACAAGTCTGAAGTGAAAGCCCGGGGC
TCAACCCCGGGACTGCTTTGGAAACTGCCTTGCTAGAGTGCTGGAGAGGTAAGTGGAATT
CCTAGT------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
---------------------------------------------------
>EF405039.1.1494_consensus
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
-------------------------------------TACGTAGGTGGCAAGCGTTATCC
GGATTTACTGGGTGTAAAGGGCGTGCAGCCGGGTCTGCAAGTCAGATGTGAAATCCATGG
GCTCAACCCATGAACTGCATTTGAAACTGTAGATCTTGAGTGTCGGAGGGGCAATCGGAA
TTCCTAGT----------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------
```


### Examples ###


Extract the consensus from all the references from a sorted and indexed BAM file, in FASTA format:


```
consensus.py mybam.sorted.bam
```

Extract the consensus from all the references from an unsorted BAM file, in FASTA format:


```
consensus.py --sortindex mybam.sorted.bam 
```

Extract the consensus of genome_1 and genome_2 from a BAM file. Positions with coverage lower than 5 are ignored (- is reported instead of base-call):


```
consensus.py --mincov 5 -c genome_1,genome_2 mybam.sorted.bam
```

Extract the consensus of genome_1 and genome_2 from a BAM file. Positions with coverage lower than 5 "high quality" bases are ignored (- is reported instead of base-call). Additionally, positions with less than 50% majoritary-letters will be substituted by a "-":


```
consensus.py --mincov 5 --minqual 30 -c genome_1,genome_2 --dominant_frq_thrsh 0.5 mybam.sorted.bam
```

Same as above, but a FASTA file is used to filter references instead:


```
consensus.py --mincov 5 --minqual 30 -c FILTER_FASTA.fasta --dominant_frq_thrsh 0.5 mybam.sorted.bam
```