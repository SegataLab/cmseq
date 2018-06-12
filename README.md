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

## Usage as Python Program ##

### `breadth_depth.py` - Breadth and Depth of coverage

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
 
#### Examples: ####

** Extract breadth and depth of coverage for all the references within a sorted and indexed `BAM` file **

```
breadth_depth.py mybam.sorted.bam
```

** Extract breadth and depth of coverage for all the references within an unsorted `BAM` file **

```
breadth_depth.py --sortindex mybam.sorted.bam 
```

** Extract breadth and depth of coverage for all the references within a sorted `BAM` file, count only the reads with minimum quality of 25 and positions with a minimum coverage of 10 ***

```
breadth_depth.py --mincov 10 --minqual 20 mybam.bam
```

** Extract breadth and depth of coverage for the references: genome_1 and genome_2 within a sorted `BAM` file **

```
breadth_depth.py -c genome_1,genome_2 mybam.bam
```

** Extract breadth and depth of coverage for the references present in MYFASTA.fasta, within a sorted `BAM` file **

```
breadth_depth.py -c MYFASTA.fasta mybam.sorted.bam
```

### `poly.py` - Polymorphic Rate

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

* referenceID
* dominant_allele_distr_mean
* dominant_allele_distr_percentile_10
* dominant_allele_distr_percentile_20
* dominant_allele_distr_percentile_30
* dominant_allele_distr_percentile_40
* dominant_allele_distr_percentile_50
* dominant_allele_distr_percentile_60
* dominant_allele_distr_percentile_70
* dominant_allele_distr_percentile_80
* dominant_allele_distr_percentile_90
* dominant_allele_distr_percentile_95
* dominant_allele_distr_percentile_98
* dominant_allele_distr_percentile_99
* dominant_allele_distr_sd
* total_covered_bases
* total_polymorphic_bases
* total_polymorphic_rate

As for ``breadh_depth.py``, also the polymorphic rate analyisis is subjected to ``mincov``, ``minqual``, and ``minlen``. Additionally, two parameters can be set to decide when a site is polymorphic:

* ``dominant_frq_thrsh`` is a percentage: if the majoritary allele frequency at position x is greater than the threshold, x is considered non-polymorphic. Otherwise, a binomial test is performed to assure that x is polymorpfic (polymorphic if p < ``pvalue``)

#### Examples ####

** extract polymorphic rate from a sorted and indexed bam file **
```
poly.py mybam.sorted.bam
```
** extract polymorphic rate from an unsorted bam file **
```
poly.py --sortindex mybam.sorted.bam 
```
** extract polymorphic rate from an unsorted bam file, counting only bases with minimum quality of 30 and minimum position-coverage of 10 **
```
poly.py --sortindex --mincov 10 --minqual 30 mybam.unsorted.bam
```

** extract polymorphic rate from an unsorted bam file, only for reads aligning against genome_1 or genome_2. Consider polymorphic only sites with majoritary-allele-freq < 70% **
```
poly.py --sortindex -c genome_1,genome_2 --dominant_frq_thrsh 0.7 mybam.unsorted.bam
```
 
 ### Subcommand `consensus` ###

Provides the Reference Free consensus for the references in a BAM alignment file, reconstructing the sequence from the raw reads, in FASTA format to standard output. The file must be indexed and sorted (alternatively, --sortindex can be used)


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

#### Examples ####


** extract the consensus from all the references from a sorted and indexed BAM file, in FASTA format **
```
consensus.py mybam.sorted.bam
```
** extract the consensus from all the references from an unsorted BAM file, in FASTA format **
```
consensus.py --sortindex mybam.sorted.bam 
```
** extract the consensus of genome_1 and genome_2 from a BAM file. Positions with coverage lower than 5 are ignored (- is reported instead of base-call).**
```
consensus.py --mincov 5 -c genome_1,genome_2 mybam.sorted.bam
```
** extract the consensus of genome_1 and genome_2 from a BAM file. Positions with coverage lower than 5 "high quality" bases are ignored (- is reported instead of base-call). Additionally, positions with less than 50% majoritary-letters will be substituted by a "-"**
```
consensus.py --mincov 5 --minqual 30 -c genome_1,genome_2 --dominant_frq_thrsh 0.5 mybam.sorted.bam
```
**Same as above, but a FASTA file is used to filter references instead**
```
consensus.py --mincov 5 --minqual 30 -c FILTER_FASTA.fasta --dominant_frq_thrsh 0.5 mybam.sorted.bam
```

 