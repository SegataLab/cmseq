#!/usr/bin/env python3
from __future__ import print_function

import logging
import multiprocessing as mp
import os
import sys
from collections import defaultdict
from typing import Dict, List, Tuple, Union

import numpy as np
import pysam
from BCBio import GFF
from Bio import Seq, SeqFeature, SeqIO, SeqRecord
from scipy import stats

__author__ = "Moreno Zolfo (moreno.zolfo@unitn.it), Nicolai Karcher (nicolai.karcher@unitn.it), Kun Huang (kun.huang@unitn.it)"
__version__ = "1.0.4"
__date__ = "21 August 2021"
logger = logging.getLogger(__name__)


def _initt(terminating_, _consensus_bamFile, _consensus_args):
    global terminating
    global consensus_args
    global consensus_bamFile
    terminating = terminating_
    consensus_args = _consensus_args
    consensus_bamFile = _consensus_bamFile


class CMSEQ_DEFAULTS:
    minqual = 30
    mincov = 1
    minlen = 0
    poly_error_rate = 0.001
    poly_pvalue_threshold = 0.01
    poly_dominant_frq_thrsh = 0.8
    trimReads = None


class BamFile:
    bam_handle = None
    bamFile = None
    contigs = {}

    def __init__(self, bamFile,
                 sort=False, index=False, stepper="nofilter", minlen=CMSEQ_DEFAULTS.minlen,
                 filtRefGenomes: Union[str, List[str]] = None, minimumReadsAligning=0):
        """
         * @description:
         * @param {bamFile} name of bamfile
         * @param {sort}    If true, sort the bamFile
         * @param {index}   If true, index the bamFile
                            or if false *,bam.bai exists, it will recognized automatically
         * @param {stepper} param in pysam
                            > The stepper controls how the iterator advances. Possible options for the stepper are
                                - all: skip reads in which any of the following flags are set:
                                        BAM_FUNMAP, BAM_FSECONDARY, BAM_FQCFAIL, BAM_FDUP
                                - nofilter: uses every single read turning off any filtering.
                                - samtools: same filter and read processing as in csamtools pileup.
         * @param {refGenomes} Filter of genome, it can in 3 different formats:
                                    1.  list of contigId
                                    2.  str of genome file of FASTA file
                                    3.  string of contigIds seperated by commas (",")
                                    3.  If not specified, all contigIds will be used
         * @return {*}
        """
        if not os.path.isfile(bamFile):
            raise Exception(bamFile + " is not accessible, or is not a file")

        fp = bamFile
        if sort:
            import subprocess
            fp += ".sorted"
            subprocess.call(["samtools", "sort", bamFile, "-o", fp])

        if index:
            pysam.index(fp)

        self.bamFile = fp
        self.bam_handle: pysam.AlignmentFile = pysam.AlignmentFile(fp, "rb")

        bamHandle = self.bam_handle
        toList = []
        if filtRefGenomes is not None:
            if isinstance(filtRefGenomes, list):
                toList = filtRefGenomes
            elif os.path.isfile(filtRefGenomes):
                toList = [record.id for record in SeqIO.parse(filtRefGenomes, "fasta")]
            else:
                toList = [element for element in filtRefGenomes.split(",")]

        self.contigs = {reference_seq_name: BamContig(self.bam_handle, reference_seq_name, length, stepper)
                        for reference_seq_name, length in zip(bamHandle.references, bamHandle.lengths)
                        if (length > minlen and
                            (not toList or reference_seq_name in toList) and
                            (not minimumReadsAligning or bamHandle.count(contig=reference_seq_name, read_callback=stepper) >= minimumReadsAligning))}

    def get_contigs(self):
        return iter(self.contigs.keys())

    def get_contigs_obj(self):
        return iter(self.contigs.values())

    def get_contig_by_label(self, contigID):
        return (self.contigs[contigID] if contigID in self.contigs else None)

    def parse_gff(self, inputGFF, inputGenome: str = ""):
        """
        get a list of contigs plus 0-indexed gene-coordinates and sense-ness of protein coding regions from a gff file.
        Only tested with prokka GFF files.
        """
        genome: Dict[str, Seq.Seq] = {}
        if inputGenome:
            try:
                genome = SeqIO.index(inputGenome, "fasta")
            except Exception:
                genome = SeqIO.to_dict(inputGenome, "fasta")
        with open(inputGFF) as in_handle:
            rec: SeqRecord.SeqRecord = None
            allGenesCount = 0
            allGoodGenes = set()
            allBadStartGene = set()
            allBadEndGene = set()
            for rec in GFF.parse(in_handle):
                if str(rec.id) not in self.contigs:
                    logger.warn(f"{rec.id} is not tracked by the BAMFile.")
                    continue

                contigFeatures: Dict[str, Tuple[Tuple[int, int], int]] = {}
                r: SeqFeature.SeqFeature = None
                for r in rec.features:
                    if "minced" in r.qualifiers["source"][0] or "Minced" in r.qualifiers["source"][0]:
                        # This catches CRISPR repeats.
                        continue
                    if r.sub_features:
                        prodigal_bool = ("Prodigal" in r.sub_features[0].qualifiers["source"][0] or
                                         "prodigal" in r.sub_features[0].qualifiers["source"][0])
                    else:
                        prodigal_bool = ("Prodigal" in r.qualifiers["source"][0] or
                                         "prodigal" in r.qualifiers["source"][0])

                    if prodigal_bool:
                        # Prokka not only finds protein sequences, but also t-/r-RNA sequences. In order to only parse protein coding sequences,
                        # I search for Prodigal/Prodigal in the source entry of the sub_features attribute.

                        # the sub_features attribute of a seq_record object is apparently deprecated. I couldn"t find any other way to access
                        # the required information, though. Should probably be fixed when I can.
                        indices = (r.location.start.position, r.location.end.position)
                        sense = r.location.strand

                        if genome:
                            gene_seq: Seq.Seq = r.extract(genome[rec.id]).seq
                            geneId = f"{rec.id}_{r.id.split('_')[1]}"
                            contigLen = len(genome[rec.id])

                            if not (str(gene_seq[0:3]) == "ATG" or str(gene_seq[0:3]) == "GTG" or str(gene_seq[0:3]) == "TTG"):
                                logger.info(f"{geneId} start abnormally at {indices[0] if sense == 1 else indices[1]} of {contigLen}")
                                allBadStartGene.add(geneId)
                            if not (str(gene_seq[-3:]) == "TAG" or str(gene_seq[-3:]) == "TAA" or str(gene_seq[-3:]) == "TGA"):
                                logger.info(f"{geneId} stop abnormally at {indices[1] if sense == 1 else indices[0]} of {contigLen}")
                                allBadEndGene.add(geneId)

                            if geneId not in allBadStartGene and geneId not in allBadEndGene:
                                allGoodGenes.add(geneId)

                        contigFeatures[geneId] = (indices, sense)
                        allGenesCount += 1

                self.contigs[str(rec.id)].annotations.update(contigFeatures)
            logger.warn(f'record {allGenesCount} genes')
            if genome:
                logger.warn(f'Validation: {len(allGoodGenes)} complete genes, '
                            f'{len(allBadStartGene)} genes without start codon, {len(allGoodGenes)} genes without end codon')

    def parallel_reference_free_consensus(self, ncores=4, **kwargs):
        terminating = mp.Event()

        with mp.Pool(initializer=_initt, initargs=(terminating, self.bamFile, kwargs), processes=ncores) as pool:
            res = [x for x in pool.imap_unordered(BamFile._parallel_consensus_worker, self.contigs.keys())]
        return res

    @staticmethod
    def _parallel_consensus_worker(contigName):
        if not terminating.is_set():
            try:
                t = BamFile(consensus_bamFile, filtRefGenomes=contigName)
                return (contigName, t.get_contig_by_label(contigName).reference_free_consensus(**consensus_args))
            except Exception:
                terminating.set()
                raise
        else:
            terminating.set()


class BamContig:
    def __init__(self, bamHandle, contigName, contigLength, stepper="nofilter"):
        self.coverage = None
        self.consensus = ""
        self.annotations: Dict[str, Tuple[Tuple[int, int], int]] = {}

        self.bam_handle: pysam.AlignmentFile = bamHandle
        self.name = contigName
        self.length = contigLength
        self.stepper = self.set_stepper(stepper)

    def set_stepper(self, ns):
        if ns in ["all", "nofilter"]:
            self.stepper = ns

    def majority_rule(data_array):
        freq_array = data_array["base_freq"]

        if any(v > 0 for v in freq_array.values()):
            return max(sorted(freq_array), key=freq_array.get)
        else:
            return "N"

    def majority_rule_polymorphicLoci(data_array):
        # Masks the consensus sequence with "*" when a polymorphic locus is found according
        # to dominant_frq_thrsh defined p-value

        freq_array = data_array["base_freq"]
        poly_pvalue = data_array["p"]

        if poly_pvalue <= 0.05:
            return "*"
        elif any([v > 0 for k, v in freq_array.items() if k != "N"]):
            return max(sorted(freq_array), key=freq_array.get)
        else:
            return "N"

    def reference_free_consensus(self, consensus_rule=majority_rule,
                                 mincov=CMSEQ_DEFAULTS.mincov, minqual=CMSEQ_DEFAULTS.minqual,
                                 dominant_frq_thrsh=CMSEQ_DEFAULTS.poly_dominant_frq_thrsh,
                                 noneCharacter="-", BAM_tagFilter=None, trimReads=None):
        consensus_positions = {}

        #print("A",mincov,minqual,dominant_frq_thrsh)
        for pileupcolumn, position_data in self.get_base_stats(min_read_depth=mincov, min_base_quality=minqual,
                                                               dominant_frq_thrsh=dominant_frq_thrsh,
                                                               BAM_tagFilter=BAM_tagFilter,
                                                               trimReads=trimReads,
                                                               error_rate=CMSEQ_DEFAULTS.poly_error_rate).items():
            consensus_positions[pileupcolumn] = consensus_rule(position_data)

        if len(consensus_positions) > 0 :
            self.consensus = "".join(((consensus_positions[position]
                                      if position in consensus_positions else noneCharacter)
                                      for position in range(1, self.length + 1)))
        else:
            self.consensus = noneCharacter * self.length

        #del consensus_positions
        return self.consensus

    def baseline_PSR(self, mincov=10, minqual=30, pvalue=0.01, error_rate=0.001, dominant_frq_thrsh=0.8, binom=None):
        """
            This function estimates the polymorphic site rate over the input contig assuming that there are no truely polymorphic sites
            (Meaning that all observed polymorphisms are due to random sequencing error). The test also puts a threshold on the "dominance"
            of the allele, meaning that it only reports a polymorphic base if the binomial test indicates significance AND the base is NOT sufficiently
            dominated by the dominant base. Defaults to 0.8 dominance (dominant / all).
        """
        polymorphic_empirical_loci = 0

        # Get coverage as well as values of contig
        depthsList = self.get_all_base_values("base_cov", min_base_qualit=minqual, min_read_depth=mincov)
        # Also get dominant allele frequency of contig
        dominantFreq = self.get_all_base_values("ratio_max2all", min_base_quality=minqual, min_read_depth=mincov)

        # For each position, draw depth-times from bernoulli with success rate 1-error_rate.
        # Determine significance based on a binomial test, as you would in the regular test for polymorphism.
        for depth, da_freq in zip(depthsList, dominantFreq):
            base_max = sum(stats.bernoulli.rvs(1 - error_rate, size=depth))
            if binom and base_max in binom and depth in binom[base_max]:
                p = binom[base_max][depth]
            else:
                p = stats.binom.cdf(base_max, depth, 1.0 - error_rate)
            if p < pvalue and da_freq < dominant_frq_thrsh:
                polymorphic_empirical_loci += 1
        PSR = float(polymorphic_empirical_loci) / float(len(depthsList))
        return PSR

    def get_base_stats_for_poly(self, minqual=CMSEQ_DEFAULTS.minqual
                                ) -> List[Tuple[Tuple[int, int, int, int], int]] :
        ATCG = ("A", "C", "G", "T")

        def rev_pos(cur_pos, gene_start, gene_end):
            """
            >>> rev_pos(0, 0, 101)
            100
            >>> rev_pos(49, 0, 101)
            51
            >>> rev_pos(50, 0, 101)
            50
            """
            return gene_start - 1 - cur_pos + gene_end

        base_freq = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0}
        base_stats: List[Tuple[Tuple[int, int, int, int], int]] = []
        if not self.annotations:
            base_stats = [None] * self.length
            start, end = 0, self.length

            base_pileup: pysam.PileupColumn = None
            for base_pileup in self.bam_handle.pileup(self.name, start=start, end=end, stepper=self.stepper):
                for quality, base in zip(base_pileup.get_query_qualities(),
                                         base_pileup.get_query_sequences(add_indels=True)):
                    if quality >= minqual and base in ATCG:
                        base_freq[base] += 1

                if sum(base_freq.values()) > 0:
                    base_stats[base_pileup.pos] = ((base_freq["A"], base_freq["C"], base_freq["G"], base_freq["T"]), base_pileup.pos)
        else:
            # Generate pileups gene-wise
            # I use the "truncate" parameter to only obtain the parsed start and stop positions. Without truncate, all positions with reads covering the parsed positions are returned.
            # I wrote a function that reverses a given gene position, which is used to effectively revert genes on the anti-sense strand.
            # Furthermore, for each read"s nucleotide over a given position I write out the complement
            for (start, end), strand in self.annotations.values():
                gene_stats = [None] * (end - start)
                pos_on_gene = 0
                bam_pileup = self.bam_handle.pileup(self.name, start=start, end=end, stepper=self.stepper, truncate=True)
                if strand == 1:
                    # If the gene is on the sense-strand, do the same as before.
                    for base_pileup in bam_pileup:
                        for quality, base in zip(base_pileup.get_query_qualities(),
                                                 base_pileup.get_query_sequences(add_indels=True)):
                            if quality >= minqual and base in ATCG:
                                base_freq[base] += 1
                        if sum(base_freq.values()) > 0:
                            gene_stats[pos_on_gene] = ((base_freq["A"], base_freq["C"], base_freq["G"], base_freq["T"]), base_pileup.pos)
                        pos_on_gene += 1
                else:
                    # If the gene is on the anti-sense strand, effectively return the reverse complement
                    # by mapping positions on a gene to it"s mirrored position (using rev_pos)
                    # and then also converting each nucleotide to it"s complement.
                    for base_pileup in bam_pileup:
                        for quality, base in zip(base_pileup.get_query_qualities(),
                                                 base_pileup.get_query_sequences(add_indels=True)):
                            if quality >= minqual and base in ATCG:
                                base_freq[base] += 1
                        if sum(base_freq.values()) > 0:
                            out_pos = rev_pos(cur_pos=int(pos_on_gene), gene_start=0, gene_end=len(gene_stats))
                            contig_pos = rev_pos(cur_pos=int(base_pileup.pos), gene_start=start, gene_end=end)
                            gene_stats[out_pos] = ((base_freq["A"], base_freq["C"], base_freq["G"], base_freq["T"]), contig_pos)
                        pos_on_gene += 1

                if len(gene_stats) % 3 != 0:
                    print("One of your genes' length is not a multiple of three. Check your gff file / gene calls.")
                    print("Contig name", self.name)
                    print("Gene position", (start, end), strand)
                    sys.exit()
                base_stats.extend(gene_stats)

        return base_stats

    def easy_polymorphism_rate(self, mincov=CMSEQ_DEFAULTS.mincov, minqual=CMSEQ_DEFAULTS.minqual,
                               dominant_frq_thrsh=CMSEQ_DEFAULTS.poly_dominant_frq_thrsh):
        #list N-long where N is the number of covered bases (N <= L(contig))
        dominanceList = []
        # DN: non-synonymous mutations
        # DS: synonymous mutations
        # D?: unknown aa
        mutationStats = {"DN": 0, "DS": 0, "D?": 0}

        codon_f1 = []
        codon_f2 = []

        for positionData in self.get_base_stats_for_poly(minqual=minqual):
            # positionData= ((A,C,G,T),position) if covered, None if not.
            bases = ["N"]

            if positionData:
                nuclAbundance, position = positionData
                base_depth = sum(nuclAbundance)
                dominance = float(max(nuclAbundance)) / float(base_depth)

                if base_depth > mincov:
                    dominanceList.append(dominance)
                    # TODO: should threshold be >20%
                    bases = [k for k, v in sorted(zip(["A", "C", "G", "T"], nuclAbundance),
                                                  key=lambda x: x[1], reverse=True) if v > 0]
                else:
                    dominanceList.append(np.nan)
            else:
                dominanceList.append(np.nan)

            first_base = bases[0]
            second_base = bases[1] if (len(bases) > 1 and dominance < dominant_frq_thrsh) else bases[0]

            codon_f1.append(first_base)
            codon_f2.append(second_base)

            # BUG: what if some bases are lost?
            if len(codon_f1) == 3 and len(codon_f2) == 3:

                codon_s1 = Seq.Seq("".join(codon_f1))
                codon_s2 = Seq.Seq("".join(codon_f2))
                codon_t1 = codon_s1.translate()
                codon_t2 = codon_s2.translate()

                #positionLabel = positionData[1] if positionData else "ND"
                RD = None
                if codon_t1 != codon_t2:
                    RD = "DN"
                elif codon_t1 == codon_t2 and (codon_s1 != codon_s2):
                    RD = "DS"
                else:  #if codon_t1 == "X" or codon_t2 == "X":
                    RD = "D?"
                mutationStats[RD] += 1

                codon_f1 = []
                codon_f2 = []

        return (dominanceList, mutationStats)

    def polymorphism_rate(self, mincov=CMSEQ_DEFAULTS.mincov, minqual=CMSEQ_DEFAULTS.minqual,
                          pvalue=CMSEQ_DEFAULTS.poly_pvalue_threshold,
                          error_rate=CMSEQ_DEFAULTS.poly_error_rate,
                          dominant_frq_thrsh=CMSEQ_DEFAULTS.poly_dominant_frq_thrsh):

        base_values = self.get_base_stats(min_read_depth=mincov, min_base_quality=minqual,
                                          error_rate=error_rate, dominant_frq_thrsh=dominant_frq_thrsh)

        rv = {}
        rv["total_covered_bases"] = len(base_values)
        rv["total_polymorphic_bases"] = 0

        if len(base_values) > 0:
            pb = sum((1 if (info["p"] < pvalue and info["ratio_max2all"] < dominant_frq_thrsh) else 0)
                     for pox, info in base_values.items())

            rv["total_polymorphic_bases"] = pb
            rv["total_polymorphic_rate"] = float(pb) / float(len(base_values))

            # If we have at least one polymorphic site
            if pb > 0:

                rv["ratios"] = [info["ratio_max2all"] for pox, info in base_values.items() if (info["p"] < pvalue and info["ratio_max2all"] < dominant_frq_thrsh)]
                rv["dominant_allele_distr_mean"] = np.mean(rv["ratios"])
                rv["dominant_allele_distr_sd"] = np.std(rv["ratios"])

                for i in [10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 98, 99]:
                    rv["dominant_allele_distr_perc_" + str(i)] = np.percentile(rv["ratios"], i)

        return rv

    def breadth_and_depth_of_coverage(self, mincov=10, minqual=30, trunc=0):
        coverage_positions = {}
        if self.length > trunc * 2:
            # Check if the contig is long enough to be truncated
            consid_r = range(int(trunc), int(self.length - trunc))
        else:
            # If a contig is too short to be truncated, ignore the truncation.
            # This is not nice and should be improved.
            consid_r = range(0, int(self.length))

        for pileupcolumn in self.bam_handle.pileup(self.name, stepper=self.stepper):
            #for each position
            if pileupcolumn.pos in consid_r:
                tCoverage = 0
                for pileupread in pileupcolumn.pileups:
                    #for each base at the position
                    if (not pileupread.is_del and
                            not pileupread.is_refskip and
                            pileupread.alignment.query_qualities[pileupread.query_position] >= minqual and
                            pileupread.alignment.query_sequence[pileupread.query_position].upper() in ("A", "T", "C", "G")):
                        tCoverage += 1

                    if tCoverage >= mincov:
                        coverage_positions[pileupcolumn.pos] = tCoverage

        if (len(coverage_positions.keys())) > 0:
            breadth = float(len(coverage_positions.keys())) / len(consid_r)
            vals = list(coverage_positions.values())
            avgdepth = np.mean(vals)
            mediandepth = np.median(vals)

            return (breadth, avgdepth, mediandepth, coverage_positions.values())
        else:
            return (np.nan, np.nan, np.nan, [np.nan])

    def depth_of_coverage(self, mincov=10, minqual=30):
        return self.breadth_and_depth_of_coverage(mincov, minqual)[1]
        #coverage_positions = {}
        #for pileupcolumn in self.bam_handle.pileup(self.name,stepper=self.stepper):
        #    if pileupcolumn.n >= mincov: coverage_positions[pileupcolumn.pos] = len([1 for pileupread in pileupcolumn.pileups if not pileupread.is_del and not pileupread.is_refskip and pileupread.alignment.query_qualities[pileupread.query_position] >= args.minqual and pileupread.alignment.query_sequence[pileupread.query_position].upper() in ("A","T","C","G") ])
        #
        #return (np.mean(coverage_positions.values()),np.median(coverage_positions.values()))

    def breadth_of_coverage(self, mincov=10, minqual=30):
        return self.breadth_and_depth_of_coverage(mincov, minqual)[0]

#------------------------------------------------------------------------------

    def get_base_stats(self, min_read_depth=CMSEQ_DEFAULTS.mincov, min_base_quality=CMSEQ_DEFAULTS.minqual,
                       error_rate=CMSEQ_DEFAULTS.poly_error_rate,
                       dominant_frq_thrsh=CMSEQ_DEFAULTS.poly_dominant_frq_thrsh,
                       BAM_tagFilter=None, trimReads=None):
        """
        get base frequencies and quality stats,
        to use in get_all_base_values() and other functions
        """
        if trimReads:
            mask_head_until = int(trimReads[0]) if (trimReads[0] is not None and trimReads[0] != "") else 0
            mask_tail_before = int(trimReads[1]) if (trimReads[1] is not None and trimReads[1] != "") else 0

        base_stats = defaultdict(dict)

        ATCG = ("A", "T", "C", "G")

        for base_pileup in self.bam_handle.pileup(self.name, stepper=self.stepper, min_base_quality=min_base_quality):
            base_freq = {"A": 0, "T": 0, "C": 0, "G": 0, "N": 0}

            pos = base_pileup.pos + 1  # 1-based

            #for each read composing the pile
            for matched_read in base_pileup.pileups:
                if not matched_read.is_del and not matched_read.is_refskip:
                    b = matched_read.alignment.query_sequence[matched_read.query_position].upper()
                    #q = matched_read.alignment.query_qualities[matched_read.query_position]
                    #print("I am get_base_stats and this is read", matched_read.alignment.query_name, " (L=",matched_read.alignment.query_length," ) at position", matched_read.query_position, "it"s a ", b)

                    thisPositionBase = "N"

                    if not trimReads or (trimReads and ((matched_read.query_position >= mask_head_until) and
                                         (matched_read.query_position <= (matched_read.alignment.query_length - mask_tail_before)))):
                        if b in ATCG:
                            if BAM_tagFilter is None or all(globals()[func](matched_read.alignment.get_tag(tag), limitValue) is True
                                                            for (tag, func, limitValue) in BAM_tagFilter):
                                thisPositionBase = b

                    base_freq[thisPositionBase] += 1

            # calculate quality stats, ignoring N"s
            base_sum = sum([base_freq[b] for b in ATCG])
            base_max = float(max([base_freq[b] for b in ATCG]))

            if base_sum >= min_read_depth:
                r = base_max / base_sum
                #print r, dominant_frq_thrsh
                if r < dominant_frq_thrsh:
                    #it makes sense to calculate pvalue
                    p = stats.binom.cdf(base_max, base_sum, 1.0 - error_rate)
                else:
                    p = 1.0

                base_stats[pos]["p"] = p                  # quality measure
                base_stats[pos]["ratio_max2all"] = r      # dominant base versus others
                base_stats[pos]["base_cov"] = base_sum    # number of reads covering the base, not counting N"s
                base_stats[pos]["base_freq"] = base_freq  # dict: {"A":4,"T":1,"C":2,"G":0,"N":0}

        return base_stats

    def get_all_base_values(self, stats_value,  *f_args, **f_kwargs):
        """
        get list of p values (or "ratio_max2all" etc) for all bases that pass argument thresholds
        p_all = a.get_contig_by_label("CONTIGNAME").get_all_base_values("p", min_base_quality=30)
        """
        base_stats = self.get_base_stats(*f_args, **f_kwargs)
        return [base_stats[k].get(stats_value, "NaN") for k in base_stats]


def loc_gte(a, b):
    return a >= b


def loc_lte(a, b):
    return a <= b


def loc_gt(a, b):
    return a > b


def loc_lt(a, b):
    return a < b


def loc_leq(a, b):
    return a == b


class bcolors:
    HEADER = "\033[95m"
    OKBLUE = "\033[94m"
    OKGREEN = "\033[92m"
    WARNING = "\033[93m"
    FAIL = "\033[91m"
    ENDC = "\033[0m"
    OKGREEN2 = "\033[42m\033[30m"
    RED = "\033[1;91m"
    CYAN = "\033[0;37m"
