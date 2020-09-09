import pandas as pd
import numpy as np
import argparse
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .cmseq import CMSEQ_DEFAULTS
from .cmseq import BamFile

def consensus_from_file():
	parser = argparse.ArgumentParser(description="outputs the consensus in FASTA format. Non covered positions (or quality-trimmed positions) are reported as a dashes: -")
	parser.add_argument('BAMFILE', help='The file on which to operate')
	parser.add_argument('-c','--contig', help='Focus on a subset of references in the BAM file. Can be a list of references separated by commas or a FASTA file (the IDs are used to subset)', metavar="REFERENCE ID" ,default=None)
	parser.add_argument('-f', help='If set unmapped (FUNMAP), secondary (FSECONDARY), qc-fail (FQCFAIL) and duplicate (FDUP) are excluded. If unset ALL reads are considered (bedtools genomecov style). Default: unset',action='store_true')
	parser.add_argument('--sortindex', help='Sort and index the file',action='store_true')
	parser.add_argument('--minqual', help='Minimum base quality. Bases with quality score lower than this will be discarded. This is performed BEFORE --mincov. Default: '+str(CMSEQ_DEFAULTS.minqual), type=int, default=CMSEQ_DEFAULTS.minqual)
	parser.add_argument('--mincov', help='Minimum position coverage to perform the polymorphism calculation. Position with a lower depth of coverage will be discarded (i.e. considered as zero-coverage positions). This is calculated AFTER --minqual. Default: '+str(CMSEQ_DEFAULTS.minlen), type=int, default=CMSEQ_DEFAULTS.mincov)
	parser.add_argument('--dominant_frq_thrsh', help='Cutoff for degree of `allele dominance` for a position to be considered polymorphic. Default: '+str(CMSEQ_DEFAULTS.poly_dominant_frq_thrsh), type=float, default=CMSEQ_DEFAULTS.poly_dominant_frq_thrsh)
	parser.add_argument('--minlen', help='Minimum Reference Length for a reference to be considered. Default: '+str(CMSEQ_DEFAULTS.minlen),default=CMSEQ_DEFAULTS.minlen, type=int)
	parser.add_argument('--trim', help='Trim the reads before computing the consensus. A value of 10:10 means that the first and last 10 positions of each read will be ignored. Default: '+str(CMSEQ_DEFAULTS.trimReads),default=CMSEQ_DEFAULTS.trimReads, type=str)
	
	args = parser.parse_args()
	si = True if args.sortindex else False
	mode = 'all' if args.f else 'nofilter'

	bf = BamFile(args.BAMFILE,sort=si,index=si,stepper=mode,minlen=args.minlen,filterInputList=args.contig)
	#tl = [bf.get_contig_by_label(contig) for contig in args.contig.split(',')] if args.contig is not None else list(bf.get_contigs_obj())
 
	lst = []

	for i in bf.get_contigs_obj():

		
		trimParam = tuple(args.trim.strip().split(':')) if args.trim else args.trim

		sq = i.reference_free_consensus(mincov=args.mincov,minqual=args.minqual,dominant_frq_thrsh=args.dominant_frq_thrsh,noneCharacter='N',trimReads=trimParam)
		
		if sq is not None:
			lst.append(SeqRecord(Seq(sq), id=i.name+"_consensus", description=''))
	SeqIO.write(lst,sys.stdout,'fasta')


if __name__ == "__main__":
	consensus_from_file()