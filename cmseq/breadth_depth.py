from __future__ import print_function

from .cmseq import CMSEQ_DEFAULTS
from .cmseq import BamFile

import pandas as pd
import numpy as np
import argparse



def bd_from_file():
	parser = argparse.ArgumentParser(description="calculate the Breadth and Depth of coverage of BAMFILE.")

	parser.add_argument('BAMFILE', help='The file on which to operate')
	parser.add_argument('-c','--contig', help='Focus on a subset of references in the BAM file. Can be a list of references separated by commas or a FASTA file (the IDs are used to subset)', metavar="REFERENCE ID" ,default=None)
	parser.add_argument('-f', help='If set unmapped (FUNMAP), secondary (FSECONDARY), qc-fail (FQCFAIL) and duplicate (FDUP) are excluded. If unset ALL reads are considered (bedtools genomecov style). Default: unset',action='store_true')
	parser.add_argument('--sortindex', help='Sort and index the file',action='store_true')
	parser.add_argument('--minlen', help='Minimum Reference Length for a reference to be considered',default=CMSEQ_DEFAULTS.minlen, type=int)
	parser.add_argument('--minqual', help='Minimum base quality. Bases with quality score lower than this will be discarded. This is performed BEFORE --mincov. Default: 30', type=int, default=CMSEQ_DEFAULTS.minqual)
	parser.add_argument('--mincov', help='Minimum position coverage to perform the polymorphism calculation. Position with a lower depth of coverage will be discarded (i.e. considered as zero-coverage positions). This is calculated AFTER --minqual. Default: 1', type=int, default=CMSEQ_DEFAULTS.mincov)
	parser.add_argument('--truncate', help='Number of nucleotides that are truncated at either contigs end before calculating coverage values.', type=float, default=0)
	parser.add_argument('--combine', help='Combine all contigs into one giant contig and report it at the end', action='store_true')

	#print vars(args)
	args = parser.parse_args()
	si = True if args.sortindex else False
	mode = 'all' if args.f else 'nofilter'

	bf = BamFile(args.BAMFILE,sort=si,index=si,stepper=mode,minlen=args.minlen,filterInputList=args.contig,minimumReadsAligning=args.mincov)

	print('Contig\tBreadth\tDepth avg\tDepth median')

	all_coverage_values = []
	for i in bf.get_contigs_obj():
		bd_result = i.breadth_and_depth_of_coverage(minqual=args.minqual,mincov=args.mincov,trunc=args.truncate)

		if not all(np.isnan(x) for x in [bd_result[0],bd_result[1],bd_result[2]]):
			print (i.name+'\t'+str(bd_result[0])+'\t'+str(bd_result[1])+'\t'+str(bd_result[2]))
			
			if args.combine:
				all_coverage_values.extend(bd_result[3])

	if args.combine:
		if np.all(np.isnan(all_coverage_values)):
			print ("all_contigs"+'\t-\t'+str("NaN")+'\t'+str("NaN"))
		else:
			print ("all_contigs"+'\t-\t'+str(np.nanmean(all_coverage_values)) + '\t'+str(np.nanmedian(all_coverage_values)))


if __name__ == "__main__":
	bd_from_file()