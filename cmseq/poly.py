import pandas as pd
import numpy as np
import argparse
import sys

from .cmseq import CMSEQ_DEFAULTS
from .cmseq import BamFile

def poly_from_file():
	parser = argparse.ArgumentParser(description="Reports the polymorpgic rate of each reference (polymorphic bases / total bases). Focuses only on covered regions (i.e. depth >= 1)")
	parser.add_argument('BAMFILE', help='The file on which to operate')
	parser.add_argument('-c','--contig', help='Focus on a subset of references in the BAM file. Can be a list of references separated by commas or a FASTA file (the IDs are used to subset)', metavar="REFERENCE ID" ,default=None)
	parser.add_argument('-f', help='If set unmapped (FUNMAP), secondary (FSECONDARY), qc-fail (FQCFAIL) and duplicate (FDUP) are excluded. If unset ALL reads are considered (bedtools genomecov style). Default: unset',action='store_true')
	parser.add_argument('--sortindex', help='Sort and index the file',action='store_true')
	parser.add_argument('--minlen', help='Minimum Reference Length for a reference to be considered. Default: '+str(CMSEQ_DEFAULTS.minlen),default=CMSEQ_DEFAULTS.minlen, type=int)
	parser.add_argument('--minqual', help='Minimum base quality. Bases with quality score lower than this will be discarded. This is performed BEFORE --mincov. Default: 30', type=int, default=CMSEQ_DEFAULTS.minqual)
	parser.add_argument('--mincov', help='Minimum position coverage to perform the polymorphism calculation. Position with a lower depth of coverage will be discarded (i.e. considered as zero-coverage positions). This is calculated AFTER --minqual. Default:'+str(CMSEQ_DEFAULTS.mincov), type=int, default=CMSEQ_DEFAULTS.mincov)
	parser.add_argument('--pvalue', help='Binomial p-value threshold for the binomal-polymorphic test. Default: '+str(CMSEQ_DEFAULTS.poly_pvalue_threshold), type=float, default=CMSEQ_DEFAULTS.poly_pvalue_threshold)
	parser.add_argument('--seq_err', help='Sequencing error rate. Default: '+str(CMSEQ_DEFAULTS.poly_error_rate), type=float, default=CMSEQ_DEFAULTS.poly_error_rate)
	parser.add_argument('--dominant_frq_thrsh', help='Cutoff for degree of `allele dominance` for a position to be considered polymorphic. Default: '+str(CMSEQ_DEFAULTS.poly_dominant_frq_thrsh), type=float, default=CMSEQ_DEFAULTS.poly_dominant_frq_thrsh)
	args = parser.parse_args()

	import pandas as pd

	si = True if args.sortindex else False
	mode = 'all' if args.f else 'nofilter'
	
	bf = BamFile(args.BAMFILE,sort=si,index=si,stepper=mode,minlen=args.minlen,filterInputList=args.contig)

	outputDF = []
	allRatios = [] 
	allGenomeCol = {'referenceID': '-GENOME-','total_covered_bases':0,'total_polymorphic_bases':0,'total_polymorphic_rate':np.nan}

	for element in bf.get_contigs_obj():
		
		tld = element.polymorphism_rate(minqual=args.minqual,mincov=args.mincov,error_rate=args.seq_err,dominant_frq_thrsh=args.dominant_frq_thrsh)
		tld['referenceID'] = element.name
	
		allGenomeCol['total_covered_bases'] += tld['total_covered_bases']
		allGenomeCol['total_polymorphic_bases'] += tld['total_polymorphic_bases'] 
		if 'ratios' in tld:
			allRatios = allRatios + tld['ratios']
			del tld['ratios']

		outputDF.append(tld)
		del tld


	if float(allGenomeCol['total_covered_bases']) and float(allGenomeCol['total_polymorphic_bases']) > 0:

		allGenomeCol['total_polymorphic_rate'] = float(allGenomeCol['total_polymorphic_bases']) / float(allGenomeCol['total_covered_bases'])
		allGenomeCol['dominant_allele_distr_mean'] = np.mean(allRatios)
		allGenomeCol['dominant_allele_distr_sd'] = np.std(allRatios)
		
		for i in [10,20,30,40,50,60,70,80,90,95,98,99]:
			allGenomeCol['dominant_allele_distr_perc_'+str(i)] = np.percentile(allRatios,i)


	outputDF.append(allGenomeCol)

	pd.DataFrame.from_dict(outputDF).set_index('referenceID').to_csv(sys.stdout,sep='#')

if __name__ == "__main__":
	poly_from_file()