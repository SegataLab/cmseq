#!/usr/bin/env python 

import os,pysam
import numpy as np
import math
import sys


class BamFile:
	bam_handle = None
	contigs = {}

	def __init__(self,bamFile,sort=False,index=False,stepper='nofilter',minlen=0):
		if not os.path.isfile(bamFile):
			raise Exception(bamFile+' is not accessible, or is not a file')

		if sort:
			import subprocess
			fp = bamFile+'.sorted'
			subprocess.call(['samtools','sort',bamFile,'-o',bamFile+'.sorted'])
		else: fp = bamFile

		if index: pysam.index(fp)

		bamHandle = pysam.AlignmentFile(fp, "rb")
		self.bam_handle = bamHandle
		self.contigs = dict((r,BamContig(self.bam_handle,r,l,stepper)) for r,l in zip(bamHandle.references,bamHandle.lengths) if l > minlen)

	def get_contigs(self): return iter(self.contigs.keys())
	def get_contigs_obj(self): return iter(self.contigs.values())
	def get_contig_by_label(self,contigID): return (self.contigs[contigID] if contigID in self.contigs else None)


class BamContig:

	coverage = None
	consensus = ''
	name = None
	length = None
	stepper = 'nofilter'

	def __init__(self,bamHandle,contigName,contigLength,stepper='nofilter'):

		self.name = contigName
		self.length = contigLength
		self.bam_handle = bamHandle 
		self.stepper=stepper


	def set_stepper(self,ns):
		if ns in ['all','nofilter']: self.stepper = ns

	def reference_free_consensus(self,consensus_rule=lambda array: max(array, key=array.get),mincov=1,minqual=30,fast=False):

		if (fast) : return self.fast_reference_free_consensus(consensus_rule)

		consensus_positions = {}

		for pileupcolumn,position_data in self.get_base_stats(min_read_depth=mincov, min_base_quality=minqual, error_rate=0.01).items():
			consensus_positions[pileupcolumn] = consensus_rule(position_data['base_freq'])
			#print 'GAMMA'+str((pileupcolumn)-1)+'_'+repr(position_data['base_freq'])+'_'+repr(consensus_positions[pileupcolumn])

		if len(consensus_positions) > 0 :
			self.consensus = ''.join([(consensus_positions[position] if position in consensus_positions else '-') for position in range(1,self.length+1)])
		else:
			self.consensus = None

		del consensus_positions
		return self.consensus


	def fast_reference_free_consensus(self,consensus_rule=lambda array: max(array, key=array.get)):

		consensus_positions = {}

		for pileupcolumn in self.bam_handle.pileup(self.name,stepper=self.stepper):
			consensus_positions[pileupcolumn.pos] = {'A':0,'T':0,'C':0,'G':0,'N':0}
			for pileupread in pileupcolumn.pileups:
				if not pileupread.is_del and not pileupread.is_refskip:
					base = pileupread.alignment.query_sequence[pileupread.query_position].upper()
					if base in ['A','T','C','G']: consensus_positions[pileupcolumn.pos][base]+=1
					else: consensus_positions[pileupcolumn.pos]['N']+=1

			consensus_positions[pileupcolumn.pos] = consensus_rule(consensus_positions[pileupcolumn.pos])

		if len(consensus_positions) > 0 : 
			self.consensus = ''.join([(consensus_positions[position] if position in consensus_positions else 'N') for position in range(0,self.length)])
		else:
			self.consensus = None

		del consensus_positions
		return self.consensus

	def polymorphism_rate(self,mincov=1,minqual=30,pvalue=0.05,precomputedBinomial=None):

		base_values = self.get_base_stats(min_read_depth=mincov, min_base_quality=minqual, error_rate=math.pow(10, (-minqual/10)),precomputedBinomial=precomputedBinomial)

		rv={}
		rv['covered_bases'] = len(base_values)
		rv['polymorphic_bases'] = 0

		if len(base_values) > 0:
			pb=sum([(1 if info['p'] < pvalue else 0) for pox,info in base_values.items()])

			rv['polymorphic_bases']= pb
			rv['polymorphic_rate'] = float(pb)/float(len(base_values))

			# If we have at least one polymorphic site
			if pb > 0 :
				rv['ratios'] = [info['ratio_max2all'] for pox,info in base_values.items() if info['p'] < pvalue]
				rv['dominant_allele_distr_mean'] = np.mean(rv['ratios'])
				rv['dominant_allele_distr_sd'] = np.std(rv['ratios'])

				for i in [10,20,30,40,50,60,70,80,90,100]:
					rv['dominant_allele_distr_perc_'+str(i)] = np.percentile(rv['ratios'],i)

		return rv


	def breadth_and_depth_of_coverage(self,mincov=1,minqual=30):
		coverage_positions = {}
		ptl=0

		for pileupcolumn in self.bam_handle.pileup(self.name,stepper=self.stepper):
			#for each position
			tCoverage = 0
			for pileupread in pileupcolumn.pileups:
				#for each base at the position

				if not pileupread.is_del and not pileupread.is_refskip and pileupread.alignment.query_qualities[pileupread.query_position] >= minqual and pileupread.alignment.query_sequence[pileupread.query_position].upper() in ('A','T','C','G'):
					tCoverage +=1


				if tCoverage >= mincov:
					coverage_positions[pileupcolumn.pos] = tCoverage

		if (len(coverage_positions.keys())) > 0:
			breadth = float(len(coverage_positions.keys()))/self.length
			avgdepth = np.mean(coverage_positions.values())
			mediandepth = np.median(coverage_positions.values())

			return (breadth,avgdepth,mediandepth)

		else: return (np.nan,np.nan,np.nan)


	def depth_of_coverage(self,mincov=1,minqual=30):
		return self.breadth_and_depth_of_coverage(mincov,minqual)[1]
		#coverage_positions = {}
		#for pileupcolumn in self.bam_handle.pileup(self.name,stepper=self.stepper):
		#	if pileupcolumn.n >= mincov: coverage_positions[pileupcolumn.pos] = len([1 for pileupread in pileupcolumn.pileups if not pileupread.is_del and not pileupread.is_refskip and pileupread.alignment.query_qualities[pileupread.query_position] >= args.minqual and pileupread.alignment.query_sequence[pileupread.query_position].upper() in ('A','T','C','G') ])
		#
		#return (np.mean(coverage_positions.values()),np.median(coverage_positions.values()))


	def breadth_of_coverage(self,mincov=1,minqual=30):
		return self.breadth_and_depth_of_coverage(mincov,minqual)[0]
		#coverage_positions = {}
		#ptl=0
		#
		#if minqual == 0:
		#	for pileupcolumn in self.bam_handle.pileup(self.name,stepper=self.stepper):
		#		if len([1 for pileupread in pileupcolumn.pileups if not pileupread.is_del and not pileupread.is_refskip and pileupread.alignment.query_sequence[pileupread.query_position].upper() in ('A','T','C','G') ]) >= mincov:
		#			ptl += 1

		#else:
		#	for pileupcolumn in self.bam_handle.pileup(self.name,stepper=self.stepper):
		#		#print minqual,len([1 for pileupread in pileupcolumn.pileups if not pileupread.is_del and not pileupread.is_refskip and pileupread.alignment.query_qualities[pileupread.query_position] >= args.minqual and pileupread.alignment.query_sequence[pileupread.query_position].upper() in ('A','T','C','G') ]), len([1 for pileupread in pileupcolumn.pileups if not pileupread.is_del and not pileupread.is_refskip and pileupread.alignment.query_sequence[pileupread.query_position].upper() in ('A','T','C','G') ])
		#		if len([1 for pileupread in pileupcolumn.pileups if not pileupread.is_del and not pileupread.is_refskip and pileupread.alignment.query_qualities[pileupread.query_position] >= args.minqual and pileupread.alignment.query_sequence[pileupread.query_position].upper() in ('A','T','C','G') ]) >= mincov:
		#			ptl += 1
		#
		#return float(ptl)/self.length

	def plot_coverage(self,flavour='polar',path='./out.pdf',smooth=0,l_avoid=False,s_avoid=False,l_color='#000000',s_color='#000000'):
		
		import matplotlib
		matplotlib.use('Agg')
		import matplotlib.pyplot as plt 
		import numpy as np
		import seaborn as sns
		
		iret={}
		for column in self.bam_handle.pileup(self.name,stepper=self.stepper):
			iret[column.pos] = column.n

		for i in range(0,self.length):
			if i not in iret: iret[i] = 0

		
		if smooth == 0:
			y_smooth = iret.values()
		else: 
			box_pts = smooth
			y_smooth = np.convolve(iret.values(), (np.ones(box_pts)/box_pts), mode='same')

		if flavour == 'polar':
			ax = plt.subplot(111, projection='polar')
			if not l_avoid: ax.plot([float(k)*2*np.pi/float(self.length) for k in iret.keys()],y_smooth,label='Coverage',linewidth=0.5,c=l_color,zorder=1)
			if not s_avoid: ax.scatter([float(k)*2*np.pi/float(self.length) for k in iret.keys()],y_smooth,label='Coverage',s=1,c=s_color, zorder=2)
			ax.grid(True)
			ax.set_xticklabels(['0', str(int(self.length/8)), str(int(self.length/4)), str(int(self.length*3/8)), str(int(self.length/2)), str(int(self.length*5/8)), str(int(self.length/4*3)), str(int(self.length*7/8))])
		else:

			if not l_avoid: plt.plot(iret.keys(),y_smooth,label='Coverage',linewidth=0.5,c=l_color,zorder=1)
			if not s_avoid: plt.scatter(iret.keys(),y_smooth,label='Coverage',s=1,c=s_color, zorder=2)
			plt.grid(True)
			plt.axis('tight')
			

		plt.savefig(path)
		plt.clf()
		plt.close()
		
#------------------------------------------------------------------------------	
	def get_base_stats(self, min_read_depth=1, min_base_quality=30, error_rate=0.01,precomputedBinomial=None):
		'''
		get base frequencies and quality stats,
		to use in get_all_base_values() and other functions
		'''
		from scipy import stats
		from collections import defaultdict
		import pickle,os


		base_stats = defaultdict(dict)



		ATCG=('A','T','C','G')
		for base_pileup in self.bam_handle.pileup(self.name,stepper=self.stepper):
			base_freq = {'A':0,'T':0,'C':0,'G':0,'N':0}
			for matched_read in base_pileup.pileups:
				
				
				if not matched_read.is_del and not matched_read.is_refskip:
					b = matched_read.alignment.query_sequence[matched_read.query_position].upper()
					q = matched_read.alignment.query_qualities[matched_read.query_position]	

					# if b=='N': print(matched_read.alignment.query_sequence)
					if q >= min_base_quality:
						if b in ATCG:
							base_freq[b] += 1
						else:
							base_freq['N']+=1
			# calculate quality stats, ignoring N's 

			base_sum=sum([base_freq[b] for b in ATCG]) 
			base_max=float(max([base_freq[b] for b in ATCG]))
			if base_sum >= min_read_depth:
				r = base_max / base_sum
				

				if precomputedBinomial and os.path.isfile(precomputedBinomial) and base_max in precomputedBinomial and base_sum in precomputedBinomial[base_max]:
					p = precomputedBinomial[base_max][base_sum]
				else:
					p = stats.binom.cdf(base_max, base_sum, 1.0 - error_rate)

				pos=base_pileup.pos+1 # 1-based
				base_stats[pos]['p']=p                 # quality measure
				base_stats[pos]['ratio_max2all']=r     # dominant base versus others
				base_stats[pos]['base_cov'] =base_sum  # number of reads covering the base, not counting N's
				base_stats[pos]['base_freq']=base_freq # dict: {'A':4,'T':1,'C':2,'G':0,'N':0}
				# base_stats[pos]['ref_base']='?' # in case of reference sequence
		return base_stats
	
	def get_all_base_values(self, stats_value,  *f_args, **f_kwargs):
		'''
		get list of p values (or 'ratio_max2all' etc) for all bases that pass argument thresholds
		p_all = a.get_contig_by_label('CONTIGNAME').get_all_base_values('p', min_base_quality=30)
		'''
		base_stats = self.get_base_stats(*f_args, **f_kwargs)
		return [ base_stats[k].get(stats_value, 'NaN') for k in base_stats]
#------------------------------------------------------------------------------


def get_contig_list(inputList):
	'''
	get a contig list from file (if argument is a file, or list of comma-separated-elements) 
	'''
	import sys,os
	from Bio import SeqIO

	toList=[]

	if os.path.isfile(inputList):
		with open(inputList, "r") as infile:
			for record in SeqIO.parse(infile, "fasta"):
				toList.append(record.id)
	else:
		toList = inputList.split(',')
	
	return toList

if __name__ == "__main__":

	def poly_from_file(args):

		import pandas as pd
		if args.precomputed and os.path.isfile(args.precomputed):
			import pickle
			binomPrecomputed = pickle.load( open( args.precomputed, "rb" ))
		else: binomPrecomputed=None

		si = True if args.sortindex else False
		mode = 'all' if args.f else 'nofilter'
		bf = BamFile(args.BAMFILE,sort=si,index=si,stepper=mode,minlen=args.minlen)

		
		get_contig_list(args.contig)

		outputDF = []
		allRatios = []
		allGenomeCol = {'referenceID': '-GENOME-','covered_bases':0,'polymorphic_bases':0,'polymorphic_rate':np.nan}
		for element in [bf.get_contig_by_label(contig) for contig in get_contig_list(args.contig)] if args.contig is not None else list(bf.get_contigs_obj()):
			tld = element.polymorphism_rate(minqual=args.minqual,mincov=args.mincov,precomputedBinomial=binomPrecomputed)
			tld['referenceID'] = element.name
		
			allGenomeCol['covered_bases'] += tld['covered_bases']
			allGenomeCol['polymorphic_bases'] += tld['polymorphic_bases'] 
			if 'ratios' in tld:
				allRatios = allRatios + tld['ratios']
				del tld['ratios']

			outputDF.append(tld)
			del tld

		if float(allGenomeCol['covered_bases']) > 0: 

			allGenomeCol['polymorphic_rate'] = float(allGenomeCol['polymorphic_bases']) / float(allGenomeCol['covered_bases'])

			allGenomeCol['dominant_allele_distr_mean'] = np.mean(allRatios)
			allGenomeCol['dominant_allele_distr_sd'] = np.std(allRatios)
			for i in [10,20,30,40,50,60,70,80,90,100]:
				allGenomeCol['dominant_allele_distr_perc_'+str(i)] = np.percentile(allRatios,i)

		outputDF.append(allGenomeCol)

		pd.DataFrame.from_dict(outputDF).set_index('referenceID').to_csv(sys.stdout,sep='\t')

		#print np.mean(k['polymorphic_rate'])




	def bd_from_file(args):

		si = True if args.sortindex else False
		mode = 'all' if args.f else 'nofilter'

		bf = BamFile(args.BAMFILE,sort=si,index=si,stepper=mode,minlen=args.minlen)

		tl = [bf.get_contig_by_label(contig) for contig in args.contig.split(',')] if args.contig is not None else list(bf.get_contigs_obj())

		print 'Contig\tBreadth\tDepth (avg)\tDepth (median)'

		for i in tl:
			bd_result = i.breadth_and_depth_of_coverage(minqual=args.minqual,mincov=args.mincov)
			print i.name+'\t'+str(bd_result[0])+'\t'+str(bd_result[1])+'\t'+str(bd_result[2])
		

	def consensus_from_file(args):
		from Bio import SeqIO
		from Bio.Seq import Seq
		from Bio.SeqRecord import SeqRecord
		from Bio.Alphabet import IUPAC
		import sys

		si = True if args.sortindex else False
		mode = 'all' if args.f else 'nofilter'



		bf = BamFile(args.BAMFILE,sort=si,index=si,stepper=mode,minlen=args.minlen)
		tl = [bf.get_contig_by_label(contig) for contig in args.contig.split(',')] if args.contig is not None else list(bf.get_contigs_obj())

		lst = []

		for i in tl:
			sq = i.reference_free_consensus(mincov=args.mincov,minqual=args.minqual,fast=(True if args.fast else False))
			if sq is not None:
				lst.append(SeqRecord(Seq(sq, IUPAC.IUPACAmbiguousDNA), id=i.name+"_consensus", description=''))
		SeqIO.write(lst,sys.stdout,'fasta')
		
				


	def plot_coverage_from_file(args):

		si = True if args.sortindex else False
		mode = 'all' if args.f else 'nofilter'
		bf = BamFile(args.BAMFILE,sort=si,index=si,stepper=mode,minlen=args.minlen)

		tl = [bf.get_contig_by_label(contig) for contig in args.contig.split(',')] if args.contig is not None else list(bf.get_contigs_obj())

		for i in tl:
			i.plot_coverage(path='./'+i.name+'.'+args.format,smooth=args.smooth,l_avoid=args.l_avoid,s_avoid=args.s_avoid,l_color=args.l_color,s_color=args.s_color,flavour=args.flavour)


	import argparse
	parser = argparse.ArgumentParser()


	

	subparsers = parser.add_subparsers(title='subcommands',description='valid subcommands')
	parser_breadth = subparsers.add_parser('bd',description="calculate the Breadth and Depth of coverage of BAMFILE.")
	parser_breadth.add_argument('BAMFILE', help='The file on which to operate')
	parser_breadth.add_argument('-c','--contig', help='Gets the breadth and depth of a specific reference within a BAM Can be a string or a list of strings separated by comma.', metavar="REFERENCE ID" ,default=None)
	parser_breadth.add_argument('-f', help='If set unmapped (FUNMAP), secondary (FSECONDARY), qc-fail (FQCFAIL) and duplicate (FDUP) are excluded. If unset ALL reads are considered (bedtools genomecov style). Default: unset',action='store_true')
	parser_breadth.add_argument('--sortindex', help='Sort and index the file',action='store_true')
	parser_breadth.add_argument('--minlen', help='Minimum Reference Length for a reference to be considered',default=0, type=int)
	parser_breadth.add_argument('--minqual', help='Minimum base quality. Bases with quality score lower than this will be discarded. This is performed BEFORE --mincov. Default: 30', type=int, default=30)
	parser_breadth.add_argument('--mincov', help='Minimum position coverage to perform the polymorphism calculation. Position with a lower depth of coverage will be discarded (i.e. considered as zero-coverage positions). This is calculated AFTER --minqual. Default: 1', type=int, default=1)
	
	parser_breadth.set_defaults(func=bd_from_file)

	parser_poly = subparsers.add_parser('poly',description="Reports the polymorpgic rate of each reference (polymorphic bases / total bases). Focuses only on covered regions (i.e. depth >= 1)")
	parser_poly.add_argument('BAMFILE', help='The file on which to operate')
	parser_poly.add_argument('-c','--contig', help='Gets the polymorphic rate _of a specific reference within a BAM Can be a string or a list of strings separated by comma.', metavar="REFERENCE ID" ,default=None)
	parser_poly.add_argument('-f', help='If set unmapped (FUNMAP), secondary (FSECONDARY), qc-fail (FQCFAIL) and duplicate (FDUP) are excluded. If unset ALL reads are considered (bedtools genomecov style). Default: unset',action='store_true')
	parser_poly.add_argument('--sortindex', help='Sort and index the file',action='store_true')
	
	parser_poly.add_argument('--minlen', help='Minimum Reference Length for a reference to be considered',default=0, type=int)
	parser_poly.add_argument('--minqual', help='Minimum base quality. Bases with quality score lower than this will be discarded. This is performed BEFORE --mincov. Default: 30', type=int, default=30)
	parser_poly.add_argument('--mincov', help='Minimum position coverage to perform the polymorphism calculation. Position with a lower depth of coverage will be discarded (i.e. considered as zero-coverage positions). This is calculated AFTER --minqual. Default: 1', type=int, default=1)
	parser_poly.add_argument('--precomputed', help='Path to pickled dictionary for precomputed probabilities of BINOM function')

 
	
	parser_poly.set_defaults(func=poly_from_file)

	parser_consensus = subparsers.add_parser('consensus',description="outputs the consensus in FASTA format. Non covered positions (or quality-trimmed positions) are reported as a dashes: -")
	parser_consensus.add_argument('BAMFILE', help='The file on which to operate')
	parser_consensus.add_argument('-c','--contig', help='Gets the breadth and depth of a specific reference within a BAM Can be a string or a list of strings separated by comma.',default=None)
	parser_consensus.add_argument('-f', help='If set unmapped (FUNMAP), secondary (FSECONDARY), qc-fail (FQCFAIL) and duplicate (FDUP) are excluded. If unset ALL reads are considered (bedtools genomecov style). Default: unset',action='store_true')
	parser_consensus.add_argument('--sortindex', help='Sort and index the file',action='store_true')
	parser_consensus.add_argument('--minqual', help='Minimum base quality. Bases with quality score lower than this will be discarded. This is performed BEFORE --mincov. Default: 0', type=int, default=0)
	parser_consensus.add_argument('--mincov', help='Minimum position coverage to perform the polymorphism calculation. Position with a lower depth of coverage will be discarded (i.e. considered as zero-coverage positions). This is calculated AFTER --minqual. Default: 1', type=int, default=1)
	parser_consensus.add_argument('--minlen', help='Minimum Reference Length for a reference to be considered',default=0, type=int)
	parser_consensus.add_argument('--fast', help='Fast version (does not account for polymorphism)',action='store_true')


	parser_consensus.set_defaults(func=consensus_from_file)


	parser_coverageplot = subparsers.add_parser('coverageplot',description="Plot the coverage of a contig")
	parser_coverageplot.add_argument('BAMFILE', help='The file on which to operate')

	
	parser_coverageplot.add_argument('-c','--contig', help='Get the breadth of a specific contig',default=None)
	parser_coverageplot.add_argument('-f', help='If set unmapped (FUNMAP), secondary (FSECONDARY), qc-fail (FQCFAIL) and duplicate (FDUP) are excluded. If unset ALL reads are considered (bedtools genomecov style). Default: unset',action='store_true')
	parser_coverageplot.add_argument('--flavour', choices=['polar','linear'], help='choose from linear or polar plot', default='polar')
	parser_coverageplot.add_argument('--sortindex', help='Sort and index the file',action='store_true')
	parser_coverageplot.add_argument('--format',help='Output format',choices=['svg','pdf','svg'],default='png')
	parser_coverageplot.add_argument('--smooth', help='Smooth factor. Default: 0 (no smooth)',default=0, type=int)
	parser_coverageplot.add_argument('--l_avoid', help='Suppresses line-plot',action='store_true')
	parser_coverageplot.add_argument('--l_color', help='Line color for matplotlib, in HTML format (#XXXXXX)',default='#000000')
	parser_coverageplot.add_argument('--s_color', help='Scatter color in HTML format (#XXXXXX)',default='#000000')
	parser_coverageplot.add_argument('--s_avoid', help='Suppresses scatter-plot',action='store_true')
	parser_coverageplot.add_argument('--minlen', help='Minimum Reference Length for a reference to be considered',default=0, type=int)
	parser_coverageplot.set_defaults(func=plot_coverage_from_file)


	args = parser.parse_args()
	args.func(args)
