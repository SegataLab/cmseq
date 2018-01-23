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
		
		self.contigs = dict((r,BamContig(self.bam_handle,r,l,stepper)) for r,l in zip(bamHandle.references,bamHandle.lengths) if l > minlen )

	def get_contigs(self): return iter(self.contigs.keys())
	def get_contigs_obj(self): return iter(self.contigs.values())
	def get_contig_by_label(self,contigID): return (self.contigs[contigID] if contigID in self.contigs else None)

	def parse_gff(self, inputGFF):
		'''
		get a list of contigs plus 0-indexed gene-coordinates and sense-ness of protein coding regions from a gff file.
		Only tested with prokka GFF files.
		'''
		from BCBio import GFF
		import Bio
		import re
		import warnings

		def rev_comp(string):
			string = string.upper()
			complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
			bases = list(string) 
			bases = [complement[base] for base in bases]
			bases.reverse()
			return ''.join(bases)

		in_file = inputGFF
		in_handle = open(in_file)
		#gene_locations = {}
		try:
			parsed_gff = GFF.parse(in_handle)
		except AttributeError:
			print 'Parsing of GFF failed. This is probably because your biopython version is too new. Try downgrading to 1.67.'
		for rec in GFF.parse(in_handle):
			tmp = []
			for r in rec.features:
				if "minced" in r.qualifiers['source'][0] or "Minced" in r.qualifiers['source'][0]:
					# This catches CRISPR repeats.
					continue
				if 'Prodigal' in r.sub_features[0].qualifiers['source'][0] or 'prodigal' in r.sub_features[0].qualifiers['source'][0]:
					# Prokka not only finds protein sequences, but also t-/r-RNA sequences. In order to only parse protein coding sequences,
					# I search for Prodigal/Prodigal in the source entry of the sub_features attribute.
					
					# the sub_features attribute of a seq_record object is apparently deprecated. I couldn't find any other way to access
					# the required information, though. Should probably be fixed when I can.
					indices = str(r.location).split('[')[1].split(']')[0].split(':')
					indices = [int(x) for x in indices]
					sense = str(r.location).split('(')[1].split(')')[0]
					if sense == "-":
						gene_seq = rev_comp(rec.seq[indices[0]:indices[1]])
					else:
						gene_seq = rec.seq[indices[0]:indices[1]]

					if (str(gene_seq[0:3]) == "ATG" or str(gene_seq[0:3]) == "GTG" or str(gene_seq[0:3]) == "TTG"):
						pass
					else:
						warnings.warn(str(r.id) + " doesn't start with a common start codon. Beware. Continuing.")

					if (str(gene_seq[-3:]) == "TAG" or str(gene_seq[-3:]) == "TAA" or str(gene_seq[-3:]) == "TGA"):
						pass
					else:
						warnings.warn(str(r.id) + " doesn't stop with a usual stop codon. Beware. Continuing.")
					tmp.append((indices, sense))
			
			if str(rec.id) in self.contigs:
				self.contigs[str(rec.id)].annotations.append(tmp)
			else:
				warnings.warn(str(rec.id) + " is not tracked by the BAMFile.")
		in_handle.close()
		


class BamContig:

	coverage = None
	consensus = ''
	name = None
	length = None
	stepper = 'nofilter'
	annotations = None

	def __init__(self,bamHandle,contigName,contigLength,stepper='nofilter'):

		self.name = contigName
		self.length = contigLength
		self.bam_handle = bamHandle 
		self.stepper=stepper
		self.annotations = []


	def set_stepper(self,ns):
		if ns in ['all','nofilter']: self.stepper = ns

	def reference_free_consensus(self,consensus_rule=lambda array: max(array, key=array.get),mincov=10,minqual=30,fast=False):

		if (fast) : return self.fast_reference_free_consensus(consensus_rule)

		consensus_positions = {}

		for pileupcolumn,position_data in self.get_base_stats(min_read_depth=mincov, min_base_quality=minqual, error_rate=0.001).items():
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

	def baseline_PSR(self,mincov=10,minqual=30,pvalue=0.01,error_rate=0.001,dominant_frq_thrsh=0.8,binom=None):
		# This function estimates the polymorphic site rate over the input contig assuming that there are no truely polymorphic sites
		# (Meaning that all observed polymorphisms are due to random sequencing error). The test also puts a threshold on the "dominance"
		# of the allele, meaning that it only reports a polymorphic base if the binomial test indicates significance AND the base is NOT sufficiently
		# dominated by the dominant base. Defaults to 0.8 dominance (dominant / all).
		from scipy import stats 

		polymorphic_empirical_loci = 0

		# Get coverage as well as values of contig
		depthsList = self.get_all_base_values('base_cov', min_base_quality=minqual,min_read_depth=mincov)

		# Also get dominant allele frequency of contig
		dominantFreq = self.get_all_base_values('ratio_max2all', min_base_quality=minqual,min_read_depth=mincov)

		# For each position, draw depth-times from bernoulli with success rate 1-error_rate. 
		# Determine significance based on a binomial test, as you would in the regular test for polymorphism.
		for depth, da_freq in zip(depthsList, dominantFreq):
			base_max = sum(stats.bernoulli.rvs(1-error_rate, size=depth))
			if binom and base_max in binom and depth in binom[base_max]:
					p = binom[base_max][depth]
			else:
					p = stats.binom.cdf(base_max, depth,1.0-error_rate)
			if p < pvalue and da_freq < dominant_frq_thrsh:
				polymorphic_empirical_loci+=1
		PSR = float(polymorphic_empirical_loci) / float(len(depthsList))
		return PSR

	def get_base_stats_for_poly(self,minqual=30):
		
		from scipy import stats
		import numpy
		from collections import defaultdict
		import pickle,os
		from itertools import chain
		import sys
		import pandas as pd

		ATCG=('A','C','G','T')
		rev_dict={'A':'T', 'T':'A', 'G':'C', 'C':'G'}
		def rev_pos(cur_pos, gene_start, gene_end):
			# This function mirrors nucleotide positions in a gene on a gene's 'mid-part'
			# (So Nucleotide 1 in a gene is mapped to the last nucleotide in the gene, nucleotide 2 is mapped to one before the last nucleotide on the gene and so forth)
			# The gene_length variable is misnamed. It should be called gene_length_minus_one..
			gene_length = ((gene_end - 1) - gene_start)
			distance_from_start = cur_pos - gene_start
			return(cur_pos + (gene_length - 2 * distance_from_start))

		if not self.annotations:
			base_stats = [None] * self.length
			for base_pileup in self.bam_handle.pileup(self.name,stepper=self.stepper):
				
				base_freq = {'A':0,'C':0,'G':0,'T':0,'N':0}
				for matched_read in base_pileup.pileups:

					if not matched_read.is_del and not matched_read.is_refskip:
						b = matched_read.alignment.query_sequence[matched_read.query_position].upper()
						q = matched_read.alignment.query_qualities[matched_read.query_position]	

						#print self.name,matched_read.query_position, b,q

						if q >= minqual and b in ATCG: base_freq[b] += 1
					 	#else: print "Q",q,"B",b
					#print "Filling",base_pileup.pos,"with", base_freq
					if sum(base_freq.values()) > 0:
						base_stats[base_pileup.pos] = ((base_freq['A'],base_freq['C'],base_freq['G'],base_freq['T']),base_pileup.pos)
		else:
			base_stats = []
			# Generate pileups gene-wise
			# I use the 'truncate' parameter to only obtain the parsed start and stop positions. Without truncate, all positions with reads covering the parsed positions are returned.
			# I wrote a function that reverses a given gene position, which is used to effectively revert genes on the anti-sense strand.
			# Furthermore, for each read's nucleotide over a given position I write out the complement
			genes_and_positions = dict()
			for gene_idx in range(0, len(self.annotations[0])):
				genes_and_positions[gene_idx] = self.annotations[0][gene_idx]
			
			for gene_idx in genes_and_positions:
				gene_stats = [None] * (genes_and_positions[gene_idx][0][1] - genes_and_positions[gene_idx][0][0])
				pos_on_gene = 0
				bam_pileup = self.bam_handle.pileup(self.name, int(genes_and_positions[gene_idx][0][0]), int(genes_and_positions[gene_idx][0][1]), stepper=self.stepper, truncate = True)
				if genes_and_positions[gene_idx][1] == "+":
					# If the gene is on the sense-strand, do the same as before.
					for base_pileup in bam_pileup:
						base_freq = {'A':0,'C':0,'G':0,'T':0,'N':0}
						for matched_read in base_pileup.pileups:
							if not matched_read.is_del and not matched_read.is_refskip:
								b = matched_read.alignment.query_sequence[matched_read.query_position].upper()
								q = matched_read.alignment.query_qualities[matched_read.query_position]	
								if q >= minqual and b in ATCG: base_freq[b] += 1
							if sum(base_freq.values()) > 0:
								gene_stats[pos_on_gene] = ((base_freq['A'],base_freq['C'],base_freq['G'],base_freq['T']), base_pileup.pos)
						pos_on_gene += 1
					base_stats.extend(gene_stats)
				else:
					# If the gene is on the anti-sense strand, effectively return the reverse complement by mapping positions on a gene to it's mirrored position (using rev_pos)
					# and then also converting each nucleotide to it's complement.
					for base_pileup in bam_pileup:
						base_freq = {'A':0,'C':0,'G':0,'T':0,'N':0}
						for matched_read in base_pileup.pileups:
							if not matched_read.is_del and not matched_read.is_refskip:
								b = matched_read.alignment.query_sequence[matched_read.query_position].upper()
								q = matched_read.alignment.query_qualities[matched_read.query_position]	
								# We have to increment the COMPLEMENT of each base when gene calls are on the reverse strand.
								if q >= minqual and b in ATCG: base_freq[rev_dict[b]] += 1
							if sum(base_freq.values()) > 0:
								out_pos = rev_pos(cur_pos = int(pos_on_gene), gene_start = 0, gene_end = len(gene_stats))
								contig_pos = rev_pos(cur_pos = int(base_pileup.pos), gene_start = genes_and_positions[gene_idx][0][0], gene_end = genes_and_positions[gene_idx][0][1])
								gene_stats[out_pos] = ((base_freq['A'],base_freq['C'],base_freq['G'],base_freq['T']), contig_pos)
						pos_on_gene += 1
					if len(gene_stats) % 3 != 0:
						print("One of your genes' length is not a multiple of three. Check your gff file / gene calls. Exiting")
						sys.exit()
					base_stats.extend(gene_stats)

		return base_stats

	def easy_polymorphism_rate(self,mincov=10,minqual=30,dominant_frq_thrsh=0.8):

		from Bio.Seq import Seq
		from Bio.Alphabet import IUPAC

		bases = self.get_base_stats_for_poly(minqual=minqual)
		
		#list N-long where N is the number of covered bases (N <= L(contig))
		dominanceList = []
 		mutationStats={'DN':0,'DS':0,'D?':0}
 		
 		explainList=[]

		codon_f1 = []
		codon_f2 = []

		for positionData in bases:
			# positionData= ((A,C,G,T),position) if covered, None if not.
			bases = ['N']

			if positionData:
				nuclAbundance,position = positionData
				base_sum=sum(nuclAbundance)
				base_max=float(max(nuclAbundance))

				if base_sum > mincov:
					dominance = float(base_max) / float(base_sum)
					dominanceList.append(dominance)
		
					tmpDict = dict((k,v) for k,v in zip(['A','C','G','T'],nuclAbundance))
					bases = [k for k,v in sorted(tmpDict.items(),key= lambda x: x[1], reverse=True) if v>0]	
			else:
				dominanceList.append(np.nan)
			
			first_base = bases[0]
			second_base = bases[1] if (len(bases) > 1 and dominance < dominant_frq_thrsh) else bases[0]
 
			codon_f1.append(first_base)
			codon_f2.append(second_base)

			if len(codon_f1) == 3 and len(codon_f2) == 3:

				codon_s1 = Seq(''.join(codon_f1),IUPAC.ambiguous_dna)
				codon_s2 = Seq(''.join(codon_f2),IUPAC.ambiguous_dna)
				codon_t1 = codon_s1.translate()
				codon_t2 = codon_s2.translate()

				positionLabel = positionData[1] if positionData else 'ND'
				
				RD=None
				if codon_t1 == "X" or codon_t2 == "X":
					mutationStats['D?'] +=1
					RD="D?"
				elif codon_t1 != codon_t2:
					mutationStats['DN'] +=1
					RD="DN"

				elif (codon_t1 == codon_t2) and (codon_s1 != codon_s2):
					mutationStats['DS'] +=1
					RD="DS"

				#if we have a mutation, save it
				if RD and positionData: explainList.append((positionLabel,RD,codon_s1,codon_s2, codon_t1,codon_t2))

				#print positionData[1] if positionData else 'ND',codon_s1,codon_s2,codon_t1,codon_t2, "RD:",RD
				codon_f1 = []
				codon_f2 = []

		return (dominanceList,mutationStats,explainList)


	def polymorphism_rate(self,mincov=10,minqual=30,pvalue=0.01,error_rate=0.001,dominant_frq_thrsh=0.8,precomputedBinomial=None):

		base_values = self.get_base_stats(min_read_depth=mincov, min_base_quality=minqual,error_rate=error_rate,precomputedBinomial=precomputedBinomial)

		rv={}
		rv['total_covered_bases'] = len(base_values)
		rv['total_polymorphic_bases'] = 0

		if len(base_values) > 0:
			pb=sum([(1 if (info['p'] < pvalue and info['ratio_max2all'] < dominant_frq_thrsh) else 0) for pox, info in base_values.items()])

			rv['total_polymorphic_bases']= pb
			rv['total_polymorphic_rate'] = float(pb)/float(len(base_values))

			# If we have at least one polymorphic site
			if pb > 0 :
				rv['ratios'] = [info['ratio_max2all'] for pox,info in base_values.items() if (info['p'] < pvalue and info['ratio_max2all'] < dominant_frq_thrsh)]
				rv['dominant_allele_distr_mean'] = np.mean(rv['ratios'])
				rv['dominant_allele_distr_sd'] = np.std(rv['ratios'])

				for i in [10,20,30,40,50,60,70,80,90,95,98,99]:
					rv['dominant_allele_distr_perc_'+str(i)] = np.percentile(rv['ratios'],i)

		return rv


	def breadth_and_depth_of_coverage(self,mincov=10,minqual=30):
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


	def depth_of_coverage(self,mincov=10,minqual=30):
		return self.breadth_and_depth_of_coverage(mincov,minqual)[1]
		#coverage_positions = {}
		#for pileupcolumn in self.bam_handle.pileup(self.name,stepper=self.stepper):
		#	if pileupcolumn.n >= mincov: coverage_positions[pileupcolumn.pos] = len([1 for pileupread in pileupcolumn.pileups if not pileupread.is_del and not pileupread.is_refskip and pileupread.alignment.query_qualities[pileupread.query_position] >= args.minqual and pileupread.alignment.query_sequence[pileupread.query_position].upper() in ('A','T','C','G') ])
		#
		#return (np.mean(coverage_positions.values()),np.median(coverage_positions.values()))


	def breadth_of_coverage(self,mincov=10,minqual=30):
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
	def get_base_stats(self, min_read_depth=10, min_base_quality=30, error_rate=0.001,dominant_frq_thrsh=0.8,precomputedBinomial=None):
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
				if precomputedBinomial and base_max in precomputedBinomial and base_sum in precomputedBinomial[base_max]:
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
		return [base_stats[k].get(stats_value, 'NaN') for k in base_stats]

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




class bcolors:
	HEADER = '\033[95m'
	OKBLUE = '\033[94m'
	OKGREEN = '\033[92m'
	WARNING = '\033[93m'
	FAIL = '\033[91m'
	ENDC = '\033[0m'        
	OKGREEN2 = '\033[42m\033[30m'
	RED = '\033[1;91m'
	CYAN = '\033[0;37m'

if __name__ == "__main__":

	def polymut_from_file(args):
		
		import pandas as pd
		import numpy as np

		outputDicts=[]

		si = True if args.sortindex else False
		mode = 'all' if args.f else 'nofilter'

		bf = BamFile(args.BAMFILE,sort=si,index=si,stepper=mode,minlen=args.minlen)

		if (args.gff_file):
			bf.parse_gff(args.gff_file)

		tl = [bf.get_contig_by_label(contig) for contig in get_contig_list(args.contig)] if args.contig is not None else list(bf.get_contigs_obj())
		for i in tl:
			dominanceArray,mutationStats,explainList = i.easy_polymorphism_rate(minqual=args.minqual,mincov=args.mincov,dominant_frq_thrsh=args.dominant_frq_thrsh)

			explain = [str(positionLabel)+":"+(bcolors.OKGREEN2 if codon_t1 == codon_t2 else bcolors.FAIL)+str(codon_t1)+'>'+str(codon_t2)+bcolors.ENDC+' ' for (positionLabel,RD,codon_s1,codon_s2, codon_t1,codon_t2) in explainList]  
			if not all(np.isnan(dominanceArray)):
				outputDicts.append({'Ref':i.name,'Len':i.length,'DN':mutationStats['DN'],'DS':mutationStats['DS'],'D?':mutationStats['D?'],'Dominance Mean': np.nanmean(dominanceArray) ,'Dominance STD':  np.nanstd(dominanceArray),'info':' '.join(explain)})
			else:
				outputDicts.append({'Ref':i.name,'Len':i.length,'DN':mutationStats['DN'],'DS':mutationStats['DS'],'D?':mutationStats['D?'],'Dominance Mean': np.nan ,'Dominance STD':  np.nan,'info':' '.join(explain)})

		out_df = pd.DataFrame.from_dict(outputDicts).set_index('Ref')
		# Check the number of non-nan entries in dominanceArray.
		considered_positions = [x for x in dominanceArray if x  != np.nan]
		print float(np.sum(out_df["DN"])), float(np.sum(out_df["DS"])), len(considered_positions)

	def poly_from_file(args):

		import pandas as pd
		if args.precomputed and os.path.isfile(args.precomputed):
			import pickle
			binomPrecomputed = pickle.load( open( args.precomputed, "rb" ))
		else: binomPrecomputed=None

		si = True if args.sortindex else False
		mode = 'all' if args.f else 'nofilter'
		bf = BamFile(args.BAMFILE,sort=si,index=si,stepper=mode,minlen=args.minlen)

		outputDF = []
		allRatios = []
		PSR_estimates = []
		allGenomeCol = {'referenceID': '-GENOME-','total_covered_bases':0,'total_polymorphic_bases':0,'total_polymorphic_rate':np.nan}
		for element in [bf.get_contig_by_label(contig) for contig in get_contig_list(args.contig)] if args.contig is not None else list(bf.get_contigs_obj()):
			
			# PSR_LIST = []
			# for k in range(0, 10):
			# 	ee=element.baseline_PSR(minqual=args.minqual,mincov=args.mincov,error_rate=args.seq_err,dominant_frq_thrsh=args.dominant_frq_thrsh,binom=binomPrecomputed)
			# 	PSR_LIST.append(ee)

			# PSR_estimates is a list of lists, with each internal list containing the monte-carlo estimates of PSR for each contig.
			# PSR_estimates.append(PSR_LIST)
			#print '1',contig
			#print '2',element
			#print '3',element.name

			tld = element.polymorphism_rate(minqual=args.minqual,mincov=args.mincov,error_rate=args.seq_err,dominant_frq_thrsh=args.dominant_frq_thrsh, precomputedBinomial=binomPrecomputed)
			tld['referenceID'] = element.name
		
			allGenomeCol['total_covered_bases'] += tld['total_covered_bases']
			allGenomeCol['total_polymorphic_bases'] += tld['total_polymorphic_bases'] 
			if 'ratios' in tld:
				allRatios = allRatios + tld['ratios']
				del tld['ratios']

			# for i in [10,20,30,40,50,60,70,80,90,95,98,99]:
			# 	tld['estimated_baseline_PSR_distr_perc'+str(i)] = np.percentile(PSR_LIST,i)

			outputDF.append(tld)
			del tld

		# In order to compute the genome-wide baseline PSR distribution, for each iteration compute the weighted averages of the PSR estimates over all contigs.
		# genome_PSR_distr = []
		# for iteration in range(len(PSR_estimates[0])):
		# 	tmp = []
		# 	for contig in range(len(PSR_estimates)):
		# 		tmp.append(PSR_estimates[contig][iteration])
		# 	w_av = np.average(tmp, weights = [x['total_covered_bases'] for x in outputDF])
		# 	genome_PSR_distr.append(w_av)


		if float(allGenomeCol['total_covered_bases']) > 0: 

			allGenomeCol['total_polymorphic_rate'] = float(allGenomeCol['total_polymorphic_bases']) / float(allGenomeCol['total_covered_bases'])

			allGenomeCol['dominant_allele_distr_mean'] = np.mean(allRatios)
			allGenomeCol['dominant_allele_distr_sd'] = np.std(allRatios)
			for i in [10,20,30,40,50,60,70,80,90,95,98,99]:
				allGenomeCol['dominant_allele_distr_perc_'+str(i)] = np.percentile(allRatios,i)
			# for i in [10,20,30,40,50,60,70,80,90,95,98,99]:
			# 	allGenomeCol['estimated_baseline_PSR_distr_perc'+str(i)] = np.percentile(genome_PSR_distr,i)


		outputDF.append(allGenomeCol)

		pd.DataFrame.from_dict(outputDF).set_index('referenceID').to_csv(sys.stdout,sep='\t')

		#print np.mean(k['polymorphic_rate'])




	def bd_from_file(args):

		si = True if args.sortindex else False
		mode = 'all' if args.f else 'nofilter'

		bf = BamFile(args.BAMFILE,sort=si,index=si,stepper=mode,minlen=args.minlen)

		tl = [bf.get_contig_by_label(contig) for contig in get_contig_list(args.contig)] if args.contig is not None else list(bf.get_contigs_obj())

		print('Contig\tBreadth\tDepth (avg)\tDepth (median)')

		for i in tl:
			bd_result = i.breadth_and_depth_of_coverage(minqual=args.minqual,mincov=args.mincov)
			print(i.name+'\t'+str(bd_result[0])+'\t'+str(bd_result[1])+'\t'+str(bd_result[2]))
		

	def consensus_from_file(args):
		from Bio import SeqIO
		from Bio.Seq import Seq
		from Bio.SeqRecord import SeqRecord
		from Bio.Alphabet import IUPAC
		import sys

		si = True if args.sortindex else False
		mode = 'all' if args.f else 'nofilter'


		bf = BamFile(args.BAMFILE,sort=si,index=si,stepper=mode,minlen=args.minlen)
		#tl = [bf.get_contig_by_label(contig) for contig in args.contig.split(',')] if args.contig is not None else list(bf.get_contigs_obj())
		tl = [bf.get_contig_by_label(contig) for contig in get_contig_list(args.contig)] if args.contig is not None else list(bf.get_contigs_obj())

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
	parser_breadth.add_argument('-c','--contig', help='Focus on a subset of references in the BAM file. Can be a list of references separated by commas or a FASTA file (the IDs are used to subset)', metavar="REFERENCE ID" ,default=None)
	parser_breadth.add_argument('-f', help='If set unmapped (FUNMAP), secondary (FSECONDARY), qc-fail (FQCFAIL) and duplicate (FDUP) are excluded. If unset ALL reads are considered (bedtools genomecov style). Default: unset',action='store_true')
	parser_breadth.add_argument('--sortindex', help='Sort and index the file',action='store_true')
	parser_breadth.add_argument('--minlen', help='Minimum Reference Length for a reference to be considered',default=0, type=int)
	parser_breadth.add_argument('--minqual', help='Minimum base quality. Bases with quality score lower than this will be discarded. This is performed BEFORE --mincov. Default: 30', type=int, default=30)
	parser_breadth.add_argument('--mincov', help='Minimum position coverage to perform the polymorphism calculation. Position with a lower depth of coverage will be discarded (i.e. considered as zero-coverage positions). This is calculated AFTER --minqual. Default: 1', type=int, default=1)
	
	parser_breadth.set_defaults(func=bd_from_file)

	parser_poly = subparsers.add_parser('poly',description="Reports the polymorpgic rate of each reference (polymorphic bases / total bases). Focuses only on covered regions (i.e. depth >= 1)")
	parser_poly.add_argument('BAMFILE', help='The file on which to operate')
	parser_poly.add_argument('-c','--contig', help='Focus on a subset of references in the BAM file. Can be a list of references separated by commas or a FASTA file (the IDs are used to subset)', metavar="REFERENCE ID" ,default=None)
	parser_poly.add_argument('-f', help='If set unmapped (FUNMAP), secondary (FSECONDARY), qc-fail (FQCFAIL) and duplicate (FDUP) are excluded. If unset ALL reads are considered (bedtools genomecov style). Default: unset',action='store_true')
	parser_poly.add_argument('--sortindex', help='Sort and index the file',action='store_true')
	
	parser_poly.add_argument('--minlen', help='Minimum Reference Length for a reference to be considered',default=0, type=int)
	parser_poly.add_argument('--minqual', help='Minimum base quality. Bases with quality score lower than this will be discarded. This is performed BEFORE --mincov. Default: 30', type=int, default=30)
	parser_poly.add_argument('--mincov', help='Minimum position coverage to perform the polymorphism calculation. Position with a lower depth of coverage will be discarded (i.e. considered as zero-coverage positions). This is calculated AFTER --minqual. Default: 10', type=int, default=10)
	parser_poly.add_argument('--seq_err', help='Sequencing error rate.', type=float, default=0.001)
	parser_poly.add_argument('--dominant_frq_thrsh', help='Cutoff for degree of `allele dominance` for a position to be considered polymorphic.', type=float, default=0.8)


	parser_poly.add_argument('--precomputed', help='Path to a pickled dictionary containing the precomputed probabilities of scipy.stats.binom function.')

	parser_poly.set_defaults(func=poly_from_file)


	###
	parser_polymut = subparsers.add_parser('polymut',description="Reports the polymorpgic rate of each reference (polymorphic bases / total bases). Focuses only on covered regions (i.e. depth >= 1)")
	parser_polymut.add_argument('BAMFILE', help='The file on which to operate')
	parser_polymut.add_argument('-c','--contig', help='Focus on a subset of references in the BAM file. Can be a list of references separated by commas or a FASTA file (the IDs are used to subset)', metavar="REFERENCE ID" ,default=None)
	parser_polymut.add_argument('--gff_file', help='Focus on protein coding gene calls on your contigs', metavar="GFF" ,default=None)
	parser_polymut.add_argument('-f', help='If set unmapped (FUNMAP), secondary (FSECONDARY), qc-fail (FQCFAIL) and duplicate (FDUP) are excluded. If unset ALL reads are considered (bedtools genomecov style). Default: unset',action='store_true')
	parser_polymut.add_argument('--sortindex', help='Sort and index the file',action='store_true')
	parser_polymut.add_argument('--minlen', help='Minimum Reference Length for a reference to be considered',default=0, type=int)
	parser_polymut.add_argument('--minqual', help='Minimum base quality. Bases with quality score lower than this will be discarded. This is performed BEFORE --mincov. Default: 30', type=int, default=30)
	parser_polymut.add_argument('--mincov', help='Minimum position coverage to perform the polymorphism calculation. Position with a lower depth of coverage will be discarded (i.e. considered as zero-coverage positions). This is calculated AFTER --minqual. Default: 10', type=int, default=10)

	parser_polymut.add_argument('--dominant_frq_thrsh', help='Cutoff for degree of `allele dominance` for a position to be considered polymorphic.', type=float, default=0.8)
	parser_polymut.set_defaults(func=polymut_from_file)
	###


	parser_consensus = subparsers.add_parser('consensus',description="outputs the consensus in FASTA format. Non covered positions (or quality-trimmed positions) are reported as a dashes: -")
	parser_consensus.add_argument('BAMFILE', help='The file on which to operate')
	parser_consensus.add_argument('-c','--contig', help='Focus on a subset of references in the BAM file. Can be a list of references separated by commas or a FASTA file (the IDs are used to subset)', metavar="REFERENCE ID" ,default=None)
	parser_consensus.add_argument('-f', help='If set unmapped (FUNMAP), secondary (FSECONDARY), qc-fail (FQCFAIL) and duplicate (FDUP) are excluded. If unset ALL reads are considered (bedtools genomecov style). Default: unset',action='store_true')
	parser_consensus.add_argument('--sortindex', help='Sort and index the file',action='store_true')
	parser_consensus.add_argument('--minqual', help='Minimum base quality. Bases with quality score lower than this will be discarded. This is performed BEFORE --mincov. Default: 0', type=int, default=30)
	parser_consensus.add_argument('--mincov', help='Minimum position coverage to perform the polymorphism calculation. Position with a lower depth of coverage will be discarded (i.e. considered as zero-coverage positions). This is calculated AFTER --minqual. Default: 1', type=int, default=10)
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
