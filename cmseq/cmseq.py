#!/usr/bin/env python3
from __future__ import print_function
import os
import pysam
import numpy as np
import math
import sys
from scipy import stats
from collections import defaultdict
import pickle,os

__author__ = 'Moreno Zolfo (moreno.zolfo@unitn.it), Nicolai Karcher (nicolai.karcher@unitn.it), Kun Huang (kun.huang@unitn.it)'
__version__ = '1.0.3'
__date__ = '23 September 2020'

def _initt(terminating_,_consensus_bamFile,_consensus_args):
	global terminating
	global consensus_args
	global consensus_bamFile
	terminating = terminating_
	consensus_args = _consensus_args
	consensus_bamFile = _consensus_bamFile
	

class CMSEQ_DEFAULTS:
	minqual = 30
	mincov  = 1
	minlen  = 0
	poly_error_rate = 0.001
	poly_pvalue_threshold = 0.01
	poly_dominant_frq_thrsh = 0.8
	trimReads = None


class BamFile:
	bam_handle = None
	bamFile = None
	contigs = {}

	def __init__(self,bamFile,sort=False,index=False,stepper='nofilter',minlen=CMSEQ_DEFAULTS.minlen,filterInputList=None,minimumReadsAligning=None):
		if not os.path.isfile(bamFile):
			raise Exception(bamFile+' is not accessible, or is not a file')

		if sort:
			import subprocess
			fp = bamFile+'.sorted'
			subprocess.call(['samtools','sort',bamFile,'-o',bamFile+'.sorted'])
		else: fp = bamFile

		if index: pysam.index(fp)

		self.bamFile = fp
		
		bamHandle = pysam.AlignmentFile(fp, "rb")
		
		self.bam_handle = bamHandle
		
		if filterInputList is not None:
			
			toList=[]
			if isinstance(filterInputList, list):
				toList = filterInputList
			
			elif os.path.isfile(filterInputList):
				from Bio import SeqIO
				
				with open(filterInputList, "r") as infile:

					for record in SeqIO.parse(infile, "fasta"):
						
						toList.append(record.id)
			else:
				toList = [element for element in filterInputList.split(',')]

			if minimumReadsAligning:
				self.contigs = dict((r,BamContig(self.bam_handle,r,l,stepper)) for r,l in zip(bamHandle.references,bamHandle.lengths) if (l > minlen and r in toList and bamHandle.count(contig=r,read_callback=stepper) >= minimumReadsAligning))
			else:
				self.contigs = dict((r,BamContig(self.bam_handle,r,l,stepper)) for r,l in zip(bamHandle.references,bamHandle.lengths) if (l > minlen and r in toList))
			
		else:
			if minimumReadsAligning:
				self.contigs = dict((r,BamContig(self.bam_handle,r,l,stepper)) for r,l in zip(bamHandle.references,bamHandle.lengths) if (l > minlen and bamHandle.count(contig=r,read_callback=stepper) >= minimumReadsAligning))
			else: 
				self.contigs = dict((r,BamContig(self.bam_handle,r,l,stepper)) for r,l in zip(bamHandle.references,bamHandle.lengths) if (l > minlen))
			
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
			complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N' : 'N'} 
			bases = list(string) 
			bases = [complement[base] for base in bases]
			bases.reverse()
			return ''.join(bases)
 
		try:
			with open(inputGFF) as in_handle:
				_ = next(GFF.parse(in_handle))
		except:
			print ('Parsing of GFF failed. This is probably because your biopython version is too new. Try downgrading to 1.76 or older')
			sys.exit(1)

		with open(inputGFF) as in_handle:
		
			for rec in GFF.parse(in_handle):
				tmp = []
				for r in rec.features:
					if "minced" in r.qualifiers['source'][0] or "Minced" in r.qualifiers['source'][0]:
						# This catches CRISPR repeats.
						continue
					if r.sub_features:
						prodigal_bool = 'Prodigal' in r.sub_features[0].qualifiers['source'][0] or 'prodigal' in r.sub_features[0].qualifiers['source'][0]
					else:
						prodigal_bool = 'Prodigal' in r.qualifiers['source'][0] or 'prodigal' in r.qualifiers['source'][0]
					
					if prodigal_bool:
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

		
		

	def parallel_reference_free_consensus(self,ncores=4,**kwargs):
		import multiprocessing as mp
 
		terminating = mp.Event()
		
		with mp.Pool(initializer=_initt, initargs=(terminating,self.bamFile,kwargs),processes=ncores) as pool:
			res= [x for x in pool.imap_unordered(BamFile._parallel_consensus_worker, self.contigs.keys())]
		return res

	@staticmethod
	def _parallel_consensus_worker(contigName):

		if not terminating.is_set():
			try:
				t=BamFile(consensus_bamFile,filterInputList=[contigName])
				return (contigName,t.get_contig_by_label(contigName).reference_free_consensus(**consensus_args))
			except Exception as e:
				terminating.set()
				raise
		else:
			terminating.set()
	
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


	def majority_rule(data_array):
		freq_array= data_array['base_freq']
		

		if any([v>0 for v in freq_array.values()]):
			return max(sorted(freq_array), key=freq_array.get)
		else: 
			return 'N'

	def majority_rule_polymorphicLoci(data_array):

		# Masks the consensus sequence with "*" when a polymorphic locus is found according
		# to dominant_frq_thrsh defined p-value
		
		freq_array= data_array['base_freq']
		poly_pvalue= data_array['p']

		if poly_pvalue <= 0.05: 
			return "*"
		elif any([v>0 for k,v in freq_array.items() if k != 'N']):
			return max(sorted(freq_array), key=freq_array.get)
		else: 
			return 'N'

	def reference_free_consensus(self,consensus_rule=majority_rule,mincov=CMSEQ_DEFAULTS.mincov,minqual=CMSEQ_DEFAULTS.minqual,dominant_frq_thrsh=CMSEQ_DEFAULTS.poly_dominant_frq_thrsh,noneCharacter='-',BAM_tagFilter=None, trimReads=None):

		consensus_positions = {}

		#print("A",mincov,minqual,dominant_frq_thrsh)
		for pileupcolumn,position_data in self.get_base_stats(min_read_depth=mincov, min_base_quality=minqual,dominant_frq_thrsh=dominant_frq_thrsh,BAM_tagFilter=BAM_tagFilter,trimReads=trimReads,error_rate=CMSEQ_DEFAULTS.poly_error_rate).items():
			consensus_positions[pileupcolumn] = consensus_rule(position_data)

		if len(consensus_positions) > 0 :
			self.consensus = ''.join([(consensus_positions[position] if position in consensus_positions else noneCharacter) for position in range(1,self.length+1)])
		else:
			self.consensus = noneCharacter*self.length

		#del consensus_positions

		return self.consensus



	def baseline_PSR(self,mincov=10,minqual=30,pvalue=0.01,error_rate=0.001,dominant_frq_thrsh=0.8,binom=None):
		# This function estimates the polymorphic site rate over the input contig assuming that there are no truely polymorphic sites
		# (Meaning that all observed polymorphisms are due to random sequencing error). The test also puts a threshold on the "dominance"
		# of the allele, meaning that it only reports a polymorphic base if the binomial test indicates significance AND the base is NOT sufficiently
		# dominated by the dominant base. Defaults to 0.8 dominance (dominant / all).
		from scipy import stats 

		polymorphic_empirical_loci = 0

		# Get coverage as well as values of contig
		depthsList = self.get_all_base_values('base_cov', min_base_qualit=yminqual,min_read_depth=mincov)

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

	def get_base_stats_for_poly(self,minqual=CMSEQ_DEFAULTS.minqual):
		
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
						print("One of your genes' length is not a multiple of three. Check your gff file / gene calls.")
						print("Contig name", self.name)
						print("Gene position", genes_and_positions[gene_idx])
						sys.exit()
					base_stats.extend(gene_stats)

		return base_stats

	def easy_polymorphism_rate(self,mincov=CMSEQ_DEFAULTS.mincov,minqual=CMSEQ_DEFAULTS.minqual,dominant_frq_thrsh=CMSEQ_DEFAULTS.poly_dominant_frq_thrsh):

		from Bio.Seq import Seq
		#from Bio.Alphabet import IUPAC

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
				dominance = float(base_max) / float(base_sum)
				
				if base_sum > mincov:	
					
					dominanceList.append(dominance)
					tmpDict = dict((k,v) for k,v in zip(['A','C','G','T'],nuclAbundance))
					bases = [k for k,v in sorted(tmpDict.items(), key = lambda x: x[1], reverse=True) if v>0]	
				else:
					dominanceList.append(np.nan)
			else:
				dominanceList.append(np.nan)
			
			first_base = bases[0]
			second_base = bases[1] if (len(bases) > 1 and dominance < dominant_frq_thrsh) else bases[0]
 
			codon_f1.append(first_base)
			codon_f2.append(second_base)

			if len(codon_f1) == 3 and len(codon_f2) == 3:

				codon_s1 = Seq(''.join(codon_f1))
				codon_s2 = Seq(''.join(codon_f2))
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

				codon_f1 = []
				codon_f2 = []

		return (dominanceList,mutationStats)


	def polymorphism_rate(self,mincov=CMSEQ_DEFAULTS.mincov,minqual=CMSEQ_DEFAULTS.minqual,pvalue=CMSEQ_DEFAULTS.poly_pvalue_threshold,error_rate=CMSEQ_DEFAULTS.poly_error_rate,dominant_frq_thrsh=CMSEQ_DEFAULTS.poly_dominant_frq_thrsh):

		base_values = self.get_base_stats(min_read_depth=mincov, min_base_quality=minqual,error_rate=error_rate,dominant_frq_thrsh=dominant_frq_thrsh)
		

		rv={}
		rv['total_covered_bases'] = len(base_values)
		rv['total_polymorphic_bases'] = 0

		if len(base_values) > 0:
			pb=sum([(1 if (info['p'] < pvalue and info['ratio_max2all'] < dominant_frq_thrsh) else 0) for pox, info in base_values.items()])
			

			rv['total_polymorphic_bases']= pb
			rv['total_polymorphic_rate'] = float(pb)/float(len(base_values))

			# If we have at least one polymorphic site
			if pb > 0:

				rv['ratios'] = [info['ratio_max2all'] for pox,info in base_values.items() if (info['p'] < pvalue and info['ratio_max2all'] < dominant_frq_thrsh)]
				rv['dominant_allele_distr_mean'] = np.mean(rv['ratios'])
				rv['dominant_allele_distr_sd'] = np.std(rv['ratios'])

				for i in [10,20,30,40,50,60,70,80,90,95,98,99]:
					rv['dominant_allele_distr_perc_'+str(i)] = np.percentile(rv['ratios'],i)

		return rv


	def breadth_and_depth_of_coverage(self,mincov=10,minqual=30,trunc=0):
		coverage_positions = {}
		if self.length > trunc*2: 
			# Check if the contig is long enough to be truncated
			consid_r = range(int(trunc), int(self.length - trunc))
		else:
			# If a contig is too short to be truncated, ignore the truncation. 
			# This is not nice and should be improved.
			consid_r = range(0, int(self.length))

		for pileupcolumn in self.bam_handle.pileup(self.name,stepper=self.stepper):
			#for each position
			if pileupcolumn.pos in consid_r:
				tCoverage = 0
				for pileupread in pileupcolumn.pileups:
					#for each base at the position	
					if not pileupread.is_del and not pileupread.is_refskip and pileupread.alignment.query_qualities[pileupread.query_position] >= minqual and pileupread.alignment.query_sequence[pileupread.query_position].upper() in ('A','T','C','G'):
							tCoverage +=1

					if tCoverage >= mincov:
						coverage_positions[pileupcolumn.pos] = tCoverage


		if (len(coverage_positions.keys())) > 0:
			breadth = float(len(coverage_positions.keys()))/len(consid_r)
			vals = list(coverage_positions.values())
			avgdepth = np.mean(vals)
			mediandepth = np.median(vals)
			
			return (breadth,avgdepth,mediandepth,coverage_positions.values())
		else: 
			return (np.nan,np.nan,np.nan,[np.nan])

	def depth_of_coverage(self,mincov=10,minqual=30):
		return self.breadth_and_depth_of_coverage(mincov,minqual)[1]
		#coverage_positions = {}
		#for pileupcolumn in self.bam_handle.pileup(self.name,stepper=self.stepper):
		#	if pileupcolumn.n >= mincov: coverage_positions[pileupcolumn.pos] = len([1 for pileupread in pileupcolumn.pileups if not pileupread.is_del and not pileupread.is_refskip and pileupread.alignment.query_qualities[pileupread.query_position] >= args.minqual and pileupread.alignment.query_sequence[pileupread.query_position].upper() in ('A','T','C','G') ])
		#
		#return (np.mean(coverage_positions.values()),np.median(coverage_positions.values()))


	def breadth_of_coverage(self,mincov=10,minqual=30):
		return self.breadth_and_depth_of_coverage(mincov,minqual)[0]
		
#------------------------------------------------------------------------------	

	
	def get_base_stats(self, min_read_depth=CMSEQ_DEFAULTS.mincov, min_base_quality=CMSEQ_DEFAULTS.minqual, error_rate=CMSEQ_DEFAULTS.poly_error_rate,dominant_frq_thrsh=CMSEQ_DEFAULTS.poly_dominant_frq_thrsh,BAM_tagFilter=None,trimReads=None):
		'''
		get base frequencies and quality stats,
		to use in get_all_base_values() and other functions
		'''


		if trimReads:
			mask_head_until = int(trimReads[0]) if (trimReads[0] is not None and trimReads[0] != '') else 0
			mask_tail_before = int(trimReads[1]) if (trimReads[1] is not None and trimReads[1] != '') else 0

		base_stats = defaultdict(dict) 

		ATCG=('A','T','C','G')


		#for each position (column)

		
		
		for base_pileup in self.bam_handle.pileup(self.name,stepper=self.stepper,min_base_quality=min_base_quality):
			base_freq = {'A':0,'T':0,'C':0,'G':0,'N':0}
			
			pos=base_pileup.pos+1 # 1-based


			#for each read composing the pile
			for matched_read in base_pileup.pileups:
				if not matched_read.is_del and not matched_read.is_refskip:

					b = matched_read.alignment.query_sequence[matched_read.query_position].upper()
					q = matched_read.alignment.query_qualities[matched_read.query_position]	
					#print("I am get_base_stats and this is read", matched_read.alignment.query_name, " (L=",matched_read.alignment.query_length," ) at position", matched_read.query_position, "it's a ", b)

					thisPositionBase = 'N'
					
					if not trimReads or (trimReads and ((matched_read.query_position >= mask_head_until) and (matched_read.query_position <= (matched_read.alignment.query_length-mask_tail_before) ) ) ):
						if b in ATCG: 
							if BAM_tagFilter is None or all(globals()[func](matched_read.alignment.get_tag(tag),limitValue) == True for (tag,func,limitValue) in BAM_tagFilter ):
								thisPositionBase = b								
					
					base_freq[thisPositionBase] += 1

			# calculate quality stats, ignoring N's 
			base_sum=sum([base_freq[b] for b in ATCG]) 
			base_max=float(max([base_freq[b] for b in ATCG]))

			if base_sum >= min_read_depth:
				r = base_max / base_sum
				#print r, dominant_frq_thrsh
				if r < dominant_frq_thrsh:
					#it makes sense to calculate pvalue
					p = stats.binom.cdf(base_max, base_sum, 1.0 - error_rate)
				else:
					p = 1.0

				
				base_stats[pos]['p']=p                 # quality measure
				base_stats[pos]['ratio_max2all']=r     # dominant base versus others
				base_stats[pos]['base_cov'] =base_sum  # number of reads covering the base, not counting N's
				base_stats[pos]['base_freq']=base_freq # dict: {'A':4,'T':1,'C':2,'G':0,'N':0}
			
		return base_stats

	
	def get_all_base_values(self, stats_value,  *f_args, **f_kwargs):
		'''
		get list of p values (or 'ratio_max2all' etc) for all bases that pass argument thresholds
		p_all = a.get_contig_by_label('CONTIGNAME').get_all_base_values('p', min_base_quality=30)
		'''
		base_stats = self.get_base_stats(*f_args, **f_kwargs)
		return [base_stats[k].get(stats_value, 'NaN') for k in base_stats]

	
		
def loc_gte(a,b):
	return a>=b

def loc_lte(a,b):
	return a<=b

def loc_gt(a,b):
	return a>b

def loc_lt(a,b):
	return a<b

def loc_leq(a,b):
	return a==b

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
