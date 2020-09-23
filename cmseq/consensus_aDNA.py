from .cmseq import CMSEQ_DEFAULTS
from .cmseq import BamFile
from .cmseq import BamContig
import os
import pysam
import math

import pandas as pd
import numpy as np
import argparse
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord






__author__ = 'Kun D. Huang (kun.huang@unitn.it), Moreno Zolfo (moreno.zolfo@unitn.it)'
__version__ = '1.0.3'
__date__ = '09 September 2020'

class CMSEQ_DEFAULTS_Ancient(CMSEQ_DEFAULTS):
	position_specific_prob = None
	position_specific_prob_thrsh = None


class BamFileAncient(BamFile):

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
				self.contigs = dict((r,BamContigAncient(self.bam_handle,r,l,stepper)) for r,l in zip(bamHandle.references,bamHandle.lengths) if (l > minlen and r in toList and bamHandle.count(contig=r,read_callback=stepper) >= minimumReadsAligning))
			else:
				self.contigs = dict((r,BamContigAncient(self.bam_handle,r,l,stepper)) for r,l in zip(bamHandle.references,bamHandle.lengths) if (l > minlen and r in toList))
			
		else:
			if minimumReadsAligning:
				self.contigs = dict((r,BamContigAncient(self.bam_handle,r,l,stepper)) for r,l in zip(bamHandle.references,bamHandle.lengths) if (l > minlen and bamHandle.count(contig=r,read_callback=stepper) >= minimumReadsAligning))
			else: 
				self.contigs = dict((r,BamContigAncient(self.bam_handle,r,l,stepper)) for r,l in zip(bamHandle.references,bamHandle.lengths) if (l > minlen))

	def get_contigs_obj(self): return iter(self.contigs.values())
			

class BamContigAncient(BamContig):
	## Get base stats code goes here

	def reference_free_consensus(self,consensus_rule=BamContig.majority_rule,mincov=CMSEQ_DEFAULTS.mincov,
		minqual=CMSEQ_DEFAULTS.minqual,dominant_frq_thrsh=CMSEQ_DEFAULTS.poly_dominant_frq_thrsh,
		noneCharacter='-',BAM_tagFilter=None,trimReads=None,post_damage_prob=None,
		pos_prob_db=CMSEQ_DEFAULTS_Ancient.position_specific_prob,refseq_idx=None):

		consensus_positions = {}

		for pileupcolumn,position_data in self.get_base_stats(min_read_depth=mincov, min_base_quality=minqual,
			dominant_frq_thrsh=dominant_frq_thrsh,BAM_tagFilter=BAM_tagFilter,trimReads=trimReads,
			post_damage_prob=post_damage_prob,pos_prob_db=pos_prob_db, refseq_idx=refseq_idx).items():
			ref_base_idx = self.name + '__' + str(pileupcolumn)
			

			if float(position_data['ratio_max2all']) >= float(dominant_frq_thrsh):
				consensus_positions[pileupcolumn] = consensus_rule(dict((k,v) for k,v in position_data['base_freq'].items() if k != 'N'))

		if len(consensus_positions) > 0 :
			self.consensus = ''.join([(consensus_positions[position] if position in consensus_positions else noneCharacter) for position in range(1,self.length+1)])
		else:
			self.consensus = noneCharacter*self.length

		del consensus_positions
		return self.consensus

	def get_base_stats(self, min_read_depth=CMSEQ_DEFAULTS.mincov, min_base_quality=CMSEQ_DEFAULTS.minqual,
	 error_rate=CMSEQ_DEFAULTS.poly_error_rate,dominant_frq_thrsh=CMSEQ_DEFAULTS.poly_dominant_frq_thrsh,
	 BAM_tagFilter=None,trimReads=None,post_damage_prob=CMSEQ_DEFAULTS_Ancient.position_specific_prob_thrsh,
	 pos_prob_db=CMSEQ_DEFAULTS_Ancient.position_specific_prob,refseq_idx=None):
		
		'''
		get base frequencies and quality stats,
		to use in get_all_base_values() and other functions
		'''

		from scipy import stats
		from collections import defaultdict
		import pickle,os


		base_stats = defaultdict(dict) 

		ATCG=('A','T','C','G')

		#for each position (column)
		for base_pileup in self.bam_handle.pileup(self.name,stepper=self.stepper):
			base_freq = {'A':0,'T':0,'C':0,'G':0,'N':0}
			

			pos=base_pileup.pos+1 # 1-based

			#for each read composing the pile

			for matched_read in base_pileup.pileups:
				if not matched_read.is_del and not matched_read.is_refskip:


					b = matched_read.alignment.query_sequence[matched_read.query_position].upper()
					q = matched_read.alignment.query_qualities[matched_read.query_position]	

					thisPositionBase = 'N'
					
					if post_damage_prob and pos_prob_db and refseq_idx: # Enter position-specific mode
						ref_base_key = self.name + '__' + str(pos)
						ref_base = refseq_idx[ref_base_key]	
						if matched_read.query_position <= 11: # Check position if on the left end 
							left_pos = matched_read.query_position + 1 # 1-based
							sub = ref_base+b
							if sub == 'CT' or sub == 'GA': # Check if the RefSeq-Read base if CT or GA
								prob = pos_prob_db[left_pos][sub]
								# print(sub, prob, matched_read.alignment.query_length, matched_read.query_position, left_pos)
								if (q >= min_base_quality) and (b in ATCG) and (prob <= post_damage_prob):
									if BAM_tagFilter is None:
										thisPositionBase = b
									elif BAM_tagFilter and all(globals()[func](matched_read.alignment.get_tag(tag),limitValue) == True for (tag,func,limitValue) in BAM_tagFilter):
										thisPositionBase = b
							else:
								if (q >= min_base_quality) and (b in ATCG):
									if BAM_tagFilter is None:
										thisPositionBase = b
									elif BAM_tagFilter and all(globals()[func](matched_read.alignment.get_tag(tag),limitValue) == True for (tag,func,limitValue) in BAM_tagFilter):
										thisPositionBase = b

						elif (matched_read.alignment.query_length - matched_read.query_position) <= 11: # check position if on the right end
							right_pos = matched_read.query_position - matched_read.alignment.query_length 
							sub = ref_base+b
							if sub == 'CT' or sub == 'GA':
								prob = pos_prob_db[right_pos][sub]
								# print(sub, prob, matched_read.alignment.query_length, matched_read.query_position, right_pos)
								if (q >= min_base_quality) and (b in ATCG) and (prob <= post_damage_prob):
									if BAM_tagFilter is None:
										thisPositionBase = b
									elif BAM_tagFilter and all(globals()[func](matched_read.alignment.get_tag(tag),limitValue) == True for (tag,func,limitValue) in BAM_tagFilter):
										thisPositionBase = b
							else:
								if (q >= min_base_quality) and (b in ATCG):
									if BAM_tagFilter is None:
										thisPositionBase = b
									elif BAM_tagFilter and all(globals()[func](matched_read.alignment.get_tag(tag),limitValue) == True for (tag,func,limitValue) in BAM_tagFilter):
										thisPositionBase = b
						else:
						# 	print(sub, matched_read.alignment.query_length, matched_read.query_position, "X")
							if (q >= min_base_quality) and (b in ATCG):
								if BAM_tagFilter is None:
									thisPositionBase = b
								elif BAM_tagFilter and all(globals()[func](matched_read.alignment.get_tag(tag),limitValue) == True for (tag,func,limitValue) in BAM_tagFilter):
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

def consensus_from_file():

	parser = argparse.ArgumentParser(description="outputs the consensus in FASTA format. Non covered positions (or quality-trimmed positions) are reported as a dashes: -")
	parser.add_argument('BAMFILE', help='The file on which to operate')
	parser.add_argument('-c','--contig', help='Focus on a subset of references in the BAM file. Can be a list of references separated by commas or a FASTA file (the IDs are used to subset)', metavar="REFERENCE ID" ,default=None)
	parser.add_argument('-f', help='If set unmapped (FUNMAP), secondary (FSECONDARY), qc-fail (FQCFAIL) and duplicate (FDUP) are excluded. If unset ALL reads are considered (bedtools genomecov style). Default: unset',action='store_true')
	parser.add_argument('-r', '--refseq', help='Input the refrence genome sequence', type=str)
	parser.add_argument('--sortindex', help='Sort and index the file',action='store_true')
	parser.add_argument('--minqual', help='Minimum base quality. Bases with quality score lower than this will be discarded. This is performed BEFORE --mincov. Default: '+str(CMSEQ_DEFAULTS.minqual), type=int, default=CMSEQ_DEFAULTS.minqual)
	parser.add_argument('--mincov', help='Minimum position coverage to perform the polymorphism calculation. Position with a lower depth of coverage will be discarded (i.e. considered as zero-coverage positions). This is calculated AFTER --minqual. Default: '+str(CMSEQ_DEFAULTS.minlen), type=int, default=CMSEQ_DEFAULTS.mincov)
	parser.add_argument('--dominant_frq_thrsh', help='Cutoff for degree of `allele dominance` for a position to be considered polymorphic. Default: '+str(CMSEQ_DEFAULTS.poly_dominant_frq_thrsh), type=float, default=CMSEQ_DEFAULTS.poly_dominant_frq_thrsh)
	parser.add_argument('--minlen', help='Minimum Reference Length for a reference to be considered. Default: '+str(CMSEQ_DEFAULTS.minlen),default=CMSEQ_DEFAULTS.minlen, type=int)
	parser.add_argument('--pos_specific_prob_tab', help='Stats_out_MCMC_correct_prob table produced from mapdamage2. It contains the position specific probability of observing a C->T or G->A due to a post-mortem damage.',default=CMSEQ_DEFAULTS_Ancient.position_specific_prob, type=str)
	parser.add_argument('--pos_damage_prob_thrsh', help = 'Maximum post-mortem damage probability for a nucletide on a read to be considered when building consensus.', default=CMSEQ_DEFAULTS_Ancient.position_specific_prob_thrsh, type = float)

	args = parser.parse_args()

	si = True if args.sortindex else False
	mode = 'all' if args.f else 'nofilter'

	bf = BamFileAncient(args.BAMFILE,sort=si,index=si,stepper=mode,minlen=args.minlen,filterInputList=args.contig)
	#tl = [bf.get_contig_by_label(contig) for contig in args.contig.split(',')] if args.contig is not None else list(bf.get_contigs_obj())
 
	lst = []
	if args.pos_specific_prob_tab and args.pos_damage_prob_thrsh and args.refseq: 
		pos_specific_prob_db = [i.rstrip().split(',') for i in open(args.pos_specific_prob_tab).readlines()][1:]
		stats_db = {}
		for i in pos_specific_prob_db:
			pos_ = int(i[1])
			CT_ = float(i[2])
			GA_ = float(i[3])
			stats_db[pos_] = {'CT': CT_, 'GA': GA_}
		pos_stats_db = stats_db
		pos_prob_thrsh = args.pos_damage_prob_thrsh

		RefSeq_dict = SeqIO.to_dict(SeqIO.parse(open(args.refseq), "fasta"))
		RefSeq_idx = {}
		for i in RefSeq_dict:
			seq = RefSeq_dict[i].seq
			for b_idx in range(len(seq)):
				RefSeq_idx[i+'__'+str(b_idx+1)]=seq[b_idx]

	else:
		pos_stats_db, pos_prob_thrsh, RefSeq_idx = None, None, None
		sys.exit("Please input position-specific probability table from mapdamage2, reference sequence, and damage probability cap!")


	for i in bf.get_contigs_obj():

		
		sq = i.reference_free_consensus(mincov=args.mincov,minqual=args.minqual,
			dominant_frq_thrsh=args.dominant_frq_thrsh,noneCharacter='N',
			trimReads=None,post_damage_prob=pos_prob_thrsh,pos_prob_db=pos_stats_db, refseq_idx=RefSeq_idx)
		
		if sq is not None:
			lst.append(SeqRecord(Seq(sq), id=i.name+"_consensus", description=''))
	SeqIO.write(lst,sys.stdout,'fasta')


if __name__ == "__main__":
	consensus_from_file()
