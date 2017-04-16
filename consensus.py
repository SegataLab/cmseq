import os,pysam
import numpy as np

# EXAMPLE OF USAGE (temporary)
## import consensus
## import numpy as np
## 
## a = consensus.BamFile('CF_TNFC005IS_t2Q15.s1.bam.sorted')
## 
## for i in a.get_contigs():
## 	print i,a.get_contig_by_label(i).reference_free_consensus()
## 	print a.get_contig_by_label(i).depth_of_coverage()  #(mean,median)
## 	print a.get_contig_by_label(i).breadth_of_coverage()


class BamFile:
	bam_handle = None
	contigs = {}

	def __init__(self,bamFile,sort=False,index=False):
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
		self.contigs = dict((r,BamContig(self.bam_handle,r,l)) for r,l in zip(bamHandle.references,bamHandle.lengths))

	def get_contigs(self): return iter(self.contigs.keys())
	def get_contig_by_label(self,contigID): return (self.contigs[contigID] if contigID in self.contigs else None)


class BamContig:

	coverage = None

	consensus = ''
	
	name = None
	length = None

	
	def __init__(self,bamHandle,contigName,contigLength):
		
		self.name = contigName
		self.length = contigLength
		self.bam_handle = bamHandle 
 

	def reference_free_consensus(self,consensus_rule=lambda array: max(array, key=array.get)):

		consensus_positions = {}
		for pileupcolumn in self.bam_handle.pileup(self.name,stepper='nofilter'):
			consensus_positions[pileupcolumn.pos] = {'A':0,'T':0,'C':0,'G':0,'N':0}
			for pileupread in pileupcolumn.pileups:
				if not pileupread.is_del and not pileupread.is_refskip:
					base = pileupread.alignment.query_sequence[pileupread.query_position].upper()
							
					if base in ['A','T','C','G']: consensus_positions[pileupcolumn.pos][base]+=1
					else: consensus_positions[pileupcolumn.pos]['N']+=1

			consensus_positions[pileupcolumn.pos] = consensus_rule(consensus_positions[pileupcolumn.pos])


		self.consensus = ''.join([(consensus_positions[position] if position in consensus_positions else 'N') for position in range(0,self.length)])
		
		del consensus_positions
		return self.consensus
		
	def depth_of_coverage(self):
		coverage_positions = {}
		for pileupcolumn in self.bam_handle.pileup(self.name,stepper='nofilter'):
			coverage_positions[pileupcolumn.pos] = pileupcolumn.n 
		
		return (np.mean(coverage_positions.values()),np.median(coverage_positions.values()))

	def breadth_of_coverage(self):
		coverage_positions = {}
		ptl=0
		
		for pileupcolumn in self.bam_handle.pileup(self.name,stepper='nofilter'):			
			if len([1 for pileupread in pileupcolumn.pileups if not pileupread.is_del and not pileupread.is_refskip and pileupread.alignment.query_sequence[pileupread.query_position].upper() in ['A','T','C','G'] ]) > 0:
				ptl += 1

		
		return float(ptl)/self.length
