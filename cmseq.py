#!/usr/bin/env python 

import os,pysam
import numpy as np


class BamFile:
	bam_handle = None
	contigs = {}

	def __init__(self,bamFile,sort=False,index=False,stepper='nofilter'):
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
		self.contigs = dict((r,BamContig(self.bam_handle,r,l,stepper)) for r,l in zip(bamHandle.references,bamHandle.lengths))

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

	def reference_free_consensus(self,consensus_rule=lambda array: max(array, key=array.get)):

		consensus_positions = {}
		for pileupcolumn in self.bam_handle.pileup(self.name,stepper=self.stepper):
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
		for pileupcolumn in self.bam_handle.pileup(self.name,stepper=self.stepper):
			coverage_positions[pileupcolumn.pos] = pileupcolumn.n 
		
		return (np.mean(coverage_positions.values()),np.median(coverage_positions.values()))

	def breadth_of_coverage(self):
		coverage_positions = {}
		ptl=0
		
		for pileupcolumn in self.bam_handle.pileup(self.name,stepper=self.stepper):			
			if len([1 for pileupread in pileupcolumn.pileups if not pileupread.is_del and not pileupread.is_refskip and pileupread.alignment.query_sequence[pileupread.query_position].upper() in ['A','T','C','G'] ]) > 0:
				ptl += 1

		
		return float(ptl)/self.length

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


if __name__ == "__main__":

	def bd_from_file(args):

		si = True if args.sortindex else False
		mode = 'all' if args.f else 'nofilter'
		bf = BamFile(args.BAMFILE,sort=si,index=si,stepper=mode)

		print 'Contig\tBreadth\tDepth (avg)\tDepth (median)'

		if args.contig is None:
			for i in bf.get_contigs_obj():
				print i.name+'\t'+str(i.breadth_of_coverage())+'\t'+str(i.depth_of_coverage()[0])+'\t'+str(i.depth_of_coverage()[1])
		else:
			cn= bf.get_contig_by_label(args.contig)
			if cn is not None:
				print cn.name+'\t'+str(cn.breadth_of_coverage())+'\t'+str(cn.depth_of_coverage()[0])+'\t'+str(cn.depth_of_coverage()[1])


	def consensus_from_file(args):
		from Bio import SeqIO
		from Bio.Seq import Seq
		from Bio.SeqRecord import SeqRecord
		from Bio.Alphabet import IUPAC

		si = True if args.sortindex else False
		mode = 'all' if args.f else 'nofilter'

		bf = BamFile(args.BAMFILE,sort=si,index=si,stepper=mode)

		if args.contig is None:
			lst=[]
			for i in bf.get_contigs_obj():
				lst.append(SeqRecord(Seq(i.reference_free_consensus(), IUPAC.IUPACAmbiguousDNA), id=i.name+"_consensus", description=''))
			SeqIO.write(lst,sys.stdin,'fasta')
		else:
			cn=bf.get_contig_by_label(args.contig)
			if cn is not None:
				print '>'+cn.name+'_consensus'
				print str(cn.reference_free_consensus())


	def plot_coverage_from_file(args):

		si = True if args.sortindex else False
		mode = 'all' if args.f else 'nofilter'
		bf = BamFile(args.BAMFILE,sort=si,index=si,stepper=mode)


		if args.contig is None: l = bf.get_contigs_obj()
		else:
			cn= bf.get_contig_by_label(args.contig)
			if cn is not None:
				l = [cn]

		for i in l:
			i.plot_coverage(path='./'+i.name+'.'+args.format,smooth=args.smooth,l_avoid=args.l_avoid,s_avoid=args.s_avoid,l_color=args.l_color,s_color=args.s_color,flavour=args.flavour)


	import argparse
	parser = argparse.ArgumentParser()
	subparsers = parser.add_subparsers(title='subcommands',description='valid subcommands')
	parser_breadth = subparsers.add_parser('bd',description="calculate the Breadth and Depth of coverage of BAMFILE")
	parser_breadth.add_argument('BAMFILE', help='The file on which to operate')
	parser_breadth.add_argument('-c','--contig', help='Get the breadth of a specific contig',default=None)
	parser_breadth.add_argument('-f', help='If set unmapped (FUNMAP), secondary (FSECONDARY), qc-fail (FQCFAIL) and duplicate (FDUP) are excluded. If unset ALL reads are considered (bedtools genomecov style). Default: unset',action='store_true')
	parser_breadth.add_argument('--sortindex', help='Sort and index the file',action='store_true')
	parser_breadth.set_defaults(func=bd_from_file)


	parser_consensus = subparsers.add_parser('consensus',description="outputs the consensus in FASTA format")
	parser_consensus.add_argument('BAMFILE', help='The file on which to operate')
	parser_consensus.add_argument('-c','--contig', help='Get the consensus of a specific contig',default=None)
	parser_consensus.add_argument('-f', help='If set unmapped (FUNMAP), secondary (FSECONDARY), qc-fail (FQCFAIL) and duplicate (FDUP) are excluded. If unset ALL reads are considered (bedtools genomecov style). Default: unset',action='store_true')
	parser_consensus.add_argument('--sortindex', help='Sort and index the file',action='store_true')
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
	
	parser_coverageplot.set_defaults(func=plot_coverage_from_file)


	args = parser.parse_args()
	args.func(args)
