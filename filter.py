#!/usr/bin/env python3

import pysam,sys
import argparse
import numpy as np

parser = argparse.ArgumentParser()

parser.add_argument('--minlen', help='Minimum length of alignment for a read to pass', type=int, default=70)
parser.add_argument('--minqual', help='Minimum average quality for a read to pass. It is computed over the fastq Phred-scores of each read', type=int, default=30)
parser.add_argument('--maxsnps', help='Maximum edit distance on the alignment for a read to pass. It is computed over NM value from the BAM)', type=float, default=1.0)
parser.add_argument('--exclude_targets', help='Exclude these entries (FASTA file to filter out)')

args = parser.parse_args()
if args.exclude_targets:
	to_exclude = list(set([rec.strip() for rec in open(args.exclude_targets)]))
else:
	to_exclude = []

samfile = pysam.AlignmentFile("-", "rb")
passingReads = pysam.AlignmentFile("-", "wb", template=samfile)

#if args.logfile:
#	logFile = open(args.logfile,'w')

for read in samfile.fetch():
	#TODOTODO
	alignment_len = int(read.query_alignment_length)
	snps = read.get_tag('NM')
	qualities = read.query_qualities
	refname = read.reference_name

	snps_rate =float(snps) / float(read.query_alignment_length)
	meanqualities =np.mean(read.query_qualities)
 	
#	print('--\nrefname: '+str(refname)+' alignment_len: '+str(alignment_len)+' snps: '+str(snps)+' snps_rate: '+str(snps_rate)+' qualities: '+str(qualities)+' meanqualities: '+str(meanqualities))

	if (not read.is_secondary) and (alignment_len >= args.minlen) and (snps_rate <= args.maxsnps) and (meanqualities >= args.minqual) and (refname not in to_exclude):

#	if not read.is_secondary and \
#	   (int(read.query_alignment_length) >= args.minlen) and \
#	   ((float(read.get_tag('NM')) / float(read.query_length)) <= args.maxsnps) and \
#	   (np.mean(read.query_qualities) >= args.minqual) and \
#	   (read.reference_name not in to_exclude):

#		#print (read.query_name+'\t'+'PASS')
		passingReads.write(read)
#	else:
		#print (read.query_name+'\t'+'EXCLUDE')

		

##	else:
##		if args.logfile: 
##			if read.is_secondary: reason = "IS SECONDARY"
##			elif (int(read.query_alignment_length) < args.minlen): reason = "SHORT ALIGNMENT ("+str(int(read.query_alignment_length))
##			elif ((float(read.get_tag('NM')) / float(read.query_length)) > args.maxsnps): reason = "IS TOO SNIPPY ("+str( float(read.get_tag('NM')) / float(read.query_length) )+")"
##			elif (np.mean(read.query_qualities) < args.minqual): reason = "IS TOO LOW QUALITY ("+str(np.mean(read.query_qualities))+")"
##			elif (read.reference_name in to_exclude): reason = "IS IN THE EXCLUDE LIST"
##
##			logFile.write(str(read.query_name)+"\tREMOVE\t"+reason+'\n')
#sys.exit(0)
passingReads.close() 
samfile.close()
