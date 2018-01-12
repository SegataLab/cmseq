#!/usr/bin/env python 

import pysam,sys
import argparse
import numpy as np

parser = argparse.ArgumentParser()

parser.add_argument('--minlen', help='Minimum breadth of alignment for a read to pass (i.e. 0.7 means 70% of the read-length must be aligning) ', type=float, default=0.70)
parser.add_argument('--minqual', help='Minimum average quality for a read to pass. It is computed over the fastq Phred-scores of each read', type=int, default=30)
parser.add_argument('--maxsnps', help='Maximum SNPS+GAPS for a read to pass. It is computed over XM and XO values (Snps+gaps)', type=int, default=5)
parser.add_argument('--logfile', help='output a logfile reporting the results')

args = parser.parse_args()
samfile = pysam.AlignmentFile("-", "rb")
passingReads = pysam.AlignmentFile("-", "wb", template=samfile)

if args.logfile:
	logFile = open(args.logfile,'w')

for read in samfile.fetch():
	if not read.is_secondary and float(read.query_alignment_length) / float(read.query_length) >= args.minlen and int(read.get_tag('XM'))+int(read.get_tag('XO')) <= args.maxsnps and np.mean(read.query_qualities[1]) >= args.minqual:
		if logFile:
			logFile.write(str(read.query_name)+"\tPASS\n")
			

		passingReads.write(read)
	else:
		if logFile:
			if read.is_secondary: reason = "IS SECONDARY"
			elif float(read.query_alignment_length) / float(read.query_length) < args.minlen: reason = "SHORT ALIGNMENT ("+str(float(read.query_alignment_length) / float(read.query_length)*100)+"% of the read-length)"
			elif int(read.get_tag('XM'))+int(read.get_tag('XO')) > args.maxsnps: reason = "IS TOO SNIPPY ("+str( int(read.get_tag('XM'))+int(read.get_tag('XO')) )+")"
			elif np.mean(read.query_qualities[1]) < args.minqual: reason = "IS TOO LOW QUALITY ("+str(np.mean(read.query_qualities[1]))+")"
			logFile.write(str(read.query_name)+"\tREMOVE\t"+reason+'\n')

passingReads.close()
if logFile:
	logFile.close()
samfile.close()
