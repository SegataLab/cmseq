#!/usr/bin/env python3
import pysam,sys
import argparse
import numpy as np

parser = argparse.ArgumentParser()

parser.add_argument('--minlen', help='Minimum length of alignment for a read to pass', type=int, default=70)
parser.add_argument('--minqual', help='Minimum average quality for a read to pass. It is computed over the fastq Phred-scores of each read', type=int, default=30)
parser.add_argument('--maxsnps', help='Maximum edit distance on the alignment for a read to pass. It is computed over NM value from the BAM)', type=float, default=1.0)
parser.add_argument('--exclude_targets', help='Exclude these entries (FASTA file to filter out)')
parser.add_argument('--exclude_reads_bam', help='Exclude these entries (BAM file to filter out)')

args = parser.parse_args()
if args.exclude_targets:
	to_exclude = list(set([rec.strip() for rec in open(args.exclude_targets)]))
else:
	to_exclude = []

if args.exclude_reads_bam:
	ex_samfile = pysam.AlignmentFile(args.exclude_reads_bam, "rb")
	reads_to_exclude= list(set([''.join(i.query_name.split('_')[:-1]) for i in ex_samfile.fetch(until_eof=True)]))

else:
	reads_to_exclude = []

samfile = pysam.AlignmentFile("-", "rb")
passingReads = pysam.AlignmentFile("-", "wb", template=samfile)

for read in samfile.fetch():
	alignment_len = int(read.query_alignment_length)
	snps = read.get_tag('NM')

	qualities = read.query_qualities
	refname = read.reference_name
	readname = read.query_name
	snps_rate =float(snps) / float(read.query_alignment_length)
	meanqualities =np.mean(read.query_qualities)

	if (not read.is_secondary) and (alignment_len >= args.minlen) and (snps_rate <= args.maxsnps) and (meanqualities >= args.minqual) and (refname not in to_exclude) and (readname not in reads_to_exclude):
		passingReads.write(read)


passingReads.close() 
samfile.close()
