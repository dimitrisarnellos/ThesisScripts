#!/usr/bin/python3
#
# Program that outputs selected sequences from a FASTA file.
# User supplies the FASTA file as first argument and the file with the sequence IDs 
# whose sequences we want to output as the second argument.
import sys, re

idList = []

with open(sys.argv[2], 'r') as f:
	for line in f:
		line = line.rstrip()
		idList.append(line)
found = False
with open(sys.argv[1], 'r') as fu:
	for line in fu:
		if re.match('\>', line):
			seqid = re.match('\>(\S+)', line)
			if seqid.group(1) in idList:
				found = True
				line = line.rstrip()
				print(line)
			else:
				found = False
		elif found is True:
			line = line.rstrip()
			print(line)