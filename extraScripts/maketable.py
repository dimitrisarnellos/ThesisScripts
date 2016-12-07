#! /usr/bin/python3
#
# Author: Dimitris Arnellos
# Usage: Counts number of mapped reads to contigs

import sys, re

#Get SAM file to memory
contigsreadsdict = {}
with open(sys.argv[1]) as f:
	for line in f:
		if re.match('\@', line):
			continue
		theline = re.match('(\S+)\s+\S+\s+(\S+)', line)
		if theline.group(2) not in contigsreadsdict:
			contigsreadsdict[theline.group(2)] = []
			contigsreadsdict[theline.group(2)].append(theline.group(1))
		else:
			contigsreadsdict[theline.group(2)].append(theline.group(1))

#Get IDs of reads from FASTA files to memory
IDdict = {}
for i in range(len(sys.argv) - 2):
	filename = sys.argv[i + 2].split('.')
	IDdict[filename[0]] = []
	with open(sys.argv[i + 2]) as file:
		for line in file:
			theline = re.match('\>(\S*)', line)			
			if theline:
				IDdict[filename[0]].append(theline.group(1))

######Make this dictionary much faster to parse
setIDs = {}
for key, val in IDdict.items():
	setIDs[key] = set(val)

####Count the reads on the contigs
score = {}
counter = 0
for i, j in contigsreadsdict.items():
	counter += 1
	score[i] = {}
	###
	'''
	if counter % 1000 == 0 or counter/len(contigsreadsdict) == 1:
		print(int(100*counter/len(contigsreadsdict)),'%', sep = "")
	'''
	###
	for k, v in setIDs.items():
		score[i][k] = 0
		for z in j:
			if z in v:
				score[i][k] += 1

###Output the table
for i in score:
	for k in IDdict:
		print(i, k, score[i][k])