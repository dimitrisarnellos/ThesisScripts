#! /usr/bin/python3
import sys, re, subprocess

contigsreadsdict = {}

with open(sys.argv[1]) as f:
	for line in f:
		theline = re.match('(\S+)\s+\S+\s+(\S+)', line)
		if theline.group(2) not in contigsreadsdict:
			contigsreadsdict[theline.group(2)] = {}
			contigsreadsdict[theline.group(2)][theline.group(1)] = []
		else:
			contigsreadsdict[theline.group(2)][theline.group(1)] = []



#This file contains all the hits against fungi even the spurious
diamoutput = {}
lines = int(subprocess.check_output(["wc", "-l", sys.argv[2]]).split()[0])

with open(sys.argv[2]) as f:
	for line in f:
		theline = re.match('(\S+)\s+(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)', line)
		if theline.group(1) not in diamoutput:
			diamoutput[theline.group(1)] = []
			diamoutput[theline.group(1)].append(theline.group(2))
			diamoutput[theline.group(1)].append(theline.group(3))
		elif diamoutput[theline.group(1)][1] < theline.group(3):
			diamoutput[theline.group(1)][0] = theline.group(2)
			diamoutput[theline.group(1)][1] = theline.group(3)

for i in contigsreadsdict:
	switch = True
	for j in contigsreadsdict[i]:
		if j in diamoutput:
			if switch:
				print(i, end="\t")
				print(diamoutput[j][0])
				switch = False
			else:
				print("\t\t"+diamoutput[j][0])
