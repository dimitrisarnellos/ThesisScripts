#! /usr/bin/python3

#Usage: make sense of kegg taxonomy output: translated read version
import sys, re

	

def readorContig():
	if len(theline.group(1).split("_")[0]) > 20:
		if len(theline.group(1).split("_")) == 3:
			orfs[2].append(theline.group(1).split("_")[1] + "_" + theline.group(1).split("_")[2])
		elif len(theline.group(1).split("_")) == 4:
			orfs[2].append(theline.group(1).split("_")[1] + "_" + theline.group(1).split("_")[2] + "_" + theline.group(1).split("_")[3])
	else:
		if len(theline.group(1).split("_")) == 4:
			orfs[2].append(theline.group(1).split("_")[2] + "_" + theline.group(1).split("_")[3])
		elif len(theline.group(1).split("_")) == 5:
			orfs[2].append(theline.group(1).split("_")[2] + "_" + theline.group(1).split("_")[3] + "_" + theline.group(1).split("_")[4])

def pushLineInOrfs():
	if re.match('K\d+', theline.group(2)):
		orfs[0].append(theline.group(3))
		orfs[1].append(float(theline.group(4)))
		readorContig()
	else:
		orfs[0].append(theline.group(2))
		orfs[1].append(float(theline.group(4)))
		readorContig()

contig = ""
orfs = [[], [], []]
with open(sys.argv[1], 'r') as f:
	for line in f:
		theline = re.match('(\S+)\s+(\S+)\s+(\S+).*\s(\S+)$', line)
		#print(theline.group(1).split("_"))
		theContig = theline.group(1).split("_")[0]
		### For translated contigs
		if len(theContig) < 20:
			theContig = theContig + "_" + theline.group(1).split("_")[1]
		###

		if contig != "":
			if theContig == contig:
				pushLineInOrfs()
			else:

				##
				if orfs[0][orfs[1].index(max(orfs[1]))] == "Fungi":
					print(contig.lstrip("user:") + "_" + orfs[2][orfs[1].index(max(orfs[1]))])
				##
				#print(orfs[0][orfs[1].index(max(orfs[1]))])
				contig = theContig
				orfs = [[], [], []]
				pushLineInOrfs()
		else:
			contig = theContig
			pushLineInOrfs()
	if orfs[0][orfs[1].index(max(orfs[1]))] == "Fungi":
		print(contig.lstrip("user:") + "_" + orfs[2][orfs[1].index(max(orfs[1]))])#.lstrip("user:").rstrip("_Frame").rstrip("_reverse"))
	#print(orfs[0][orfs[1].index(max(orfs[1]))])