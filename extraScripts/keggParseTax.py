#! /usr/bin/python3

#Usage: make sense of kegg taxonomy output
import sys, re

#readorcontig = ""
'''
def checkInput():
	readorcontig = input("Does the file contain reads OR contigs? (r / c) \n")
	if readorcontig != "r" and readorcontig != "c":
		checkInput()
	else:
		return readorcontig

readorcontig = checkInput()
#print(readorcontig)
'''
def pushLineInOrfs():
	if re.match('K\d+', theline.group(2)):
		orfs[0].append(theline.group(3))
		orfs[1].append(float(theline.group(4)))
	else:
		orfs[0].append(theline.group(2))
		orfs[1].append(float(theline.group(4)))

contig = ""
orfs = [[], []]
with open(sys.argv[1], 'r') as f:
	for line in f:
		theline = re.match('(\S+)\s+(\S+)\s+(\S+).*\s+(\S+)$', line)
		contParts = ""
		theContig = ""
		contParts = theline.group(1).split("_")
		theContig = contParts[0] + "_" + contParts[1]
		'''
		if readorcontig == "c":
			#print("toplel")
			contParts = theline.group(1).split("_")
			theContig = contParts[0] + "_" + contParts[1]
		elif readorcontig == "r":
			theContig == theline.group(1).split("_")[0].lstrip("user:")
		'''

		if contig != "":
			if theContig == contig:
				pushLineInOrfs()
			else:
				#print("toplel")
				##
				if orfs[0][orfs[1].index(max(orfs[1]))] == "Fungi":
					print(contig.lstrip("user:"))
				##
				#print(line.rstrip())
				#print(orfs[0][orfs[1].index(max(orfs[1]))])
				contig = theContig
				orfs = [[], []]
				pushLineInOrfs()
		else:
			contig = theContig
			pushLineInOrfs()
	#print(orfs[0][orfs[1].index(max(orfs[1]))])
	##
	if orfs[0][orfs[1].index(max(orfs[1]))] == "Fungi":
		print(contig.lstrip("user:"))
		#print("Fungi")
	##
