#! /usr/bin/python3
# Akin to trueFungalUniprot.py edited for a special occation of translated contigs

import sys, re

#parse input
fungdict = {}
with open(sys.argv[1]) as f:
	for line in f:
		theline = re.match('(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+(\S+)', line)
		contig = theline.group(1).split('_')
		contigName = contig.pop(0) + '_' + contig.pop(0)
		frameName = "_".join(contig)
		if (contigName not in fungdict):
			fungdict[contigName] = {}
			temparray = []
			temparray.append(float(theline.group(2)))
			temparray.append(float(theline.group(3)))
			fungdict[contigName][frameName] = temparray
		else:
			temparray = []
			temparray.append(float(theline.group(2)))
			temparray.append(float(theline.group(3)))
			fungdict[contigName][frameName] = temparray


#Keeping the contig frame with the highest score
fungOneFrame = {}
for i, j in fungdict.items():
	evalue = []
	score = []
	frame = []
	for k, l in j.items():
		evalue.append(l[0])
		score.append(l[1])
		frame.append(k)
	fungOneFrame[i] = [frame[score.index(max(score))], evalue[score.index(max(score))], max(score)]


#for i in fungOneFrame:
#	print(i)

#parse second input
bacdict = {}
with open(sys.argv[2]) as f:
	for line in f:
		theline = re.match('(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+(\S+)', line)
		contig = theline.group(1).split('_')
		contigName = contig.pop(0) + '_' + contig.pop(0)
		frameName = "_".join(contig)
		if (contigName not in bacdict):
			bacdict[contigName] = {}
			temparray = []
			temparray.append(float(theline.group(2)))
			temparray.append(float(theline.group(3)))
			bacdict[contigName][frameName] = temparray
		else:
			temparray = []
			temparray.append(float(theline.group(2)))
			temparray.append(float(theline.group(3)))
			bacdict[contigName][frameName] = temparray
#print(bacdict)
#Keeping the contig frame with the highest score
bacOneFrame = {}
for i, j in bacdict.items():
	evalue = []
	score = []
	frame = []
	for k, l in j.items():
		evalue.append(l[0])
		score.append(l[1])
		frame.append(k)
	bacOneFrame[i] = [frame[score.index(max(score))], evalue[score.index(max(score))], max(score)]


#See if contig is fungal or bacterial
for i, j in fungOneFrame.items():
	if i in bacOneFrame:
		if j[2] >= bacOneFrame[i][2]:
			#print(i, j[0])
			continue
		'''		
		if j[2] > bacOneFrame[i][2]:
			print(i, j[1]*4.33452520, bacOneFrame[i][1])
		'''
		

