#! /usr/bin/python3

#preparing data for sqlite database
import sys, re

speclist = []
with open(sys.argv[1], 'r') as f:
	for line in f:
		theline = line.split()
		lineType = len(theline)
		if lineType >= 3:
			print(theline[0], theline[1], theline[2].rstrip(":"), sep="\t")
			#speclist.append(theline)
#speclist = set(tuple(x) for x in speclist)