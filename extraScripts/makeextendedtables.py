#! /usr/bin/python3

# Usage: ./makeextendedtables.py [matriKo.txt|matrixpath.txt] ../../AnalysesR/theTable.txt [ko|path]
import sys, re

dictio = {}

if sys.argv[3] == "ko".rstrip() or sys.argv[3] == "path".rstrip():
	with open(sys.argv[1]) as f:
		for line in f:
			if sys.argv[3] == "path".rstrip():
				theline = re.match('(\S+)\s+\S+\s+(\S+.*)', line)
				dictio[theline.group(1)] = theline.group(2)
			if sys.argv[3] == "ko".rstrip():
				theline = re.match('(\S+)\s+(\S+.*)', line)
				dictio[theline.group(1)] = theline.group(2)			

	with open(sys.argv[2]) as f:
		for line in f:
			theline = re.match('(\S+)\s+(\S+)\s+(\S+)', line)
			if theline.group(1) in dictio:
				print(theline.group(1), theline.group(2), theline.group(3), dictio[theline.group(1)], sep="\t")
else:
	print("Invalid mode")
