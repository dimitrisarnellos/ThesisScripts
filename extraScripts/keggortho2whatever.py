#!/usr/bin/python3
#
# Usage: ./keggortho2whatever.py annContigs.txt [ko|path]
import sys, re
from bioservices import KEGG

k = KEGG()

if sys.argv[2] == "ko".rstrip() or sys.argv[2] == "path".rstrip():
	with open(sys.argv[1], 'r') as f:
		for line in f:
			theline = re.match('(\S+)\s+(\S+)', line)
			res = k.get("ko:"+theline.group(2))
			#print(res)
			#break
			
			les = res.split('\n')
			for l in les:
				if sys.argv[2] == "ko".rstrip():
					if re.match('DEFINITION', l):
						print(theline.group(1).split("_")[0] + "_" + theline.group(1).split("_")[1], l.lstrip('DEFINITION\t '))
				elif sys.argv[2] == "path".rstrip():
					if re.match('PATHWAY', l):
						print(theline.group(1).split("_")[0] + "_" + theline.group(1).split("_")[1], l.lstrip('PATHWAY\s+\S+\s(\S+)'))
else:
	print("Invalid mode")