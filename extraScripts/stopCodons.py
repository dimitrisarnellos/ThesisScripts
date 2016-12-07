#! /usr/bin/python3

import sys, re
'''
with open(sys.argv[1], 'r') as f:
	for line in f:
		if re.match('\>', line):
			print(line.strip())
		else:
			print(line.count("*"))
'''
##SHOULD I DO AN ORF SPLITER INSTEAD???!?!?

with open(sys.argv[1], 'r') as f:
	for line in f:
		if re.match('\>', line):
			header = line.strip()
		else:
			orfs = line.strip().split("*")
			#print(orfs)
			if len(orfs) <= 2:
				properOrfs = []
				for i in range(len(orfs)):
					if len(orfs[i]) > 10:
						properOrfs.append(orfs[i])
				
				for i in range(len(properOrfs)):
					#if len(orfs[i]) > 10:
					print(header + "_ORF" + str(i + 1))
					print(properOrfs[i])
				