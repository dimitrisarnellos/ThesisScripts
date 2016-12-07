#! /usr/bin/python3

import sys, re

read = ""
with open(sys.argv[1], 'r') as f:
	for line in f:
		theline = re.match('(\S+)\s+', line)
		if read == "":
			print(line.rstrip())
			read = theline.group(0)
		elif read != theline.group(0):
			print(line.rstrip())
			read = theline.group(0)