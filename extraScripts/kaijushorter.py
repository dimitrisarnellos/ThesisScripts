#! /usr/bin/python3

import sys, re

tax = []
with open(sys.argv[1], 'r') as f:
	for line in f:
		tax.append(line.rstrip())

tax = set(tax)

with open(sys.argv[2], 'r') as f:
	for line in f:
		theline = re.match('\S+\s+(\S+)\s+[^\t]*\t(\S+)', line)
		if theline.group(2) in tax:
			print(theline.group(1))