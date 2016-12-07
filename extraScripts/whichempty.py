#! /usr/bin/python3
#
from Bio import SeqIO
import sys

thefile = list(SeqIO.parse(sys.argv[1], "fasta"))
thefile.sort(key = len, reverse=True)
for i in thefile:
	print(i)