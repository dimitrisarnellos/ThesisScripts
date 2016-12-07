#! /usr/bin/python3
#
from Bio import SeqIO
from Bio.Seq import Seq
import sys

thefile = list(SeqIO.parse(sys.argv[1], "fasta"))
theProtein = []
for i in thefile:
	revcom = i.seq.reverse_complement()
	for j in range(3):
		print(">", i.id, "_Frame_", j + 1, sep="")
		print(i.seq[j:].translate())
		print(">", i.id, "_reverse_Frame_", j + 1, sep="")
		print(revcom[j:].translate())
		'''
		if len(i.seq[j:]) % 3 == 2:
			newseq = i.seq[j:-2]
			print(newseq.translate())
		'''