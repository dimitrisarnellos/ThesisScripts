#! /usr/bin/python3
from Bio import SeqIO
import sys

max_len = 0
min_len = 0
counter = 0
for record in SeqIO.parse(sys.argv[1], "fasta"):
	counter += 1
	if len(record) > max_len:
		max_len = len(record)
	if counter == 1:
		min_len = len(record)
	elif len(record) < min_len:
		min_len = len(record)
print("Max length:", max_len)
print("Min length:", min_len)


#####N50 trial
thefile = list(SeqIO.parse(sys.argv[1], "fasta"))
thefile.sort(key = len, reverse=True)
lensum = 0
for i in thefile:
	lensum += len(i)
n50num = 0
for i in thefile:
	n50num += len(i)
	if n50num >= lensum/2:
		print("N50:", len(i))
		break
print("Number of contigs:", len(thefile))
print("Total contig length:", lensum)

####N50 for contigs >=500 bp

# largecontigs = []
# for record in SeqIO.parse(sys.argv[1], "fasta"):
# 	if len(record) >= 500:
# 		largecontigs.append(record)
# largecontigs.sort(key = len, reverse=True)

# print("For contigs >= 500 bp:")
# lensum2 = 0
# for i in largecontigs:
# 	lensum2 += len(i)
# n50num2 = 0
# for i in largecontigs:
# 	n50num2 += len(i)
# 	if n50num2 >= lensum2/2:
# 		print("N50:", len(i))
# 		break
# print("Total contig length:", lensum2)