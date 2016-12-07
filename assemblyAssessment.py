#! /usr/bin/python3
import sys, re, argparse

usage = "Author: Dimitris Arnellos\n\nThis program assigns Uniprot taxonomic IDs from Blast 6 type outpout files to contigs that have being assembled from the reads that were blasted."
parser = argparse.ArgumentParser(description=usage,
	formatter_class=argparse.RawDescriptionHelpFormatter,
	)


parser.add_argument(
	'a',
	help="Alignment sam file",
	metavar='SAM',
	type=argparse.FileType('r')
	)

parser.add_argument(
	'b',
	help="Blast 6 output file",
	metavar='Blast',
	type=argparse.FileType('r')
	)

parser.add_argument(
	'-s',
	action="store_true",
	help="Gene - contig summary"
	)

parser.add_argument(
	'-t',
	action="store_true",
	help="Taxon - contig summary"
	)

parser.add_argument(
	'-u',
	action="store_true",
	help="Unique gene ID per contig"
	)

parser.add_argument(
	'-v',
	action="store_true",
	help="All genes ID per contig"
	)

parser.add_argument(
	'-p',
	action="store_true",
	help="Per contig read hit"
	)

args = parser.parse_args()
parser.parse_args()


#functions
def base():
	countdict = {}
	taxcountdict = {}
	for i in contigsreadsdict:
		countdict[i] = []
		taxcountdict[i] = []
		for j in contigsreadsdict[i]:
			if j in diamoutput:
				countdict[i].append(diamoutput[j][0])
				taxcountdict[i].append(diamoutput[j][0])
	for i in list(countdict):
		if countdict[i] == []:
			del countdict[i]
		else:
			countdict[i] = list(set(countdict[i]))
	for i in list(taxcountdict):
		if taxcountdict[i] == []:
			del taxcountdict[i]
		else:
			taxcountdict[i] = list(set(taxcountdict[i]))
			onlytaxa = []
			for j in taxcountdict[i]:
				idsplit = j.split('_')
				onlytaxa.append(idsplit[1])
			taxcountdict[i] = list(set(onlytaxa))
	countarr = []
	taxcountarr = []
	for i in countdict:
		countarr.append(len(countdict[i]))
		taxcountarr.append(len(taxcountdict[i]))
		if args.u:
			uniqgene(countdict, i)
	if args.s:
	 	summary(countarr)
	if args.t:
	 	summary(taxcountarr)
	if args.p:
		makeHitTable()

def summary(arrey):
	for i in range(1, 21):
		print(i, ":", arrey.count(i))
	print("Total:", sum(arrey))

def uniqgene(dictc, nickeln):
	switch = True
	for j in dictc[nickeln]:
		if switch:
			print(nickeln, end="\t")
			print(j)
			switch = False
		else:
 			print("\t\t"+j)

def makeHitTable():
	tempdict = dict(contigsreadsdict)
	sum = 0
	for i in tempdict:
		sum += len(list(set(tempdict[i])))
		print(i, len(list(set(tempdict[i]))))
	print("Total:", sum)





#Main part

contigsreadsdict = {}

with args.a as f:
	for line in f:
		if re.match('\@', line):
			continue
		theline = re.match('(\S+)\s+\S+\s+(\S+)', line)
		if theline.group(2) not in contigsreadsdict:
			contigsreadsdict[theline.group(2)] = {}
			contigsreadsdict[theline.group(2)][theline.group(1)] = []
		else:
			contigsreadsdict[theline.group(2)][theline.group(1)] = []

#This file contains all the hits against fungi even the spurious
diamoutput = {}
with args.b as f:
	for line in f:
		theline = re.match('(\S+)\s+(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)', line)
		if theline.group(1) not in diamoutput:
			diamoutput[theline.group(1)] = []
			diamoutput[theline.group(1)].append(theline.group(2))
			diamoutput[theline.group(1)].append(theline.group(3))
		elif diamoutput[theline.group(1)][1] < theline.group(3):
			diamoutput[theline.group(1)][0] = theline.group(2)
			diamoutput[theline.group(1)][1] = theline.group(3)


##Call the functions


if args.u or args.s or args.t or args.p:
	base()

if args.v:
	for i in contigsreadsdict:
		switch = True
		for j in contigsreadsdict[i]:
			if j in diamoutput:
				if switch:
					print(i, end="\t")
					print(diamoutput[j][0])
					switch = False
				else:
					print("\t\t"+diamoutput[j][0])

