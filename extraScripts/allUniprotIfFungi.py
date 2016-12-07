#! /usr/bin/python3

#Usage: Take result.m8 from diamonding contigs to whole of uniprot and find out if fungi.
# Use as allUniprotFungi.py speclist.txt taxonomy-ancestor%3A4751.list result.m8 | sort | uniq -c

import sys, re, sqlite3

###
conn = sqlite3.connect(sys.argv[1])
cursorObject = conn.cursor()
###

'''
speclist = []
with open(sys.argv[1], 'r') as f:
	for line in f:
		theline = line.split()
		lineType = len(theline)
		if lineType >= 3:
			speclist.append(theline)
speclist = set(tuple(x) for x in speclist)
'''

taxaID = []
with open(sys.argv[2], 'r') as f:
	for line in f:
		taxaID.append(line.rstrip())
taxaID = set(taxaID)


contig = ""
hits = [[], []]
with open(sys.argv[3], 'r') as f:
	count = 0 ###
	for line in f:
		theline = re.match('(\S+)\s+(\S+).*\s(\S+)$', line)
		theContig = theline.group(1)
		if contig != "":
			if theContig == contig:
				hits[0].append(theline.group(2))
				hits[1].append(float(theline.group(3)))
			else:
				#taking the second part of the x_y ID of uniprot which is taxonomy, with the highest score
				contigsTaxon = hits[0][hits[1].index(max(hits[1]))]
				uniprotTaxon = contigsTaxon.split("_")[1]
				#
				query = 'SELECT taxid FROM tab WHERE uniId=?'
				cursorObject.execute(query, (uniprotTaxon,))
				for i in cursorObject:
					if i[0] in taxaID:
						print(contig)
						#print("Fungus")
				'''
				for i in speclist:
					if uniprotTaxon == i[0]:
						taxonId = i[2].rstrip(":")						
						if taxonId in taxaID:
							print("Fungus")
						elif i[1] == "E":
							print("Other Eucaryot")
						elif i[1] == "B":
							print("Bacterium")
						elif i[1] == "A":
							print("Archeum")
				'''
				contig = theContig
				hits = [[], []]
				hits[0].append(theline.group(2))
				hits[1].append(float(theline.group(3)))						
		else:
			contig = theContig
			hits[0].append(theline.group(2))
			hits[1].append(float(theline.group(3)))
			
	#taking the second part of the x_y ID of uniprot which is taxonomy, with the highest score
	contigsTaxon = hits[0][hits[1].index(max(hits[1]))]
	uniprotTaxon = contigsTaxon.split("_")[1]
	#
	query = 'SELECT taxid FROM tab WHERE uniId=?'
	cursorObject.execute(query, (uniprotTaxon,))
	for i in cursorObject:
		if i[0] in taxaID:
			print(contig)
		#uniprotTaxon = theline.group(2).split("_")[1]
		#[i[2] for i in speclist if i[0] == uniprotTaxon]
		#print(uniprotTaxon)