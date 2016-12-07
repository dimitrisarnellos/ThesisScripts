#! /usr/bin/python3

#Usage: Take blast result from blasting NT and find out if fungi.
# Use as allNCBIFungi.py fungi.ncbi acces2taxidGB.sqlite(nucl_gb.accession2taxid) G1_merd_hiQual.fna.blastn | sort | uniq -c

import sys, re, sqlite3

conn = sqlite3.connect(sys.argv[2])
cursorObject = conn.cursor()
##storing IDs from fungi.ncbi
#print("Reading fungi IDs")
taxaID = []
with open(sys.argv[1], 'r') as f:
	for line in f:
		taxaID.append(line.rstrip())
taxaID = set(taxaID)
'''
#getting accession number to taxon ID mapping
speclist = {}
with open(sys.argv[2], 'r') as f:
	counter = 0
	print("Reading accession number reference. It might take a few minutes")
	for line in f:
		counter += 1
		theline = re.match('\S+\s+(\S+)\s(\S+)', line)
		#print(theline.group(1), theline.group(2))
		speclist[theline.group(1)] = theline.group(2)
		#if counter == 50:
		#	break
#print(speclist)
#speclist = set(tuple(x) for x in speclist)
'''
#print("Reading blast result file")
read = ""
hits = [[], []]
with open(sys.argv[3], 'r') as f:
	count = 0 ###
	for line in f:
		theline = re.match('(\S+)\s+(\S+).*\s(\S+)$', line)
		theRead = theline.group(1)
		if read != "":
			if theRead == read:
				hits[0].append(theline.group(2))
				hits[1].append(float(theline.group(3)))
			else:
				#getting the accession number with highest score
				readsTaxon = hits[0][hits[1].index(max(hits[1]))]
				ncbiTaxon = readsTaxon.split("|")[3]
				#
				query = 'SELECT taxid FROM tab WHERE accnr=?'
				cursorObject.execute(query, (ncbiTaxon,))
				for i in cursorObject:
					#print(i[0])
					
					if i[0] in taxaID:
						#print("Fungus")
						print(read)
					
				'''
				if ncbiTaxon in speclist:
					#print(ncbiTaxon)
					taxonId = speclist[ncbiTaxon]
					if taxonId in taxaID:
						print("Fungus")
				'''
				read = theRead
				hits = [[], []]
				hits[0].append(theline.group(2))
				hits[1].append(float(theline.group(3)))						
		else:
			read = theRead
			hits[0].append(theline.group(2))
			hits[1].append(float(theline.group(3)))

	query = 'SELECT taxid FROM tab WHERE accnr=?'
	cursorObject.execute(query, (ncbiTaxon,))
	for i in cursorObject:
		#print(i[0])
		if i[0] in taxaID:
			#print("Fungus")
			print(read)