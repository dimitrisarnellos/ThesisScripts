#!/home/dimitris/Documents/FinalProject/forVenv/bin/python3
#
# Author: Dimitrios Arnellos
#
# Description: The program makes a kegg search of pathways based on a list of genes that is provided and returns
# lists of genes that interact/relate to the genes on the list, each list being a file of the interacting gene,
# one file for each pathway.
#
# Functions: Two functions are defined to reduce repetitive code and add recursiveness to the script. Function
# getGeneNameFromEntries(entryName) accepts a kegg gene code id and returns the common gene id. Function 
# getTargetsFromRelations(rels, theEntry) accepts a 'relations' dictionary as outputed from the module, and based
# on that it builds the list of the interacting/related genes that we are looking for.
#
# Usage: The user runs the program from the command line, providing a file that contains the genes for which
# interactions are being investigated. The genes could be one in each line or tab - delimited. Example:
# ./keggPP.py list.txt
# Note: The program must be run having the package bioservices, as modified by the author of this program.

import sys, os
from bioservices import KEGG
k = KEGG()

############### The two needed functions
def getGeneNameFromEntries(entryName):
	genes = []
	geneIDs = entryName.split(' ')
	for oneGene in geneIDs:
		geneIDsEntry = k.get(oneGene)
		geneEntryDict = k.parse(geneIDsEntry)
		if 'NAME' in geneEntryDict:
			genes.append(geneEntryDict['NAME'][0].split(', ')[0])
	return genes


def getTargetsFromRelations(rels, theEntry):
	for relation in rels['relations']:
		
		if relation['entry1'] == theEntry['id']:
			for ent in rels['entries']:
				if ent['id'] == relation['entry2']:
					if ent['type'] == 'gene':
						for i in getGeneNameFromEntries(ent['name']):
							genes_result.append(i)
					elif ent['type'] != 'map': #to be sure we don't get errors
						getTargetsFromRelations(rels, ent)
		
		#for the opposite direction as well:	
		if relation['entry2'] == theEntry['id']:
			for ent in rels['entries']:
				if ent['id'] == relation['entry1']:
					if ent['type'] == 'gene':
						for i in getGeneNameFromEntries(ent['name']):
							genes_result.append(i)

#############################################


############################# Handle input file

with open(sys.argv[1], 'r') as fh:
	data = fh.read()
data = data.replace('\t', '\n')
inputGenes = data.split('\n')
if '' in inputGenes:
	inputGenes.remove('')

#############################


############## Main part
result = {}

for inputGene in inputGenes:
	print('Processing gene ' + inputGene + ':')
	try:
		pathways = k.get_pathway_by_gene(inputGene, 'hsa')
	except (AttributeError):
		print ('Invalid gene identifier')
		continue

	else:
		if pathways:
			for pathway in pathways:
				if pathway != 'hsa01100':
					print('\tProcessing pathway ' + pathway)
					rel = k.parse_kgml_pathway(pathway)
					genes_result = []
		
					####### Part where the functions are called
					for entry in rel['entries']:
						if entry['type'] == 'gene':
							for aGene in getGeneNameFromEntries(entry['name']):
								if aGene == inputGene:
									getTargetsFromRelations(rel, entry)
					############################

					unique_genes = list(set(genes_result))

					if pathway in result:
						result[pathway] += genes_result
					else:
						if (genes_result != []):
							result[pathway] = genes_result

#############################



########### Outputing the data
if result != {}:
	print('Outputing files')
	for i in result:
		result[i] = list(set(result[i]))
	#print(result)

	outfile = sys.argv[1].split('.')[0] + 'OUTPUT'
	if not os.path.exists(outfile):
		os.makedirs(outfile)
		for i in result:
			with open(outfile + '/' + i, 'w') as out:
				for j in result[i]:
					out.write(j + '\n')
	else:
		for i in result:
			with open(outfile + '/' + i, 'w') as out:
				for j in result[i]:
					out.write(j + '\n')

###########

