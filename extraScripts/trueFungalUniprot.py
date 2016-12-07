#
# I want this script to take the common IDs from bacterial and fungal results and keep the fungal ones based on e value
# Use: fungal DIAMOND output as first argument, bacterial output as second.
# Update: that was the old thing. Now its for every taxon so first input argument is the one you want to keep, so first
# 	  fungal if you want that, and then all the rest (6 more input files, bacteria plants viruses vertebrates invertebrates archaea

import sys, re


#parsing the fungal output
fungdict = {}
fungkey = ''
with open(sys.argv[1]) as f:
	for line in f:
		theline = re.match('(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+(\S+)', line)
		if (theline.group(1) not in fungdict):
			fungdict[theline.group(1)] = []
			temparray = []
			temparray.append(float(theline.group(2)))
			temparray.append(float(theline.group(3)))
			fungdict[theline.group(1)].append(temparray)
		else:
			temparray = []
			temparray.append(float(theline.group(2)))
			temparray.append(float(theline.group(3)))
			fungdict[theline.group(1)].append(temparray)

#parsing the bacterial output
bacdict = {}
backey = ''
with open(sys.argv[2]) as f:
	for line in f:
		theline = re.match('(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+(\S+)', line)
		if (theline.group(1) not in bacdict):
			bacdict[theline.group(1)] = []
			temparray = []
			temparray.append(float(theline.group(2)))
			temparray.append(float(theline.group(3)))
			bacdict[theline.group(1)].append(temparray)
		else:
			temparray = []
			temparray.append(float(theline.group(2)))
			temparray.append(float(theline.group(3)))
			bacdict[theline.group(1)].append(temparray)

#parsing the plant output
plantdict = {}
plantkey = ''
with open(sys.argv[3]) as f:
	for line in f:
		theline = re.match('(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+(\S+)', line)
		if (theline.group(1) not in plantdict):
			plantdict[theline.group(1)] = []
			temparray = []
			temparray.append(float(theline.group(2)))
			temparray.append(float(theline.group(3)))
			plantdict[theline.group(1)].append(temparray)
		else:
			temparray = []
			temparray.append(float(theline.group(2)))
			temparray.append(float(theline.group(3)))
			plantdict[theline.group(1)].append(temparray)
'''
#parsing the virus output
virdict = {}
virkey = ''
with open(sys.argv[4]) as f:
	for line in f:
		theline = re.match('(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+(\S+)', line)
		if (theline.group(1) not in virdict):
			virdict[theline.group(1)] = []
			temparray = []
			temparray.append(float(theline.group(2)))
			temparray.append(float(theline.group(3)))
			virdict[theline.group(1)].append(temparray)
		else:
			temparray = []
			temparray.append(float(theline.group(2)))
			temparray.append(float(theline.group(3)))
			virdict[theline.group(1)].append(temparray)

#parsing the vertebrate output
vertdict = {}
vertkey = ''
with open(sys.argv[5]) as f:
	for line in f:
		theline = re.match('(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+(\S+)', line)
		if (theline.group(1) not in vertdict):
			vertdict[theline.group(1)] = []
			temparray = []
			temparray.append(float(theline.group(2)))
			temparray.append(float(theline.group(3)))
			vertdict[theline.group(1)].append(temparray)
		else:
			temparray = []
			temparray.append(float(theline.group(2)))
			temparray.append(float(theline.group(3)))
			vertdict[theline.group(1)].append(temparray)

'''
#parsing the invertebrate output
invertdict = {}
invertkey = ''
with open(sys.argv[4]) as f:
	for line in f:
		theline = re.match('(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+(\S+)', line)
		if (theline.group(1) not in invertdict):
			invertdict[theline.group(1)] = []
			temparray = []
			temparray.append(float(theline.group(2)))
			temparray.append(float(theline.group(3)))
			invertdict[theline.group(1)].append(temparray)
		else:
			temparray = []
			temparray.append(float(theline.group(2)))
			temparray.append(float(theline.group(3)))
			invertdict[theline.group(1)].append(temparray)

#parsing the archaea output
archdict = {}
archkey = ''
with open(sys.argv[5]) as f:
	for line in f:
		theline = re.match('(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+(\S+)', line)
		if (theline.group(1) not in archdict):
			archdict[theline.group(1)] = []
			temparray = []
			temparray.append(float(theline.group(2)))
			temparray.append(float(theline.group(3)))
			archdict[theline.group(1)].append(temparray)
		else:
			temparray = []
			temparray.append(float(theline.group(2)))
			temparray.append(float(theline.group(3)))
			archdict[theline.group(1)].append(temparray)

allkeys = []
for i in fungdict:
	allkeys.append(i)
for i in bacdict:
	allkeys.append(i)
for i in plantdict:
	allkeys.append(i)
'''
for i in virdict:
	allkeys.append(i)
for i in vertdict:
	allkeys.append(i)
'''
for i in invertdict:
	allkeys.append(i)

for i in archdict:
	allkeys.append(i)

allkeys = list(set(allkeys))
#print(len(allkeys)) this command is nice, it gives the number of reads from the sample for which we got hits


for i in allkeys:
	dictarray = {}
	if i in fungdict:
		dictarray["Fungus"] = fungdict[i]
	if i in bacdict:
		dictarray["Bacterium"] = bacdict[i]
	if i in plantdict:
		dictarray["Plant"] = plantdict[i]
	'''
	if i in virdict:
		dictarray["Virus"] = virdict[i]
	if i in vertdict:
		dictarray["Vertebrate"] = vertdict[i]
	'''
	if i in invertdict:
		dictarray["Invertebrate"] = invertdict[i]
	
	if i in archdict:
		dictarray["Archaium"] = archdict[i]
	#print(i)
	#print(dictarray)

	allevalue = {}
	allscore = {}
	for k in dictarray:
		evalue = []
		score = []
		for j in dictarray[k]:
			evalue.append(j[0])
			score.append(j[1])
		if k == "Fungus":
			allevalue["Fungus"] = min(evalue)
			allscore["Fungus"] = max(score)
		if k == "Bacterium":
			allevalue["Bacterium"] = 0.23138452*min(evalue)
			allscore["Bacterium"] = max(score)
		if k == "Plant":
			allevalue["Plant"] = 1.85285762*min(evalue)
			allscore["Plant"] = max(score)
		'''
		if k == "Virus":
			allevalue["Virus"] = 2.05879689*min(evalue)
			allscore["Virus"] = max(score)
		if k == "Vertebrate":
			allevalue["Vertebrate"] = 2.47001966*min(evalue)
			allscore["Vertebrate"] = max(score)
		'''
		if k == "Invertebrate":
			allevalue["Invertebrate"] = 1.12400184*min(evalue)
			allscore["Invertebrate"] = max(score)
		
		if k == "Archaium":
			allevalue["Archaium"] = 7.94784927*min(evalue)
			allscore["Archaium"] = max(score)

	#print(i)

	smallestevalue = {}
	biggestscore = {}
	for key, value in allevalue.items():
		if value == min(allevalue.values()):
			smallestevalue[key] = value

	for key, value in allscore.items():
		if value == max(allscore.values()):
			biggestscore[key] = value

	if len(smallestevalue) > 1 or len(biggestscore) > 1:
		if 'Fungus' in biggestscore:
			#print(biggestscore, smallestevalue)
			#if next (iter (smallestevalue.values())) < 1.0e-10:
			with open('Fungus', 'a') as out:
				out.write(i + '\n')
		#edw gia evalue htan elif:
		#elif min(smallestevalue.values()) < 1.0e-10:
		else:
			with open('Spurious', 'a') as out:
				out.write(i + '\n')
	elif len(smallestevalue) == 1 and len(biggestscore) == 1:
		if next (iter (smallestevalue.keys())) == next (iter (biggestscore.keys())):
			#if next (iter (smallestevalue.values())) < 1.0e-10:
			with open(next (iter (smallestevalue.keys())), 'a') as out:
				out.write(i + '\n')
				#print(next (iter (smallestevalue.keys())), next (iter (smallestevalue.values())))
		else:
			print('Wrong', i)


