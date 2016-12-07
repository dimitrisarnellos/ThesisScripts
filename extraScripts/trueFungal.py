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


#parsing the invertebrate output
invertdict = {}
invertkey = ''
with open(sys.argv[6]) as f:
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

#parsing the invertebrate output
archdict = {}
archkey = ''
with open(sys.argv[7]) as f:
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
for i in virdict:
	allkeys.append(i)
for i in vertdict:
	allkeys.append(i)
for i in invertdict:
	allkeys.append(i)
for i in archdict:
	allkeys.append(i)

allkeys = list(set(allkeys))
#print(len(allkeys)) this command is nice, it gives the number of reads from the sample for which we got hits


for i in allkeys:
	dictarray = []
	if i in fungdict:
		dictarray.append(fungdict)
	if i in bacdict:
		dictarray.append(bacdict)
	if i in plantdict:
		dictarray.append(plantdict)
	if i in virdict:
		dictarray.append(virdict)
	if i in vertdict:
		dictarray.append(vertdict)
	if i in invertdict:
		dictarray.append(invertdict)
	if i in archdict:
		dictarray.append(archdict)

	allevalue = {}
	allscore = {}
	for k in dictarray:
		for j in k:
			if j == i:
				evalue = []
				score = []
				for l in k[j]:
					evalue.append(l[0])
					score.append(l[1])
				if k == fungdict:
					allevalue["Fungus"] = min(evalue)
					allscore["Fungus"] = max(score)
				if k == bacdict:
					allevalue["Bacterium"] = 0.14842337*min(evalue)
					allscore["Bacterium"] = max(score)
				if k == plantdict:
					allevalue["Plant"] = 1.03201316*min(evalue)
					allscore["Plant"] = max(score)
				if k == virdict:
					allevalue["Virus"] = 2.05879689*min(evalue)
					allscore["Virus"] = max(score)
				if k == vertdict:
					allevalue["Vertebrate"] = 2.47001966*min(evalue)
					allscore["Vertebrate"] = max(score)
				if k == invertdict:
					allevalue["Invertebrate"] = 1.47997259*min(evalue)
					allscore["Invertebrate"] = max(score)
				if k == archdict:
					allevalue["Archaium"] = 2.76149162*min(evalue)
					allscore["Archaium"] = max(score)

	smallestevalue = []
	biggestscore = []
	for key, value in allevalue.items():
		if value == min(allevalue.values()):
			smallestevalue.append(key)
			#print(key, end=' ')

	for key, value in allscore.items():
		if value == max(allscore.values()):
			biggestscore.append(key)

	if len(smallestevalue) > 1 or len(biggestscore) > 1:
		with open('Spurious', 'a') as out:
			out.write(i + '\n')
	elif len(smallestevalue) == 1 and len(biggestscore) == 1:
		if smallestevalue[0] == biggestscore[0]:
			with open(smallestevalue[0], 'a') as out:
				out.write(i + '\n')
		else:
			print('Wrong', i)


'''
for i in fungdict:
	if i in bacdict:
		fungevaltemp = []
		bacevaltemp = []
		fungscoretemp = []
		bacscoretemp = []
		for j in range(len(fungdict[i])):
			fungevaltemp.append(fungdict[i][j][0])
			fungscoretemp.append(fungdict[i][j][1])
		for h in range(len(bacdict[i])):
			bacevaltemp.append(bacdict[i][h][0])
			bacscoretemp.append(bacdict[i][h][1])
		#print(6.737*min(fungevaltemp), max(fungscoretemp), min(bacevaltemp), max(bacscoretemp), end=' ')
		if 6.737*min(fungevaltemp) < min(bacevaltemp) and max(fungscoretemp) > max(bacscoretemp):
			print(i)
'''
'''
			print('Fungal')
		elif 6.737*min(fungevaltemp) < min(bacevaltemp):
			print('Spurious')
		else:
			print('Bacterial')
'''
