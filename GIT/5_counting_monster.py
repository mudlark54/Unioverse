#!/usr/bin/env python
#Script written by Jesse Breinholt (jessebreinholt@gmail.com) and Chandra Earl (sunray1@ufl.edu) to count the loci per taxa and put into a tab separated matrix
#To use the script, the processed fasta file must contain sequences data on a single line. The script will not work on a multiple line fasta file. 
#The script needs a single line formatted fasta file containing sequences from many loci and taxa and a character as a delimiter to split the loci from the taxa name. 
#I suggest you use an underscore. The name of the loci must be in the first part of the sequence name, then the delimiter, then what ever else you want
#(do not use dashes"-" or spaces in sequences names) for example: L001_Bombyx_mori. 
# Other delimiters have not been tested but should work.
#The script also expects each sequences to end in either _R being a reference sequences or contain the word "comp" which is added in the IBA assembly for each sequence.
#example: ./counting_monster.py fastafile _ 

import sys

arguments = sys.argv
if "-h" in arguments:
	sys.exit("counting_monster.py fastafile delimiter")
elif len(arguments) < 3:
	sys.exit("counting_monster.py fastafile delimiter")

filein = arguments[1]
delim = arguments[2]
locuslist = []

#open input file and parse out a set of loci
inputopen = open(filein, "r")
line = inputopen.readline()
while line:
	if line[0] == ">":
		linesplit = line.strip(">").strip().split(delim)
		locus = linesplit[0]
		locuslist.append(locus)
		line = inputopen.readline()
	else:
		line = inputopen.readline()
inputopen.close()
locuslist = set(locuslist)
locuslist = list(locuslist)

identifier_lists = [[] for i in range(0, len(locuslist))]
count_list = [0]*len(locuslist)

#open input file again and count the number of loci/taxa
inputopen = open(filein, "r")
line = inputopen.readline()
while line:
	if line[0] == ">":
		linesplit = line.strip(">").strip().split(delim)
		locus = linesplit[0]
		if "comp" in line:
			identifier_number = locuslist.index(locus)
			half1, half2 = line.strip().split(delim + delim)
			halflist = half1.split(delim)
			del(halflist[0])
			half1 = delim.join(halflist)
			identifier_lists[identifier_number].append(half1)
			line = inputopen.readline()
		if delim + "R" in line:
			identifier_number = locuslist.index(locus)
			line = line.strip()
			linelist = line.split(delim+delim)
			del(linelist[0])
			linestring = ''.join(linelist)
			identifier_lists[identifier_number].append(linestring)
			line = inputopen.readline()
		else:
			line = inputopen.readline()
	else:
		line = inputopen.readline()
inputopen.close()

#Write everything to a file

tableout = open("tableout.txt", "w")
for locus in locuslist:
	tableout.write("\t"+locus)
tableout.write("\n")

identifier_list_all = sum(identifier_lists, [])

for identifier in set(identifier_list_all):
	tableout.write(identifier)
	count_number = 0
	for list in identifier_lists:
		count = list.count(identifier)
		tableout.write("\t" + str(count))
		current = count_list[count_number] + count
		count_list[count_number] = current
		count_number = count_number + 1
	count2 = identifier_list_all.count(identifier)
	tableout.write("\t" + str(count2))
	tableout.write("\n")

for number in count_list:
	tableout.write("\t" + str(number))	
tableout.close()

