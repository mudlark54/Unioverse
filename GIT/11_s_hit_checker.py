#!/usr/bin/env python

#Script written by Jesse Breinholt (jessebreinholt@gmail.com) and Chandra Earl (sunray1@ufl.edu) to
#process the output of BLAST to find sequences that fit the single hit criteria. This script processes a blast
#formatted table in blast -outfmt 6 from a blast analysis to find the best hits on 3 different contigs and
#allows 3 of the best hits per contig.
#(for example blast output from this: blastn -task blastn -query infile -db genomedatabasename -out outfile -outfmt 6 -max_target_seqs 3 -max_hsps 3)
#To run you need the blast table and a percent(0-1) of the second best hit bitscore/ best hit bitscore that you
#consider to be too close to determine if the sequence is single copy. For example if the best hit bit score is 100 and
#the second best hit is 90 then the second best hit bit score is >= to .90 and too close to confidently consider
#it as single hit.  We have found .90 to be a good indication that the sequences is single copy but setting this
#lower will increases the probability the sequences is single copy. Example: ./s_hit_checker.py inblasttable .90


import sys
import os

arguments = sys.argv
if '-h' in arguments:
	sys.exit("./s_hit_checker.py blast_table percent(.90)")
elif len(arguments) < 3:
	sys.exit("./s_hit_checker.py blast_table percent(.90)")

tablefile = arguments[1]
percent = arguments[2]
lname, crap1 = tablefile.split(".")
bitcheck=0
seqcheck = ''
namelist = set()
outfile = open(lname + '_del_list' + percent + ".txt", "w")
table = open(tablefile,"r")
line = table.readline()


#open blast table and parse to find sequences that have hits with bit score above or equal to the set percent of the best hit bit score
while line:
	seqname, scaf, id, alen, missm, gap, qs, qend, ts, tend, eval, bit =line.split()
	bit = float(bit)
#	print(seqname, seqcheck, bit, bitcheck)
	if seqname == seqcheck:
		if bit/bitcheck >= float(percent):
			namelist.add(seqname)
			line = table.readline()
		if bit/bitcheck < float(percent):
			line = table.readline()
	if seqname != seqcheck:
		seqcheck = seqname
		bitcheck = bit
		line = table.readline()
table.close()


#write the list of sequences to delete
for x in namelist:
	outfile.write(x + "\n")
#print the number of sequnces that have more then 1 hit fitting the percent threshold
print str(len(namelist)) + " seqs with > 1 signficant hits"
outfile.close()
