#!/usr/bin/env python

# Script written by Jesse Breinholt (jessebreinholt@gmail.com) and Chandra Earl (sunray1@ufl.edu) to identify and remove sequences for each taxon that had more than one sequence per locus. 
#This script processes a fasta file and deletes sequences when there remains two copies for the same loci and taxon. It requires Biopython.
#example: ./remove_duplicates.py fastafile

import sys
import os
from Bio import SeqIO

arguments = sys.argv
fastain = arguments[1]

if "-h" in sys.argv:
	sys.exit("./remove_duplicates.py fastafile")
elif len(arguments) < 1:
	sys.exit("./remove_duplicates.py fastafile")

lname, crap1 = fastain.split(".")
outfile1 =lname + "_single.fas"
outfile2 =lname + "_dups.list"
outfile3 =lname + "_keep.list"
DL=open(outfile2, "w")
KL=open(outfile3, "w")

seqdic= {}
dup=set()
keep=set()
goodcount=0
dupcount=0

fasta_sequences = SeqIO.parse(open(fastain),'fasta')
end = False

for seq in fasta_sequences:
	nparts=seq.id.split("comp")
	if nparts[0] in seqdic:
		seqdic[nparts[0]] += 1
	else:
		seqdic[nparts[0]] = 1


fasta_sequences2 = SeqIO.parse(open(fastain),'fasta')
end = False
with open(outfile1, "w") as f:		
	for seq in fasta_sequences2:
		nparts=seq.id.split("comp")
		if seqdic[nparts[0]] is 1:		
			SeqIO.write([seq], f, "fasta")
			keep.add(seq.id)
			goodcount += 1
		if seqdic[nparts[0]] > 1:
			dup.add(seq.id)
			dupcount +=1
			print seq.id + "\t" + str(seqdic[nparts[0]])
for x in dup:
	DL.write(x + "\n" )
	
for y in keep:
	KL.write(y + "\n" )
	
print "Total Seqs processed "+ str(goodcount+dupcount)
print "Seqs to Keep " + str(goodcount)
print "Seqs to delete " + str(dupcount)
per = 100* float(dupcount)/float(goodcount)
print "% delete " + str(per)
	
DL.close()
KL.close()
sys.exit()