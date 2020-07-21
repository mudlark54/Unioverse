#!/usr/bin/env python
#Script written by Jesse Breinholt (jessebreinholt@gmail.com) and Chandra Earl (sunray1@ufl.edu)
#to process the output of BLAST to find if the location of the best hit on the genome is the same
#location as the probe target from that genome. This script processes a blast formatted table in
#blast -outfmt 6 from a blast analysis to find the single best hit location of a sequence
#(for example blast output from this: blastn -task blastn -query infile -db genomedatabasename -out outfile -outfmt 6 -max_target_seqs 1 -max_hsps 1).
#To run you need the blast table and the name of reference taxa that must be in the name of
#sequences of each loci from the reference genome. The fasta file you use to blast must contain
#sequence data from the reference genome for each loci you what to test. It will then compare the
#location where the data from the reference hit its genome to where the other sequences hit the
#genome to check for orthology. Example: ./ortholog_filter.py blasttable BMORI

import sys

arguments = sys.argv
if "-h" in arguments:
	sys.exit("./ortholog_filter blasttable BMORI(reference taxa name)")
elif len(arguments) < 3:
	sys.exit("./ortholog_filter blasttable BMORI(reference taxa name)")

tablefile = arguments[1]
probename = arguments[2]
lname, extension = tablefile.split(".")
outfile1 =open(lname + "_del_list.txt", "w")
outfile2 =open(lname + "_keep_list", "w")

hit=[]
loci=set([])
count=0
seqset=set([])


#open table and parse for loci in the blast table
with open(tablefile, "r") as table:
	makein=table.readlines()
	for i in makein:
		loci.add(i.split("_",-1)[0])
		ALL_loci=list(loci)
table.close()

#process by loci recording the references scaffold and coordinates
for x in ALL_loci:
	print "Processing " + x + " .............\n"
	with open(tablefile, "r") as table2:
		makein2=table2.readlines()
		for i in makein2:
			taxa, scaf, id, length, mismatch, gaps, qs, qend, ts, tend, evalue, bit=i.split()
			if taxa.startswith(str(x) + "_"):
				if taxa == str(x) + "_" + probename + "_R":
					hit.append(scaf)
					print taxa + " scaffold : " + scaf
					leftcoord=int(ts)
					rightcoord=int(tend)
					if int(ts) < int(tend):
						direction= int(1)
					if int(ts) > int(tend):
						direction = int(0)
#					print taxa + " direction(fwd=1, rev=0) : " + str(direction)
	table2.close()		


# open table again to check to see if each sequence hit scaf and coordinates of reference
	with open(tablefile, "r") as table3:
		makein3=table3.readlines()
		for i in makein3:
			taxa3, scaf3, id3, length3, mismatch3, gaps3, qs3, qend3, ts3, tend3, evalue3, bit3=i.split()
			seqset.add(taxa3)
			if int(ts3) < int(tend3):
				seqdirection= int(1)
			if int(ts3) > int(tend3):
				seqdirection = int(0)
#			print "seq direction(fwd=1, rev=0) : " + str(seqdirection)
			if taxa3.startswith(str(x) + "_"):
				if scaf3 not in hit:
					print "diffent scaffold " + taxa3 + " scaffold : " + scaf3
					outfile1.write(taxa3 + "\n")
					count +=1
				if scaf3 in hit: 
					if direction is 1 and seqdirection is 1:
						if int(ts3) < rightcoord and int(tend3) > leftcoord:
							outfile2.write(taxa3 + "\n")
						else:
							outfile1.write(taxa3 + "\n")
							print "Same scaffold, different location Direction (ref fwd : seq fwd) " + str(direction) + ":"+ str(seqdirection)
							print str(leftcoord) + " " + str(rightcoord) + "|" + i
							count +=1
					
					if direction is 1 and seqdirection is 0:
						if int(tend3) < rightcoord and int(ts3) > leftcoord:
							outfile2.write(taxa3 + "\n")
						else:
							outfile1.write(taxa3 + "\n")
							print "Same scaffold, different location Direction(ref fwd: seq rev) " + str(direction) + ":"+ str(seqdirection)
							print str(leftcoord) + " " + str(rightcoord) + "|" + i
							count +=1					
					if direction is 0 and seqdirection is 0:
						if int(tend3) < leftcoord and int(ts3) > rightcoord:
							outfile2.write(taxa3 + "\n")
						else:
							outfile1.write(taxa3 + "\n")
							print "Same scaffold, different location Direction(ref rev: seq rev) " + str(direction) + ":"+ str(seqdirection)
							print str(leftcoord) + " " + str(rightcoord) + "|" + i
							count +=1
					if direction is 0 and seqdirection is 1:
						if int(ts3) < leftcoord and int(tend3) > rightcoord:
							outfile2.write(taxa3 + "\n")
						else:
							outfile1.write(taxa3 + "\n")
							print "Same scaffold, different location Direction(ref rev: seq fwds) " + str(direction) + ":"+ str(seqdirection)
							print str(leftcoord) + " " + str(rightcoord) + "|" + i
							count +=1







print str(count) + "/" + str(len(seqset)) + " (deleted sequences/total sequences)"


table3.close()
outfile1.close()
outfile2.close()