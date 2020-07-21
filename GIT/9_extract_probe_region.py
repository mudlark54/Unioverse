#!/usr/bin/env python

#Script written by Jesse Breinholt (jessebreinholt@gmail.com) and Chandra Earl (sunray1@ufl.edu) to split alignments into head, probe, and tail regions based on the beginning and end of a reference sequence in the alignment for a list of alignment fasta files where the data is on a single line.
#To use the script the aligned fasta file must contain sequences data on a single line. The script will not work on a multiple line alignment fasta file. You must give the script a list of files you want to process, name of the reference sequence you want to use to trim the alignment, and a out directory (./extract_probe_region.py list.txt refname outdir). 
#The list must contain the name of a single fasta file per line and the files must be located in the directory where you submit the script unless the full path to the alignments are provided in the list.
#The name of the reference sequences will define the probe region with the first base (ATGC) and last BASE (ATGC) in the alignment. 
#example: ./extract_probe_region.py inlist BMORI outdir



import os, sys

arguments = sys.argv

inputfiles = []
seqlist = []

if "-h" in sys.argv:
	sys.exit("./extract_probe_region.py list.txt refname outdir")
elif len(arguments) < 3:
	sys.exit("./extract_probe_region.py list.txt refname outdir")

listinput = arguments[1]
refname = arguments[2]
outdir = arguments[3]

listin = open(listinput, "r")
line = listin.readline()
while line:
	inputfiles.append(line.strip())
	line = listin.readline()
listin.close()

for file in inputfiles:
	fileopen = open(file, "r")
	line = fileopen.readline()
	while line:
		if line[0] == ">":
			if refname in line:
				line = fileopen.readline()
				seq = line.strip()
				seqlist = list(seq)
				seqstrip = seq.strip("-")
				try:
					beg = seqlist.index(seqstrip[0])
				except:
					beg = "-"
				revlist = list(reversed(seqlist))
				revstr = seq[::-1]
				revseqstrip = revstr.strip("-")
				try:
					revend = revstr.index(revseqstrip[0])
					end = len(seq) - revend
				except:
					end = "-"
				line = fileopen.readline()
			else:
				line = fileopen.readline()
		else:
				line = fileopen.readline()
	fileopen.close()

	try:
		os.mkdir(outdir)
	except:
		pass
	if beg == "-" and end == "-":
		fileout = open(outdir + "/" + file + ".nodata", "w")
		fileopen = open(file, "r")
		line = fileopen.readline()
		while line:
			fileout.write(line)
			line = fileopen.readline()
		fileout.close()

	else:
		fileout = open(outdir + "/" + file + ".trimmed", "w")
		if beg == 0:
			pass
		else:
			filehead = open(outdir + "/" + file + ".header", "w")
		if end == len(seq):
			pass
		else:
			filetail = open(outdir + "/" + file + ".tail", "w")
		fileopen = open(file, "r")
		line = fileopen.readline()
		while line:
			if line[0] == ">":
				fileout.write(line)
				if beg == 0:
					pass
				else:
					filehead.write(line)
				if end == len(seq):
					pass
				else:
					filetail.write(line)
			else:
				line = fileopen.readline()
			line = fileopen.readline()
			lineseq = line
			fileout.write(lineseq[beg:end] + "\n")
			if beg == 0:
				pass
			else:
				filehead.write(lineseq[0:beg] + "\n")
			if end == len(seq):
				pass
			else:
				filetail.write(lineseq[end:])
			line = fileopen.readline()
		fileout.close()
		if beg == 0:
			pass
		else:
			filehead.close()
		if end == len(seq):
			pass
		else:
			filetail.close()
		fileopen.close()
