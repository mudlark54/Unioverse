#! /usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import AlignInfo, MultipleSeqAlignment
from Bio.Alphabet.IUPAC import ambiguous_dna
from Bio.SubsMat import FreqTable
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from statistics import mean, stdev as STD
from math import log2 as log2
import math
import sys

if "-h" in sys.argv:
	sys.exit("Usage: python3 alignment_DE_trim.py INfile OUTfile %density %saturation graph-1/0 taxaforprobe\nExample: python3 alignment_DE_trim.py L1.fa L1_DEtrim.fa 60 1.5 1 Bmori\n\"taxaforprobe\" is optional if included the sequence with the name provided in the alignment will be used to define a probe region and density and entrpy filter will not trim inside of this region")

try:
	inFile = sys.argv[1]
	outFile = sys.argv[2]
	pctDen = int(sys.argv[3])/100
	pctSat = float(sys.argv[4])#/100
	graph = int(sys.argv[5])
	alignIN = AlignIO.read(inFile,'fasta',alphabet=ambiguous_dna)
except Exception: 
	print("Usage: python3 alignment_DE_trim.py INfile OUTfile %density %saturation graph-1/0 taxaforprobe\nExample: python3 alignment_DE_trim.py L1.fa L1_DEtrim.fa 60 1.5 1 Bmori\n\"taxaforprobe\" is optional if included the sequence with the name provided in the alignment will be used to define a probe region and density and entrpy filter will not trim inside of this region ")
	sys.exit(1)
try:
	PROBEname = sys.argv[6]
except Exception: 
	print("No PROBE sequence processing")
	PROBE = False
else:
	PROBE = True
	

colsL = alignIN.get_alignment_length()
seqsL = len(alignIN)
colsX = []
Density = []
Cons = []
Is = []
EN = []
recid=[]
start = 0
end = 0
belowSTCUT = 0
passD = 0

for rec in alignIN:
	ID = rec.id
	if (PROBE and not ID.find(PROBEname) >= 0) :
		recid.append(0)
	elif (PROBE and ID.find(PROBEname) >= 0):
		recid.append(1)
	else:
		recid.append(0)


if (1 not in recid and PROBE):
	print("\n#############################\nNo sequence with probe name given found in alignment:\n\tWill not save columns in probe region\n#############################\n\n")
	sys.exit()
	
for rec in alignIN:
	ID = rec.id
	#print(ID)
	Seq = rec.seq
	SEQ = Seq.upper()
	#print(ID.find(PROBEname))
	if (PROBE and ID.find(PROBEname) >= 0) :
#		print(ID)
		start = min( SEQ.find("A"), SEQ.find("C"), SEQ.find("G"), SEQ.find("T"))
		end = max( SEQ.rfind("A"), SEQ.rfind("C"), SEQ.rfind("G"), SEQ.rfind("T"))
		print("Probe region for ", ID, start, end)
		continue




for colIDX in range(colsL):
	col = alignIN[:,colIDX]
	col = col.upper()
	cnts = [
	col.count('A'), #0
	col.count('C'), #1
	col.count('G'), #2
	col.count('T'), #3
	col.count('-'), #4
	col.count('R'), #5
	col.count('Y'), #6
	col.count('W'), #7
	col.count('S'), #8
	col.count('K'), #9
	col.count('M'), #10
	col.count('D'), #11
	col.count('V'), #12
	col.count('H'), #13
	col.count('B'), #14
	col.count('X'), #15
	col.count('N'), #16
	]
	density = (
		cnts[0] + cnts[1] + cnts[2]  + cnts[3] + 
		(cnts[5] + cnts[6]  + cnts[7] + cnts[8] + cnts[9] + cnts[10])/2 +
		(cnts[11] + cnts[12] + cnts[13] + cnts[14])/3 +
		(cnts[15] +cnts[16])/4
		)
	Density.append(density)
	totACGT = cnts[0] +cnts[1] + cnts[2] +cnts[3]
	Pi = [cnts[i]/totACGT for i in range(4) if totACGT != 0]
	entropy = -sum([(f*log2(f)) for f in Pi if f != 0])
	EN.append(entropy)

# new alignment for selecting columns to save
alignN = alignIN[:,1:0]
if PROBE:
	for colIDX in range(colsL):
		if(start <= colIDX <= end):
			alignN = alignN[:] + alignIN[:,colIDX:colIDX+1]
		
		elif (Density[colIDX]/seqsL >= pctDen) :
			passD += 1
			if (EN[colIDX] <= pctSat) :  #save the column by saturation
				belowSTCUT += 1
				alignN = alignN[:] + alignIN[:,colIDX:colIDX+1]
		continue


if not PROBE:
	for colIDX in range(colsL):
		if (Density[colIDX]/seqsL >= pctDen) :
			passD += 1
			if (EN[colIDX] <= pctSat) :  #save the column by saturation
				belowSTCUT += 1
				alignN = alignN[:] + alignIN[:,colIDX:colIDX+1]
		continue
		
ABOVEIS=100*(belowSTCUT/passD)
print("Number of columns that pass Density\t", str(passD))
print("% of columns that pass Entropy cutoff " + str(ABOVEIS))
AlignIO.write(alignN[:,1:], outFile, "fasta")


print("Old alignment ", len(alignIN), "sequences in alignment of length ", colsL , ".")
colsN = alignN.get_alignment_length()
print("New alignment ", len(alignN), "sequences in alignment of length ", colsN-1 , ".")
print("New alignment", outFile, "written.")



#PLOT	
plt.figure(num=sys.argv[0])
plt.suptitle(inFile, fontsize=12, x=.40)
ix = [x for x in range(colsL)]
a = [(m/seqsL)*100 for m in Density]
c = Cons
e= EN
plt.subplot(311)
plt.title('Density')
plt.ylabel('Nuc. Density')
plt.plot(ix,a,'k-')
plt.legend(["Red Shade Excluded"], fontsize='small')
plt.axhline(y=pctDen*100, color='r', lw=2, ls='--')
if PROBE:
	plt.axvspan(start,end, color='#87CEEB', alpha=.5)
plt.axhspan(pctDen*100, 0, color='r', alpha=.1)



plt.subplot(312)
plt.title('entropy')
plt.ylabel('entropy')
plt.ylim(0,2)
plt.plot(ix,EN,'k-')
plt.axhline(y=pctSat, color='#808080', lw=2, ls='--')
if PROBE:
	plt.axvspan(start,end, color='#87CEEB', alpha=.5)
plt.axhspan(pctSat, 2, color='r', alpha=.1)


if PROBE:
	plt.axvspan(start,end, color='#87CEEB', alpha=.5)
plt.tight_layout(pad=1.0, h_pad=1.2, w_pad=1.2)
if graph == 1:
	plt.savefig(outFile + ".pdf")
sys.exit()
