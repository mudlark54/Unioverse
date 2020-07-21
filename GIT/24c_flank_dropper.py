#! /usr/bin/env python3

# Script was written by Jesse Breinholt (jessebreinholt@gamil.com) and Chandra Earl (sunray1@ufl.edu)
#This script takes a Anchored hybrid alignment and defines the probe region based on a seq in the file.
#For data outside of the probe (Head and Tail), it calculates a simple distance score from a 50% consensus
#sequence by counting mismatches (ignores gaps) and scales the distance by the number of nongap postions
#for each sequences. Then using the mean and the number of standard deviation of the scores it removes
#leading and tailing data that have scores above the given standard deviation above the mean score.
#Both head and tail are allowed different stdev cut offs. If the head or tail region do not have
#any large outliers, the seqs with the most differences are cut off, even though there is likely valid
#data. Therefore I suggest you look at the data before and after a cut off to see
#if you need to set the tail or head region to a much higher cut off to avoid removing good data.
#Example: ./flank_dropper.py INfile OUTfile BMORI 2 2




from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import AlignInfo, MultipleSeqAlignment
from Bio.Alphabet.IUPAC import ambiguous_dna
from statistics import mean, stdev as STD
import sys

if "-h" in sys.argv:
	sys.exit("Usage: flank_dropper.py INfile OUTfile taxaNameForProbe StdevHead StdevTail\nExample: ./flank_dropper.py INfile OUTfile BMORI 2 2")

try:
	inFile = sys.argv[1]
	outFile = sys.argv[2]
	PROBEname = sys.argv[3]	
	cutoff = float(sys.argv[4])
	cutoff_tail = float(sys.argv[5])
	alignIN = AlignIO.read(inFile,'fasta',alphabet=ambiguous_dna)
except Exception: 
	print("Usage: flank_dropper.py INfile OUTfile taxaNameForProbe StdevHead StdevTail\nExample: ./flank_dropper.py INfile OUTfile BMORI 2 2")
	sys.exit(1)
	


colsL = alignIN.get_alignment_length()
seqsL = len(alignIN)
start = 0
end = 0
count=0
score=0
datac=0
score2=0
datac2=0
hfix=0
tfix=0
#empty dictionary
sdict={}
sdict2={}


#go through seqs, find probe taxa, record start and end
for rec in alignIN:
	ID = rec.id
	Seq = rec.seq
	SEQ = Seq.upper()
	if PROBEname in ID:
		#print(ID)
		start = min( SEQ.find("A"), SEQ.find("C"), SEQ.find("G"), SEQ.find("T"))
		end = max( SEQ.rfind("A"), SEQ.rfind("C"), SEQ.rfind("G"), SEQ.rfind("T"))
		break
#print("Probe region- ID:", ID, "start:",start,"end:", end,"\n")

#make a 50% consesus seq
info = AlignInfo.SummaryInfo(alignIN)
consensus = info.dumb_consensus(threshold=.50,consensus_alpha=ambiguous_dna)
#print("Consensus:",consensus,"\n")
#consensus for the head and tail region used to get number of gaps needed to insert below
headcon=consensus[:start+1]
tailcon=consensus[end:]
colsH = len(headcon)
colsT = len(tailcon)

#loop through each column by record, record scores for head and tail regions and add to dictionary
for rec in alignIN:
	Seq = rec.seq
	ID = rec.id
	for colIDX in range(colsL):
		col = alignIN[count:count+1,colIDX]
		if colIDX < start:
			if col != consensus[colIDX] and col != "-":
				score+=1
			if col != "-":
				datac+=1
			if score==0:
				sdict[ID]=score
			if score>=1:
				sdict[ID]=score/datac
		if colIDX < start and colIDX > end:
			pass
		if colIDX > end:
			if col != consensus[colIDX] and col != "-":	
				score2+=1
			if col != "-":
				datac2+=1
			if score2==0:
				sdict2[ID]=score2
			if score2>=1:
				sdict2[ID]=score2/datac2
#reset scores and data counts and increases count to loop through the taxa
	score=0
	datac=0
	score2=0
	datac2=0
	count+=1


#calculate mean and stdev and cut off for head and tail
#do not calculate the mean or cut off for head or tail if there is no data outside of probe region
if int(colsH) > 1:
	hmean = mean(sdict.values())
	hCUT = hmean + (cutoff * STD(sdict.values()))
	Hstdev=STD(sdict.values())
	#print("Head scores-  mean:",hmean,"cut off given input stdev:",hCUT,"Taxa scores: ", sdict,"\n\n")

if int(colsT) > 1:
	tmean = mean(sdict2.values())
	tCUT = tmean + (cutoff_tail * STD(sdict2.values()))
	Tstdev=STD(sdict2.values())
	#print("Tail scores-  mean:",tmean,"cut off given input stdev:",tCUT,"Taxa scores: ", sdict2,"\n\n")



# edit the alignment by dropping and replacing seqs for each record

for rec in alignIN:
	Seq = rec.seq
	ID = rec.id
	
	if int(colsH) >= 25 and int(colsT) >= 25: # if head and tail 25 or above process
		if sdict[ID] > hCUT and sdict2[ID] > tCUT:              #if the tail and head are ugly
			alignIN = alignIN[1:, :]                            #drop first seq
			fixedseq = SeqRecord("-"*(int(colsH)-1)+Seq[start:end+1]+"-"*(int(colsT)-1), id=ID, description='')
			# make a fixed seq by replacing ugly head and tail data with gaps
			alignIN.append(fixedseq)
			#append fixed seq to aligment
			hfix+=1
			tfix+=1
		if sdict[ID] > hCUT and sdict2[ID] <= tCUT:             #if head is ugly and tail is fine
			alignIN = alignIN[1:, :] 
			fixedseq = SeqRecord("-"*(int(colsH)-1)+Seq[start:], id=ID, description='')
			alignIN.append(fixedseq)	
			hfix+=1
		if sdict[ID] <= hCUT and sdict2[ID] > tCUT:             #if tail is ugly and head is fine
			alignIN = alignIN[1:, :] 
			fixedseq = SeqRecord(Seq[:end+1]+"-"*(int(colsT)-1), id=ID, description='')
			alignIN.append(fixedseq)	
			tfix+=1
		if sdict[ID] <= hCUT and sdict2[ID] <= tCUT:             #if both tail and head are fine
			alignIN = alignIN[1:, :]
			fixedseq = SeqRecord(Seq[:], id=ID, description='')
			alignIN.append(fixedseq)


	if colsH >= 25 and colsT < 25: # if head above 25 and tail below process head only
		if sdict[ID] > hCUT:             #if head is ugly and tail is fine
			alignIN = alignIN[1:, :] 
			fixedseq = SeqRecord("-"*(int(colsH)-1)+Seq[start:], id=ID, description='')
			alignIN.append(fixedseq)	
			hfix+=1
		if sdict[ID] <= hCUT:             #if both tail and head are fine
			alignIN = alignIN[1:, :]
			fixedseq = SeqRecord(Seq[:], id=ID, description='')
			alignIN.append(fixedseq)
	
	if colsH < 25 and colsT >= 25: # if tail above 25 and head below process head only
		if sdict2[ID] > tCUT:             #if tail is ugly and head is fine
			alignIN = alignIN[1:, :] 
			fixedseq = SeqRecord(Seq[:end+1]+"-"*(int(colsT)-1), id=ID, description='')
			alignIN.append(fixedseq)	
			tfix+=1
		if sdict2[ID] <= tCUT:             #if both tail and head are fine
			alignIN = alignIN[1:, :]
			fixedseq = SeqRecord(Seq[:], id=ID, description='')
			alignIN.append(fixedseq)
	if colsH < 25 and colsT < 25: # if both below 25 output original alignment
		alignIN = alignIN[1:, :]
		fixedseq = SeqRecord(Seq[:], id=ID, description='')
		alignIN.append(fixedseq)




#print("Number of taxa the leading data was replaced with gaps: ",hfix) 
#print("Number of taxa the taling data was replaced with gaps: ",tfix)

# account for when there is no head or tail; and Hstdev or Tstdev have not been made yet
try:
    Hstdev
except NameError:
    Hstdev = str(0)
try:
    Tstdev
except NameError:
   Tstdev = str(0)
#print("FileName\tHeadLength\tTailLength\tHeadStdev\tTailStdev\t#HeadRemoved\t#TailRemoved")
print(inFile, "\t", str(int(colsH)-1), "\t", str(int(colsT)-1), "\t", Hstdev, "\t", Tstdev, "\t", hfix, "\t", tfix)
#write aligmnet
AlignIO.write(alignIN, outFile, "fasta")
