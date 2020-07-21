#!/usr/bin/env python
# -*- coding: utf-8 -*-

#script follows suggestion of Eric Normandeau on Biostars.org https://www.biostars.org/p/2822/ 

import sys
from Bio import SeqIO
arguments = sys.argv
if '-h' in arguments:
    sys.exit("./removelist.py input listtoget output")
elif len(arguments) < 4:
	sys.exit("./removelist.py input listtoget output")

fasta_file = sys.argv[1]  # Input fasta file
number_file = sys.argv[2] # Input interesting numbers file, one per line
result_file = sys.argv[3] # Output fasta file

remove = set()
with open(number_file) as f:
    for line in f:
        line = line.strip()
        if line != "":
            remove.add(line)
fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
end = False
with open(result_file, "w") as f:
    for seq in fasta_sequences:
        if seq.id in remove:
            SeqIO.write([seq], f, "fasta")