### modified from https://stackoverflow.com/questions/22886699/how-to-remove-all-n-sequence-entries-from-fasta-files

import sys
import Bio
from Bio import SeqIO

INPUT  = sys.argv[1]  # Input fasta file
OUTPUT = sys.argv[2] # output fasta file

def main():
    records = Bio.SeqIO.parse(INPUT, 'fasta')
    filtered = (rec for rec in records if any(ch != '-' for ch in rec.seq))
    Bio.SeqIO.write(filtered, OUTPUT, 'fasta')

if __name__=="__main__":
    main()



