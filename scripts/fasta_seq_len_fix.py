import sys
from Bio import SeqIO

infile=sys.argv[1]
outfile=sys.argv[2]

SeqIO.convert(infile, "fasta", outfile, "fasta")
