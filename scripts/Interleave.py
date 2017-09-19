#This Python script requires Biopython 1.51 or later
from Bio import SeqIO
import itertools
import argparse

parser = argparse.ArgumentParser()
    
parser.add_argument("file1", help="fastq R1 file")
    
parser.add_argument("file2", help="fastq R2 file")

parser.add_argument("file_out", help="out file")
    
args = parser.parse_args()

format = "fastq" #or "fastq-illumina", or "fasta", or ...

def interleave(iter1, iter2) :
    for (forward, reverse) in itertools.izip(iter1,iter2):
        assert forward.id == reverse.id
        forward.id += "/1"
        reverse.id += "/2"
        yield forward
        yield reverse

records_f = SeqIO.parse(open(args.file1,"rU"), format)
records_r = SeqIO.parse(open(args.file2,"rU"), format)

handle = open(args.file_out, "w")
count = SeqIO.write(interleave(records_f, records_r), handle, format)
handle.close()
print "%i records written to %s" % (count, args.file_out)
