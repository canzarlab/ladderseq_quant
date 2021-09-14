from Bio import SeqIO
import sys

fqs = [SeqIO.parse(f, "fastq") for f in sys.argv[1:]]
while True:
    for fq in fqs:
        try:
            print(next(fq).format("fastq"))
        except StopIteration:
            fqs.remove(fq)
    if len(fqs) == 0:
        break
