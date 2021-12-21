import sys

args = sys.argv
genomeseq = args[1]

sequence = open(genomeseq)
for contig in sequence:
    seqcontent = contig.split()
    print(">" + seqcontent[0])
    print(seqcontent[1])
    print("//")


