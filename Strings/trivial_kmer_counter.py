import sys
import math

'''
Trivial k-mer counter for fasta files.
Intended output validation and time comparisons.
Usage: (for k=2)
    cat my.fasta | python3 trivial_kmer_counter.py 2
Input: stdin, fasta, exactly one line per sequence.
Output: stdout, csv, one line per input sequence,
    columns: seqnum,seqlen,counts...
The counts are in lexical order: AA, AC, AG, ...
'''

if len(sys.argv)>1:
    k=int(sys.argv[1])
else:
    k=2

values={}
values['A']=0
values['C']=1
values['G']=2
values['T']=3
seqnum=0
totalmers=pow(4,k)
for line in sys.stdin:
    line = line.rstrip()
    if line[0] != '>':
        kmercounts=[0]*totalmers
        seqnum += 1
        strlen = len(line)
        sys.stdout.write("%d,%d"%(seqnum,strlen))
        for i in range(strlen-k+1):
            value = 0
            power=pow(4,k-1)
            for j in range(k):
                chr = line[i+j]
                val = values[chr]
                value += power*val
                power = power // 4
            kmercounts[value]+=1
        for count in kmercounts:
            sys.stdout.write(",%d"%count)
        sys.stdout.write('\n')
