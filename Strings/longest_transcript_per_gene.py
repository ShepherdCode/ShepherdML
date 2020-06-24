import sys

saved={}
lengths={}
for rawline in sys.stdin:
    line=rawline.rstrip('\n')
    fields = line.split('\t')
    gene=fields[0]
    seq=fields[3]
    length=len(seq)
    if gene not in lengths or lengths[gene] < length:
        saved[gene]=line
        lengths[gene]=length

for line in saved.values():
    fields = line.split('\t')
    gene=fields[0][0:15]
    transcript=fields[1]
    biotype=fields[2]
    seq=fields[3]
    print (">"+gene+" "+transcript+" "+biotype)
    print (seq)





