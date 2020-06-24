import sys

genes={}

idfile=sys.argv[1]
with open(idfile,"r") as f:
    for rawline in f:
        line=rawline.rstrip('\n')
        fields = line.split('\t')
        gene = fields[0]
        genes[gene]=0

good = False
for rawline in sys.stdin:
    line=rawline.rstrip('\n')
    if line[0]=='>':
        good = False
        gene = line[1:16]
        if gene in genes and genes[gene]==0:
            print(line)
            good = True
            genes[gene] = 1
    else:
        if good:
            print(line)

with open("notfound","w") as f:
    for gene in genes:
        if genes[gene]==0:
            print(gene,file=f)
