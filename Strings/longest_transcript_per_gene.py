import sys

'''
Input FASTA on STDIN.
Assume input defline formatted like this:
>geneID otherstuff...
Output FASTA on STDOUT.
Output contains only the longest sequence per geneID.

This is designed for parsing ensemble lncRNA fasta with deflines like
>ENSG-ID ENST-ID otherstuff

Usage:
# cat infasta | longest_transcript_per_gene.py > outfasta
'''

deflines={}
sequences={}
lengths={}
gene=""
for rawline in sys.stdin:
    line=rawline.rstrip('\n')
    if line[0]=='>':
        defline=line
        fields=line.split(' ')
        gene=fields[0]
    else:
        sequence=line
        length=len(line)
        if gene not in lengths or lengths[gene] < length:
            deflines[gene]=defline
            sequences[gene]=sequence
            lengths[gene]=length
        gene=""

for gene in deflines.keys():
    defline = deflines[gene]
    sequence = sequences[gene]
    print(defline)
    print(sequence)





