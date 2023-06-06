import sys
for line in sys.stdin:
    fields=line.strip().split('\t')
    ID = fields[2].split('.')[0]
    POS = fields[3]
    CIGAR = fields[5]
    MATE_POS = fields[7]
    INSERT_SIZE = fields[8]
    SCORE = '-1000'
    MISMATCHES = '0'
    for j in range(11,len(fields)):
        if fields[j].startswith('AS:i:'):
            SCORE = fields[j][5:] # align score
        if fields[j].startswith('XM:i:'):
            MISMATCHES = fields[j][5:] # mismatch count
    print(ID,CIGAR,SCORE,MISMATCHES,POS)
