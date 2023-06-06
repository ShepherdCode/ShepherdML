infile='alignments.csv'
lenfile='CDS.length'
id_to_len=dict()
with open (lenfile,'r') as fin:
    for line in fin:
        line=line.strip()
        ID,LEN = line.split(' ')
        id_to_len[ID]=LEN
with open (infile,'r') as fin:
    for line in fin:
        line=line.strip()
        ID,CIGAR,SCORE,MISMATCH,POS = line.split(' ')
        LEN = int(id_to_len[ID])
        POS = int(POS)
        SCORE = int(SCORE)
        MISMATCH = int(MISMATCH)
        end_pos = max(0,LEN-POS-150)
        min_dist = min(POS,end_pos)
        print(ID,CIGAR,SCORE,MISMATCH,POS,end_pos,min_dist)
