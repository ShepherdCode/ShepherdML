import sys
# parse outputs of samview_to_csv.sh
for line in sys.stdin:
    fields=line.strip().split('\t')
    AS = '-1000'
    XM = '0'
    for j in range(8,len(fields)):
        if fields[j].startswith('AS:i:'):
            AS = fields[j][5:] # align score
        if fields[j].startswith('XM:i:'):
            XM = fields[j][5:] # mismatch count
    print(fields[1],fields[2],fields[4],fields[6],AS,XM)
# 1 read name    N
# 2 bit flags    N
# 3 target name  Y
# 4 left pos     Y
# 5 mapq         N always 0
# 6 cigar        Y
# 7 mate or =    N always =
# 8 mate pos     Y
# 9 insert size  N
# 10 sequence    N
# 11 qualities   N
# AS align score Y
# XM mismatches  Y
