import sys

sums=[]
for i in range(1,7):
    sums.append(0)

lines=[]
for line in sys.stdin:
    lines.append(line)
    fields=line.split(' ')
    for i in range(1,6):
        sums[i] = sums[i] + int(fields[i])

pct=[]
for i in range(1,7):
    pct.append(0)
for line in lines:
    fields=line.split(' ')
    outline=fields[0]
    for i in range(1,6):
        pct[i]=1.0*int(fields[i])/sums[i]
        spct="%.5f" % pct[i]
        outline=outline+" "+spct
    print(outline)



