#!/usr/bin/env python
# coding: utf-8

'''
Combine K-mer counts for multiple K.
Yeah, our program should generate these.
'''
import sys
import pandas as pd

SUFFIX='.features.csv'
prefix2='pcRNA.2mer'
prefix3='pcRNA.3mer'
outfix='pcRNA.2-3'
if len(sys.argv)==4:
    prefix2=sys.argv[1]
    prefix3=sys.argv[2]
    outfix=sys.argv[3]
else:
    print("Incorrect number of command-line arguments.")
    print("Using defaults.")
infile2= prefix2 + SUFFIX
infile3= prefix3 + SUFFIX
outfile= outfix  + SUFFIX
print("inputs: "+infile2+", "+infile3)  

df2 = pd.read_csv (infile2)
df3 = pd.read_csv (infile3)
dfm=df3.drop(columns=['seqnum','seqlen'])
dfc=pd.concat([df2,dfm],axis='columns')
dfc.to_csv(outfile)

print("outputs: "+outfile)
