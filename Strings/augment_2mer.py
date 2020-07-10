#!/usr/bin/env python
# coding: utf-8

'''
Input: 2mer count csv
Output: input plus added columns for co-occurence
'''
import sys
import pandas as pd

SUFFIX='.features.csv'
inprefix='ncRNA.2mer'
outprefix='ncRNA.2mer_co'
if len(sys.argv)==3:
    inprefix=sys.argv[1]
    outprefix=sys.argv[2]
else:
    print("Incorrect number of command-line arguments.")
    print("Using defaults.")
infile= inprefix + SUFFIX
outfile= outprefix + SUFFIX
print("in/out: "+infile+" / "+outfile)

pairs_of_interest=[]
# indicators of nc
pairs_of_interest.append('CA-CG')
pairs_of_interest.append('CC-CG')
pairs_of_interest.append('CG-AC')
pairs_of_interest.append('CG-AG')
pairs_of_interest.append('CG-GC')
pairs_of_interest.append('CG-CT')
pairs_of_interest.append('CG-GA')
pairs_of_interest.append('CG-GG')
pairs_of_interest.append('CG-GT')
pairs_of_interest.append('CG-TG')
pairs_of_interest.append('CG-AA')
pairs_of_interest.append('CG-AT')
pairs_of_interest.append('CG-TA')
pairs_of_interest.append('CG-TC')
pairs_of_interest.append('CG-TT')
pairs_of_interest.append('GC-GG')
# indicators of pc
pairs_of_interest.append('AA-AT')
pairs_of_interest.append('AA-CA')
pairs_of_interest.append('AA-CT')
pairs_of_interest.append('AA-TT')
pairs_of_interest.append('AC-TA')
pairs_of_interest.append('CC-TA')
pairs_of_interest.append('CT-TT')
pairs_of_interest.append('GA-TA')
pairs_of_interest.append('TA-AA')
pairs_of_interest.append('TA-AG')
pairs_of_interest.append('TA-AT')
pairs_of_interest.append('TA-CA')
pairs_of_interest.append('TA-CT')
pairs_of_interest.append('TA-TC')
pairs_of_interest.append('TA-TG')
pairs_of_interest.append('TA-TT')

all_seqs=[]
df2 = pd.read_csv (infile)
rows=df2.shape[0]
for r in range(rows):
    features_per_seq=[]
    for pair in pairs_of_interest:
        mer0=pair[:2]
        mer1=pair[3:]
        val0=df2.iloc[r].loc[mer0]
        val1=df2.iloc[r].loc[mer1]
        minval=min(val0,val1)
        features_per_seq.append(minval)
    all_seqs.append(features_per_seq)
df3=pd.DataFrame(all_seqscolumns=pairs_of_interest)
#dfm=df3.drop(columns=['seqnum','seqlen'])
#dfc=pd.concat([df2,dfm],axis='columns')
#dfc.to_csv(outfile)
