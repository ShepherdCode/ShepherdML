#!/usr/bin/env python
# coding: utf-8

# # Coding vs non-coding
# K-mer count features extracted from GenCode pc and nc RNA.
# Use GenCode 34.
# Use one RNA per gene; the one transcript with median length (use floor where count is even).

# Same process for protein coding and non-coding.
# * Start with GenCode 34 fasta file.
# * Run gencode_preprocess.py to make all caps and remove long seqs, short seqs, seqs with N. For each gene with multiple transcripts, choose the one transcript with median length. Among the remaining, remove any sequences that are duplicates of previous ones (We find up to 7 exact duplicate sequences per sequence identifier).
# * Run spot_dupes to make sure the dupes are gone.
# * Run fasta_to_feature.py to generate CSV file of K-mer counts.
# * Run trivial_kmer_counter.py on subset of sequences to verify the K-mer counting.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ncfile='ncRNA.2mer.features.csv'
nc_features = pd.read_csv(ncfile,header=0)

pcfile='pcRNA.2mer.features.csv'
pc_features = pd.read_csv(pcfile,header=0)

# ## Generate train set, test set
# Introduce the labels 0=non-coding, 1=protein-coding.
# We are worried that the longest sequences are a special case.
# Use stratified split to ensure an even split of train/test by sequence length.

# Manufacture labels for the two datasets
nc_labels_temp=[0]*nc_features.shape[0]
pc_labels_temp=[1]*pc_features.shape[0]
nc_labels=pd.core.frame.DataFrame(nc_labels_temp,columns=['label'])
pc_labels=pd.core.frame.DataFrame(pc_labels_temp,columns=['label'])
nc_all=pd.concat([nc_labels,nc_features],axis='columns')
pc_all=pd.concat([pc_labels,pc_features],axis='columns')
# And combine non-coding + protein-coding into one data structure.
all_instances=pd.concat([nc_all,pc_all])

def sizebin(df):
    return pd.cut(df["seqlen"],
                              bins=[0,1000,2000,4000,8000,16000,np.inf],
                              labels=[0,1,2,3,4,5])
bin_labels= sizebin(all_instances)
from sklearn.model_selection import StratifiedShuffleSplit
splitter = StratifiedShuffleSplit(n_splits=1, test_size=0.2, random_state=37863)
# split(x,y) expects that y is the labels.
# Trick: Instead of y, give it it the bin labels that we generated.
for train_index,test_index in splitter.split(all_instances,bin_labels):
    train_set = all_instances.iloc[train_index]
    test_set = all_instances.iloc[test_index]
train_set.shape

all_labels= sizebin(all_instances)
train_labels= sizebin(train_set)
tot_all=len(all_labels)
tot_train=len(train_labels)
print("cat all_num   all_pct     train_num   train_pct")
for i in range(6):
    # Using value_counts returns unique count, and anyway hits a recursion limit here.
    truth=all_labels.apply(lambda x: True if x==i else False)
    cnt1=len(truth[truth==True].index)
    truth=train_labels.apply(lambda x: True if x==i else False)
    cnt2=len(truth[truth==True].index)
    print("%3d %7d %10f vs %7d %10f"%(i,cnt1,cnt1/tot_all,cnt2,cnt2/tot_train))

# Move the seqnum and seqlen columns to a separate matrix.
# Move the labels column to a separate matrix.
X_train_ids=train_set[['seqnum','seqlen']].copy()
y_train=    train_set[['label']].copy()
X_train=    train_set.drop(columns=['label','seqnum','seqlen'])

X_test_ids= test_set[['seqnum','seqlen']].copy()
y_test=     test_set[['label']].copy()
X_test=     test_set.drop(columns=['label','seqnum','seqlen'])

X_train.to_pickle("ncRNA.pcRNA.X_train.pkl")
X_train_ids.to_pickle("ncRNA.pcRNA.X_train_ids.pkl")
y_train.to_pickle("ncRNA.pcRNA.y_train.pkl")

X_test.to_pickle("ncRNA.pcRNA.X_test.pkl")
X_test_ids.to_pickle("ncRNA.pcRNA.X_test_ids.pkl")
y_test.to_pickle("ncRNA.pcRNA.y_test.pkl")
