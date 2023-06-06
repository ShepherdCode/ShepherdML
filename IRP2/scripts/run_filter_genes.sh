#!/bin/sh

# Purpose:
# Use homozygous sequencing to decide which genes are non-performant.
# For example, a gene that generates no reads should not be included in adjusted P-value computation.
# Or, a gene whose reads map to the wrong parent should not be used at all.

# Input:
# The model.tsv assigns each file of counts to a sample and replicate.
# This will help us filter by replicates per sample.

# Ouput:
# Three counts for each gene:
#    min reads in any one replicate
#    min reads in any one sequencing sample (some replicates were sequenced twice)
#    fold change defined as (true-false)/false where 1 is the minimum false value
#    i.e. 500 reads to true parent and 1 to false parent is about 500 fold
#    (but 5 reads to true parent and 0 to false parent is (5-1)/1= 4 fold)

# Set up your python environment first.
#SRC="/cluster/home/jasonrm/Source/MOLBAR/src"
SRC='.'

echo "Write counts to *.tsv ..."

# python3 ${SRC}/filter_homozygous.py model.tsv > min_and_fold.tsv
python3 ${SRC}/filter_genes.py model.tsv > min_and_fold.tsv

echo "Write genes to *.genes_pass_filter ..."
echo "Filter for min fold change 5"
# The counter c is used to skip the input header line

cat min_and_fold.tsv | awk '{if ((c++ >0) && ($4>=5)) print $1;}' > genes_pass_filter.txt

