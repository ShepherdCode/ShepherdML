#!/bin/sh

# Produced a balanced set of ALL transcripts per gene (minus the 1000 test set).

echo "How many protein coding sequences?"
cp ../long_and_short/pcRNA.gc36.all.fasta pcRNA.gc36.balance.fasta
total=$(grep -c '^>' pcRNA.gc36.balance.fasta)
echo $total
echo

echo "How many non-coding sequences initially?"
cp ../long_and_short/ncRNA.gc36.all.fasta ncRNA.tmp
grep -c '^>' ncRNA.tmp
echo

echo "How many sequences finally?"
cat ncRNA.tmp | \
    awk '{if(substr($1,1,1)==">") x=$1; else {print x "&" $1;}}' \
	> concat.tmp
shuf concat.tmp | head -n $total > random.tmp
cat random.tmp | tr '&' '\n' > ncRNA.gc36.balance.fasta
grep -c '^>' *.fasta
echo

echo "Cleanup..."
#rm -v *.tmp
