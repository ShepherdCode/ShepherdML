#!/bin/sh

date
echo "Total read pairs in bam:"
samtools view Primary.bam | cut -f 1 | uniq | wc -l

date
echo "processing ..."
samtools view Primary.bam | python bam_to_ml_csv.py > ml_stats.csv
echo -n $?
echo " exit status"

date
echo "Total read pairs in cvs:"
wc -l ml_stats.csv

date
