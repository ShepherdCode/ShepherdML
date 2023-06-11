#!/bin/sh

echo "If necessary, load the samtools module."
echo "If necessary, load the python module."

SRC="/cluster/projects/nn9525k/hybrids/jasonrm/Arenosa/scripts"

date
echo "Processing ..."
samtools view Primary.bam | python ${SRC}/bam_to_ml_csv.py > ml_stats.csv
echo -n $?
echo " exit status"

date
echo "Total read pairs in bam:"
samtools view Primary.bam | cut -f 1 | uniq | wc -l

date
echo "Total read pairs in cvs:"
wc -l ml_stats.csv

gzip -v ml_stats.csv

date
