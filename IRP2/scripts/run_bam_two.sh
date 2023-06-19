#!/bin/sh

SRC=/cluster/projects/nn9525k/hybrids/jasonrm/Arenosa/scripts

# Use --irp because we used FASTA files containing _MxM or _SxS in the defline.

date
python ${SRC}/bam_two_targets.py --irp \
       ../map_MxM_BR4_to_M/Sorted.bam \
       ../map_MxM_BR4_to_S/Sorted.bam \
       > ml_stats.csv
echo -n $?
echo " exit status"

gzip -v ml_stats.csv

date
ls -l
