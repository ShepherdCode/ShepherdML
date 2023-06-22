#!/bin/sh

SRC=/cluster/projects/nn9525k/hybrids/jasonrm/Arenosa/scripts

module --force purge
module load StdEnv 
module load GCC/11.3.0
module load Python/3.10.4-GCCcore-11.3.0
module load SAMtools/1.16.1-GCC-11.3.0

date
echo 'Extract lyrata read stats from its two BAM files'
python ${SRC}/bam_two_targets.py \
       map_lyrata_to_lyrata/Sorted.bam \
       map_lyrata_to_halleri/Sorted.bam \
       > lyrata_read_stats.csv
echo -n $?
echo " exit status"

date
echo 'Extract halleri read stats from its two BAM files'
python ${SRC}/bam_two_targets.py \
       map_halleri_to_lyrata/Sorted.bam \
       map_halleri_to_halleri/Sorted.bam \
       > halleri_read_stats.csv
echo -n $?
echo " exit status"

date
gzip -v *stats.csv

date
ls -l
