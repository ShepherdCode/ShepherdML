#!/bin/sh

echo load module

module --force purge
module load StdEnv 

module load GCC/11.3.0
module load Bowtie2/2.4.5-GCC-11.3.0
module load Python/3.10.4-GCCcore-11.3.0
module load SAMtools/1.16.1-GCC-11.3.0

echo merge fasta

python make_diploid_fasta.py MxM.fasta --name1 MxM SxS.fasta --name2 SxS
# default output filename is diploid.fasta

echo build index

bowtie2-build diploid.fasta diploid

date
echo done
