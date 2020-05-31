#!/bin/bash
###PBS -q debug
#PBS -q standby
###PBS -q testqueue
###PBS -q training
####PBS -l nodes=4:ppn=1
####PBS -l procs=4
#PBS -l nodes=1:ppn=1
#PBS -l walltime=4:00:00
###PBS -m ae   ## email onabort and exit
#PBS -m n   ## no email
#PBS -M jrm0122@mix.wvu.edu
#PBS -Nhum1
cd /users/jrm0122/myscratch/count_kmers
module list
echo -n "START AT "
date
python3 KmerCount.py ../GRCh38.oneline.fasta human --k 6 --fasta
echo -n "END AT "
date

