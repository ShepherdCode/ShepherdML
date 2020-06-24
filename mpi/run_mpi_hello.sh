#!/bin/bash
#PBS -q training
###PBS -q debug
###PBS -q standby
####PBS -q testqueue
#PBS -l nodes=1:ppn=8
#PBS -m ae   ##mail onabort and exit
#PBS -M jrm0122@mix.wvu.edu
#PBS -Njob1
cd /users/jrm0122/CS560
module load mpi/intel/5.1.3
module list
echo -n "START AT "
date
mpirun -np 8  ./mpi_hello
echo -n "FINISH AT "
date