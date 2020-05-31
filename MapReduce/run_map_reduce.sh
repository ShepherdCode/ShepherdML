#!/bin/bash
###PBS -q debug
#PBS -q standby
###PBS -q testqueue
###PBS -q training
####PBS -l nodes=4:ppn=1
####PBS -l procs=4
#PBS -l nodes=1:ppn=4
#PBS -l walltime=00:02:00
###PBS -m ae   ## email onabort and exit
#PBS -m n   ## no email
#PBS -M jrm0122@mix.wvu.edu
#PBS -Njob6
cd /users/jrm0122/CS560
module load mpi/intel/5.1.3
module list
echo -n "START AT "
date
mpirun -np 4  ./map_reduce
echo -n "FINISH AT "
date