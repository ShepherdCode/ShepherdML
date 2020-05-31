helpdesk@hpc.wvu.edu


#include <stdio.h>
#include "mpi.h"
int main( argc, argv )
int  argc;
char **argv;
{
    int rank, size;
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &size );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    printf( "Hello world from process %d of %d\n", rank, size );
    MPI_Finalize();
    return 0;
}

module load mpi/intel/5.1.3
mpicc -o mpi_hello  mpi_hello.c

cat runner.sh
#!/bin/bash
#PBS -q training
#PBS -l nodes=1:ppn=8
#PBS -m ae   ## mail on abort and exit
#PBS -M jrm0122@mix.wvu.edu
#PBS -N job1
cd /users/jrm0122/CS560
mpirun -np 8  ./mpi_hello

./runner.sh
Hello world from process 1 of 8
Hello world from process 2 of 8
Hello world from process 3 of 8
Hello world from process 4 of 8
Hello world from process 5 of 8
Hello world from process 6 of 8
Hello world from process 7 of 8
Hello world from process 0 of 8

qsub ./runner.sh
qsub: submit error (Unauthorized Request  MSG=group ACL is not satisfied: user jrm0122@srih0001.hpc.wvu.edu, queue training)

It seems we don't have access to the training queue.

Try -q standby ...
qstat -u jrm0122
srih0001.hpc.wvu.edu: 
                                                                                  Req'd       Req'd       Elap
Job ID                  Username    Queue    Jobname          SessID  NDS   TSK   Memory      Time    S   Time
----------------------- ----------- -------- ---------------- ------ ----- ------ --------- --------- - ---------
4761771.srih0001.hpc.w  jrm0122     standby  job1                --      1      8       3gb  04:00:00 Q       -- 

qstat -u jrm0122  # this will show status
      	     C - Completed (my guess)
      	     E - Job is exiting after having run.
             H - Job is held.
             Q - job is queued, eligable to run or routed.
	     R - job is	running.
       	     T - job is being moved to new location.
	     W - job is waiting	for its execution time (-a option) to be reached.
	     S - (Unicos only) job is suspend.

checkjob 4761772  # only valid once job is running

# cancel is deprecated. Use this instead.
mjobctl -c    4782305
