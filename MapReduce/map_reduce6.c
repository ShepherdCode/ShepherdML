// Use map/reduce to compute the max in an array of positive random integers.
// Use MPI to divide up the work. Use MPI scatter and MPI gather.
//
// This tutorial was helpful for parameters to MPI_Scatter & MPI_Gather:
//  https://mpitutorial.com/tutorials/mpi-scatter-gather-and-allgather/
//
// My program computes the global array without MPI before starting.
// This allows us to validate the result computed with MPI.
//
// My program handles 3 cases:
// Case 1: mpi_size > array_size
//    The extra workers waste time on uninitialized memory but
//    we only harvest the results from the first array_size workers.
// Case 2: array_size % mpi_size = 0
//    Nothing special needs to be done.
// Case 3: array_size % mpi_size = mod > 0
//    The last mod numbers are left out of the MPI computation.
//    The master adds them to the end of the results array.
//    This is not ideal because the master runs longer than the others.
//    We have not investigate MPI_Scatterv which might offer a better way.
//
// The user can select either of two data sets at compile time.
//    The sequential data is handy for debugging. Every array[i]=i.
//    The random data is required for the homework assignment.

#include <limits.h>
#include <stdio.h>
#include <time.h>
#include <malloc.h>
#include "mpi.h"
const int ARRAY_SIZE = 11;
const int MPI_MASTER = 0;
int* input_data_buffer;
int* input_subset;
int* results_data_buffer;
int mpi_rank, mpi_size;
//int DATA_SET = 0;  // 0=sequential
int DATA_SET = 1;  // 1=random

int computeMax (int* arry, int num, int start) {
  int i;
  int m = INT_MIN;
  for (i=0; i<num; i++) {
    if (i >= start) {
      //printf("Compare [%d]=%d to %d.\n",i,arry[i],m);
      if (arry[i] > m) {
	m = arry[i];
      }
    }
  }
  return m;
}

int describe_array ( int* arry, int num, int start) {
  char* values = malloc (10000);
  int i,written;
  char* offset = values;
  for (i=0; i<num; i++) {
    if (i<start) {
      written=sprintf (offset, " (%d)=%d", i, arry[i]);
    } else {
      written=sprintf (offset, " [%d]=%d", i, arry[i]);
    }
    offset += written;
  }
  printf ("Worker %d: Array size %d %s\n", mpi_rank, num, values);
  free ( values );
}

int main( int argc, char** argv ) {
  int hidden_global_max;
  int computed_global_max;
  int computed_local_max;
  int NUM_VALUES_RETURNED=1;  // each worker returns one number
  int chunk_size = 0;
  int start_leftovers = 0;
  int num_leftovers = 0;
  int included;
  //input_data_buffer = NULL;
  //input_subset = NULL;
  //results_data_buffer = NULL;

  // MPI SETUP
  MPI_Init( &argc, &argv );
  MPI_Comm_size( MPI_COMM_WORLD, &mpi_size );
  MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );
  if (mpi_size > ARRAY_SIZE) {
    // pathological case
    chunk_size = 1;
    included = ARRAY_SIZE;
    start_leftovers = ARRAY_SIZE;
  } else {
    chunk_size = ARRAY_SIZE / mpi_size;  
    num_leftovers = ARRAY_SIZE % mpi_size;
    start_leftovers = chunk_size * mpi_size;
    included = mpi_size + num_leftovers;
  }
  printf("Numbers=%d Workers=%d ChunkSize=%d NumLeftovers=%d StartLeftovers=%d\n",
	 ARRAY_SIZE,mpi_size,chunk_size,num_leftovers,start_leftovers);

  // MASTER SETUP
  input_subset = malloc ( chunk_size * sizeof(MPI_INT));
  if (input_subset == NULL) {
    printf ("No memory\n");
    return (2);
  }
  if (mpi_rank==MPI_MASTER) {
    input_data_buffer = malloc 
      ( ARRAY_SIZE * sizeof(MPI_INT));
    results_data_buffer = malloc 
      ( (included) * sizeof(MPI_INT) * NUM_VALUES_RETURNED);
    if (input_data_buffer == NULL || 
	results_data_buffer == NULL) {
      printf ("No memory\n");
      return (3);
    }
    array_setup(ARRAY_SIZE, DATA_SET);
    hidden_global_max = computeMax 
      ( input_data_buffer, ARRAY_SIZE, 0);
    describe_array 
      ( input_data_buffer, ARRAY_SIZE, 0 );
    printf("MASTER SAYS: ArraySize %d, True max = %d.\n", 
	   ARRAY_SIZE, hidden_global_max);
  }

  // WORKER COMPUTATION
  MPI_Scatter(input_data_buffer,chunk_size,MPI_INT,
	      input_subset,chunk_size,MPI_INT,
	      MPI_MASTER,MPI_COMM_WORLD);
  describe_array 
    ( input_subset, chunk_size, 0 );
  computed_local_max = computeMax 
    ( input_subset, chunk_size, 0 );
  printf("Worker %d says: local max = %d.\n", mpi_rank, computed_local_max);   
  MPI_Gather(&computed_local_max,NUM_VALUES_RETURNED,MPI_INT,
	     results_data_buffer,NUM_VALUES_RETURNED,MPI_INT,
	     MPI_MASTER,MPI_COMM_WORLD);

  // MASTER WRAP UP
  // TO DO: include leftovers
  int result = 0;
  if (mpi_rank==MPI_MASTER) {
    int i;
    for (i=0; i<num_leftovers; i++) {
      results_data_buffer[i+mpi_size]=input_data_buffer[i+start_leftovers];
    }
    describe_array
      ( results_data_buffer, included, 0 );
    computed_global_max = computeMax 
      ( results_data_buffer, included, 0 );
    printf("MASTER says: global max = %d.\n", computed_global_max); 
    if (computed_global_max == hidden_global_max) {
      printf("SUCCESS!\n");
      result = 0;
    } else {
      printf("FAILURE!\n");
      result = 1;
    }
  }
  MPI_Finalize();
  return result;
}

int array_setup(int size, int data_set) {
  int i;
  srand(time(0));
  for (i=0; i< size; i++) {
    if (data_set == 0) {
      input_data_buffer[i] = i;   // good for debugging
    } else {
      input_data_buffer[i] = rand();   // values 0 to RAND_MAX
    }
  }
}


