#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>


int main(int argc, char *argv[])
{
  int me, nprocs;
  int t;
  int n= t = 0;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
  
  printf("Hi from node %d of %d\n", me, nprocs);
  
  MPI_Finalize();
  return 0;
}
