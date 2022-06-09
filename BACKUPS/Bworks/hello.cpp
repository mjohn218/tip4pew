#include <iostream>
#include <mpi++.h>
#include <fstream>


using namespace std;
int main(int argc, char **argv)
{
  int me, nprocs;
  int t;
  int n= t = 0;
  MPI::Init(argc, argv);
  me=MPI::COMM_WORLD.Get_rank();
  nprocs=MPI::COMM_WORLD.Get_size();
  
  cout << "Hi from node "<< me << " of " <<  nprocs<< '\n';
  
  MPI::Finalize();
  return 0;
}
