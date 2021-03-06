/*
  Program: Tip4p Ewald in cpp
  Programmer: Margaret Johnson
  Start: 1/10/2005
  This is a program to run molecular dynamics simulations given initial 
congfigurations of solution volumes, i.e. water. To run in NPT, NVT, NVE 
ensembles. Using Frenkel and Smit Molecular simulation book, and fortran
code from THG updated by Rajesh  Murarka. Also code from Elizabeth Verschell



*/

#include <mpi++.h>
#include <fstream>
#include <ctime>

#include <math.h>
#include <cmath>
#include <complex>
#include "Ensemble.h"
#include <iomanip>
#include "Parameters.h"

#include <unistd.h>

using namespace std;
/*Define global variables Natom is 3 for water, nsite is 4 [added fictitious 
site]*/

//Functional routines, integration, force, energy calculations, 
//specific code is in Ensemble.cpp

void NVEintegrateMotion(Ensemble &NVE, Constants cons, Parameters parms, WaterSend *mols, WaterSend *local_mols, MPI::Datatype &Watertype, int t);

/*Global Variables*/

int N_Water=512;
int num_local;

int main(int argc, char* argv[])
{
  double mstart, mend;
  mstart=clock();

  /*Build datatype describing class*/
  MPI::Datatype Watertype;
  MPI::Datatype types[3]={MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
  MPI::Aint disp[3];
  MPI::Aint lb, extent;

  
  MPI::Init(argc, argv);
  
  int rank=MPI::COMM_WORLD.Get_rank();
  int nprocs=MPI::COMM_WORLD.Get_size();
    
  cout << "My rank: "<< rank <<endl;
  Ensemble NVE;
  Parameters parms;
  Constants cons;
  Averages avg;
   
  /*Create a datatype to pass around that is a class. Includes
   both positions of all atoms and all Kvectors.*/
  
  int blocklen[3]={4, 4, 4};
  
  /*Number of mols per processor*/
  cons.NperP=N_Water/nprocs;
  WaterSend * mols= new WaterSend[N_Water];
  WaterSend * local_mols = new WaterSend[cons.NperP];
  cout << "nperP: "<<cons.NperP <<endl;
  

  /*Compute displacements of class components*/
  /*
  disp[0]=MPI::Get_address(&mols[0].x);
  disp[1]=MPI::Get_address(&mols[0].y);
  disp[2]=MPI::Get_address(&mols[0].z);
  */
  disp[0]=0;
  disp[1]=32;
  disp[2]=62;
  Watertype=Watertype.Create_struct( 3, blocklen, disp, types);
  int size=Watertype.Get_size();
  cout <<" Size is : " << size <<endl;
  Watertype.Get_extent(lb, extent);
  cout <<"lb is: " << lb <<endl;
  Watertype.Commit();

  MPI::Datatype Kv;
    
  srand(static_cast<unsigned>(time(0)));

  /*All initial positions read in on Processor 0 as well as all parameters*/
  //Thermostat therm;
  //Barostat baro;
  

  /*Timing parms*/

  int t;


  cout.precision(16);

  //Read in the data coordinates and simulation parameters.
  /*All processors read in initial coordinates and parameters
   coordinates go into WaterSend mols object.Then split them up into a local
  set*/
  char * PARMFILE=argv[1];
  char * COORDFILE=argv[2];
  if(!argv[1]) cerr << "no parameter file \n"; 
  if(!argv[2]) cerr << "no coord file \n";


  NVE.readParms(PARMFILE, cons, parms);
  NVE.readCoords(COORDFILE, mols);
  
  cout << "Read in molecules: " <<mols[0].x[0]<<endl;
  cout <<"rhomm: " <<parms.Rhomm<<endl;

  cons.N=N_Water;
  cons.rank=rank;

  if(cons.N != 512)
    {
      cout <<"arrays not initialiazed to correct size!" <<endl;
      return(1);
    }

  /*
    Initialize variables, water molecule constants and Tip4p-Ewald
    parameters 
  */
  
  NVE.InitialParms(parms, cons, avg);

  NVE.RigidMol(parms, cons, mols);

  
  /*
    Make local molecules point to a part of the main array
  */
  
  local_mols=&mols[rank*cons.NperP];
  cout << "Local molecule start: "<<local_mols[0].x[0]<<endl;
  cout << "Local molecule fourth: "<<local_mols[0].x[3]<<endl;

  if(rank==0)
    {
      char *filename="StartCoord4.out";

      NVE.printAllCoords(filename, cons, mols);
    }

  /*Pull all molecules assigned locally into mols again*/
    MPI::COMM_WORLD.Allgather(local_mols, cons.NperP, Watertype, mols, cons.NperP, Watertype);

    
  if(rank==0)
    {
      char *filename="GatherCoord4.out";

      NVE.printAllCoords(filename, cons, mols);
    }

  
  /*
    Calculate initial force and energy 
  */
  
  NVE.calcEForce(NVE, cons, parms, mols, local_mols);

  if(rank==0)
    {
      char * ffile="StartForce.out";
      NVE.writeAllForce(ffile, cons);
      
    }
  /*
    For new coordinates set new velocities from Maxwell Boltzman, then 
    remove components along normal to bonds
  */

  WaterSend * InitVel=new WaterSend[N_Water];
    
  if(cons.Equil==1)
    {
      if(rank==0)
	{
	  NVE.Maxwell(cons, parms);
	  NVE.VelocityNormal(cons, parms, mols);
	}    
    }
  
  /*
    Initial Kinetic Energy and Temp, Scale Velocity scales  velocities
    to set temperature and evaluates Kinetic Energy
  */

  /*Read in velocities from file so same for all tests, only P0 reads
   in velocities, and then distributes values to all other processors
  */

  char * velfile="velocity.inp";
  NVE.readVelocity(velfile, InitVel);
  
  int i;
  for(i=0;i<cons.NperP;i++)
    {
      NVE.vx[i][0]=InitVel[i+rank*cons.NperP].x[0];
      NVE.vx[i][1]=InitVel[i+rank*cons.NperP].x[1];
      NVE.vx[i][2]=InitVel[i+rank*cons.NperP].x[2];
      NVE.vy[i][0]=InitVel[i+rank*cons.NperP].y[0];
      NVE.vy[i][1]=InitVel[i+rank*cons.NperP].y[1];
      NVE.vy[i][2]=InitVel[i+rank*cons.NperP].y[2];
      NVE.vz[i][0]=InitVel[i+rank*cons.NperP].z[0];
      NVE.vz[i][1]=InitVel[i+rank*cons.NperP].z[1];
      NVE.vz[i][2]=InitVel[i+rank*cons.NperP].z[2];

    }
  delete[] InitVel;
  //  NVE.ScaleVelocity(NVE, cons, parms);

  /*NVE.printVel(velout, cons);
  NVE.checkIdeal(geofile, cons, parms);
  */
  /*
    Evaluate the total Energy: Potential plus Kinetic plus Long range
    LJ correction
  */
  
  NVE.LongRangeCorrect(cons, parms);

  NVE.ComputeAllEnergy(NVE, cons, parms);

  /*Set initial Energy*/
  
  NVE.Einitial=NVE.Etotal;
  cout << "Einitial: "<<NVE.Einitial<<endl;

  /*****************************************************************/
  double Akin, InstTemp;
  /*
    Loop over desired time steps
  */
  
  cout <<"start integration:\n";
  for(t=0;t<cons.Tsteps;t++)
    {


      /*
	Velocity Verlet move with Rattle and Shake constraints on molecule
	inside NVEintegrateMotion is an AllGather MPI call.
       */
      
      NVEintegrateMotion(NVE,cons,parms, mols, local_mols, Watertype, t);
  
       
      /*Evaluate instantaneous Temperature and Kinetic Energy*/
      
      Akin=NVE.KineticA(cons);
      InstTemp=Akin/parms.ConstM;
    
      /*
      if(t%1000==0)
	NVE.ScaleVelocity(NVE, cons, parms);
      */
      /*mflops=N^2/1.0E6 and mflop_per_s=mflops/secs*/
      
      NVE.ComputeAllEnergy(NVE, cons, parms);
      
      
      
    }
  /*
    NVE.checkIdeal(endfile, cons, parms);
  
    NVE.printCoords(filename, cons);

  */
  mend=clock();
  cout << (mend-mstart)/CLOCKS_PER_SEC <<endl;
  cout << "done. "<<endl;
  MPI::Finalize();


}



void NVEintegrateMotion(Ensemble &NVE, Constants cons, Parameters parms, WaterSend *mols, WaterSend *local_mols, MPI::Datatype &Watertype, int t)
{
  /*NVE Velocity Verlet*/
  int i;
  for(i=0;i<cons.NperP;i++)
    {
      //update particle NVE.velocities (O, H, H)
      NVE.vx[i][0]+=cons.dt2*parms.InvMass[0]*NVE.Fx[i][0];
      NVE.vy[i][0]+=cons.dt2*parms.InvMass[0]*NVE.Fy[i][0];
      NVE.vz[i][0]+=cons.dt2*parms.InvMass[0]*NVE.Fz[i][0];

      NVE.vx[i][1]+=cons.dt2*parms.InvMass[1]*NVE.Fx[i][1];
      NVE.vy[i][1]+=cons.dt2*parms.InvMass[1]*NVE.Fy[i][1];
      NVE.vz[i][1]+=cons.dt2*parms.InvMass[1]*NVE.Fz[i][1];
      
      NVE.vx[i][2]+=cons.dt2*parms.InvMass[2]*NVE.Fx[i][2];
      NVE.vy[i][2]+=cons.dt2*parms.InvMass[2]*NVE.Fy[i][2];
      NVE.vz[i][2]+=cons.dt2*parms.InvMass[2]*NVE.Fz[i][2];
      

    }
  cout <<"integrate...moved vel.\n";

  //Apply PBC

  /*
    apply constraints to position and velocity and then update particles
    Move the Particles!!!
  */
  
  NVE.RattlePos(cons, parms, local_mols);
  
  /*
  char * vel1="vel_moved.out";
  char * pos1="pos_moved.out";
  NVE.printCoords(pos1, cons, mols);
  NVE.printVel(vel1, cons);
  NVE.writeParms(parms, cons);
  */
  cout <<"done RattlePos \n";
  
  /*Now do Mpi alltoall in order to get all the same positions per proc*/

  /*Send all updated positions to all processors*/
  MPI::COMM_WORLD.Allgather(local_mols, cons.NperP, Watertype, mols, cons.NperP, Watertype);
  cout << "Alltoall done: "<<endl;
  MPI::COMM_WORLD.Barrier();
  if(cons.rank==0)
    {
      char file[80];
      sprintf(file,"step.%d",t);
      char * filestr=file;
      NVE.printCoords(filestr, cons, mols);
    }

  //calculate force
  NVE.calcEForce(NVE, cons, parms, mols, local_mols);
  if(cons.rank==0)
    {
      char * forfile="myforce_0.out";
      NVE.writeAllForce(forfile, cons);
    }
  MPI::COMM_WORLD.Barrier();    
 
  for(i=0;i<cons.NperP;i++)
    {
      //update particle velocities
      NVE.vx[i][0]+=cons.dt2*parms.InvMass[0]*NVE.Fx[i][0];
      NVE.vy[i][0]+=cons.dt2*parms.InvMass[0]*NVE.Fy[i][0];
      NVE.vz[i][0]+=cons.dt2*parms.InvMass[0]*NVE.Fz[i][0];
      
      NVE.vx[i][1]+=cons.dt2*parms.InvMass[1]*NVE.Fx[i][1];
      NVE.vy[i][1]+=cons.dt2*parms.InvMass[1]*NVE.Fy[i][1];
      NVE.vz[i][1]+=cons.dt2*parms.InvMass[1]*NVE.Fz[i][1];
      
      NVE.vx[i][2]+=cons.dt2*parms.InvMass[2]*NVE.Fx[i][2];
      NVE.vy[i][2]+=cons.dt2*parms.InvMass[2]*NVE.Fy[i][2];
      NVE.vz[i][2]+=cons.dt2*parms.InvMass[2]*NVE.Fz[i][2];
    } 
  //apply constraints
  cout <<"start rattlevel:\n ";
  NVE.RattleVel(cons, parms, local_mols);
    

}
