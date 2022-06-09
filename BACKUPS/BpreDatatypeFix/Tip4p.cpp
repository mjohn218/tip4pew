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

void NVEintegrateMotion(Ensemble &NVE, Constants cons, Parameters parms, WaterSend *mols, WaterSend *local_mols, MPI::Datatype Watertype, int t, double *Kreal, double *Kimag, double *SKreal, double *SKimag);
void mKWALD(Ensemble &NNN, Constants cons, Parameters parms, WaterSend *mols, WaterSend *local_mols, double *Kreal, double *Kimag, double *SKreal, double *SKimag);
void mcalcEForce(Ensemble &NNN, Constants cons, Parameters parms, WaterSend *mols, WaterSend *local_mols, double *Kreal, double *Kimag, double *SKreal, double *SKimag);

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

  /*Define Kvectors, each processor has a K vector*/
  double *Kreal=new double[2442];
  double *Kimag=new double[2442];
  double *SKreal=new double[2442];
  double *SKimag=new double[2442];
  
  if(rank==0)
    {
      char *filename="StartCoord4.out";

      NVE.printAllCoords(filename, cons, mols);
    }

  /*Pull all molecules assigned locally into mols again*/
    MPI::COMM_WORLD.Allgather(local_mols, cons.NperP, Watertype, mols, cons.NperP, Watertype);
    int i;
    /*    for(i=0;i<2442;i++)
     {
       Kreal[i]=0.0;
       Kimag[i]=0.0;
       SKreal[i]=0.0;
       SKimag[i]=0.0;
       }*/
  if(rank==0)
    {
      char *filename="GatherCoord4.out";

      NVE.printAllCoords(filename, cons, mols);
    }

  
  /*
    Calculate initial force and energy 
  */
  
  mcalcEForce(NVE, cons, parms, mols, local_mols, Kreal, Kimag, SKreal, SKimag);

  if(rank==0)
    {
      ofstream kfile("kreal.out");
      ofstream k2file("kimag.out");
      ofstream mkfile("mykreal.out");
      ofstream mk2file("mykimag.out");
      for(i=0;i<2442;i++)
	{
	  kfile << SKreal[i] <<endl;
	  k2file << SKimag[i] <<endl;
	  mkfile <<Kreal[i]<<endl;
	  mk2file <<Kimag[i]<<endl;
	}
    }
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
      
      NVEintegrateMotion(NVE,cons,parms, mols, local_mols, Watertype, t, Kreal,  Kimag, SKreal, SKimag);
  
       
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



void NVEintegrateMotion(Ensemble &NVE, Constants cons, Parameters parms, WaterSend *mols, WaterSend *local_mols, MPI::Datatype Watertype, int t, double *Kreal, double *Kimag, double *SKreal, double *SKimag)
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
      NVE.printAllCoords(filestr, cons, mols);
      char *file2="local_0.out";
      NVE.printAllMyCoords(file2, cons, local_mols);
    }

  //calculate force
 
  mcalcEForce(NVE, cons, parms, mols, local_mols, Kreal, Kimag, SKreal, SKimag);
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


void mcalcEForce(Ensemble &NNN, Constants cons, Parameters parms, WaterSend *mols, WaterSend *local_mols, double *Kreal, double *Kimag, double *SKreal, double *SKimag)
{
  //calculate energy: derivative of the force
  //forces are short range (lennard jones) and long range (coulombic)
  //also external forces, in this case no external forces
  // reads in also AMASS


 

  int i;

  double CONSTANT=418.4; //joules conversion



  for(i=0;i<cons.NperP;i++)
    {
      NNN.Fx[i][0]=0.0;
      NNN.Fy[i][0]=0.0;
      NNN.Fz[i][0]=0.0;
      
      NNN.Fx[i][1]=0.0;
      NNN.Fy[i][1]=0.0;
      NNN.Fz[i][1]=0.0;
      
      NNN.Fx[i][2]=0.0;
      NNN.Fy[i][2]=0.0;
      NNN.Fz[i][2]=0.0;

      NNN.Fx[i][3]=0.0;
      NNN.Fy[i][3]=0.0;
      NNN.Fz[i][3]=0.0;


    }

  //Center of mass of the molecules
  //assumes you have 3 atoms per molecule (water) and evidently
  //each molecule per atom can have a different mass?
  NNN.CenterMassWater(cons, mols);
  //find n bonded energy and derivatives in real space
  /*
    if(cons.rank==0)
    {
      cout << "PRE RWALD\n"<<endl;
      NNN.writeParms(parms, cons);
    }
  */

  NNN.RWALD(cons, parms, mols, local_mols);

  char * myforcefile="rwaldForce.out";
  NNN.writeAllForce(myforcefile, cons);

  //find n bonded energy and derivatives in reciprocal space
  
  
  
  mKWALD(NNN, cons, parms, mols, local_mols, Kreal, Kimag, SKreal, SKimag);
  
  

  /* cout <<"fx0: " <<NNN.Fx[0][1]<<endl;
  cout <<"fy0: " <<NNN.Fy[0][1]<<endl;
  cout <<"fz0: " <<NNN.Fz[0][1]<<endl;
  */


  //evaluate forces on primary atoms due to secondary atom
  //and update the forces on all real atoms


  NNN.Primder(cons, parms, local_mols);
  
  /*  cout <<"primder\n";
  cout <<"fx0: " <<NNN.Fx[0][1]<<endl;
  cout <<"fy0: " <<NNN.Fy[0][1]<<endl;
  cout <<"fz0: " <<NNN.Fz[0][1]<<endl;
  */

  /*get energy sum: E0 from Rwald, EC from Rwald, Vk from Kwald*/


  NNN.EPOT=NNN.E0+NNN.EC+NNN.VK;
  cout << "E0: "<< NNN.E0 << "EC: "<<NNN.EC<<"VK: "<<NNN.VK<<endl;
  cout <<"actual total potential energy: "<< NNN.EPOT <<endl;

  //forces in units for the verlet equation
  double Fxtotal=0.0;
  double Fytotal=0.0;
  double Fztotal=0.0;


  for(i=0;i<cons.NperP;i++)
    { 
      NNN.Fx[i][0]*=CONSTANT;
      NNN.Fy[i][0]*=CONSTANT;
      NNN.Fz[i][0]*=CONSTANT;

      NNN.Fx[i][1]*=CONSTANT;
      NNN.Fy[i][1]*=CONSTANT;
      NNN.Fz[i][1]*=CONSTANT;

      NNN.Fx[i][2]*=CONSTANT;
      NNN.Fy[i][2]*=CONSTANT;
      NNN.Fz[i][2]*=CONSTANT;
      
      Fxtotal+=NNN.Fx[i][0]+NNN.Fx[i][1]+NNN.Fx[i][2];
      Fytotal+=NNN.Fy[i][0]+NNN.Fy[i][1]+NNN.Fy[i][2];
      Fztotal+=NNN.Fz[i][0]+NNN.Fz[i][1]+NNN.Fz[i][2];
    }
  /*char * forceout="forcefile.out";
  NNN.writeForce(forceout, cons);
  */
  //end of function calcEForce
  /*  cout << "Force of Fx[0][1]: "<< Fx[0][1]<<endl; 
  cout << "Force of Fy[0][1]: "<< Fy[0][1]<<endl; 
  cout << "Force of Fz[0][1]: "<< Fz[0][1]<<endl; 
  */
  cout <<"\nFxtotal: "<<Fxtotal<<" proc: "<<cons.rank<<endl;
  cout << "\nFytotal: "<<Fytotal<<" proc "<<cons.rank <<endl;
  cout << "\nFztotal: "<<Fztotal<<"proc " <<cons.rank <<endl;
  
}

void mKWALD(Ensemble & NNN, Constants cons, Parameters parms, WaterSend *mols, WaterSend *local_mols, double *Kreal, double *Kimag, double *SKreal, double *SKimag)
{
  //Now we are in Fourier space
  //what do we need from real space?
  /*
    This function needs to be adjusted if N changes or Nsite or
    *****Kmax******** array allocations and loops depend Kmax=10.
  */
  
  //Energy and derivatives in Fourier space (reciprocal space frequency)
  //Want VK and KVir_LJ

  /*Return to main*/
  NNN.VK=0.0;
  
  double FACT2=2.0*cons.Pi*parms.InvBOXL;
  double A22I=1.0/(2.0*cons.Alpha2);
  double FACTVK=2.0*FACT2*parms.InvBOXL*parms.InvBOXL;
  double VIRKxx=0.0;
  double VIRKyy=0.0;
  double VIRKzz=0.0;
  double VIRMxx=0.0;
  double VIRMyy=0.0;
  double VIRMzz=0.0;
  double VIRxx;   //VIRKxx-VIRMxx
  double VIRyy;
  double VIRzz;
  int i, j, k;
  int zi, yi;

  /*for this program Thetx from 0 through 10 are k=0:10
    Thetx from 11 through 20 corresponds to k=-10:-1 i.e. Thetx[i][j][20] is kx=-1
    Replaced complex expix, expiy etc. with just their argument, all have magnitude 1 and grossly
    simplifies computing expixk and exp(irk)
  */
  double SumKRe[3];
  double SumKIm[3];
  double SKRe, SKIm;
  double Thetx[11][cons.NperP][3];
  double Thety[21][cons.NperP][3];
  double Thetz[21][cons.NperP][3];

  double Thetrki[3];
  
  double CFact;
  double magnitude;
  double mag2;
  double argSK;
  double argx, argy, argz;
  double FCT;
  double Veck;
  int Ksq;
  int Kscut=105;
  double Xk, Yk, Zk;
  double RFact;
  double Vex, Vs;
  double Rx1, Ry1, Rz1;
  double TERM;
  double Xp, Yp, Zp;
  double RKsq, InvRKsq;
  double alphPi=cons.ALPHA/sqrt(cons.Pi);
  double qq;

  double Fxik, Fyik, Fzik;
  double R2, R;

  double tSam;
  ofstream kf("compRe.out");
  ofstream k2f("compIm.out");


  /*initialize fourier terms for each particle for k vectors of [0,0,0], [1,1,1], and [n,-1,-1]
   Indexing is strange in ky and kz, thetz[0][][]=kz=-10, thetz[20][][]=kz=10*/
  for(i=0;i<cons.NperP;i++)
    {
      //H1 index 10 is k=0
      Thetx[0][i][0]=0.0;
      Thety[10][i][0]=0.0;
      Thetz[10][i][0]=0.0;
      
      argx=FACT2*local_mols[i].x[1];
      argy=FACT2*local_mols[i].y[1];
      argz=FACT2*local_mols[i].z[1];
      //11 is k=+1
      Thetx[1][i][0]=argx;
      Thety[11][i][0]=argy;
      Thetz[11][i][0]=argz;
      
      //9 is k=-1
      //not looping over negative kx vectors
      Thety[9][i][0]=-argy;
      Thetz[9][i][0]=-argz;
	
      //H2
      Thetx[0][i][1]=0.0;
      Thety[10][i][1]=0.0;
      Thetz[10][i][1]=0.0;
      
      argx=FACT2*local_mols[i].x[2];
      argy=FACT2*local_mols[i].y[2];
      argz=FACT2*local_mols[i].z[2];
      Thetx[1][i][1]=argx;
      Thety[11][i][1]=argy;
      Thetz[11][i][1]=argz;
      
      Thety[9][i][1]=-argy;
      Thetz[9][i][1]=-argz;

      //Msite      
      Thetx[0][i][2]=0.0;
      Thety[10][i][2]=0.0;
      Thetz[10][i][2]=0.0;
      
      argx=FACT2*local_mols[i].x[3];
      argy=FACT2*local_mols[i].y[3];
      argz=FACT2*local_mols[i].z[3];
      Thetx[1][i][2]=argx;
      Thety[11][i][2]=argy;
      Thetz[11][i][2]=argz;

      Thety[9][i][2]=-argy;
      Thetz[9][i][2]=-argz;
	
    }
  int ki;
  int kx, ky, kz;

  //12 is +2
  /*Only loop over your molecules, but all loop over all k-vectors*/     
  for(ki=2;ki<11;ki++)
    {
      for(i=0;i<cons.NperP;i++)
	{
	  for(j=0;j<3;j++)
	    {
	      Thetx[ki][i][j]=ki*Thetx[1][i][j];
	      
	      Thety[ki+10][i][j]=ki*Thety[11][i][j];
	      Thety[10-ki][i][j]=-ki*Thety[11][i][j];
	      
	      Thetz[ki+10][i][j]=ki*Thetz[11][i][j];
	      Thetz[10-ki][i][j]=-ki*Thetz[11][i][j];
	    }
	  
	}
    }
  
  /*Have just computed exp(irk) for all k values for all particles(each particle has one r value)
    now need to compute sum over all r values for each k
*/
  double start, end;
  int vec=0;
  start=clock();
  
  for(kx=0;kx<11;kx++)
    {
      if(kx==0) FCT=1.0;
      else FCT=2.0;
      yi=0;
      for(ky=-10;ky<11;ky++)
	{
	  zi=0;
	  for(kz=-10;kz<11;kz++)
	    {

	      Ksq=kx*kx+ky*ky+kz*kz;
	      if(Ksq!=0 && Ksq <= Kscut)
		{
		  Xk=kx*FACT2;
		  Yk=ky*FACT2;
		  Zk=kz*FACT2;
		  /*Rksq is for one k vector k^2 * (2pi/L)^2*/
		  RKsq=Xk*Xk+Yk*Yk+Zk*Zk;
		  InvRKsq=1.0/RKsq;
		  
		  /*sumK is sum over all atoms q*exp(-irk) for one k vector loop over i and j*/
		  
		  SumKRe[0]=0.0;
		  SumKIm[0]=0.0;
		  SumKRe[1]=0.0;
		  SumKIm[1]=0.0;
		  SumKRe[2]=0.0;
		  SumKIm[2]=0.0;
		  for(i=0;i<cons.NperP;i++)
		    {
		      Thetrki[0]=Thetx[kx][i][0]+Thety[yi][i][0]+Thetz[zi][i][0];
		      SumKRe[0]+=cos(Thetrki[0]);
		      SumKIm[0]+=sin(Thetrki[0]);
		      
		      Thetrki[1]=Thetx[kx][i][1]+Thety[yi][i][1]+Thetz[zi][i][1];
		      SumKRe[1]+=cos(Thetrki[1]);
		      SumKIm[1]+=sin(Thetrki[1]);
		      
		      Thetrki[2]=Thetx[kx][i][2]+Thety[yi][i][2]+Thetz[zi][i][2];
		      SumKRe[2]+=cos(Thetrki[2]);
		      SumKIm[2]+=sin(Thetrki[2]);
		      
		    }
		  SumKRe[0]*=parms.Q[1];
		  SumKIm[0]*=parms.Q[1];
		  SumKRe[1]*=parms.Q[2];
		  SumKIm[1]*=parms.Q[2];
		  SumKRe[2]*=parms.Q[3];
		  SumKIm[2]*=parms.Q[3];
		  
		  SKRe=SumKRe[0]+SumKRe[1]+SumKRe[2];
		  SKIm=SumKIm[0]+SumKIm[1]+SumKIm[2];
		  kf << SKRe << "     ";
		  Kreal[vec]=SKRe;
		  
		  kf << SKRe << "    " << Kreal[vec] << endl;
		  k2f << SKIm <<endl;
		  Kimag[vec]=SKIm;
		  vec++;
		  // cout << "hello world\n";
		}//End if in Ksq limits
	      zi++;
	      
	    }
	  yi++;
	  
	}
    }
  /*End of loop over k vectors*/
  cout << "numK vecs: "<<vec <<endl;
  cout << "should be -313: "<<Kreal[0]<<endl;

		  /*   StartSend************************************************************************************************************************************************************************
Need to do a global send here to turn the SKRe below over all K vectors, so
end K loop here and restart another afterwards*/
  MPI::COMM_WORLD.Allreduce(Kreal, SKreal, 2442, MPI::DOUBLE, MPI::SUM);
  MPI::COMM_WORLD.Allreduce(Kimag, SKimag, 2442, MPI::DOUBLE, MPI::SUM);

  /*Now that everyone has the k-vector, they can update the forces on their atoms*/
  vec=0;
  /*
    after the alltoall, the kvector should be summer over processors and now
    inside of kvector.Re[vec] should be total over all N particles for each Kvector
  */
  for(kx=0;kx<11;kx++)
    {
      if(kx==0) FCT=1.0;
      else FCT=2.0;
      yi=0;
      for(ky=-10;ky<11;ky++)
	{
	  zi=0;
	  for(kz=-10;kz<11;kz++)
	    {

	      Ksq=kx*kx+ky*ky+kz*kz;
	      if(Ksq!=0 && Ksq <= Kscut)
		{
		  Xk=kx*FACT2;
		  Yk=ky*FACT2;
		  Zk=kz*FACT2;
		  /*Rksq is for one k vector k^2 * (2pi/L)^2*/
		  RKsq=Xk*Xk+Yk*Yk+Zk*Zk;
		  InvRKsq=1.0/RKsq;


		  Veck=FACTVK*exp(-cons.SIGMA*RKsq)/RKsq;
		  //magni^2 of complex number SK=Sum_r(exp(irk))
		  SKRe=SKreal[vec];
		  SKIm=SKimag[vec];
		  vec++;
		  
		  mag2=SKRe*SKRe+SKIm*SKIm;
		  magnitude=sqrt(mag2);
		  argSK=atan(SKIm/SKRe);
		  /*atan is only -pi/2 to pi/2 need whole range*/
		  if(SKRe<0)
		    argSK+=cons.Pi;
		  /*Potential term in K-space to electrostatic potential*/
		  tSam=0.5*FCT*Veck*mag2;

		  NNN.VK+=tSam;

		  VIRKxx+=tSam*(1.0-2.0*Xk*Xk*InvRKsq-Xk*Xk*A22I);
		  VIRKyy+=tSam*(1.0-2.0*Yk*Yk*InvRKsq-Yk*Yk*A22I);
		  VIRKzz+=tSam*(1.0-2.0*Zk*Zk*InvRKsq-Zk*Zk*A22I);
		  for(i=0;i<cons.NperP;i++)
		    {
		      Thetrki[0]=Thetx[kx][i][0]+Thety[yi][i][0]+Thetz[zi][i][0];
		      Thetrki[1]=Thetx[kx][i][1]+Thety[yi][i][1]+Thetz[zi][i][1];
		      Thetrki[2]=Thetx[kx][i][2]+Thety[yi][i][2]+Thetz[zi][i][2];

		      for(j=0;j<3;j++)
			{
			  RFact=-FCT*parms.Q[j+1]*Veck;
			  CFact=magnitude*sin(argSK-Thetrki[j]);
			  /*Force is derivative of Vk*/
			  Fxik=RFact*Xk*CFact;
			  Fyik=RFact*Yk*CFact;
			  Fzik=RFact*Zk*CFact;
			  NNN.Fx[i][j+1]+=Fxik;
			  NNN.Fy[i][j+1]+=Fyik;
			  NNN.Fz[i][j+1]+=Fzik;
			  Xp=local_mols[i].x[j+1]-NNN.XCOM[i];
			  Yp=local_mols[i].y[j+1]-NNN.YCOM[i];
			  Zp=local_mols[i].z[j+1]-NNN.ZCOM[i];
			  VIRMxx+=Xp*Fxik;
			  VIRMyy+=Yp*Fyik;
			  VIRMzz+=Zp*Fzik;
			}
		    }
		}//End if in Ksq limits
	      zi++;
	    }
	  yi++;
	}
    }
  //end loop over kx, ky, kz
  end=clock();
  cout << "kloop: "<<(end-start)/CLOCKS_PER_SEC<<endl;

  VIRxx=VIRKxx-VIRMxx;
  VIRyy=VIRKyy-VIRMyy;
  VIRzz=VIRKzz-VIRMzz;
  /*Return to main*/
  NNN.KVir_LJ=VIRxx+VIRyy+VIRzz;


  /* Last two pieces of electorstatic potential Vex is electorstatic coulomb sum from charges interacting within the same molecule.  Vs is spurious self interaction term of continuos gaussian with the point charge at the
center of it.*/

  Vex=0.0;
  Vs=0.0;

  for(j=1;j<cons.NSITE;j++)
    {
      for(k=1;k<cons.NSITE;k++)
	{
	  for(i=0;i<cons.NperP;i++)
	    {
	      if(j!=k)
		{
		  Rx1=local_mols[i].x[j]-local_mols[i].x[k];
		  Ry1=local_mols[i].y[j]-local_mols[i].y[k];
		  Rz1=local_mols[i].z[j]-local_mols[i].z[k];
		  R2=Rx1*Rx1+Ry1*Ry1+Rz1*Rz1;
		  R=sqrt(R2);
		  TERM=0.5*parms.Q[j]*parms.Q[k];
		  if(R!=0.0){
		    
		    Vex+=TERM*erf(cons.ALPHA*R)/R;  
		  }
		  		  
		}
	      	      
	    }
	}
    }//end for loops over j, k, i
  qq=parms.Q[1]*parms.Q[1]+parms.Q[2]*parms.Q[2]+parms.Q[3]*parms.Q[3];
  Vs=cons.N*alphPi*qq;
  cout <<"VS: "<<Vs<<" Vex: "<<Vex<<" VK: "<<NNN.VK<<endl;
  NNN.VK=NNN.VK-Vs-Vex;
  /*  cout <<"in kwald\n";
  cout <<"fx0: " <<Fx[0][1]<<endl;
  cout <<"fy0: " <<Fy[0][1]<<endl;
  cout <<"fz0: " <<Fz[0][1]<<endl;
  */
			     

}
