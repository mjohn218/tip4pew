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


//Functional routines, integration, force, energy calculations, 
//specific code is in Ensemble.cpp

void NVEintegrateMotion(Ensemble &NVE, Constants cons, Parameters parms, WaterSend *mols, WaterSend *local_mols, MPI::Datatype Watertype, int t, double *Kreal, double *Kimag, double *SKreal, double *SKimag, Kindex *kxyz, DiffSend *myopos);

void mKWALD(Ensemble &NNN, Constants cons, Parameters parms, WaterSend *mols, WaterSend *local_mols, double *Kreal, double *Kimag, double *SKreal, double *SKimag, Kindex *kxyz);

void mcalcEForce(Ensemble &NNN, Constants cons, Parameters parms, WaterSend *mols, WaterSend *local_mols, double *Kreal, double *Kimag, double *SKreal, double *SKimag, Kindex *kxyz);

void Ksetup1(Constants &cons);
void Ksetup2(Constants cons, Kindex *kxyz, double &Vself, Parameters parms);
void ComputeEnergyInit(Ensemble & NNN, double *Efinish, Constants cons, Parameters parms);
void ComputeEnergy(Ensemble & NNN, double *Efinish, Constants cons, Parameters parms);
void ScaleVelocity(Ensemble &NNN, Constants cons, Parameters parms);

/*Global Variables*/

int N_Water=512;


int main(int argc, char* argv[])
{
  double mstart, mend;
  int i;
  mstart=clock();

  /*Build datatype describing class*/
  MPI::Datatype Watertype;
  MPI::Datatype types[3]={MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
  MPI::Aint disp[3];
  MPI::Aint lb, extent;
  MPI::Datatype Difftype;

  MPI::Aint ddp[3];
  
  MPI::Init(argc, argv);
  
  int rank=MPI::COMM_WORLD.Get_rank();
  int nprocs=MPI::COMM_WORLD.Get_size();
    
  Ensemble NVE;
  Parameters parms;
  Constants cons;
  Averages a;
   
  /*Create a datatype to pass around that is a class. Includes
   positions of all atoms .Also a datatype for keeping track of 
   positions not put back inside the periodic box (for Diffusion
   Calculations)*/
  
  int blocklen[3]={4, 4, 4};
  int dblock[3]={3, 3, 3};
  /*Number of mols per processor*/

  cons.NperP=N_Water/nprocs;
  WaterSend * mols= new WaterSend[N_Water];
  WaterSend * local_mols = new WaterSend[cons.NperP];
  cout << "nperP: "<<cons.NperP <<endl;
  DiffSend * opos=new DiffSend[N_Water];
  DiffSend * myopos=new DiffSend[cons.NperP];

  /*Compute displacements of class components*/

  disp[0]=0;
  disp[1]=32;  //bytes=4*8. (4 doubles in each block) 
  disp[2]=64;
  Watertype=Watertype.Create_struct( 3, blocklen, disp, types);
  //int size=Watertype.Get_size();
  Watertype.Get_extent(lb, extent);
  Watertype.Commit();

  ddp[0]=0;
  ddp[1]=24;  //bytes=3*8 (3 double in each block)
  ddp[2]=48;
  Difftype=Difftype.Create_struct(3, dblock, ddp, types);
  Difftype.Commit();

  srand(static_cast<unsigned>(time(0)));

  /*All initial positions read in on Processor 0 as well as all parameters*/
  //Thermostat therm;
  //Barostat baro;

  int t;
  char *filesave;
  char *velsave;
  cout.precision(12);

  /*
    All processors read in initial coordinates and parameters
   coordinates go into WaterSend mols object.Then split them up into a local
   set
  */
  
  char * PARMFILE=argv[1];
  char * COORDFILE=argv[2];
  if(!argv[1]) cerr << "no parameter file \n"; 
  if(!argv[2]) cerr << "no coord file \n";


  NVE.readParms(PARMFILE, cons, parms);
  NVE.readCoords(COORDFILE, mols);
  
  cons.N=N_Water;
  cons.rank=rank;
  cout << " N molecules: "<< cons.N <<endl;

  if(cons.N != 512)
    {
      cout <<"arrays not initialiazed to correct size!" <<endl;
      return(1);
    }

  /*
    Initialize variables, water molecule constants and Tip4p-Ewald
    parameters 
  */
  
  NVE.InitialParms(parms, cons, a);
  NVE.RigidMol(parms, cons, mols, opos);

  /*Set up the Kvectors*/

  Ksetup1(cons);
  Kindex * kxyz=new Kindex[cons.Kleng];
  Ksetup2(cons, kxyz, NVE.Vs, parms);


  /*
    Make local molecules point to a part of the main array
  */
  
  local_mols=&mols[rank*cons.NperP];
  myopos=&opos[rank*cons.NperP];

  /*Define Kvectors, each processor has a K vector, Kreal is local
   SKreal is Total vector after summing across Processors*/
  double *Kreal=new double[cons.Kleng];
  double *Kimag=new double[cons.Kleng];
  double *SKreal=new double[cons.Kleng];
  double *SKimag=new double[cons.Kleng];
  
  /*
    Calculate initial force and energy And create EnergyVector
    to be summed across processors
  */
  double *Evector=new double[7];
  NVE.Energyv=new double[7];
  
  mcalcEForce(NVE, cons, parms, mols, local_mols, Kreal, Kimag, SKreal, SKimag, kxyz);


  /*
    For new coordinates set new velocities from Maxwell Boltzman, then 
    remove components along normal to bonds
  */

  WaterSend * InitVel=new WaterSend[N_Water];
    
  if(cons.Equil==1)
    {
      cout <<"Initializing Vel from Maxwell Boltzman. \n";
      NVE.Maxwell(cons, parms, InitVel);
      NVE.VelocityNormal(cons, parms, mols, InitVel);
    }
  else
    {
      /*Read in velocities from file. Each Proc reads in velocities*/
      char * velfile="velocity.inp";
      cout <<"Reading In Vel from file: "<< velfile <<endl;
      NVE.readVelocity(velfile, InitVel);
    }

  /*Now assign each processor there chunk of velocities to accompany their
    molecules.  Thoughout remainder of simulation each Proc only knows 
    N/P velocities and Forces
  */
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

  /*
    Evaluate the total Energy: Potential plus Kinetic plus Long range
    LJ correction
  */
  
  NVE.LongRangeCorrect(cons, parms);
  ComputeEnergyInit(NVE,Evector, cons, parms);
  
  /*Set initial Energy*/
  double Einitial;
  Einitial=NVE.Etotal;
  cout << "Einitial: "<<Einitial<<endl;

  cout << "MUST CREATE DIRECTORIES TO SAVE CONFIGS INSIDE!!\n";
  cout << " Name them: DiffP and PbcP "<<endl;
  /*****************************************************************/
  /*Loop over desired time steps*/

  int Teqp1=cons.Tequil+1;
  double enint;
  double chi=0.0;
  ofstream ipropfile("inst_prop.out");
  ipropfile.precision(12);
  cout << "Equilibrating.. for: "<<cons.Tequil<<" steps."<<endl;

  for(t=1;t<Teqp1;t++)
    {
      /*Equilibrating*/

      /*
	Velocity Verlet move with Rattle and Shake constraints on molecule
	inside NVEintegrateMotion is an AllGather MPI call.
       */
      
      NVEintegrateMotion(NVE,cons,parms, mols, local_mols, Watertype, t, Kreal,  Kimag, SKreal, SKimag, kxyz, myopos);
          
      
      /*Scale Velocities for the first ~50000 steps. Measure Initial Energy
	From point after finished scaling vels*/
      
      if(t<50001)
	{
	  if(t%1000==0)
	    {
	      ScaleVelocity(NVE, cons, parms);
	    }
	}

      
      /*Compute and Write out Instantaneous properties ~500 steps*/
      if(t%parms.statwrite==0)
	{
	  ComputeEnergy(NVE,Evector, cons, parms);
	  if(rank==0)
	    {
	      ipropfile<<t<<' '<<NVE.Etotal<<' '<<NVE.EKIN<<' '<<NVE.EPOT<<' '<<NVE.InstTemp<<' '<<NVE.Pint<<endl;
	    }
	}
      if(t==50001){
	ComputeEnergy(NVE,Evector, cons, parms);
	Einitial=NVE.Etotal;
      }
    
    }
  
  
  /*Done with Equilibraton time length, Write out this equilibrated config
   And Begin Data collection part of simulation*/
  
  char file[60];
  char filec[60];
  char filev[60];
  filesave="Equil_config.out";
  if(rank==0)  NVE.printCoords(filesave, cons, mols,t);
  sprintf(filev, "Equil_vel.%d",rank);
  velsave=filev;
  NVE.printVel(velsave, cons, t); 
  cout <<" Done equilibrating..Collect statistics & trajectories\n";

  int save100=0;
  int Tstepp1=cons.Tsteps+1; //because of cpp convention, I want start at 1
  ofstream chifile("chi_deviations.out");
  ofstream posfile;
  ofstream diffile;
  diffile.precision(12);
  chifile.precision(12);
  posfile.precision(12);
    
  /*Running simulation ...Done equilibrating    */

  for(t=1;t<Tstepp1;t++)
    {
      
      NVEintegrateMotion(NVE,cons,parms, mols, local_mols, Watertype, t, Kreal,  Kimag, SKreal, SKimag, kxyz, myopos);
      
      ComputeEnergy(NVE, Evector, cons, parms);
      
      /*Keep Track of Averages*/

      a.Eavg+=NVE.Etotal;
      a.E2avg+=NVE.Etotal*NVE.Etotal;
      a.Epavg+=NVE.EPOT;
      a.Ekavg+=NVE.EKIN;
      a.Tavg+=NVE.InstTemp;
      a.Pavg+=NVE.Pint;

      /*Finished computation for this tstep*/

      if(t%parms.isave==0)
	{
	  /*Add in a Gather command so that will have the non-pbc coords*/
	  MPI::COMM_WORLD.Barrier();
	  MPI::COMM_WORLD.Gather(myopos, cons.NperP, Difftype, opos, cons.NperP, Difftype, 0);
	  
	  if(rank==0)
	    {
	      /*Save 20 configs to each file*/
	      if(save100==0)
		{
		  sprintf(file, "PbcP/save_pbc_pos.%d", t/parms.isave);
		  posfile.open(file);
		  sprintf(filec, "DiffP/ohh_diff_pos.%d", t/parms.isave);
		  diffile.open(filec);
		  //sprintf(file, "save_vel.%d", t/cons.isave);
		  //vfile.open(file);
		}
	      posfile << "Save Config: "<< t <<" Eavg "<<a.Eavg/(t*1.0)<<endl;
	      //posfile <<"\n E2avg "<<a.E2avg/(t*1.0)<<endl;
	      //posfile <<" Epavg: "<<a.Epavg/(t*1.0)<<endl;
	      //posfile <<" Ekavg: "<<a.Ekavg/(t*1.0)<<endl;
	      //posfile <<" Tavg: "<<a.Tavg/(t*1.0)<<" Pavg: "<<a.Pavg/(t*1.0)<<endl;
	      diffile << "save "<< t <<endl;
	      for(i=0;i<cons.N;i++)
		{
		  posfile <<mols[i].x[0]<<' '<<mols[i].y[0]<<' '<<mols[i].z[0]<<endl;
		  posfile <<mols[i].x[1]<<' '<<mols[i].y[1]<<' '<<mols[i].z[1]<<endl;
		  posfile <<mols[i].x[2]<<' '<<mols[i].y[2]<<' '<<mols[i].z[2]<<endl;
		  diffile <<opos[i].x[0]<<' '<<opos[i].y[0]<<' '<<opos[i].z[0]<<endl;
		  diffile <<opos[i].x[1]<<' '<<opos[i].y[1]<<' '<<opos[i].z[1]<<endl;
		  diffile <<opos[i].x[2]<<' '<<opos[i].y[2]<<' '<<opos[i].z[2]<<endl;
		}
	      save100++;
	      if(save100==100)
		{
		  save100=0;
		  posfile.close();
		  diffile.close();
		  //vfile.close();
		}
	      
	    }//End if on processor 0
	}// end if time to save a config
      
      /*Evaluate average Energy deviations*/
      enint=abs((NVE.Etotal-Einitial)/Einitial);
      chi=(chi*(t-1)+enint)/(t*1.0);
      if(rank==0)
	{
	  if(t%parms.statwrite==0)
	    {
	      chifile << t<<' '<<enint<<' '<<chi<<endl;
	      ipropfile<<t<<' '<<NVE.Etotal<<' '<<NVE.EKIN<<' '<<NVE.EPOT<<' '<<NVE.InstTemp<<' '<<NVE.Pint<<endl;
	    }
	}	  
      if(t%parms.allsave==0)
	{
	  //Do this Infrequently, opens and closes a file every time!
	  if(rank==0)
	    {
	      sprintf(filec, "config_all.%d", t/parms.allsave);
	      ofstream pout(filec);
	      pout << "Save Config: "<< t <<" Eavg "<<a.Eavg/(t*1.0);
	      pout <<"\n E2avg "<<a.E2avg/(t*1.0)<<endl;
	      pout <<" Epavg: "<<a.Epavg/(t*1.0)<<endl;
	      pout <<" Ekavg: "<<a.Ekavg/(t*1.0)<<endl;
	      pout <<" Tavg: "<<a.Tavg/(t*1.0)<<" Pavg: "<<a.Pavg/(t*1.0)<<endl;
	      for(i=0;i<cons.N;i++)
		{
		  pout <<mols[i].x[0]<<' '<<mols[i].y[0]<<mols[i].z[0]<<endl;
		  pout <<mols[i].x[1]<<' '<<mols[i].y[1]<<mols[i].z[1]<<endl;
		  pout <<mols[i].x[2]<<' '<<mols[i].y[2]<<mols[i].z[2]<<endl;
		}
	    }
	  sprintf(filev, "veloc_all.%d.%d",rank, t/parms.allsave);
	  velsave=filev;
	  NVE.printVel(velsave, cons, t); 
	}
      
    }//End of Simulation
  

  /*Finished Write out Averages*/
  ofstream afile("Averages.out");
  afile <<" Energy Avg: "<<a.Eavg/(t*1.0)<<endl;
  afile <<"\n E2avg "<<a.E2avg/(t*1.0)<<endl;
  afile <<" Epavg: "<<a.Epavg/(t*1.0)<<endl;
  afile <<" Ekavg: "<<a.Ekavg/(t*1.0)<<endl;
  afile <<" Tavg: "<<a.Tavg/(t*1.0)<<endl;
  afile <<" Pavg: "<<a.Pavg/(t*1.0)<<endl;
  afile.close();
  mend=clock();
  cout << (mend-mstart)/CLOCKS_PER_SEC <<endl;
  cout << "done. "<<endl;
  MPI::Finalize();


}



void NVEintegrateMotion(Ensemble &NVE, Constants cons, Parameters parms, WaterSend *mols, WaterSend *local_mols, MPI::Datatype Watertype, int t, double *Kreal, double *Kimag, double *SKreal, double *SKimag, Kindex *kxyz, DiffSend *myopos)
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

  /*
    apply constraints to position and velocity and then update particles
    Move the Particles!!!
  */
  
  NVE.RattlePos(cons, parms, local_mols, myopos);
  
  
  /*Now do Mpi alltoall in order to get all the same positions per proc*/

  /*Send all updated positions to all processors*/
  MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Allgather(local_mols, cons.NperP, Watertype, mols, cons.NperP, Watertype);
  

  //calculate force

  mcalcEForce(NVE, cons, parms, mols, local_mols, Kreal, Kimag, SKreal, SKimag, kxyz);
 

  /*Final step in Velocity Verlet*/

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
  
  /*apply constraints*/

  NVE.RattleVel(cons, parms, local_mols);
    

}


void mcalcEForce(Ensemble &NNN, Constants cons, Parameters parms, WaterSend *mols, WaterSend *local_mols, double *Kreal, double *Kimag, double *SKreal, double *SKimag, Kindex *kxyz)
{
  //calculate energy: derivative of the force
  //forces are short range (lennard jones) and long range (coulombic)
  //also external forces, in this case no external forces

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
  //assumes you have 3 atoms per molecule (water) 

  NNN.CenterMassWater(cons, mols);
  //find n bonded energy and derivatives in real space

  
  NNN.RWALD(cons, parms, mols, local_mols);
  
  //find n bonded energy and derivatives in reciprocal space
  
  mKWALD(NNN, cons, parms, mols, local_mols, Kreal, Kimag, SKreal, SKimag, kxyz);
  

  //evaluate forces on primary atoms due to secondary atom
  //and update the forces on all real atoms


  NNN.Primder(cons, parms, local_mols);
  

  /*Forces in units for verlet equation, From kcal/mol/Angstrom to
   g*Angstrom/ps^2/mol, then F*t/m =Anstrom/ps= velocity*/

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
      
      //Fxtotal+=NNN.Fx[i][0]+NNN.Fx[i][1]+NNN.Fx[i][2];
      // Fytotal+=NNN.Fy[i][0]+NNN.Fy[i][1]+NNN.Fy[i][2];
      // Fztotal+=NNN.Fz[i][0]+NNN.Fz[i][1]+NNN.Fz[i][2];
    }

  /*
  cout <<"\nFxtotal: "<<Fxtotal<<" proc: "<<cons.rank<<endl;
  cout << "\nFytotal: "<<Fytotal<<" proc "<<cons.rank <<endl;
  cout << "\nFztotal: "<<Fztotal<<"proc " <<cons.rank <<endl;
  */
  //end of function mcalcEForce
  
}

void mKWALD(Ensemble & NNN, Constants cons, Parameters parms, WaterSend *mols, WaterSend *local_mols, double *Kreal, double *Kimag, double *SKreal, double *SKimag, Kindex *kxyz)
{
  //Now we are in Fourier space

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

  int i, j, k;
 
  double SumKRe[3];
  double SumKIm[3];
  double SKRe, SKIm;
 
  double Thetrki[3];
  
  double CFact;
  double magnitude;
  double mag2;
  double argSK;
 
  double FCT;
  double Veck;
 
  double Xk, Yk, Zk;
  double RFact;
  
  double Rx1, Ry1, Rz1;
  double TERM;
  double Xp, Yp, Zp;
  double RKsq, InvRKsq;
  double Fxik, Fyik, Fzik;
  double R2, R;

  double tSam;
  int ki;
  double start, end;

  start=clock();
  /*
    Length of Kvector and the kx, ky, kz values of each vector are set up
    previously in Ksetup1 and Ksetup2
  */
  for(ki=0;ki<cons.Kleng;ki++)
    {
  
      if(kxyz[ki].kx==0) FCT=1.0;
      else FCT=2.0;


      Xk=kxyz[ki].kx*FACT2;
      Yk=kxyz[ki].ky*FACT2;
      Zk=kxyz[ki].kz*FACT2;
      
      /*sumK is sum over all atoms q*exp(-irk) for one k vector loop over i and j*/
      
      SumKRe[0]=0.0;
      SumKIm[0]=0.0;
      SumKRe[1]=0.0;
      SumKIm[1]=0.0;
      SumKRe[2]=0.0;
      SumKIm[2]=0.0;
      for(i=0;i<cons.NperP;i++)
	{
	  Thetrki[0]=Xk*local_mols[i].x[1]+Yk*local_mols[i].y[1]+Zk*local_mols[i].z[1];
	  SumKRe[0]+=cos(Thetrki[0]);
	  SumKIm[0]+=sin(Thetrki[0]);

	  Thetrki[1]=Xk*local_mols[i].x[2]+Yk*local_mols[i].y[2]+Zk*local_mols[i].z[2];
	  SumKRe[1]+=cos(Thetrki[1]);
	  SumKIm[1]+=sin(Thetrki[1]);

	  Thetrki[2]=Xk*local_mols[i].x[3]+Yk*local_mols[i].y[3]+Zk*local_mols[i].z[3];
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
     
      Kreal[ki]=SKRe;
      Kimag[ki]=SKIm;
  
    }

  /*End of loop over k vectors*/


  /*   StartSend************************************************************************************************************************************************************************
Need to do a global send here to turn the SKRe below over all K vectors over 
all particles, so
end K loop here and restart another afterwards*/
  MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Allreduce(Kreal, SKreal, cons.Kleng, MPI::DOUBLE, MPI::SUM);
  MPI::COMM_WORLD.Allreduce(Kimag, SKimag, cons.Kleng, MPI::DOUBLE, MPI::SUM);

  /*Now that everyone has the k-vector, they can update the forces on their atoms*/
  
  for(ki=0;ki<cons.Kleng;ki++)
    {
  
      if(kxyz[ki].kx==0) FCT=1.0;
      else FCT=2.0;


      Xk=kxyz[ki].kx*FACT2;
      Yk=kxyz[ki].ky*FACT2;
      Zk=kxyz[ki].kz*FACT2;
      /*Rksq is for one k vector k^2 * (2pi/L)^2*/
      RKsq=Xk*Xk+Yk*Yk+Zk*Zk;
      InvRKsq=1.0/RKsq;


      Veck=FACTVK*exp(-cons.SIGMA*RKsq)/RKsq;
      SKRe=SKreal[ki];
      SKIm=SKimag[ki];
		  
      //magni^2 of complex number SK=Sum_r(exp(irk))		  
      mag2=SKRe*SKRe+SKIm*SKIm;
      magnitude=sqrt(mag2);
      argSK=atan(SKIm/SKRe);
      /*atan is only -pi/2 to pi/2 need whole range*/
      if(SKRe<0)
	argSK+=cons.Pi;

      /*Potential term in K-space to electrostatic potential*/
      tSam=0.5*FCT*Veck*mag2;

      NNN.VK+=tSam;
      /*Contirbution to virial pressure (sum_r F*r) */ 
     
      VIRKxx+=tSam*(1.0-2.0*Xk*Xk*InvRKsq-Xk*Xk*A22I);
      VIRKyy+=tSam*(1.0-2.0*Yk*Yk*InvRKsq-Yk*Yk*A22I);
      VIRKzz+=tSam*(1.0-2.0*Zk*Zk*InvRKsq-Zk*Zk*A22I);
      for(i=0;i<cons.NperP;i++)
	{
	  Thetrki[0]=Xk*local_mols[i].x[1]+Yk*local_mols[i].y[1]+Zk*local_mols[i].z[1];
	  Thetrki[1]=Xk*local_mols[i].x[2]+Yk*local_mols[i].y[2]+Zk*local_mols[i].z[2];
	  Thetrki[2]=Xk*local_mols[i].x[3]+Yk*local_mols[i].y[3]+Zk*local_mols[i].z[3];

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
		
	}//end loop over myN particles

    }//End secondn loop over kvector

  end=clock();
  //  cout << "kloop: "<<(end-start)/CLOCKS_PER_SEC<<endl;

  //  VIRxx=VIRKxx-VIRMxx;
  //VIRyy=VIRKyy-VIRMyy;
  //VIRzz=VIRKzz-VIRMzz;
  NNN.VirK=VIRKxx+VIRKyy+VIRKzz;
  NNN.VirM=VIRMxx+VIRMyy+VIRMzz;

  /*Kvirial will be computed in the same way,(VirK-VirM) but first VirM must
   be summed across processors. */

  /* Last two pieces of electorstatic potential Vex is electorstatic coulomb sum from charges interacting within the same molecule.  Vs is spurious self interaction term of continuos gaussian with the point charge at the
center of it.Vs moved into Ksetup because it never changes*/

  NNN.Vex=0.0;

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
		    
		    NNN.Vex+=TERM*erf(cons.ALPHA*R)/R;  
		  }
		  		  
		}
	      	      
	    }
	}
    }//end for loops over j, k, i

  //NNN.VK=NNN.VK-Vs-Vex;
  /*Save in Vector for MPIreduce*/
  NNN.Energyv[3]=NNN.Vex;
  NNN.Energyv[4]=NNN.VirM;
			     

}

void Ksetup1(Constants &cons)
{
  /*Set up k vector to loop over, need Kscut, and Kvector and index object.*/

  int lengthK=0;
  int ix, iy, iz;
  int Ksq;


  for(ix=0;ix<cons.Kmax;ix++)
    {
      for(iy=cons.Kmin;iy<cons.Kmax;iy++)
	{
	  for(iz=cons.Kmin;iz<cons.Kmax;iz++)
	    {
	      Ksq=ix*ix+iy*iy+iz*iz;
	      if(Ksq!=0)
		{
		  if(Ksq<=cons.Kscut)
		    {
		      lengthK++;
		    }
		}
	    }
	}
    }

  cons.Kleng=lengthK;
  cout << "Size of K vector: " << cons.Kleng <<endl;
}
void Ksetup2(Constants cons, Kindex *kxyz, double &Vself, Parameters parms)
{  
  /*Allocate array to hold kvector coordinates*/
 
  int Ksq;
  int ix, iy, iz;
  int ind=0;
  for(ix=0;ix<cons.Kmax;ix++)
    {
      for(iy=cons.Kmin;iy<cons.Kmax;iy++)
	{
	  for(iz=cons.Kmin;iz<cons.Kmax;iz++)
	    {
	      Ksq=ix*ix+iy*iy+iz*iz;
	      if(Ksq!=0)
		{
		  if(Ksq<=cons.Kscut)
		    {
		      kxyz[ind].kx=ix;
		      kxyz[ind].ky=iy;
		      kxyz[ind].kz=iz;
		      ind++;
		    }
		}
	    }
	}
    }
  double qq;
  double alphPi=cons.ALPHA/sqrt(cons.Pi);
  qq=parms.Q[1]*parms.Q[1]+parms.Q[2]*parms.Q[2]+parms.Q[3]*parms.Q[3];
  Vself=cons.N*alphPi*qq;
  cout << "Vself; "<< Vself <<endl;

  /*End of Ksetup*/
}

void ComputeEnergyInit(Ensemble & NNN, double *Efinish, Constants cons, Parameters parms)
{
  /*Energyv[0]=E0
    Energyv[1]=EC
    Energyv[2]=Real_Vir
    Energyv[3]=Vex
    Energyv[4]=VirM
    Energyv[5]=kineCOM
    Energyv[6]=EKIN
  */
  double kineCOM;
  double VKtot;
  double VirTotal;
  double kine;

  kineCOM=NNN.CenterMassKinEnergy(cons);
  NNN.Energyv[5]=kineCOM;
  
  kine=NNN.KineticA(cons);
  kine=kine*0.5/418.4;
  NNN.Energyv[6]=kine;

/*Sum the energies across processors*/
  MPI::COMM_WORLD.Reduce(NNN.Energyv, Efinish, 7, MPI::DOUBLE, MPI::SUM, 0);
  /*0.5 Factor becasue everything in rwald is summed over all pairs twice
   */
  if(cons.rank==0)
    {
      Efinish[0]=Efinish[0]*0.5;
      Efinish[1]=Efinish[1]*0.5;
      Efinish[2]=Efinish[2]*0.5;
      
      VKtot=NNN.VK-NNN.Vs-Efinish[3];
      NNN.EPOT=Efinish[0]+Efinish[1]+VKtot;
      
      NNN.EKIN=Efinish[6];
      NNN.ELRCint=NNN.ELRC/parms.Volume;
      NNN.Etotal=NNN.EPOT+NNN.EKIN+NNN.ELRCint;
      
      
      NNN.KVirial=NNN.VirK-Efinish[4];
      VirTotal=Efinish[2]+NNN.KVirial;
      NNN.PLRCint=NNN.PLRC/parms.Volume;
      NNN.Pint=(Efinish[5]+NNN.PLRCint+VirTotal)/parms.pconstant;
      NNN.InstTemp=NNN.EKIN*2.0*418.4/parms.ConstM;
      
      cout << "E0 "<< Efinish[0]<<" EC " <<Efinish[1] << " VK "<< NNN.VK <<" Vex " <<Efinish[3]<<" Vs "<<NNN.Vs<<" VKall "<< VKtot<<endl;
      cout <<" EPOT: "<<NNN.EPOT <<" EKIN " <<NNN.EKIN <<" ELRCint "<<NNN.ELRCint<<endl;
      /*ELRCint should not change unless volume changes, only in NPT*/
      cout << " Etotal :"<<NNN.Etotal <<endl;
      cout <<endl;
      cout << "RealVir: "<< Efinish[2]<<" KVirial "<<NNN.KVirial<<endl;
      cout << "InstTemp: "<< NNN.InstTemp<<endl;
      cout << " Pint "<< NNN.Pint <<endl;
    }
  
  
}
void ComputeEnergy(Ensemble & NNN, double *Efinish, Constants cons, Parameters parms)
{
  /*Energyv[0]=E0
    Energyv[1]=EC
    Energyv[2]=Real_Vir
    Energyv[3]=Vex
    Energyv[4]=VirM
    Energyv[5]=kineCOM
    Energyv[6]=EKIN
  */
  double kineCOM;
  double VKtot;
  double VirTotal;
  double kine;

  kineCOM=NNN.CenterMassKinEnergy(cons);
  NNN.Energyv[5]=kineCOM;
  
  kine=NNN.KineticA(cons);
  kine=kine*0.5/418.4;
  NNN.Energyv[6]=kine;

/*Sum the energies across processors*/
  MPI::COMM_WORLD.Reduce(NNN.Energyv, Efinish, 7, MPI::DOUBLE, MPI::SUM, 0);
  /*0.5 Factor becasue everything in rwald is summed over all pairs twice
   */
  if(cons.rank==0)
    {
      Efinish[0]=Efinish[0]*0.5;
      Efinish[1]=Efinish[1]*0.5;
      Efinish[2]=Efinish[2]*0.5;
      
      VKtot=NNN.VK-NNN.Vs-Efinish[3];
      NNN.EPOT=Efinish[0]+Efinish[1]+VKtot;
      
      NNN.EKIN=Efinish[6];
      //NNN.ELRCint=NNN.ELRC/parms.Volume; Does not change.
      NNN.Etotal=NNN.EPOT+NNN.EKIN+NNN.ELRCint;
      
      
      NNN.KVirial=NNN.VirK-Efinish[4];
      VirTotal=Efinish[2]+NNN.KVirial;
      //NNN.PLRCint=NNN.PLRC/parms.Volume; Does not change. 
      NNN.Pint=(Efinish[5]+NNN.PLRCint+VirTotal)/parms.pconstant;
      NNN.InstTemp=NNN.EKIN*2.0*418.4/parms.ConstM;
      
    }
  
  
}
void ScaleVelocity(Ensemble &NNN, Constants cons, Parameters parms)
{
  /*Scale the velocities to the set temperature 
   Sets the Energy Ensemble Variable EKIN*/

  double Akin;
  double * myTemp=new double[1];
  double * InstTemp=new double[1];
  double Tf;
  int i;


  Akin=NNN.KineticA(cons);
  myTemp[0]=Akin/parms.ConstM;
  /*Must sum InstTemp variable over all particle vels use Allreduce*/
  MPI::COMM_WORLD.Allreduce(myTemp, InstTemp, 1, MPI::DOUBLE, MPI::SUM);
  //cout<<"Temp pre scale: "<<InstTemp[0]<<endl;
  Tf=sqrt(cons.Temp/InstTemp[0]);
  //scale
  for(i=0;i<cons.NperP;i++)
    {
      NNN.vx[i][0]*=Tf;
      NNN.vy[i][0]*=Tf;
      NNN.vz[i][0]*=Tf;

      NNN.vx[i][1]*=Tf;
      NNN.vy[i][1]*=Tf;
      NNN.vz[i][1]*=Tf;

      NNN.vx[i][2]*=Tf;
      NNN.vy[i][2]*=Tf;
      NNN.vz[i][2]*=Tf;

    }

}

