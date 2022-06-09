/*
  Program: Tip4p Ewald in cpp
  Programmer: Margaret Johnson
  Start: 1/10/2005
  This is a program to run molecular dynamics simulations given initial 
congfigurations of solution volumes, i.e. water. To run in NPT, NVT, NVE 
ensembles. Using Frenkel and Smit Molecular simulation book, and fortran
code from THG updated by Rajesh  Murarka. Also code from Elizabeth Verschell

Written in parallel 5/2005

Requires 
1. This main code, 2. Ensemble.cpp 3. Ensemble.h 4. Parameters.h
compile with Makefile

*/

#include <mpi++.h>
#include <fstream>
#include <ctime>
//#include "md_timer.h"
#include <math.h>
#include <cmath>
#include <complex>
#include "Ensemble.h"
#include <iomanip>
#include "Parameters.h"
#include <unistd.h>
#define nperp 72

using namespace std;

//NON MPI-specific code is in Ensemble.cpp (functions that require no communication)
class kspace
{
 public:
  double p[nperp][3];
};


void NVTintegrateMotion(Ensemble &NVE, Constants cons, Parameters parms, WaterSend *mols, WaterSend *local_mols, MPI::Datatype Watertype, int t, double *Kreal, double *Kimag, double *SKreal, double *SKimag, Kindex *kxyz, DiffSend *myopos, DiffSend *vel, WaterSend *force, Center *com, Table *terf, double *Veck, kspace *xr, kspace *xi, kspace *yr, kspace *yi, kspace *zr, kspace *zi, kspace *ekr, kspace *eki, double *Mat, double *Minv, Thermostat &th);

void mKWALD(Ensemble &NNN, Constants cons, Parameters parms, WaterSend *mols, WaterSend *local_mols, double *Kreal, double *Kimag, double *SKreal, double *SKimag, Kindex *kxyz, WaterSend *force, Center *com, double *Veck, kspace *xr, kspace *xi, kspace *yr, kspace *yi, kspace *zr, kspace *zi, kspace *ekr, kspace *eki);

void mcalcEForce(Ensemble &NNN, Constants cons, Parameters parms, WaterSend *mols, WaterSend *local_mols, double *Kreal, double *Kimag, double *SKreal, double *SKimag, Kindex *kxyz, WaterSend *force, Center *com, Table *terf, double *Veck, kspace *xr, kspace *xi, kspace *yr, kspace *yi, kspace *zr, kspace *zi, kspace *ekr, kspace *eki);

void Ksetup1(Constants &cons);
void Ksetup2(Constants cons, Kindex *kxyz, double &Vself, Parameters parms, double *Veck);
void ComputeEnergyInit(Ensemble & NNN, double *Efinish, Constants cons, Parameters parms, DiffSend *vel, Thermostat &th);
void ComputeEnergy(Ensemble & NNN, double *Efinish, Constants cons, Parameters parms, DiffSend *vel, Thermostat &th);
void ScaleVelocity(Ensemble &NNN, Constants cons, Parameters parms, DiffSend *vel);
void calcRDF(WaterSend *mols, double *gr, double boxl, int N, double delr );

int main(int argc, char* argv[])
{
  /*MPI VERSION*/

  /*Build  Two datatypes describing classes of Molecules */

  MPI::Datatype Watertype; //For Molecules with 4 positions
  MPI::Datatype types[3]={MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
  MPI::Aint disp[3];
  MPI::Aint lb, extent;
  MPI::Datatype Difftype;  //For molecules need 3 positions only
  MPI::Aint ddp[3];

  /*Initialize MPI Communication world*/
  
  MPI::Init(argc, argv);

  /* Local Variables */  
  Ensemble NVT;
  Thermostat th;
  Parameters parms;
  Constants cons;
  Averages a;
  int rank, nprocs;
  
  rank=MPI::COMM_WORLD.Get_rank();
  nprocs=MPI::COMM_WORLD.Get_size();
  
   
  /*Create a datatype to pass around that encompasses the classes. Includes
   positions of all atoms .Also a datatype for keeping track of 
   positions not put back inside the periodic box (for Diffusion
   Calculations, needs no 4th positions)*/
  
  int blocklen[3]={4, 4, 4};
  int dblock[3]={3, 3, 3};
  
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
  
  /*Finished setting up MPI datatypes Watertype and Difftype. */

  char * PARMFILE=argv[1];
  char * COORDFILE=argv[2];
  if(!argv[1]) cerr << "no parameter file \n"; 
  if(!argv[2]) cerr << "no coord file \n";

  /*Read in and set parameters and constants from input file*/

  readNVTParms(PARMFILE, cons, parms, th);

  /*Number of mols per processor, Only Set up for Perfect distribution!*/
  cons.NperP=cons.N/nprocs;
  cout << "NperP: "<<cons.NperP <<endl;
  cout << "defined nperp for kspace class: "<<nperp<<endl;
  if(cons.NperP!=nperp){
    cout << " !!!!defined nperp does not equal N/nprocs !!!!"<<endl;
    cout <<" FIX DEFINED nperp TO EQUAL: "<<cons.NperP<<" ..EXITING"<<endl;
    exit(1);
  }

  /*WaterSend and DiffSend Classes defined in Ensemble.h*/
  WaterSend * mols= new WaterSend[cons.N];
  WaterSend * local_mols = new WaterSend[cons.NperP];
  DiffSend * opos=new DiffSend[cons.N];
  DiffSend * myopos=new DiffSend[cons.NperP];
  
  Center *com=new Center[cons.N];

    /*
    All processors read in initial coordinates and parameters
   coordinates go into WaterSend mols object.Then split them up into a local
   set
  */
  
  readCoords(COORDFILE, mols);
  
  cons.rank=rank;
  cons.nproc=nprocs;
  cout << " N molecules: "<< cons.N <<endl;


  /*
    Initialize variables, water molecule constants and Tip4p-Ewald
    parameters 
  */
  
  InitialParms(parms, cons, a, NVT);
  RigidMol(parms, cons, mols, opos, com);
  setThermostat(th, cons);

  char *filesave;
  filesave="start_config.out";
  if(rank==0)  printCoords(filesave, cons, mols,0);


  int i, t;
  ofstream ipropfile("inst_prop.out");

  char *velsave;
  int Teqp1=cons.Tequil+1;
  int Tstepp1=cons.Tsteps+1; //because of cpp convention, I want start at 1
  cons.Tvscale+=1;
  double enint;
  double chi=0.0;
  /*Matrices for rattlepos and rattlevel, 3X3, read across rows, then down*/
  double *Mat=new double[9];
  double *Minv=new double[9];
  double *Evector=new double[7];
  double Einitial;
  char file[60];
  char filec[60];
  char filev[60];

  ofstream thermfile("thermostat.out");
  int save100;
  int pt;
  NVT.Energyv=new double[7];
  WaterSend * force=new WaterSend[cons.NperP];
  ipropfile.precision(12);
  thermfile.precision(12);
  cout.precision(12);
  srand(static_cast<unsigned>(time(0)));
  /*
    Make arrays to hold g(r)
  */
  double delr=.025;
  int nbins;
  int NumConfigs=0;
  nbins=int(parms.BOXLENGTH/(2.0*delr));
  double *gr=new double[nbins];
  for(i=0;i<nbins;i++)gr[i]=0.0;

  /*
    Make local molecules point to a part of the main array
  */
  
  local_mols=&mols[rank*cons.NperP];
  myopos=&opos[rank*cons.NperP];

  /*Set up the Kvectors*/
  
  Ksetup1(cons);
  Kindex * kxyz=new Kindex[cons.Kleng];
  double *Veck=new double[cons.Kleng];
  Ksetup2(cons, kxyz, NVT.Vs, parms, Veck);

  /*Define Kvectors, each processor has a K vector, Kreal is local
   SKreal is Total vector after summing across Processors*/
  double *Kreal=new double[cons.Kleng];
  double *Kimag=new double[cons.Kleng];
  double *SKreal=new double[cons.Kleng];
  double *SKimag=new double[cons.Kleng];
  
  /*kspace class defined in Ensemble.h. Speeds up Kwald summations*/
  kspace * xr=new kspace[11];
  kspace * xi=new kspace[11];
  kspace * yr=new kspace[11];
  kspace * yi=new kspace[11];
  kspace * zr=new kspace[11];
  kspace * zi=new kspace[11];
  kspace * ekr=new kspace[cons.Kleng];
  kspace * eki=new kspace[cons.Kleng];

  /*Lookup table for Erf function, does not generate much speedups*/
  Table *terf=new Table[500001];
  TableSet(terf, cons, parms);

  /*
   Evaluate Initial Forces on the molecules
  */  

  mcalcEForce(NVT, cons, parms, mols, local_mols, Kreal, Kimag, SKreal, SKimag, kxyz, force, com, terf, Veck, xr, xi, yr, yi, zr, zi, ekr, eki);


  /*
    For new coordinates set new velocities from Maxwell Boltzman, then 
    remove components along normal to bonds
  */

  DiffSend * InitVel=new DiffSend[cons.N];
  DiffSend * myvel=new DiffSend[cons.NperP];
 
  double xx, yy, zz;
  double xv, yv, zv;
  double rvdot;



  if(cons.Equil==1)
    {
      cout <<"Initializing Vel from Maxwell Boltzman. "<< cons.Temp<<endl;
      if(rank==0)
      	{
	  /*Initialize velocities only on proc 0 to ensure COM is fixed 
	    point, then scatter to all processors*/
	  Maxwell(cons, parms, InitVel);
	  VelocityNormal(cons, parms, mols, InitVel, Mat, Minv);
	}
      MPI::COMM_WORLD.Barrier();
      MPI::COMM_WORLD.Scatter(InitVel, cons.NperP, Difftype, myvel, cons.NperP, Difftype, 0);

    }
  else
    {
      /*Read in velocities from file. Proc 0 reads in velocities
	then scatters chunks to each other proc.*/
      if(rank==0)
	{
	  char * velfile="velocity.inp";
	  cout <<"Reading In Vel from file: "<< velfile <<endl;
	  ifstream VELF(velfile);
	  if(!VELF) {
	    cerr << "NO VELOCITY FILE GIVEN!\n";
	    cout <<"INITIALIZING FROM MAXWELL-BOLTZMAN..!"<<endl;
	    Maxwell(cons, parms, InitVel);
	  }
	  else
	    {
	      //VELF >> N; 
	      for(i=0;i<cons.N;i++)
		{
		  VELF>>InitVel[i].x[0]>>InitVel[i].y[0]>>InitVel[i].z[0];
		  VELF>>InitVel[i].x[1]>>InitVel[i].y[1]>>InitVel[i].z[1];
		  VELF>>InitVel[i].x[2]>>InitVel[i].y[2]>>InitVel[i].z[2];
		}
	      VELF.close();
	    }
	  //readVelocity(velfile, InitVel);
	  VelocityNormal(cons, parms, mols, InitVel, Mat, Minv);
	  
	}
      MPI::COMM_WORLD.Scatter(InitVel, cons.NperP, Difftype, myvel, cons.NperP, Difftype, 0);
    }
  MPI::COMM_WORLD.Barrier();
  int ch=cons.NperP-2;

  xx=local_mols[ch].x[0]-local_mols[ch].x[1];
  yy=local_mols[ch].y[0]-local_mols[ch].y[1];
  zz=local_mols[ch].z[0]-local_mols[ch].z[1];
  xx-=parms.BOXLENGTH*round(xx*parms.InvBOXL);
  yy-=parms.BOXLENGTH*round(yy*parms.InvBOXL);
  zz-=parms.BOXLENGTH*round(zz*parms.InvBOXL);
  xv=myvel[ch].x[0]-myvel[ch].x[1];
  yv=myvel[ch].y[0]-myvel[ch].y[1];
  zv=myvel[ch].z[0]-myvel[ch].z[1];
  rvdot=xx*xv+yy*yv+zz*zv;
  cout << "ZERO...."<<endl;
  cout << rvdot <<' '<< rank<<endl;

  if(cons.Equil==1) ScaleVelocity(NVT, cons, parms, myvel);
  
  MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Gather(myvel, cons.NperP, Difftype, InitVel, cons.NperP, Difftype, 0);
  velsave="start_vel.out";
  if(rank==0) printAllVel(velsave, cons, 0, InitVel); 

  /*  Thoughout remainder of simulation each Proc only knows 
    N/P velocities and Forces
  */
  
  /*
    Evaluate the total Energy: Potential plus Kinetic plus Long range
    LJ correction
  */
  
  LongRangeCorrect(cons, parms, NVT);
  ComputeEnergyInit(NVT,Evector, cons, parms, myvel, th);
  
  /*Set initial Energy*/
  
  Einitial=NVT.Etotal;
  cout << "Einitial: "<<Einitial<<endl;
  if(rank==0)
    {
      thermfile <<0<<endl;
      thermfile <<th.Eta[0]<<' '<<th.Eta[1]<<' '<<th.Eta[2]<<' '<<th.Eta[3]<<' '<<th.Eta[4]<<endl;
      thermfile <<th.Veta[0]<<' '<<th.Veta[1]<<' '<<th.Veta[2]<<' '<<th.Veta[3]<<' '<<th.Veta[4]<<endl;
    }
  
  cout << "\nMUST CREATE DIRECTORIES TO SAVE CONFIGS INSIDE!!\n";
  cout << " Name them: DiffP and PbcP "<<endl;
  cout <<" Number of configs saved per file: "<<parms.nperfile<<endl;

  /*****************************************************************/
  /*BEGIN SIMULATION, FIRST EQUILIBRATION ROUND, THEN DATA COLLECTION*/

  cout << "Equilibrating.. for: "<<cons.Tequil<<" steps."<<endl;
    
  for(t=1;t<Teqp1;t++)
    {
      /*Equilibrating*/

      /*
	Velocity Verlet move with Rattle and Shake constraints on molecule
	inside NVTintegrateMotion is an AllGather MPI call.
       */
      
      NVTintegrateMotion(NVT,cons,parms, mols, local_mols, Watertype, t, Kreal,  Kimag, SKreal, SKimag, kxyz, myopos, myvel, force, com, terf, Veck, xr, xi, yr, yi, zr, zi, ekr, eki, Mat, Minv, th);
      
      
      ComputeEnergy(NVT,Evector, cons, parms, myvel, th);

      /*Compute and Write out Instantaneous properties ~500 steps*/
      if(t%parms.statwrite==0)
	{
	  if(rank==0)
	    {
	      ipropfile<<t<<' '<<NVT.Etotal<<' '<<NVT.EKIN<<' '<<NVT.EPOT<<' '<<NVT.InstTemp<<' '<<NVT.Pint<<endl;
	    
	      thermfile << t<<endl;
	      thermfile <<th.Eta[0]<<' '<<th.Eta[1]<<' '<<th.Eta[2]<<' '<<th.Eta[3]<<' '<<th.Eta[4]<<endl;
	      thermfile <<th.Veta[0]<<' '<<th.Veta[1]<<' '<<th.Veta[2]<<' '<<th.Veta[3]<<' '<<th.Veta[4]<<endl;
	    }
	}


      if(t%100000==0)
	{
	  //Do this Infrequently, opens and closes a file every time!
	  if(rank==0)
	    {
	      sprintf(filec, "config_equil.%d", t/100000);
	      ofstream pout(filec);
	      pout.precision(12);
	      pout << "Save Config: "<< t <<endl;
	      for(i=0;i<cons.N;i++)
		{
		  pout <<mols[i].x[0]<<' '<<mols[i].y[0]<<' '<<mols[i].z[0]<<endl;
		  pout <<mols[i].x[1]<<' '<<mols[i].y[1]<<' '<<mols[i].z[1]<<endl;
		  pout <<mols[i].x[2]<<' '<<mols[i].y[2]<<' '<<mols[i].z[2]<<endl;
		}
	      pout.close();
	    }
	  MPI::COMM_WORLD.Barrier();
	  MPI::COMM_WORLD.Gather(myvel, cons.NperP, Difftype, InitVel, cons.NperP, Difftype, 0);
	  sprintf(filev, "veloc_equil.%d", t/100000);
	  velsave=filev;
	  if(rank==0) printAllVel(velsave, cons, t, InitVel); 

	}
    
    }
  
  /*Done with Equilibraton time length, Write out this equilibrated config
   And Begin Data collection part of simulation*/
  
  filesave="Equil_config.out";
  if(rank==0)  printCoords(filesave, cons, mols,t);
  MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Gather(myvel, cons.NperP, Difftype, InitVel, cons.NperP, Difftype, 0);
  velsave="Equil_vel.out";
  if(rank==0) printAllVel(velsave, cons, t, InitVel); 


  cout <<" Done equilibrating..Collect statistics & trajectories\n";

  save100=0;
  pt=0;
  ofstream chifile("chi_deviations.out");
  ofstream posfile;
  ofstream diffile;
  diffile.precision(12);
  chifile.precision(12);
  posfile.precision(12);
    
  /*Running simulation ...Done equilibrating    */

  for(t=1;t<Tstepp1;t++)
    {
      
      NVTintegrateMotion(NVT,cons,parms, mols, local_mols, Watertype, t, Kreal,  Kimag, SKreal, SKimag, kxyz, myopos, myvel, force, com, terf, Veck, xr, xi, yr, yi, zr, zi, ekr, eki, Mat, Minv, th);
      
      ComputeEnergy(NVT, Evector, cons, parms, myvel, th);
      
      /*Keep Track of Averages*/

      a.Eavg+=NVT.Etotal;
      a.E2avg+=NVT.Etotal*NVT.Etotal;
      a.Epavg+=NVT.EPOT;
      a.Ep2avg+=NVT.EPOT*NVT.EPOT;
      a.Ekavg+=NVT.EKIN;
      a.Tavg+=NVT.InstTemp;
      a.Pavg+=NVT.Pint;
      a.P2avg+=NVT.Pint*NVT.Pint;

      /*FINISHED COMPUTATION FOR THIS TIMESTEP*/

      if(t%parms.isave==0)
	{
	  /*Add in a Gather command so that will have the non-pbc coords*/
	  MPI::COMM_WORLD.Barrier();
	  MPI::COMM_WORLD.Gather(myopos, cons.NperP, Difftype, opos, cons.NperP, Difftype, 0);
	  
	  if(rank==0)
	    {
	      calcRDF(mols, gr, parms.BOXLENGTH, cons.N, delr);
	      NumConfigs++;
	    }
	      /*Save 100 configs to each file*/
	  //  if(save100==0)
	  //{
	  //  pt++;
	  //  sprintf(file, "PbcP/save_pbc_pos.%d", pt);
	  //  posfile.open(file);
	  //  sprintf(filec, "DiffP/ohh_diff_pos.%d", pt);
	  //  diffile.open(filec);
	  //  //sprintf(file, "save_vel.%d", t/cons.isave);
	  //  //vfile.open(file);
	  //}
	  //  posfile << "Save Config: "<< t <<" Eavg "<<a.Eavg/(t*1.0)<<endl;
	  //  diffile << "save "<< t <<endl;
	  //  for(i=0;i<cons.N;i++)
	  //{
	  //  posfile <<mols[i].x[0]<<' '<<mols[i].y[0]<<' '<<mols[i].z[0]<<endl;
	  //  posfile <<mols[i].x[1]<<' '<<mols[i].y[1]<<' '<<mols[i].z[1]<<endl;
	  //  posfile <<mols[i].x[2]<<' '<<mols[i].y[2]<<' '<<mols[i].z[2]<<endl;
	  //  diffile <<opos[i].x[0]<<' '<<opos[i].y[0]<<' '<<opos[i].z[0]<<endl;
	  //  diffile <<opos[i].x[1]<<' '<<opos[i].y[1]<<' '<<opos[i].z[1]<<endl;
	  //  diffile <<opos[i].x[2]<<' '<<opos[i].y[2]<<' '<<opos[i].z[2]<<endl;
	  //}
	  //  save100++;
	  //  if(save100==parms.nperfile)
	  //{
	  //  save100=0;
	  //  posfile.close();
	  //  diffile.close();
	  //  //vfile.close();
	  //}
	      
	  //}//End if on processor 0
	}// end if time to save a config
      
      /*Evaluate average Energy deviations*/
      enint=abs((NVT.Etotal-Einitial)/Einitial);
      chi=(chi*(t-1)+enint)/(t*1.0);
      if(t%parms.statwrite==0)
	{
	  if(rank==0)
	    {
	      chifile << t<<' '<<enint<<' '<<chi<<endl;
	      ipropfile<<t<<' '<<NVT.Etotal<<' '<<NVT.EKIN<<' '<<NVT.EPOT<<' '<<NVT.InstTemp<<' '<<NVT.Pint<<endl;
	      
	      thermfile <<t<<endl;
	      thermfile <<th.Eta[0]<<' '<<th.Eta[1]<<' '<<th.Eta[2]<<' '<<th.Eta[3]<<' '<<th.Eta[4]<<endl;
	      thermfile <<th.Veta[0]<<' '<<th.Veta[1]<<' '<<th.Veta[2]<<' '<<th.Veta[3]<<' '<<th.Veta[4]<<endl;
	    }
	}

      if(t%parms.allsave==0)
	{
	  //Do this Infrequently, opens and closes a file every time!
	  if(rank==0)
	    {
	      sprintf(filec, "config_all.%d", t/parms.allsave);
	      ofstream pout(filec);
	      pout.precision(12);
	      pout << "Save Config: "<< t <<" Eavg "<<a.Eavg/(t*1.0);
	      pout <<"\n E2avg "<<a.E2avg/(t*1.0)<<endl;
	      pout <<" Epavg: "<<a.Epavg/(t*1.0)<<" Ep2avg: "<<a.Ep2avg/(t*1.0)<<endl;
	      pout <<" Ekavg: "<<a.Ekavg/(t*1.0)<<" P2avg: "<<a.P2avg/(t*1.0)<<endl;
	      pout <<" Tavg: "<<a.Tavg/(t*1.0)<<" Pavg: "<<a.Pavg/(t*1.0)<<endl;
	      for(i=0;i<cons.N;i++)
		{
		  pout <<mols[i].x[0]<<' '<<mols[i].y[0]<<' '<<mols[i].z[0]<<endl;
		  pout <<mols[i].x[1]<<' '<<mols[i].y[1]<<' '<<mols[i].z[1]<<endl;
		  pout <<mols[i].x[2]<<' '<<mols[i].y[2]<<' '<<mols[i].z[2]<<endl;
		}
	      pout.close();
	    }
	  MPI::COMM_WORLD.Barrier();
	  MPI::COMM_WORLD.Gather(myvel, cons.NperP, Difftype, InitVel, cons.NperP, Difftype, 0);
	  sprintf(filev, "veloc_all.%d", t/parms.allsave);
	  velsave=filev;
	  if(rank==0) printAllVel(velsave, cons, t, InitVel); 

	}
      
    }//End of Simulation
  

  /*Finished Write out Averages*/
  ofstream afile;
  afile.precision(12);
  if(rank==0)
    {
      afile.open("Averages.out");
      afile <<" Energy Avg: "<<a.Eavg/(t*1.0)<<endl;
      afile <<"\n E2avg "<<a.E2avg/(t*1.0)<<endl;
      afile <<" Epavg: "<<a.Epavg/(t*1.0)<<endl;
      afile <<" Ekavg: "<<a.Ekavg/(t*1.0)<<endl;
      afile <<" Tavg: "<<a.Tavg/(t*1.0)<<endl;
      afile <<" Pavg: "<<a.Pavg/(t*1.0)<<endl;
      afile <<" P2avg: "<<a.P2avg/(t*1.0)<<endl;
      afile <<" Ep2avg: "<<a.Ep2avg/(t*1.0)<<endl;
      afile.close();
  }

  /*Average rdf over all Numconfig configurations*/
  
  double pi=acos(-1.0);
  double nideal;
  double ra, vb;
  ofstream grfile;
  if(rank==0)
    {
      sprintf(file, "grnew_tip4p.%g.%g", cons.Temp, parms.Rhomm);
      grfile.open(file);
      for(i=0;i<nbins;i++)
	{
	  ra=delr*(i+0.5);
	  vb=(pow((i+1.0),3)-pow(i*1.0,3))*pow(delr,3.0);
	  nideal=(4.0/3.0)*pi*vb*parms.Density;
	  /*Normalize gr by the value of an rdf for an ideal gas*/
	  gr[i]=gr[i]/(NumConfigs*cons.N*nideal);
	  grfile << ra<<' '<<gr[i]<<endl;
	}
     grfile.close();
    }
  




  
  filesave="Final_config.out";
  if(rank==0)  printCoords(filesave, cons, mols,t);
  MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Gather(myvel, cons.NperP, Difftype, InitVel, cons.NperP, Difftype, 0);
  velsave="Final_vel.out";
  if(rank==0) printAllVelf(velsave, cons, t, InitVel); 


  MPI::COMM_WORLD.Barrier();
  cout << "done. "<<endl;
    
  MPI::Finalize();


}



void NVTintegrateMotion(Ensemble &NVT, Constants cons, Parameters parms, WaterSend *mols, WaterSend *local_mols, MPI::Datatype Watertype, int t, double *Kreal, double *Kimag, double *SKreal, double *SKimag, Kindex *kxyz, DiffSend *myopos, DiffSend *vel, WaterSend *force, Center *com, Table *terf, double *Veck, kspace *xr, kspace *xi, kspace *yr, kspace *yi, kspace *zr, kspace *zi, kspace *ekr, kspace *eki, double *Mat, double *Minv, Thermostat &th)
{
  /*NVT Velocity Verlet*/
  int i;
  double Ekint;
  double *mykin=new double[1];
  double *Tkin=new double[1];

  Tkin[0]=NVT.EKIN;
  /*At this point KE is updated and correct*/  
  MPI::COMM_WORLD.Bcast(Tkin, 1, MPI::DOUBLE, 0);

  NHCINT(vel, Tkin[0], cons, th);

  for(i=0;i<cons.NperP;i++)
    {
      //update particle NVT.velocities (O, H, H)
      vel[i].x[0]+=cons.dt2*parms.InvMass[0]*force[i].x[0];
      vel[i].y[0]+=cons.dt2*parms.InvMass[0]*force[i].y[0];
      vel[i].z[0]+=cons.dt2*parms.InvMass[0]*force[i].z[0];

      vel[i].x[1]+=cons.dt2*parms.InvMass[1]*force[i].x[1];
      vel[i].y[1]+=cons.dt2*parms.InvMass[1]*force[i].y[1];
      vel[i].z[1]+=cons.dt2*parms.InvMass[1]*force[i].z[1];
      
      vel[i].x[2]+=cons.dt2*parms.InvMass[2]*force[i].x[2];
      vel[i].y[2]+=cons.dt2*parms.InvMass[2]*force[i].y[2];
      vel[i].z[2]+=cons.dt2*parms.InvMass[2]*force[i].z[2];
      
    }

  /*
    apply constraints to position and velocity and then update particles
    Move the Particles!!!
  */
  
  RattlePos(cons, parms, local_mols, myopos, vel, Mat, Minv);
  
  
  /*Now do Mpi alltoall in order to get all the same positions per proc*/

  /*Send all updated positions to all processors*/
  MPI::COMM_WORLD.Allgather(local_mols, cons.NperP, Watertype, mols, cons.NperP, Watertype);
  

  //calculate force

  mcalcEForce(NVT, cons, parms, mols, local_mols, Kreal, Kimag, SKreal, SKimag, kxyz, force, com, terf, Veck, xr, xi, yr, yi,zr, zi, ekr, eki );
 

  /*Final step in Velocity Verlet*/

  for(i=0;i<cons.NperP;i++)
    {
      //update particle velocities
      vel[i].x[0]+=cons.dt2*parms.InvMass[0]*force[i].x[0];
      vel[i].y[0]+=cons.dt2*parms.InvMass[0]*force[i].y[0];
      vel[i].z[0]+=cons.dt2*parms.InvMass[0]*force[i].z[0];
      
      vel[i].x[1]+=cons.dt2*parms.InvMass[1]*force[i].x[1];
      vel[i].y[1]+=cons.dt2*parms.InvMass[1]*force[i].y[1];
      vel[i].z[1]+=cons.dt2*parms.InvMass[1]*force[i].z[1];
      
      vel[i].x[2]+=cons.dt2*parms.InvMass[2]*force[i].x[2];
      vel[i].y[2]+=cons.dt2*parms.InvMass[2]*force[i].y[2];
      vel[i].z[2]+=cons.dt2*parms.InvMass[2]*force[i].z[2];
    } 
  
  /*apply constraints*/

  RattleVel(cons, parms, local_mols, vel, Mat, Minv);
    
  /*Recompute EKIN */
  mykin[0]=KineticA(cons, vel);
  
  /*Must sum Akin variable over all particle vels use Allreduce*/
  MPI::COMM_WORLD.Allreduce(mykin, Tkin, 1, MPI::DOUBLE, MPI::SUM);
  Ekint=Tkin[0]*0.5/418.4;
  
  /*Update thermostat vels and pos, and scale velocities*/
  NHCINT(vel, Ekint, cons, th);
  delete[] mykin;
  delete[] Tkin;
}


void mcalcEForce(Ensemble &NNN, Constants cons, Parameters parms, WaterSend *mols, WaterSend *local_mols, double *Kreal, double *Kimag, double *SKreal, double *SKimag, Kindex *kxyz, WaterSend *force, Center *com, Table *terf, double *Veck, kspace *xr, kspace *xi, kspace *yr, kspace *yi, kspace *zr, kspace *zi, kspace *ekr, kspace *eki)
{
  //calculate energy: derivative of the force
  //forces are short range (lennard jones) and long range (coulombic)
  //also external forces, in this case no external forces

  int i;
  double CONSTANT=418.4; //joules conversion
  
  for(i=0;i<cons.NperP;i++)
    {
      force[i].x[0]=0.0;
      force[i].y[0]=0.0;
      force[i].z[0]=0.0;
      
      force[i].x[1]=0.0;
      force[i].y[1]=0.0;
      force[i].z[1]=0.0;
      
      force[i].x[2]=0.0;
      force[i].y[2]=0.0;
      force[i].z[2]=0.0;

      force[i].x[3]=0.0;
      force[i].y[3]=0.0;
      force[i].z[3]=0.0;


    }

  //Center of mass of the molecules
  //assumes you have 3 atoms per molecule (water) 

  CenterMassWater(cons, mols, com);
  //find n bonded energy and derivatives in real space


  RWALD(cons, parms, mols, local_mols, force, com, NNN, terf);


  //find n bonded energy and derivatives in reciprocal space

  mKWALD(NNN, cons, parms, mols, local_mols, Kreal, Kimag, SKreal, SKimag, kxyz, force, com, Veck, xr, xi, yr, yi, zr, zi, ekr, eki);


  //evaluate forces on primary atoms due to secondary atom
  //and update the forces on all real atoms



  Primder(cons, parms, local_mols, force);
  

  /*Forces in units for verlet equation, From kcal/mol/Angstrom to
   g*Angstrom/ps^2/mol, then F*t/m =Anstrom/ps= velocity*/

  for(i=0;i<cons.NperP;i++)
    { 
      force[i].x[0]*=CONSTANT;
      force[i].y[0]*=CONSTANT;
      force[i].z[0]*=CONSTANT;

      force[i].x[1]*=CONSTANT;
      force[i].y[1]*=CONSTANT;
      force[i].z[1]*=CONSTANT;

      force[i].x[2]*=CONSTANT;
      force[i].y[2]*=CONSTANT;
      force[i].z[2]*=CONSTANT;
      
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

void mKWALD(Ensemble & NNN, Constants cons, Parameters parms, WaterSend *mols, WaterSend *local_mols, double *Kreal, double *Kimag, double *SKreal, double *SKimag, Kindex *kxyz, WaterSend *force, Center *com, double *Veck, kspace *xr, kspace *xi, kspace *yr, kspace *yi, kspace *zr, kspace *zi, kspace *ekr, kspace *eki)
{
  //Now we are in Fourier space

  //Energy and derivatives in Fourier space (reciprocal space frequency)
  //Want VK and KVir_LJ

  /*Return to main*/
  NNN.VK=0.0;
  
  double FACT2=2.0*cons.Pi*parms.InvBOXL;
  double f2=FACT2;
  double A22I=1.0/(2.0*cons.Alpha2);
  double VIRKxx=0.0;
  double VIRKyy=0.0;
  double VIRKzz=0.0;
  double VIRMxx=0.0;
  double VIRMyy=0.0;
  double VIRMzz=0.0;

  int i, j;
 
  double SumKRe[3];
  double SumKIm[3];
  double SKRe, SKIm;
   
  double CFact;
  double mag2;
   
  double FCT;
  //  double Veck;
  int kx, ky, kz;
  double Xk, Yk, Zk;
  double RFact;
  
  double Rx1, Ry1, Rz1;
  double TERM;
  double Xp, Yp, Zp;
  double RKsq, InvRKsq;
  double Fxik, Fyik, Fzik;
  double R2, R;
  
  double c, s;
  double tSam;
  int ki;
  int zs, ys;
  int my=cons.rank*cons.NperP;

  for(i=0;i<cons.NperP;i++)
    {
      for(j=0;j<3;j++)
	{
	  xr[0].p[i][j]=1.0;
	  xi[0].p[i][j]=0.0;
	  yr[0].p[i][j]=1.0;
	  yi[0].p[i][j]=0.0;
	  zr[0].p[i][j]=1.0;
	  zi[0].p[i][j]=0.0;
	  
	  xr[1].p[i][j]=cos(f2*local_mols[i].x[j+1]);
	  xi[1].p[i][j]=sin(f2*local_mols[i].x[j+1]);
	  yr[1].p[i][j]=cos(f2*local_mols[i].y[j+1]);
	  yi[1].p[i][j]=sin(f2*local_mols[i].y[j+1]);
	  zr[1].p[i][j]=cos(f2*local_mols[i].z[j+1]);
	  zi[1].p[i][j]=sin(f2*local_mols[i].z[j+1]);
	  
	  xr[2].p[i][j]=2.0*xr[1].p[i][j]*xr[1].p[i][j]-1.0;
	  xi[2].p[i][j]=2.0*xr[1].p[i][j]*xi[1].p[i][j];
	  yr[2].p[i][j]=2.0*yr[1].p[i][j]*yr[1].p[i][j]-1.0;
	  yi[2].p[i][j]=2.0*yr[1].p[i][j]*yi[1].p[i][j];
	  zr[2].p[i][j]=2.0*zr[1].p[i][j]*zr[1].p[i][j]-1.0;
	  zi[2].p[i][j]=2.0*zr[1].p[i][j]*zi[1].p[i][j];
      
	}
    }
  //sin(2a)=2sin(a)cos(a);  cos(2a)=2*cos(a)*cos(a)-1
  for(ki=4;ki<11;ki+=2)
    {
      for(i=0;i<cons.NperP;i++)
	{
	  for(j=0;j<3;j++)
	    {
	      //Odd modes
	      xr[ki-1].p[i][j]=xr[ki-2].p[i][j]*xr[1].p[i][j]-xi[ki-2].p[i][j]*xi[1].p[i][j];
	      xi[ki-1].p[i][j]=xi[ki-2].p[i][j]*xr[1].p[i][j]+xi[1].p[i][j]*xr[ki-2].p[i][j];
	      yr[ki-1].p[i][j]=yr[ki-2].p[i][j]*yr[1].p[i][j]-yi[ki-2].p[i][j]*yi[1].p[i][j];
	      yi[ki-1].p[i][j]=yi[ki-2].p[i][j]*yr[1].p[i][j]+yi[1].p[i][j]*yr[ki-2].p[i][j];
	      zr[ki-1].p[i][j]=zr[ki-2].p[i][j]*zr[1].p[i][j]-zi[ki-2].p[i][j]*zi[1].p[i][j];
	      zi[ki-1].p[i][j]=zi[ki-2].p[i][j]*zr[1].p[i][j]+zi[1].p[i][j]*zr[ki-2].p[i][j];
	      
	      //even modes
	      xr[ki].p[i][j]=2.0*xr[ki/2].p[i][j]*xr[ki/2].p[i][j]-1.0;
	      xi[ki].p[i][j]=2.0*xr[ki/2].p[i][j]*xi[ki/2].p[i][j];
	      yr[ki].p[i][j]=2.0*yr[ki/2].p[i][j]*yr[ki/2].p[i][j]-1.0;
	      yi[ki].p[i][j]=2.0*yr[ki/2].p[i][j]*yi[ki/2].p[i][j];
	      zr[ki].p[i][j]=2.0*zr[ki/2].p[i][j]*zr[ki/2].p[i][j]-1.0;
	      zi[ki].p[i][j]=2.0*zr[ki/2].p[i][j]*zi[ki/2].p[i][j];
	    }
	}
    }
  /*
    Length of Kvector and the kx, ky, kz values of each vector are set up
    previously in Ksetup1 and Ksetup2
  */
  for(ki=0;ki<cons.Kleng;ki++)
    {
  
      kx=kxyz[ki].kx;
      ky=kxyz[ki].ky;
      kz=kxyz[ki].kz;
      
      /*sumK is sum over all atoms q*exp(-irk) for one k vector loop over i and j*/
      ys=kxyz[ki].sy;
      zs=kxyz[ki].sz;

      SumKRe[0]=0.0;
      SumKIm[0]=0.0;
      SumKRe[1]=0.0;
      SumKIm[1]=0.0;
      SumKRe[2]=0.0;
      SumKIm[2]=0.0;
      for(i=0;i<cons.NperP;i++)
	{
	  for(j=0;j<3;j++)
	    {
	      c=xr[kx].p[i][j]*yr[ky].p[i][j]*zr[kz].p[i][j]-xi[kx].p[i][j]*yr[ky].p[i][j]*zi[kz].p[i][j]*zs-xr[kx].p[i][j]*yi[ky].p[i][j]*ys*zi[kz].p[i][j]*zs-xi[kx].p[i][j]*yi[ky].p[i][j]*ys*zr[kz].p[i][j];
	      s=xr[kx].p[i][j]*yr[ky].p[i][j]*zi[kz].p[i][j]*zs+xi[kx].p[i][j]*yr[ky].p[i][j]*zr[kz].p[i][j]+xr[kx].p[i][j]*yi[ky].p[i][j]*ys*zr[kz].p[i][j]-xi[kx].p[i][j]*yi[ky].p[i][j]*ys*zi[kz].p[i][j]*zs;
	      SumKIm[j]+=s;
	      SumKRe[j]+=c;
	      ekr[ki].p[i][j]=c;
	      eki[ki].p[i][j]=s;
	    }
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
      Yk=kxyz[ki].ky*kxyz[ki].sy*FACT2;
      Zk=kxyz[ki].kz*kxyz[ki].sz*FACT2;
      /*Rksq is for one k vector k^2 * (2pi/L)^2*/
      RKsq=Xk*Xk+Yk*Yk+Zk*Zk;
      InvRKsq=1.0/RKsq;


      //      Veck=FACTVK*exp(-cons.SIGMA*RKsq)/RKsq;
      SKRe=SKreal[ki];
      SKIm=SKimag[ki];
		  
      //magni^2 of complex number SK=Sum_r(exp(irk))		  
      mag2=SKRe*SKRe+SKIm*SKIm;
      
      /*Potential term in K-space to electrostatic potential*/
      tSam=0.5*FCT*Veck[ki]*mag2;

      NNN.VK+=tSam;
      /*Contirbution to virial pressure (sum_r F*r) */ 
     
      VIRKxx+=tSam*(1.0-2.0*Xk*Xk*InvRKsq-Xk*Xk*A22I);
      VIRKyy+=tSam*(1.0-2.0*Yk*Yk*InvRKsq-Yk*Yk*A22I);
      VIRKzz+=tSam*(1.0-2.0*Zk*Zk*InvRKsq-Zk*Zk*A22I);
      for(i=0;i<cons.NperP;i++)
	{
	  for(j=0;j<3;j++)
	    {
	      RFact=-FCT*parms.Q[j+1]*Veck[ki];
	      CFact=-eki[ki].p[i][j]*SKRe+ekr[ki].p[i][j]*SKIm; //magnitude*sin(argSK-Thetrki[j]);
	      /*Force is derivative of Vk*/
	      Fxik=RFact*Xk*CFact;
	      Fyik=RFact*Yk*CFact;
	      Fzik=RFact*Zk*CFact;
	      force[i].x[j+1]+=Fxik;
	      force[i].y[j+1]+=Fyik;
	      force[i].z[j+1]+=Fzik;
	      Xp=local_mols[i].x[j+1]-com[i+my].x;
	      Yp=local_mols[i].y[j+1]-com[i+my].y;
	      Zp=local_mols[i].z[j+1]-com[i+my].z;
	      VIRMxx+=Xp*Fxik;
	      VIRMyy+=Yp*Fyik;
	      VIRMzz+=Zp*Fzik;
	    }
		
	}//end loop over myN particles

    }//End secondn loop over kvector



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
  for(i=0;i<cons.NperP;i++)
    {
  
      Rx1=local_mols[i].x[3]-local_mols[i].x[2];
      Ry1=local_mols[i].y[3]-local_mols[i].y[2];
      Rz1=local_mols[i].z[3]-local_mols[i].z[2];
      R2=Rx1*Rx1+Ry1*Ry1+Rz1*Rz1;
      R=sqrt(R2);
      TERM=parms.Q[3]*parms.Q[2];
      NNN.Vex+=TERM*erf(cons.ALPHA*R)/R;  
		
      Rx1=local_mols[i].x[3]-local_mols[i].x[1];
      Ry1=local_mols[i].y[3]-local_mols[i].y[1];
      Rz1=local_mols[i].z[3]-local_mols[i].z[1];
      R2=Rx1*Rx1+Ry1*Ry1+Rz1*Rz1;
      R=sqrt(R2);
      TERM=parms.Q[3]*parms.Q[1];
      NNN.Vex+=TERM*erf(cons.ALPHA*R)/R;  
      	
      Rx1=local_mols[i].x[2]-local_mols[i].x[1];
      Ry1=local_mols[i].y[2]-local_mols[i].y[1];
      Rz1=local_mols[i].z[2]-local_mols[i].z[1];
      R2=Rx1*Rx1+Ry1*Ry1+Rz1*Rz1;
      R=sqrt(R2);
      TERM=parms.Q[2]*parms.Q[1];
      NNN.Vex+=TERM*erf(cons.ALPHA*R)/R;  
      
		  		  
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

void Ksetup2(Constants cons, Kindex *kxyz, double &Vself, Parameters parms, double *Veck)
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
		      if(iy<0) kxyz[ind].sy=-1;
		      else kxyz[ind].sy=1;
		      kxyz[ind].ky=abs(iy);
		      if(iz<0) kxyz[ind].sz=-1;
		      else kxyz[ind].sz=1;
		      kxyz[ind].kz=abs(iz);
		      ind++;
		    }
		}
	    }
	}
    }
  /*Compute Vself: will not change, charges on mols do not change!
    Evaluate Values for Veck, K space Vectors do not change!*/
  double FACT2=2.0*cons.Pi*parms.InvBOXL;
  double FACTVK=2.0*FACT2*parms.InvBOXL*parms.InvBOXL;
  int i;
  double Xk, Yk, Zk, RKsq, InvRKsq;
  for(i=0;i<cons.Kleng;i++)
    {

      Xk=kxyz[i].kx*FACT2;
      Yk=kxyz[i].ky*FACT2;
      Zk=kxyz[i].kz*FACT2;
      /*Rksq is for one k vector k^2 * (2pi/L)^2*/
      RKsq=Xk*Xk+Yk*Yk+Zk*Zk;
      InvRKsq=1.0/RKsq;
      Veck[i]=FACTVK*exp(-cons.SIGMA*RKsq)/RKsq;
    }
  double qq;
  double alphPi=cons.ALPHA/sqrt(cons.Pi);
  qq=parms.Q[1]*parms.Q[1]+parms.Q[2]*parms.Q[2]+parms.Q[3]*parms.Q[3];
  Vself=cons.N*alphPi*qq;
  cout << "Vself; "<< Vself <<endl;

  /*End of Ksetup*/
}

void ComputeEnergyInit(Ensemble & NNN, double *Efinish, Constants cons, Parameters parms, DiffSend *vel, Thermostat &th)
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
  int i;

  kineCOM=CenterMassKinEnergy(cons, vel);
  NNN.Energyv[5]=kineCOM;
  
  kine=KineticA(cons, vel);
  kine=kine*0.5/418.4;
  NNN.Energyv[6]=kine;

/*Sum the energies across processors*/
  MPI::COMM_WORLD.Reduce(NNN.Energyv, Efinish, 7, MPI::DOUBLE, MPI::SUM, 0);
  /*0.5 Factor becasue everything in rwald is summed over all pairs twice
   */
  if(cons.rank==0)
    {
      Efinish[0]=Efinish[0]*.5;
      Efinish[1]=Efinish[1]*.5;
      Efinish[2]=Efinish[2]*.5;
      
      VKtot=NNN.VK-NNN.Vs-Efinish[3];
      NNN.EPOT=Efinish[0]+Efinish[1]+VKtot;
      
      NNN.EKIN=Efinish[6];
      NNN.ELRCint=NNN.ELRC/parms.Volume;
      /*Thermostat total Energy*/
      th.KE=0.5*th.QMass[0]*th.Veta[0]*th.Veta[0];
      th.PE=th.gnkt*th.Eta[0];
      for(i=1;i<cons.MC;i++)
	{
	  th.KE+=0.5*th.QMass[i]*th.Veta[i]*th.Veta[i];
	  th.PE+=th.gkt*th.Eta[i];
	}
      th.TE=(th.KE+th.PE)/418.4;
 
      NNN.Etotal=NNN.EPOT+NNN.EKIN+NNN.ELRCint+th.TE;
      
      
      NNN.KVirial=NNN.VirK-Efinish[4];
      VirTotal=Efinish[2]+NNN.KVirial;
      NNN.PLRCint=NNN.PLRC/parms.Volume;
      NNN.Pint=(Efinish[5]+NNN.PLRCint+VirTotal)*parms.pconstant;
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
void ComputeEnergy(Ensemble & NNN, double *Efinish, Constants cons, Parameters parms, DiffSend *vel, Thermostat &th)
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
  int i;

  kineCOM=CenterMassKinEnergy(cons, vel);
  NNN.Energyv[5]=kineCOM;
  
  kine=KineticA(cons, vel);
  kine=kine*0.5/418.4;
  NNN.Energyv[6]=kine;

/*Sum the energies across processors*/
  MPI::COMM_WORLD.Reduce(NNN.Energyv, Efinish, 7, MPI::DOUBLE, MPI::SUM, 0);
  /*0.5 Factor becasue everything in rwald is summed over all pairs twice
   */
  if(cons.rank==0)
    {
      Efinish[0]=Efinish[0]*.5;
      Efinish[1]=Efinish[1]*.5;
      Efinish[2]=Efinish[2]*.5;
      
      VKtot=NNN.VK-NNN.Vs-Efinish[3];
      NNN.EPOT=Efinish[0]+Efinish[1]+VKtot;
      
      NNN.EKIN=Efinish[6];
      //NNN.ELRCint=NNN.ELRC/parms.Volume; Does not change.
      th.KE=0.5*th.QMass[0]*th.Veta[0]*th.Veta[0];
      th.PE=th.gnkt*th.Eta[0];
      for(i=1;i<cons.MC;i++)
	{
	  th.KE+=0.5*th.QMass[i]*th.Veta[i]*th.Veta[i];
	  th.PE+=th.gkt*th.Eta[i];
	}
      th.TE=(th.KE+th.PE)/418.4;
      NNN.Etotal=NNN.EPOT+NNN.EKIN+NNN.ELRCint+th.TE;
      
      
      NNN.KVirial=NNN.VirK-Efinish[4];
      VirTotal=Efinish[2]+NNN.KVirial;
      //NNN.PLRCint=NNN.PLRC/parms.Volume; Does not change. 
      NNN.Pint=(Efinish[5]+NNN.PLRCint+VirTotal)*parms.pconstant;
      NNN.InstTemp=NNN.EKIN*2.0*418.4/parms.ConstM;
      
    }
  
  
}
void ScaleVelocity(Ensemble &NNN, Constants cons, Parameters parms, DiffSend *vel)
{
  /*Scale the velocities to the set temperature 
   Sets the Energy Ensemble Variable EKIN*/

  double Akin;
  double * myTemp=new double[1];
  double * InstTemp=new double[1];
  double Tf;
  int i;


  Akin=KineticA(cons, vel);
  myTemp[0]=Akin/parms.ConstM;
  /*Must sum InstTemp variable over all particle vels use Allreduce*/
  MPI::COMM_WORLD.Allreduce(myTemp, InstTemp, 1, MPI::DOUBLE, MPI::SUM);
  //cout<<"Temp pre scale: "<<InstTemp[0]<<endl;
  Tf=sqrt(cons.Temp/InstTemp[0]);
  //scale
  for(i=0;i<cons.NperP;i++)
    {
      vel[i].x[0]*=Tf;
      vel[i].y[0]*=Tf;
      vel[i].z[0]*=Tf;

      vel[i].x[1]*=Tf;
      vel[i].y[1]*=Tf;
      vel[i].z[1]*=Tf;

      vel[i].x[2]*=Tf;
      vel[i].y[2]*=Tf;
      vel[i].z[2]*=Tf;

    }
  delete[] myTemp;
  delete[] InstTemp;
}


void calcRDF(WaterSend *mols, double *gr, double boxl, int N, double delr)
{
  int i, j;
  double rx, ry, rz;
  double r2, r;
  int ig;
  double invboxl=1.0/boxl;
  for(i=0;i<N-1;i++)
    {
      for(j=i+1;j<N;j++)
	{
	  rx=mols[i].x[0]-mols[j].x[0];
	  ry=mols[i].y[0]-mols[j].y[0];
	  rz=mols[i].z[0]-mols[j].z[0];
	  rx-=boxl*round(rx*invboxl);
	  ry-=boxl*round(ry*invboxl);
	  rz-=boxl*round(rz*invboxl);
	  r2=rx*rx+ry*ry+rz*rz;
	  r=sqrt(r2);
	  if(r<boxl/2.0)
	  {
	    ig=int(r/delr);
	    gr[ig]+=2.0; //contribution for i and j
	  }//end if
	}
    }//End loop over pairs

}
