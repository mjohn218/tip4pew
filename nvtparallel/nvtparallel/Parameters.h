/*Header file to create classes to hold variables, parameters, and constants
for the Tip4pEw code
classes: Parameters, Thermostat, Barostat, and Constants

Margaret Johnson
2/2/05

*/

#ifndef __PARAMETERS_H
#define __PARAMETERS_H

#include <fstream>
#include <iostream>

using namespace std;

class Constants
{
public:
   int N;    //number of molecules
   int NperP; //number of molecules per processor
   int rank; //processor i'm on
   int nproc; //number of processors
   int NUMPRIM;   //number of atoms/molecule= 3 for water
   int NSITE;   //number of sites (3+1 fictitious)
   int NTC;      //preciser convergence for each chain~5
   int NYS;     //number of evolutions expansions(split the operator)~5 
   int MC;      //number of thermostats~5
   int NCONST;  //?
   int Equil;  //Flag for equilibration 1=true do it
   int Tsteps;  //total time steps
   int Tequil; //time steps for equilibration
   int Tvscale; //Time steps for velocity rescaling
   int Kmax;   //upper limit +1 of k-space coords ~11
   int Kmin;   //lower limit of k-space coords ~-10
   int Kscut;  //diameter of spheere in k-space we calculate values within
   int Kleng;  //number of different K vectors

   double R_LOWER;  // min radius for computing LJ interactions
   double R_UPPER;  // between r lower and this value, adjust LJ
  //potential for long range interactions, above this distance, LJ truncated
   double R_LOWER2;
   double R_UPPER2;
   double Temp;
  
   double BOLTZK;   //boltzman constant
   double KbDF;
   double Pi;  //yes it's 3.1415...
   double PiDEG; //180

   double ASPC; //parameter for LJ repulsion
   double CSPC; //parameter for LJ attraction

   //For Ewald sums
   double ALPHA;  //std (width) of gaussian around point charges
  //for ewald sums
   double SIGMA;
   double Alpha2;

   double dt;
   double dt2;
  
   double AMass[3];

};
class Thermostat
{
 public:
  double Veta[5];   //velocity of thermostat
  double Gt[5];      //a function of QMass and v, m, and Veta
  double Eta[5];    //thermostat positions
  double QMass[5];   //thermostat masses (a function of frequency. kg*m^2
  
  double taup;
  double gnkt, gkt;
  
  double wdti4[5], wdti8[5], wdti2[5]; //weights for evolution operators
  double KE; //kinetic energy
  double PE; //potential enrgy
  double TE; //total energy

};


class Barostat
{
public:
  double EPS;
  double GN1KT;
  double WEPS;
  double VEPS;
  double GEPS;
  double ODNF;

  //Coeffients for the maclaurin series of sinh(x)/x for
  //modified velocity verlet position update E2=1/3!, E4=1/5!,
  double E2, E4, E6, E8;
  void setBarostat(Barostat &baro);

};

class Parameters
{
 public:
  double BOXLENGTH;  //Angstrom ~24 for 512 molecules
  double InvBOXL;
  double BOXLOLD;
  double rhoold;
  double Volume; //will change for NPT

  //LJ params for water
  double ASPC;
  double CSPC;
  double RSPC;
  double ESPC;
  
  
  double MOLWT; //=18.01528;
  double AVAGAD; //=.602214199;
  double PIDEG; //=180.0;

  /*meant to speed up by saving a computation*/

  double CoeffPi; //2/sqrt(pi) used in RWALD in rajesh code
  double ConstM; //DNF *Boltzk for Temp scaling
  double pconstant; //1/(3*volume*.00001458)
  double rincr; // for lookup table
  double irincr;

  //water molecule signature
  double THETAH;  //=104.52;
  double ROH;  //=0.9572;
  double RHH;  //=1.5139006585;

  double d01sq;//=ROH^2
  double d02sq;//same
  double d12sq;//RHH^2


  double Rhomm;  // in units g/cm^3 (i.e. water at 1atm ~1g/cm^3
  double Density; // in units molecules/Angstrom^3 (rho*AVAGAD/MolWt)
  //rho is pressure dependent therfore so is Density

  double RM;  //=0.1250;
  double MASSH;  //=1.007276466;
  double MASSO;  //=16.00072706;

  double AMass[3];  //assumes all water molecules (no other soluute molecules
  // O, H, H
  double InvMass[3]; // 1/mass for each O, H, H
  
  double Q2, Q3, Q4, Q1; //charge Q1 is Oxygen is zero for TIP4P
  double Q[4];
  
  double AS;
  double BS;
  double CS;
  
  double tol;
  
  int isave; //saves frequently, just the oxygen
  int statwrite; //saves inst properties
  int allsave; //saves full pos and vel config, do infrequently
  int iterRat; //number of iterations to rattle
  int nperfile; //number of configs saved to each file

};

class Averages
{
 public:
  double Eavg;
  double E2avg;
  double Epavg;
  double Ep2avg;
  double Ekavg;
  double Tavg;
  double T2avg;
  double Pavg;
  double P2avg;
  



};




#endif
