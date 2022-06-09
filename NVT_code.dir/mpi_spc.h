#ifndef __MPI_SPC_H
#define __MPI_SPC_H


#include <iostream>
#include "Parameters.h"
#include <complex>


using namespace std;



double GaussV();
void MatrixInv(double *Shake, double *Sinv);


class WaterSend
{
public:
  double x[4];
  double y[4];
  double z[4];
};
class DiffSend
{
 public:
  double x[3];
  double y[3];
  double z[3];
};
class Center
{
 public:
  double x;
  double y;
  double z;
};

class Kindex
{
 public:
  unsigned short int kx;
  unsigned short int ky;
  unsigned short int kz;
  short int sy;
  short int sz;
};
class Table
{
 public:
  double a;
  double b;
  double c;
  double d;
};
class Ensemble
{
 public:

  double * Energyv;
  double E0;
  double EC;
  double VK;
  double Vex;
  double Vs;
  double EKIN;
  double Etotal;
  double EPOT;
  double ELRC;
  double PLRC;
  double ELRCint;
  double PLRCint;
  double Pint;
  double InstTemp;

  double RealVir_LJ;
  double KVirial;
  double VirM;
  double VirK;
  

};  
  
//Function declarations
void readNVTParms(char* filename, Constants &cons, Parameters &parms, Thermostat&th);
void readParms(char * filename, Constants &cons, Parameters &parms);
void readCoords(char * coordFileIn, DiffSend *mols);
void readVelocity(char * velFileIn, DiffSend *vels);

//For debugging mostly
void writeForce(char * filename, Constants cons, DiffSend *force);
void writeAllForce(char * filename, Constants cons, DiffSend *force);
void printCoords(char *filename, Constants cons, DiffSend *mols, int t);
void printAllCoords(char *filename, Constants cons, DiffSend *mols);
void printAllMyCoords(char *filename, Constants cons, DiffSend *mols);
void printVel(char *filename, Constants cons, int t, DiffSend *vel);
void printAllVel(char *filename, Constants cons, int t, DiffSend *vel);
void writeParms(Parameters parms, Constants cons);
void checkIdeal(char * filename, Constants cons, Parameters parms, DiffSend *mols);  


void InitialParms(Parameters &parms, Constants &cons, Averages &avg, Ensemble &nnn);
void RigidMol(Parameters &parms, Constants cons, DiffSend *mols, DiffSend *opos, Center *com);

//void NVTintegrateMotion(Ensemble &NVT, Constants cons, Thermostat &therm);
//void NPTintegrateMotion(Ensemble &NPT, Barostat &baro, Thermostat &therm, Constants cons);
// void NHCpisoInt(Ensemble &NPT, Constants cons, Thermostat &therm, Barostat &baro);
void setThermostat(Thermostat &therm, Constants cons);
void NHCINT(DiffSend *vel, double Ekin, Constants cons, Thermostat &therm);
double KineticA(Constants cons, DiffSend *vel);
double CenterMassKinEnergy(Constants cons, DiffSend *vel);
void RWALD(Constants cons, Parameters parms, DiffSend *mols, DiffSend *local_mols, DiffSend *force, Center *com, Ensemble &nnn, Table *terf);

void Primder(Constants cons, Parameters parms, DiffSend *local_mols, DiffSend *force);

void Maxwell(Constants cons, Parameters parms, DiffSend *InitVel);
void VelocityNormal(Constants cons, Parameters parms, DiffSend *mols, DiffSend *InitVel, double *Mat, double *Minv);
void LongRangeCorrect(Constants cons, Parameters parms, Ensemble &nnn);
void CenterMassWater(Constants cons, DiffSend *mols, Center *com);

void RattlePos(Constants cons, Parameters parms, DiffSend *local_mols, DiffSend *myopos, DiffSend *vel, double *Mat, double *Minv);
void RattleVel(Constants cons, Parameters parms, DiffSend *local_mols, DiffSend *vel, double *Mat, double *Minv);
void TableSet(Table *terf, Constants cons, Parameters &parms);

#endif
