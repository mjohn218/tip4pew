#ifndef __ENSEMBLE_H
#define __ENSEMBLE_H


#include <iostream>
#include "Parameters.h"
#include <complex>
#define nperp 128

using namespace std;



double GaussV();
void MatrixInv(double ** Shake, double **Sinv);


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
class kspace
{
 public:
  double p[nperp][3];
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

void readParms(char * filename, Constants &cons, Parameters &parms);
void readCoords(char * coordFileIn, WaterSend *mols);
void readVelocity(char * velFileIn, DiffSend *vels);

//For debugging mostly
void writeForce(char * filename, Constants cons, WaterSend *force);
void writeAllForce(char * filename, Constants cons, WaterSend *force);
void printCoords(char *filename, Constants cons, WaterSend *mols, int t);
void printAllCoords(char *filename, Constants cons, WaterSend *mols);
void printAllMyCoords(char *filename, Constants cons, WaterSend *mols);
void printVel(char *filename, Constants cons, int t, DiffSend *vel);
void writeParms(Parameters parms, Constants cons);
void checkIdeal(char * filename, Constants cons, Parameters parms, WaterSend *mols);  


void InitialParms(Parameters &parms, Constants &cons, Averages &avg, Ensemble &nnn);
void RigidMol(Parameters &parms, Constants cons, WaterSend *mols, DiffSend *opos, Center *com);

//void NVTintegrateMotion(Ensemble &NVT, Constants cons, Thermostat &therm);
//void NPTintegrateMotion(Ensemble &NPT, Barostat &baro, Thermostat &therm, Constants cons);
// void NHCpisoInt(Ensemble &NPT, Constants cons, Thermostat &therm, Barostat &baro);
// void NHCInt(Ensemble &NVT, Constants cons, Thermostat &therm);
double KineticA(Constants cons, DiffSend *vel);
double CenterMassKinEnergy(Constants cons, DiffSend *vel);
void RWALD(Constants cons, Parameters parms, WaterSend *mols, WaterSend *local_mols, WaterSend *force, Center *com, Ensemble &nnn, Table *terf);

void Primder(Constants cons, Parameters parms, WaterSend *local_mols, WaterSend *force);

void Maxwell(Constants cons, Parameters parms, DiffSend *vel);
void VelocityNormal(Constants cons, Parameters parms, WaterSend *mols, DiffSend *vel);
void LongRangeCorrect(Constants cons, Parameters parms, Ensemble &nnn);
void CenterMassWater(Constants cons, WaterSend *mols, Center *com);

void RattlePos(Constants cons, Parameters parms, WaterSend *local_mols, DiffSend *myopos, DiffSend *vel);
void RattleVel(Constants cons, Parameters parms, WaterSend *local_mols, DiffSend *vel);
void TableSet(Table *terf, Constants cons, Parameters &parms);

#endif
