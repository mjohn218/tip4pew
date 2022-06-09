#ifndef __ENSEMBLE_H
#define __ENSEMBLE_H


#include <iostream>
#include "Parameters.h"
#include <complex>


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
class Kspace
{
public:
  double KvectRe[200];
  double KvectIm[200];
};

class Ensemble
{
 public:
  
  double vx[256][3];
  double vy[256][3];
  double vz[256][3];

  double Fx[256][4];
  double Fy[256][4];
  double Fz[256][4];

  double Einitial;
  double E0;
  double EC;
  double VK;
  double EKIN;
  double Etotal;
  double EPOT;
  double ELRC;
  double PLRC;
  double ELRCint;
  double PLRCint;
  double Pint;
  double kineCOM;

  double RealVir_LJ;
  double KVir_LJ;
  double VirTotal;

  double XCOM[512]; //centers of mass
  double YCOM[512];
  double ZCOM[512];
  
  
  //Function declarations

  void readParms(char * filename, Constants &cons, Parameters &parms);
  void readCoords(char * coordFileIn, WaterSend *mols);
  void readVelocity(char * velFileIn, WaterSend *vels);
  
  //For debugging mostly
  void writeForce(char * filename, Constants cons);
  void writeAllForce(char * filename, Constants cons);
  void printCoords(char *filename, Constants cons, WaterSend *mols);
  void printAllCoords(char *filename, Constants cons, WaterSend *mols);
  void printAllMyCoords(char *filename, Constants cons, WaterSend *mols);
  void printVel(char *filename, Constants cons);
  void writeParms(Parameters parms, Constants cons);
  void checkIdeal(char * filename, Constants cons, Parameters parms, WaterSend *mols);  


  void InitialParms(Parameters &parms, Constants &cons, Averages &avg);
  void RigidMol(Parameters &parms, Constants cons, WaterSend *mols);
  
  void NVEintegrateMotionC(Ensemble &NVE, Constants cons, Parameters parms, WaterSend *mols, WaterSend *local_mols);

  //void NVTintegrateMotion(Ensemble &NVT, Constants cons, Thermostat &therm);
  //void NPTintegrateMotion(Ensemble &NPT, Barostat &baro, Thermostat &therm, Constants cons);
  // void NHCpisoInt(Ensemble &NPT, Constants cons, Thermostat &therm, Barostat &baro);
  // void NHCInt(Ensemble &NVT, Constants cons, Thermostat &therm);
  
  void calcEForce(Ensemble &NNN, Constants cons, Parameters parms, WaterSend *mols, WaterSend *local_mols);
  void ComputeAllEnergy(Ensemble &NNN, Constants cons, Parameters parms);
  double KineticA(Constants cons);
  double CenterMassKinEnergy(Constants cons);
  void ScaleVelocity(Ensemble &NNN, Constants cons, Parameters parms);
  void RWALD(Constants cons, Parameters parms, WaterSend *mols, WaterSend *local_mols);
  void KWALD(Constants cons, Parameters parms, WaterSend *mols, WaterSend *local_mols);
  void Primder(Constants cons, Parameters parms, WaterSend *mols, WaterSend *local_mols);

  void Maxwell(Constants cons, Parameters parms);
  void VelocityNormal(Constants cons, Parameters parms, WaterSend *mols);
  void LongRangeCorrect(Constants cons, Parameters parms);
  void CenterMassWater(Constants cons, WaterSend *mols);
  
  void RattlePos(Constants cons, Parameters parms, WaterSend *local_mols);
  void RattleVel(Constants cons, Parameters parms, WaterSend *local_mols);


 private:


};




#endif
