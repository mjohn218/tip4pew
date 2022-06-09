/*Contains cpp code for class entries Ensemble, declared in Ensemble.h
and used in Tip4p_w.cpp

    Date:1/24/2005
    Programmer: Margaret Johnson
*/
#include "Ensemble.h"

#include <iostream>
#include <iomanip>
#include "Parameters.h"
#include <ctime>
#include <cmath>
#include <complex>

inline double dotProd(double *v,double *s)
{
  double ans;
  ans=v[0]*s[0]+v[1]*s[1]+v[2]*s[2];
  return ans;

} 
double GaussV()
{

  /*Box mueller method for gaussian distributed random number from a uniform
random number generator~ rand()*/

  double R=2.0;
  double rnum;
  double V1, V2;

  while(R>=1.0)
    {

      V1=2.0*rand()/RAND_MAX-1.0;
      V2=2.0*rand()/RAND_MAX-1.0;
      R=V1*V1+V2*V2;
    }
  rnum=V1*sqrt(-2.0*log(R)/R);
  return rnum;

}

inline double nearestInt(double vec, double boxlen)
{
  int nInt;
  double div=vec/boxlen;
  if(vec>=0)
    {
      nInt=(int)(div+.5);
      nInt=-nInt;
    }      
  else
    nInt=(int)(.5-div);
  
  
  return (double)nInt;
}


void Ensemble::readParms(char* filename, Constants &cons, Parameters &parms)
{
  /*look in  parameter.h  for a description of what each value represents.*/
  
  ifstream PARMFILE(filename);
  if(!PARMFILE) cerr << " NO parm file given\n"; 

  PARMFILE >> cons.Tsteps;// 
  PARMFILE.ignore(400, '\n');
  PARMFILE >> cons.dt;  //
  PARMFILE.ignore(400, '\n');
  PARMFILE >> cons.Temp; //
  PARMFILE.ignore(400, '\n');
  PARMFILE >> cons.N;// 
  PARMFILE.ignore(400, '\n');
  PARMFILE >> parms.Rhomm;// 
  PARMFILE.ignore(400, '\n'); 
  PARMFILE >> cons.NUMPRIM;// 
  PARMFILE.ignore(400, '\n');
  PARMFILE >> cons.Equil;// 
  PARMFILE.ignore(400, '\n');
  PARMFILE >> parms.tol;// 
  PARMFILE.ignore(400, '\n');
  PARMFILE >> cons.NCONST; // 
  PARMFILE.ignore(400, '\n');
  PARMFILE >> cons.R_LOWER; //
  PARMFILE.ignore(400, '\n');
  PARMFILE >> cons.R_UPPER; // 
  PARMFILE.ignore(400, '\n');
  PARMFILE >> parms.RSPC; //
  PARMFILE.ignore(400, '\n');
  PARMFILE >> parms.ESPC; //
  PARMFILE.ignore(400, '\n');
  PARMFILE >> parms.Q2; //
  PARMFILE.ignore(400, '\n');
  PARMFILE >> parms.Q3;// 
  PARMFILE.ignore(400, '\n');
  PARMFILE >> parms.Q4;// 
  PARMFILE.ignore(400, '\n');
  PARMFILE >> parms.BOXLOLD;// 
  PARMFILE.ignore(400, '\n');
  PARMFILE.close();
}

void Ensemble::readCoords(char *coordFileIn, WaterSend *mols)
{

  /*i.e. N=512, NATOMS=3 for water*/

  int i, N;
  ifstream COORDFILE(coordFileIn);
  if(!COORDFILE) cerr << "No coord file given\n";
  COORDFILE >> N; 

  for(i=0;i<N;i++)
    {
      COORDFILE>>mols[i].x[0]>>mols[i].y[0]>>mols[i].z[0]; //oxygen 
      COORDFILE>>mols[i].x[1]>>mols[i].y[1]>>mols[i].z[1]; //hydrogen
      COORDFILE>>mols[i].x[2]>>mols[i].y[2]>>mols[i].z[2]; //hydrogen
     }
  COORDFILE.close();
}

void Ensemble::readVelocity(char *velFileIn, WaterSend *vels)
{

  /*i.e. N=512, NATOMS=3 for water*/

  int i, N;
  ifstream VELFILE(velFileIn);
  if(!VELFILE) cerr << "No vel file given\n";
  VELFILE >> N; 

  for(i=0;i<N;i++)
    {
      VELFILE>>vels[i].x[0]>>vels[i].y[0]>>vels[i].z[0]; //ovygen 
      VELFILE>>vels[i].x[1]>>vels[i].y[1]>>vels[i].z[1]; //hydrogen
      VELFILE>>vels[i].x[2]>>vels[i].y[2]>>vels[i].z[2]; //hydrogen
     }
  VELFILE.close();
}


void Ensemble::InitialParms(Parameters &parms, Constants &cons, Averages &avg)
{
  /*Only after read in certain parameters/constants from file*/  

  double RUP_LOW;
  double RUP_LOW2;
  double RUP_LOW3;
  double RUP_LOW4;
  double RUP_LOW5;
  double RMFact=332.0637349768129609;
  //Initialize parameters to set values (i.e.atom masses, etc) (Parameters.cpp)
  

  parms.MOLWT=18.01528;
  parms.AVAGAD=.602214199;
  parms.PIDEG=180.0;

  //water molecule signature
  parms.THETAH=104.52;
  parms.ROH=0.9572;
  parms.RHH=1.5139006585;

  parms.RM=0.1250;
  parms.MASSH=1.007276466;
  parms.MASSO=16.00072706;

  parms.AMass[0]=16.00072706; //Oxygen
  parms.AMass[1]=1.007276466;  //Hydrogen
  parms.AMass[2]=1.007276466;  //Hydrogen

  parms.InvMass[0]=0.06249716005093; //1/massOxygen
  parms.InvMass[1]=0.99277609847384;
  parms.InvMass[2]=0.99277609847384;

  
  cons.BOLTZK=0.83144703385;
  cons.Pi=acos(-1.0);
  cons.dt2=cons.dt/2.0;
  cons.PiDEG=180.0;

  //3+1 fictitious
  cons.NSITE=4;

  //Same as atomic masses AMass in parms, repeated for ease in passing classes
  cons.AMass[0]=parms.AMass[0];
  cons.AMass[1]=parms.AMass[1];
  cons.AMass[2]=parms.AMass[2];

  /*Powers of rupper - rlower*/
  cons.R_LOWER2=cons.R_LOWER*cons.R_LOWER;
  cons.R_UPPER2=cons.R_UPPER*cons.R_UPPER;
  RUP_LOW=cons.R_UPPER2-cons.R_LOWER2;
  RUP_LOW2=RUP_LOW*RUP_LOW;
  RUP_LOW3=RUP_LOW2*RUP_LOW;
  RUP_LOW4=RUP_LOW2*RUP_LOW2;
  RUP_LOW5=RUP_LOW4*RUP_LOW;

  parms.AS=-10.0/RUP_LOW3;
  parms.BS=15.0/RUP_LOW4;
  parms.CS=-6.0/RUP_LOW5;

  /*Density is numMolecules/vol  is pressure dependent*/

  cout << "in initial: rho: "<<parms.Rhomm<<endl;
  parms.Density=parms.Rhomm*parms.AVAGAD/parms.MOLWT;
  parms.BOXLENGTH=pow(cons.N*1.0/parms.Density,1.0/3.0);
  parms.InvBOXL=1.0/parms.BOXLENGTH;
  parms.Volume=parms.BOXLENGTH*parms.BOXLENGTH*parms.BOXLENGTH;

  cout << "MOLECULE Density: " << parms.Density <<endl;
  cout << "BOXLENGTH: "  << parms.BOXLENGTH << endl;
  cout << "VOLUME "<<parms.Volume << endl;
  
  

  /*Set values of parameters Q is charge 2, 3 are H  4 is Oxygen*/  
  double sqrF=sqrt(RMFact);
  parms.Q2=parms.Q2*sqrF;
  parms.Q3=parms.Q3*sqrF;
  parms.Q4=parms.Q4*sqrF;

  /*Coeffieicents for LJ potential*/  
  double RSPC2, RSPC6, RSPC12;
  RSPC2=parms.RSPC*parms.RSPC;
  RSPC6=RSPC2*RSPC2*RSPC2;
  RSPC12=RSPC6*RSPC6;
  
  parms.ASPC=4.0*parms.ESPC*RSPC12;
  parms.CSPC=4.0*parms.ESPC*RSPC6;


  //charge values for each atom
  parms.Q[0]=0.0;
  parms.Q[1]=parms.Q2;
  parms.Q[2]=parms.Q3;
  parms.Q[3]=parms.Q4;


  /*Initialize enery and temperature quantities*/
  
  parms.iterRat=50;
  
  cons.ALPHA=0.35;
  double TwoPi=cons.Pi*cons.Pi;
  cons.Alpha2=cons.ALPHA*cons.ALPHA;
  
  double PFACT=0.0000145839756;
  parms.pconstant=1.0/(3.0*parms.Volume*PFACT);
  cons.SIGMA=1.0/(4.0*cons.Alpha2);
  parms.CoeffPi=2.0/sqrt(cons.Pi);
  double DNF=cons.N*(3.0*cons.NUMPRIM-1.0*cons.NCONST)-3.0;
  parms.ConstM=DNF*cons.BOLTZK;
  parms.Icheck=1;
  parms.Jcheck=1;
  parms.Kcheck=1;
  avg.ETave=0.0;
  avg.ETave2=0.0;
  avg.EPave=0.0;
  avg.EPave2=0.0;
  avg.Tave=0.0;
  avg.Tave2=0.0;
  avg.Pave=0.0;
  avg.Pave2=0.0;
  parms.Chi=0.0;
  parms.Isse=0;
  
  parms.Ichunk=0;

  /*Reciprocal space quantities*/
  
  

}//End of function Initial


void Ensemble::RigidMol(Parameters &parms, Constants cons, WaterSend *mols)
{
  /*Scale the initial configuration to a given boxlength or density, find center of mass 
and establish teh position of the imaginary 4th 'atom' of water*/

  double ThetaRad, Thetad;
  double XH20, YH20, ZH20;
  double XH10, YH10, ZH10;
  double XHH, YHH, ZHH;
  double R10sq, R10;
  double R20sq, R20;
  double RHHsq, RHH;
  double cosThetaRad;
  double Fixd[3][3];


  /*Scale Factor*/
  double SCBOXL;
  SCBOXL=parms.BOXLENGTH/parms.BOXLOLD;

  //scale the old configuration
  //find center of mass

  double sumx=0.0;
  double sumy=0.0;
  double sumz=0.0;
  double asum=0.0;
  int i;
  double am0=parms.AMass[0];
  double am1=parms.AMass[1];
  double am2=parms.AMass[2];

  double RB, XB, YB, ZB;

  for(i=0;i<cons.N;i++)
    {
      sumx=sumy=sumz=asum=0.0;
      sumx=am0*mols[i].x[0]+am1*mols[i].x[1]+am2*mols[i].x[2];
      sumy=am0*mols[i].y[0]+am1*mols[i].y[1]+am2*mols[i].y[2];
      sumz=am0*mols[i].z[0]+am1*mols[i].z[1]+am2*mols[i].z[2];
      asum=am0+am1+am2;
      
      //COM
      XCOM[i]=sumx/asum;
      YCOM[i]=sumy/asum;
      ZCOM[i]=sumz/asum;
      
      //new positions O, H, H
      mols[i].x[0]+=(SCBOXL-1.0)*XCOM[i];
      mols[i].y[0]+=(SCBOXL-1.0)*YCOM[i];
      mols[i].z[0]+=(SCBOXL-1.0)*ZCOM[i];

      mols[i].x[1]+=(SCBOXL-1.0)*XCOM[i];
      mols[i].y[1]+=(SCBOXL-1.0)*YCOM[i];
      mols[i].z[1]+=(SCBOXL-1.0)*ZCOM[i];

      mols[i].x[2]+=(SCBOXL-1.0)*XCOM[i];
      mols[i].y[2]+=(SCBOXL-1.0)*YCOM[i];
      mols[i].z[2]+=(SCBOXL-1.0)*ZCOM[i];
      
            
      /*Distance Hyd1-Oxygen*/
      XH10=mols[i].x[1]-mols[i].x[0];
      YH10=mols[i].y[1]-mols[i].y[0];
      ZH10=mols[i].z[1]-mols[i].z[0];
      
      XH10+=parms.BOXLENGTH*nearestInt(XH10,parms.BOXLENGTH);
      YH10+=parms.BOXLENGTH*nearestInt(YH10,parms.BOXLENGTH);
      ZH10+=parms.BOXLENGTH*nearestInt(ZH10,parms.BOXLENGTH);
      
      R10sq=XH10*XH10+YH10*YH10+ZH10*ZH10;
      R10=sqrt(R10sq);

      /*Distance Hyd2-Oxygen*/
      XH20=mols[i].x[2]-mols[i].x[0];
      YH20=mols[i].y[2]-mols[i].y[0];
      ZH20=mols[i].z[2]-mols[i].z[0];

      XH20+=parms.BOXLENGTH*nearestInt(XH20,parms.BOXLENGTH);
      YH20+=parms.BOXLENGTH*nearestInt(YH20,parms.BOXLENGTH);
      ZH20+=parms.BOXLENGTH*nearestInt(ZH20,parms.BOXLENGTH);
      
      R20sq=XH20*XH20+YH20*YH20+ZH20*ZH20;
      R20=sqrt(R20sq);

      /*Angle between vectors OH1 and OH2*/
      double DotProd=XH20*XH10+YH20*YH10+ZH20*ZH10;
      cosThetaRad=DotProd/R20/R10;
      ThetaRad=acos(cosThetaRad);
      Thetad=ThetaRad*(cons.PiDEG/cons.Pi);
      
      /*vector of H-H (no bond there)*/
      XHH=mols[i].x[2]-mols[i].x[1];
      YHH=mols[i].y[2]-mols[i].y[1];
      ZHH=mols[i].z[2]-mols[i].z[1];
      
      XHH+=parms.BOXLENGTH*nearestInt(XHH, parms.BOXLENGTH);
      YHH+=parms.BOXLENGTH*nearestInt(YHH, parms.BOXLENGTH);
      ZHH+=parms.BOXLENGTH*nearestInt(ZHH, parms.BOXLENGTH);
      
      RHHsq=XHH*XHH+YHH*YHH+ZHH*ZHH;
      RHH=sqrt(RHHsq);
      
      //These values of R and Thetad should remain ideal or errors! 
      
      /*Now find 4th positions*/

      XB=0.5*(XH10+XH20);
      YB=0.5*(YH10+YH20);
      ZB=0.5*(ZH10+ZH20);
      RB=sqrt(XB*XB+YB*YB+ZB*ZB);
      /*No use parallelizing this, need x[3] etc. for all processors*/
      mols[i].x[3]=mols[i].x[0]+parms.RM*XB/RB;
      mols[i].y[3]=mols[i].y[0]+parms.RM*YB/RB;
      mols[i].z[3]=mols[i].z[0]+parms.RM*ZB/RB;
      
    }
  /*Distances between O [0] H [1] and H [2]. OH are ideal, and HH is calculated above*/
  Fixd[0][1]=parms.ROH;
  Fixd[1][0]=Fixd[0][1];
  Fixd[0][2]=parms.ROH;
  Fixd[2][0]=Fixd[0][2];
  Fixd[1][2]=RHH;
  Fixd[2][1]=Fixd[1][2];

  parms.d01sq=parms.ROH*parms.ROH;
  parms.d02sq=parms.ROH*parms.ROH;
  parms.d12sq=RHH*RHH;


}//end rigid mol


void Ensemble::writeParms(Parameters parms, Constants cons)
{
  cout.precision(16);

  cout << "Total steps: " << cons.Tsteps<< endl;
  cout <<" dt: " << cons.dt<< endl;
  cout << "Temp: "<< cons.Temp<< endl;
  cout << "num molecules: "<< cons.N<< endl;
  cout << "density g/cm3: "<< parms.Rhomm<<endl;
  cout << "atoms/molecule: "<< cons.NUMPRIM<<endl;
  cout << "equilibration flag: "<< cons.Equil<<endl;
  cout << "tolerance: "<< parms.tol<< endl;
  cout << "n constant: "<< cons.NCONST<< endl;
  cout << "lj lower radius cutoff: "<< cons.R_LOWER<< endl;
  cout << "lj upper radius: "<< cons.R_UPPER<< endl;
  cout << "RSPC for LJ coeff: "<< parms.RSPC<<endl;
  cout << "ESPC for LJ: "<< parms.ESPC<<endl;
  cout << "Q2:  "<< parms.Q2<<endl;
  cout << "Q3: "<< parms.Q3<<endl;
  cout << "Q4: "<< parms.Q4<<endl;
  cout << "theta water: " <<parms.THETAH <<endl;
  cout << "OH length: " << parms.ROH <<endl;
  cout << "mass Hydrogen: " << parms.MASSH <<endl;



}

void Ensemble::printCoords(char *filename, Constants cons, WaterSend *mols)
{
  int i;
  ofstream outFile(filename);

  for(i=0;i<cons.N;i++)
    {
      outFile.precision(16);
      outFile << mols[i].x[0] << "  " <<mols[i].y[0]<< "  " <<mols[i].z[0]<<endl;
      outFile << mols[i].x[1] << "  " <<mols[i].y[1]<< "  " <<mols[i].z[1]<<endl;
      outFile << mols[i].x[2] << "  " <<mols[i].y[2]<< "  " <<mols[i].z[2]<<endl;
    }
  outFile.close();
}

void Ensemble::printAllCoords(char *filename, Constants cons, WaterSend *mols)
{
  int i;
  ofstream outFile(filename);

  for(i=0;i<cons.N;i++)
    {
      outFile.precision(16);
      outFile << mols[i].x[0] << "  " <<mols[i].y[0]<< "  " <<mols[i].z[0]<<endl;
      outFile << mols[i].x[1] << "  " <<mols[i].y[1]<< "  " <<mols[i].z[1]<<endl;
      outFile << mols[i].x[2] << "  " <<mols[i].y[2]<< "  " <<mols[i].z[2]<<endl;
      outFile << mols[i].x[3] << "  " <<mols[i].y[3]<< "  " <<mols[i].z[3]<<endl;
    }
  outFile.close();
}

void Ensemble::printAllMyCoords(char *filename, Constants cons, WaterSend *mols)
{
  int i;
  ofstream outFile(filename);

  for(i=0;i<cons.NperP;i++)
    {
      outFile.precision(16);
      outFile << mols[i].x[0] << "  " <<mols[i].y[0]<< "  " <<mols[i].z[0]<<endl;
      outFile << mols[i].x[1] << "  " <<mols[i].y[1]<< "  " <<mols[i].z[1]<<endl;
      outFile << mols[i].x[2] << "  " <<mols[i].y[2]<< "  " <<mols[i].z[2]<<endl;
      outFile << mols[i].x[3] << "  " <<mols[i].y[3]<< "  " <<mols[i].z[3]<<endl;
    }
  outFile.close();
}

void Ensemble::writeForce(char *filename, Constants cons)
{
  int i;
  ofstream outFile(filename);

  for(i=0;i<cons.NperP;i++)
    {
      outFile.precision(16);
      outFile << Fx[i][0] << "  " <<Fy[i][0]<< "  " <<Fz[i][0]<<endl;
      outFile << Fx[i][1] << "  " <<Fy[i][1]<< "  " <<Fz[i][1]<<endl;
      outFile << Fx[i][2] << "  " <<Fy[i][2]<< "  " <<Fz[i][2]<<endl;
    }
  outFile.close();
}

void Ensemble::writeAllForce(char *filename, Constants cons)
{
  int i;
  ofstream outFile(filename);

  for(i=0;i<cons.NperP;i++)
    {
      outFile.precision(16);
      outFile << Fx[i][0] << "  " <<Fy[i][0]<< "  " <<Fz[i][0]<<endl;
      outFile << Fx[i][1] << "  " <<Fy[i][1]<< "  " <<Fz[i][1]<<endl;
      outFile << Fx[i][2] << "  " <<Fy[i][2]<< "  " <<Fz[i][2]<<endl;
      outFile << Fx[i][3] << "  " <<Fy[i][3]<< "  " <<Fz[i][3]<<endl;
    }
  outFile.close();
}

void Ensemble::printVel(char *filename, Constants cons)
{
  int i;
  ofstream outFile(filename);

  for(i=0;i<cons.NperP;i++)
    {
      outFile.precision(16);
      outFile << vx[i][0] << "  " <<vy[i][0]<< "  " <<vz[i][0]<<endl;
      outFile << vx[i][1] << "  " <<vy[i][1]<< "  " <<vz[i][1]<<endl;
      outFile << vx[i][2] << "  " <<vy[i][2]<< "  " <<vz[i][2]<<endl;
    }
  outFile.close();
}


void Ensemble::NVEintegrateMotionC(Ensemble &NVE, Constants cons, Parameters parms, WaterSend *mols, WaterSend *local_mols)
{
  /*NVE Velocity Verlet*/
  int i;
  for(i=0;i<cons.NperP;i++)
    {
      //update particle velocities (O, H, H)
      vx[i][0]+=cons.dt2*parms.InvMass[0]*Fx[i][0];
      vy[i][0]+=cons.dt2*parms.InvMass[0]*Fy[i][0];
      vz[i][0]+=cons.dt2*parms.InvMass[0]*Fz[i][0];

      vx[i][1]+=cons.dt2*parms.InvMass[1]*Fx[i][1];
      vy[i][1]+=cons.dt2*parms.InvMass[1]*Fy[i][1];
      vz[i][1]+=cons.dt2*parms.InvMass[1]*Fz[i][1];
      
      vx[i][2]+=cons.dt2*parms.InvMass[2]*Fx[i][2];
      vy[i][2]+=cons.dt2*parms.InvMass[2]*Fy[i][2];
      vz[i][2]+=cons.dt2*parms.InvMass[2]*Fz[i][2];
      

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
    
  
  //calculate force
  NVE.calcEForce(NVE, cons, parms, mols, local_mols);
      
  for(i=0;i<cons.NperP;i++)
    {
      //update particle velocities
      vx[i][0]+=cons.dt2*parms.InvMass[0]*Fx[i][0];
      vy[i][0]+=cons.dt2*parms.InvMass[0]*Fy[i][0];
      vz[i][0]+=cons.dt2*parms.InvMass[0]*Fz[i][0];
      
      vx[i][1]+=cons.dt2*parms.InvMass[1]*Fx[i][1];
      vy[i][1]+=cons.dt2*parms.InvMass[1]*Fy[i][1];
      vz[i][1]+=cons.dt2*parms.InvMass[1]*Fz[i][1];
      
      vx[i][2]+=cons.dt2*parms.InvMass[2]*Fx[i][2];
      vy[i][2]+=cons.dt2*parms.InvMass[2]*Fy[i][2];
      vz[i][2]+=cons.dt2*parms.InvMass[2]*Fz[i][2];
    } 
  //apply constraints
  cout <<"start rattlevel:\n ";
  NVE.RattleVel(cons, parms, local_mols);
    

}



void Ensemble::calcEForce(Ensemble &NNN, Constants cons, Parameters parms, WaterSend *mols, WaterSend *local_mols)
{
  //calculate energy: derivative of the force
  //forces are short range (lennard jones) and long range (coulombic)
  //also external forces, in this case no external forces
  // reads in also AMASS

  double start, end;
 

  int i, j;

  double CONSTANT=418.4; //joules conversion



  for(i=0;i<cons.NperP;i++)
    {
      Fx[i][0]=0.0;
      Fy[i][0]=0.0;
      Fz[i][0]=0.0;
      
      Fx[i][1]=0.0;
      Fy[i][1]=0.0;
      Fz[i][1]=0.0;
      
      Fx[i][2]=0.0;
      Fy[i][2]=0.0;
      Fz[i][2]=0.0;

      Fx[i][3]=0.0;
      Fy[i][3]=0.0;
      Fz[i][3]=0.0;


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
  start=clock();
  NNN.RWALD(cons, parms, mols, local_mols);
  end=clock();
  cout <<"rwald: "<<(start-end)/CLOCKS_PER_SEC <<endl;
  /*  cout <<"fx0: " <<Fx[0][1]<<endl;
  cout <<"fy0: " <<Fy[0][1]<<endl;
  cout <<"fz0: " <<Fz[0][1]<<endl;
  */
  char * myforcefile="rwaldForce.out";
  NNN.writeAllForce(myforcefile, cons);

  //find n bonded energy and derivatives in reciprocal space
  
  /*
  start=clock();
  NNN.KWALD(cons, parms, mols, local_mols);
  end=clock();
  cout <<"Kwald: "<<(start-end)/CLOCKS_PER_SEC <<endl;
  */


  /* cout <<"fx0: " <<Fx[0][1]<<endl;
  cout <<"fy0: " <<Fy[0][1]<<endl;
  cout <<"fz0: " <<Fz[0][1]<<endl;
  */


  //evaluate forces on primary atoms due to secondary atom
  //and update the forces on all real atoms


  NNN.Primder(cons, parms, local_mols);
  
  /*  cout <<"primder\n";
  cout <<"fx0: " <<Fx[0][1]<<endl;
  cout <<"fy0: " <<Fy[0][1]<<endl;
  cout <<"fz0: " <<Fz[0][1]<<endl;
  */

  /*get energy sum: E0 from Rwald, EC from Rwald, Vk from Kwald*/


  EPOT=E0+EC+VK;
  cout << "E0: "<< E0 << "EC: "<<EC<<"VK: "<<VK<<endl;
  cout <<"actual total potential energy: "<< EPOT <<endl;

  //forces in units for the verlet equation
  double Fxtotal=0.0;
  double Fytotal=0.0;
  double Fztotal=0.0;


  for(i=0;i<cons.NperP;i++)
    { 
      Fx[i][0]*=CONSTANT;
      Fy[i][0]*=CONSTANT;
      Fz[i][0]*=CONSTANT;

      Fx[i][1]*=CONSTANT;
      Fy[i][1]*=CONSTANT;
      Fz[i][1]*=CONSTANT;

      Fx[i][2]*=CONSTANT;
      Fy[i][2]*=CONSTANT;
      Fz[i][2]*=CONSTANT;
      
      Fxtotal+=Fx[i][0]+Fx[i][1]+Fx[i][2];
      Fytotal+=Fy[i][0]+Fy[i][1]+Fy[i][2];
      Fztotal+=Fz[i][0]+Fz[i][1]+Fz[i][2];
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


void Ensemble::RWALD(Constants cons, Parameters parms, WaterSend *mols, WaterSend *local_mols)
{
  //Energy and derivatives in real space
  /*Want E0  EC  RealVir_LJ */
  
  //Rlower is the lower limit of the distance between the 2 molecules
  //which is 9 in this case?


  //do solvent-solvent interactions
  //this part must take into consideration that some of your neighbor 
  //molecules might be that way because they are close to you on the 
  //boundaries, i.e. if you are close to the edge of the box, and they are
  // opposite you then you will interact with them  

    
  //E0, EC, ALPHA, 
  //FOR WATER ASPC, CSPC, RSPC, ESPC, Q(N, 4)

  
  //AS, BS, CS are all set in function sim parms
  int i;
  
  for(i=0;i<cons.NperP;i++)
    {
      Fx[i][0]=0.0;
      Fy[i][0]=0.0;
      Fz[i][0]=0.0;
      
      Fx[i][1]=0.0;
      Fy[i][1]=0.0;
      Fz[i][1]=0.0;
      
      Fx[i][2]=0.0;
      Fy[i][2]=0.0;
      Fz[i][2]=0.0;

      Fx[i][3]=0.0;
      Fy[i][3]=0.0;
      Fz[i][3]=0.0;


    }
  /*
  if(cons.rank==0)
    {
      cout << "IN RWALD \n" <<endl;
      writeParms(parms, cons);
    }
  */
  int j, k, n;
  double Fx0, Fy0, Fz0;
  double FXnew, FYnew, FZnew;
  double rX0, rY0, rZ0;
  double RXnew, RYnew, RZnew;
  double xCMnow, yCMnow, zCMnow;

  double r02;
  double r2, r6, r12;
  /*Return these to main*/
  E0=0.0;
  EC=0.0;
  double Vir_LrXX=0.0;  //sum over VirijXX for all pairs
  double Vir_LrYY=0.0;
  double Vir_LrZZ=0.0;


  double VirijXX, VirijYY, VirijZZ;
  
  
  double ELJ;   //lennard jones potential

  //Temporary values to speed computing  
  double SumQQ=0.0;
  double COMFAC, COMFAC1, COMFAC2;
  double Term1, Term, Term1A, Term1AA;
  double Q1, Rij1, RijSQUARED;

  //FActors for intermediate R zone
  double ELJO;
  double ZR, ZRSquared, ZRCubed, ZRFourth, ZRFifth;
  double SZR, DSZR;

  double Tmpp;
  /*Update all forces on i, loop j over all molecules, or create neighbor
   list and loop j over those.  
   do not interact with self: hence j!=(i+rank*cons.NperP) expression 
  */
  //  ofstream atomfile("atoms.out");
  cout << "rwald boxlength: "<<parms.BOXLENGTH <<endl;
  cout << "NSITE: "<<cons.NSITE <<endl;
  for(i=0;i<cons.NperP;i++)
    {
      for(j=0;j<cons.N;j++)
	{
	   if( j!=(i+cons.rank*cons.NperP))
	     {
	       //Distance between 2 oxygen molecules
	  rX0=local_mols[i].x[0]-mols[j].x[0];
	  rY0=local_mols[i].y[0]-mols[j].y[0];
	  rZ0=local_mols[i].z[0]-mols[j].z[0];
	  
	  //Box is Periodic
	  //i.e. the molecules cannnot be further away then half the
	  //boxlength in 1 dimension, otherwise they are closer to the 
	  //periodic image of te molecule.

	  rX0=rX0+parms.BOXLENGTH*nearestInt(rX0,parms.BOXLENGTH);
	  rY0=rY0+parms.BOXLENGTH*nearestInt(rY0,parms.BOXLENGTH);
	  rZ0=rZ0+parms.BOXLENGTH*nearestInt(rZ0,parms.BOXLENGTH);
	  r02=rX0*rX0+rY0*rY0+rZ0*rZ0; //r^2=(x-x)^2+(y-y)^2+(z-z)^2

	  SumQQ=0.0;
	  VirijXX=0.0;
	  VirijYY=0.0;
	  VirijZZ=0.0;

	  //Compute the lennard jones potential between the 2 
	  //molecules if they are close enough together.(the oxygen 
	  //atoms distance) 

	  if(r02<=cons.R_LOWER2)
	    {
	      r2=1.0/r02; // 1/r^2 term for LJ potential
	      r6=r2*r2*r2;
	      r12=r6*r6;
	      
	      /*Lennard Jones Potential contribution
	       force is derivative of LJ*/

	      ELJ=(parms.ASPC*r12-parms.CSPC*r6);
	      COMFAC=(-12.0*parms.ASPC*r12+6.0*parms.CSPC*r6)*r2;
	      Fx0= -COMFAC*rX0;
	      Fy0= -COMFAC*rY0;
	      Fz0= -COMFAC*rZ0;
	      E0=E0+ELJ;

	      Fx[i][0]=Fx[i][0]+Fx0;
	      Fy[i][0]=Fy[i][0]+Fy0;
	      Fz[i][0]=Fz[i][0]+Fz0;

	       
	      //distance between 2 molecules Center Of Mass
	      xCMnow=XCOM[i]-XCOM[j];
	      yCMnow=YCOM[i]-YCOM[j];
	      zCMnow=ZCOM[i]-ZCOM[j];
	      xCMnow+=parms.BOXLENGTH*nearestInt(xCMnow,parms.BOXLENGTH);
	      yCMnow+=parms.BOXLENGTH*nearestInt(yCMnow,parms.BOXLENGTH);
	      zCMnow+=parms.BOXLENGTH*nearestInt(zCMnow,parms.BOXLENGTH);
	      
	      VirijXX+=Fx0;
	      VirijYY+=Fy0;
	      VirijZZ+=Fz0;
	      
	      //loop over other atoms in the molecules (Hydrogens) and imag
	      //potential is now a short range coulomb, hence erfc, no LJ
	      for(k=1;k<cons.NSITE;k++)
		{
		  for(n=1;n<cons.NSITE;n++)
		    {
		      RXnew=local_mols[i].x[k]-mols[j].x[n];
		      RYnew=local_mols[i].y[k]-mols[j].y[n];
		      RZnew=local_mols[i].z[k]-mols[j].z[n];
		      RXnew+=parms.BOXLENGTH*nearestInt(RXnew,parms.BOXLENGTH);
		      RYnew+=parms.BOXLENGTH*nearestInt(RYnew,parms.BOXLENGTH);
		      RZnew+=parms.BOXLENGTH*nearestInt(RZnew,parms.BOXLENGTH);
		      RijSQUARED=RXnew*RXnew+RYnew*RYnew+RZnew*RZnew;
		      Rij1=sqrt(RijSQUARED);
		      Q1=parms.Q[k]*parms.Q[n];
		      Term=erfc(cons.ALPHA*Rij1)/Rij1;
		      //complementary error function
		      Term1=Q1*Term;
		      Tmpp=parms.CoeffPi*cons.ALPHA*exp(-cons.Alpha2*RijSQUARED)+Term;
		      Term1A=Q1*Tmpp/RijSQUARED;
		      FXnew=Term1A*RXnew;
		      FYnew=Term1A*RYnew;
		      FZnew=Term1A*RZnew;
		      SumQQ+=Term1;

		      Fx[i][k]+=FXnew;
		      Fy[i][k]+=FYnew;
		      Fz[i][k]+=FZnew;
		      
		      
		      VirijXX+=FXnew;
		      VirijYY+=FYnew;
		      VirijZZ+=FZnew;

		    }//end loop over H1
		}//end loop over H2

	      VirijXX*=xCMnow;
	      VirijYY*=yCMnow;
	      VirijZZ*=zCMnow;


	    }


	  /*In this case, you are between the 2 radii, and you have to size everything by a polynomial expression.otherwise function discontinuous at boundary*/

	  else if(r02>cons.R_LOWER2 && r02 <=cons.R_UPPER2)
	    {

	      //now you are in the range where this extra S term
	      //is calculated
	      r2=1.0/r02;
	      r6=r2*r2*r2;
	      r12=r6*r6;
	      
	      //Lennard Jones Potential
	      ELJO=(parms.ASPC*r12-parms.CSPC*r6);
	      ZR=r02-cons.R_LOWER2;
	      ZRSquared=ZR*ZR;
	      ZRCubed=ZRSquared*ZR;
	      ZRFourth=ZRCubed*ZR;
	      ZRFifth=ZRFourth*ZR;

	      //some polynomial expression
	      SZR=1.0+parms.AS*ZRCubed+parms.BS*ZRFourth+parms.CS*ZRFifth;
	      DSZR=6.0*parms.AS*ZRSquared+8.0*parms.BS*ZRCubed+10.0*parms.CS*ZRFourth;
	      
	      //LJ is adjusted in this case
	      
	      ELJ=ELJO*SZR;
	      COMFAC1=(-12.0*parms.ASPC*r12+6.0*parms.CSPC*r6)*r2;
	      COMFAC2=ELJO*DSZR;
	      Fx0=-(COMFAC1*SZR+COMFAC2)*rX0;
	      Fy0=-(COMFAC1*SZR+COMFAC2)*rY0;
	      Fz0=-(COMFAC1*SZR+COMFAC2)*rZ0;

	      //energy plus LJ
	      E0=E0+ELJ;

	      //update particle forces Oxygen
	      Fx[i][0]+=Fx0;
	      Fy[i][0]+=Fy0;
	      Fz[i][0]+=Fz0;
	   
	      //update center of mass
	      
	      xCMnow=XCOM[i]-XCOM[j];
	      yCMnow=YCOM[i]-YCOM[j];
	      zCMnow=ZCOM[i]-ZCOM[j];
	      xCMnow+=parms.BOXLENGTH*nearestInt(xCMnow,parms.BOXLENGTH);
	      yCMnow+=parms.BOXLENGTH*nearestInt(yCMnow,parms.BOXLENGTH);
	      zCMnow+=parms.BOXLENGTH*nearestInt(zCMnow,parms.BOXLENGTH);
	      
	      VirijXX+=Fx0;
	      VirijYY+=Fy0;
	      VirijZZ+=Fz0;
	
	      //Put in this code for the extra hydrogen atoms

	      for(k=1;k<cons.NSITE;k++)
		{
		  for(n=1;n<cons.NSITE;n++)
		    {
		      RXnew=local_mols[i].x[k]-mols[j].x[n];
		      RYnew=local_mols[i].y[k]-mols[j].y[n];
		      RZnew=local_mols[i].z[k]-mols[j].z[n];
		      RXnew+=parms.BOXLENGTH*nearestInt(RXnew,parms.BOXLENGTH);
		      RYnew+=parms.BOXLENGTH*nearestInt(RYnew,parms.BOXLENGTH);
		      RZnew+=parms.BOXLENGTH*nearestInt(RZnew,parms.BOXLENGTH);
		      RijSQUARED=RXnew*RXnew+RYnew*RYnew+RZnew*RZnew;
		      Rij1=sqrt(RijSQUARED);
		      Q1=parms.Q[k]*parms.Q[n];
		      Term=erfc(cons.ALPHA*Rij1)/Rij1;
		      
		      //complementary error function
		      Term1=Q1*Term;
		      Tmpp=parms.CoeffPi*cons.ALPHA*exp(-cons.Alpha2*RijSQUARED)+Term;
		      Term1A=Q1*Tmpp/RijSQUARED;
		      Term1AA=Term1A*SZR;
		      FXnew=Term1AA*RXnew;
		      FYnew=Term1AA*RYnew;
		      FZnew=Term1AA*RZnew;
		      SumQQ+=Term1*SZR;
		      Fx[i][k]+=FXnew;
		      Fy[i][k]+=FYnew;
		      Fz[i][k]+=FZnew;
		      
		      VirijXX+=FXnew;
		      VirijYY+=FYnew;
		      VirijZZ+=FZnew;

		    }
		}
	      VirijXX*=xCMnow;
	      VirijYY*=yCMnow;
	      VirijZZ*=zCMnow;
	      	            
	    }// End of else if


	  //Energy term
	  EC+=SumQQ;
	  Vir_LrXX+=VirijXX;
	  Vir_LrYY+=VirijYY;
	  Vir_LrZZ+=VirijZZ;

	     }
	}//end of loop over second molecule in pair
    }// end of loop over pairs
  /*cout <<"rwald\n ";
   cout <<"energy of 0,0: "<<E0_0 <<endl;
  cout <<"fx0: " <<Fx[0][1]<<endl;
  cout <<"fy0: " <<Fy[0][1]<<endl;
  cout <<"fz0: " <<Fz[0][1]<<endl;
  */
  /*Return RealVir_LJ to main*/
  RealVir_LJ=Vir_LrXX+Vir_LrYY+Vir_LrZZ;

}// end of function rwald


void Ensemble::KWALD(Constants cons, Parameters parms, WaterSend *mols, WaterSend *local_mols)
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
  VK=0.0;
  
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
  double Thetx[11][512][3];
  double Thety[21][512][3];
  double Thetz[21][512][3];
  double Thetrk[512][3];

  
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
  int myK=11;
  double tSam;



  /*initialize fourier terms for each particle for k vectors of [0,0,0], [1,1,1], and [n,-1,-1]
   Indexing is strange in ky and kz, thetz[0][][]=kz=-10, thetz[20][][]=kz=10*/
  for(i=0;i<cons.N;i++)
    {
      //H1 index 10 is k=0
      Thetx[0][i][0]=0.0;
      Thety[10][i][0]=0.0;
      Thetz[10][i][0]=0.0;
      
      argx=FACT2*mols[i].x[1];
      argy=FACT2*mols[i].y[1];
      argz=FACT2*mols[i].z[1];
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
      
      argx=FACT2*mols[i].x[2];
      argy=FACT2*mols[i].y[2];
      argz=FACT2*mols[i].z[2];
      Thetx[1][i][1]=argx;
      Thety[11][i][1]=argy;
      Thetz[11][i][1]=argz;
      
      Thety[9][i][1]=-argy;
      Thetz[9][i][1]=-argz;

      //Msite      
      Thetx[0][i][2]=0.0;
      Thety[10][i][2]=0.0;
      Thetz[10][i][2]=0.0;
      
      argx=FACT2*mols[i].x[3];
      argy=FACT2*mols[i].y[3];
      argz=FACT2*mols[i].z[3];
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
      for(i=0;i<cons.N;i++)
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
		  for(i=0;i<cons.N;i++)
		    {
		      Thetrk[i][0]=Thetx[kx][i][0]+Thety[yi][i][0]+Thetz[zi][i][0];
		      SumKRe[0]+=cos(Thetrk[i][0]);
		      SumKIm[0]+=sin(Thetrk[i][0]);
		      
		      Thetrk[i][1]=Thetx[kx][i][1]+Thety[yi][i][1]+Thetz[zi][i][1];
		      SumKRe[1]+=cos(Thetrk[i][1]);
		      SumKIm[1]+=sin(Thetrk[i][1]);
		      
		      Thetrk[i][2]=Thetx[kx][i][2]+Thety[yi][i][2]+Thetz[zi][i][2];
		      SumKRe[2]+=cos(Thetrk[i][2]);
		      SumKIm[2]+=sin(Thetrk[i][2]);
		      
		    }
		  SumKRe[0]*=parms.Q[1];
		  SumKIm[0]*=parms.Q[1];
		  SumKRe[1]*=parms.Q[2];
		  SumKIm[1]*=parms.Q[2];
		  SumKRe[2]*=parms.Q[3];
		  SumKIm[2]*=parms.Q[3];
		  /*   StartSend************************************************************************************************************************************************************************
Need to do a global send here to turn the SKRe below over all K vectors, so
end K loop here and restart another afterwards*/

		  Veck=FACTVK*exp(-cons.SIGMA*RKsq)/RKsq;
		  
		  SKRe=SumKRe[0]+SumKRe[1]+SumKRe[2];
		  SKIm=SumKIm[0]+SumKIm[1]+SumKIm[2];

		  //magni^2 of complex number SK=Sum_r(exp(irk))
		  mag2=SKRe*SKRe+SKIm*SKIm;
		  magnitude=sqrt(mag2);
		  argSK=atan(SKIm/SKRe);
		  /*atan is only -pi/2 to pi/2 need whole range*/
		  if(SKRe<0)
		    argSK+=cons.Pi;
		  /*Potential term in K-space to electrostatic potential*/
		  tSam=0.5*FCT*Veck*mag2;

		  VK+=tSam;

		  VIRKxx+=tSam*(1.0-2.0*Xk*Xk*InvRKsq-Xk*Xk*A22I);
		  VIRKyy+=tSam*(1.0-2.0*Yk*Yk*InvRKsq-Yk*Yk*A22I);
		  VIRKzz+=tSam*(1.0-2.0*Zk*Zk*InvRKsq-Zk*Zk*A22I);
		  for(i=0;i<cons.N;i++)
		    {
		      for(j=0;j<3;j++)
			{
			  RFact=-FCT*parms.Q[j+1]*Veck;
			  CFact=magnitude*sin(argSK-Thetrk[i][j]);
			  /*Force is derivative of Vk*/
			  Fxik=RFact*Xk*CFact;
			  Fyik=RFact*Yk*CFact;
			  Fzik=RFact*Zk*CFact;
			  Fx[i][j+1]+=Fxik;
			  Fy[i][j+1]+=Fyik;
			  Fz[i][j+1]+=Fzik;
			  Xp=mols[i].x[j+1]-XCOM[i];
			  Yp=mols[i].y[j+1]-YCOM[i];
			  Zp=mols[i].z[j+1]-ZCOM[i];
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
  KVir_LJ=VIRxx+VIRyy+VIRzz;


  /* Last two pieces of electorstatic potential Vex is electorstatic coulomb sum from charges interacting within the same molecule.  Vs is spurious self interaction term of continuos gaussian with the point charge at the
center of it.*/

  Vex=0.0;
  Vs=0.0;

  for(j=1;j<cons.NSITE;j++)
    {
      for(k=1;k<cons.NSITE;k++)
	{
	  for(i=0;i<cons.N;i++)
	    {
	      if(j!=k)
		{
		  Rx1=mols[i].x[j]-mols[i].x[k];
		  Ry1=mols[i].y[j]-mols[i].y[k];
		  Rz1=mols[i].z[j]-mols[i].z[k];
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
  cout <<"VS: "<<Vs<<" Vex: "<<Vex<<" VK: "<<VK<<endl;
  VK=VK-Vs-Vex;
  /*  cout <<"in kwald\n";
  cout <<"fx0: " <<Fx[0][1]<<endl;
  cout <<"fy0: " <<Fy[0][1]<<endl;
  cout <<"fz0: " <<Fz[0][1]<<endl;
  */
			     

}

void Ensemble::Primder(Constants cons, Parameters parms, WaterSend *local_mols )
{
  /*transfer force from massless fictitious site onto real sites
  so update the forces on the real molecules 
  */

  //local vars
  double X10, Y10, Z10;//distances to oxygen atoms [0]
  double X20, Y20, Z20;
  double X30, Y30, Z30;

  double Fx0, Fy0, Fz0;
  double FxH, FyH, FzH;
  double ComFx, ComFy, ComFz;
  double Gama;
  double Xb, Yb, Zb;
  double Rb;
  double Rsq30, Fr30, Ft;
  double PRM=.125;
  
  int i;
  for(i=0;i<cons.NperP;i++)
    {
      //Distances Hydrogens [2]and [1] to Oxygen [0]
      X10=local_mols[i].x[1]-local_mols[i].x[0];
      Y10=local_mols[i].y[1]-local_mols[i].y[0];
      Z10=local_mols[i].z[1]-local_mols[i].z[0];
      X20=local_mols[i].x[2]-local_mols[i].x[0];
      Y20=local_mols[i].y[2]-local_mols[i].y[0];
      Z20=local_mols[i].z[2]-local_mols[i].z[0];
      
      Xb=0.5*(X10+X20);
      Yb=0.5*(Y10+Y20);
      Zb=0.5*(Z10+Z20);
      
      Rb=sqrt(Xb*Xb+Yb*Yb+Zb*Zb);
      Gama=PRM/Rb;
      //distances fictional [3] to Oxygen [0]
      X30=local_mols[i].x[3]-local_mols[i].x[0];
      Y30=local_mols[i].y[3]-local_mols[i].y[0];
      Z30=local_mols[i].z[3]-local_mols[i].z[0];
      
      Rsq30=X30*X30+Y30*Y30+Z30*Z30;
      Fr30=X30*Fx[i][3]+Y30*Fy[i][3]+Z30*Fz[i][3];
      Ft=Fr30/Rsq30;

      //center of mass of forces
      ComFx=Gama*(Fx[i][3]-Ft*X30);
      ComFy=Gama*(Fy[i][3]-Ft*Y30);
      ComFz=Gama*(Fz[i][3]-Ft*Z30);
      
      Fx0=Fx[i][3]-ComFx;
      Fy0=Fy[i][3]-ComFy;
      Fz0=Fz[i][3]-ComFz;
      
      FxH=0.5*ComFx;
      FyH=0.5*ComFy;
      FzH=0.5*ComFz;
      
      //Now update forces of real atoms
      Fx[i][0]+=Fx0;
      Fy[i][0]+=Fy0;
      Fz[i][0]+=Fz0;
      Fx[i][1]+=FxH;
      Fy[i][1]+=FyH;
      Fz[i][1]+=FzH;
      Fx[i][2]+=FxH;
      Fy[i][2]+=FyH;
      Fz[i][2]+=FzH;
      



    }//end loop over NperP

}//End function primder


void Ensemble::Maxwell(Constants cons, Parameters parms)
{
  int i;
  double Sigma0, Sigma1, Sigma2;

  double Rn;
  double sumx, sumy, sumz, amsum;
  Sigma0=sqrt(cons.BOLTZK*cons.Temp/parms.AMass[0]);
  Sigma1=sqrt(cons.BOLTZK*cons.Temp/parms.AMass[1]);
  Sigma2=sqrt(cons.BOLTZK*cons.Temp/parms.AMass[2]);


  for(i=0;i<cons.NperP;i++)
    {

      Rn=GaussV();
      vx[i][0]=Sigma0*Rn;
      Rn=GaussV();
      vy[i][0]=Sigma0*Rn;
      Rn=GaussV();
      vz[i][0]=Sigma0*Rn;


      Rn=GaussV();
      vx[i][1]=Sigma1*Rn;
      Rn=GaussV();
      vy[i][1]=Sigma1*Rn;
      Rn=GaussV();
      vz[i][1]=Sigma1*Rn;


      Rn=GaussV();
      vx[i][2]=Sigma2*Rn;
      Rn=GaussV();
      vy[i][2]=Sigma2*Rn;
      Rn=GaussV();
      vz[i][2]=Sigma2*Rn;
    }
  /*Make the Center of Mass of bulk water a Fixed point. Sum_x(everyAtom_x*hisx_velocity)
   amsum will be 18.01528*N  */
  
  int j;
  sumx=sumy=sumz=amsum=0.0; 
  for(i=0;i<cons.N;i++)
    {
      for(j=0;j<cons.NUMPRIM;j++)
	{
	  
	  sumx+=parms.AMass[j]*vx[i][j];
	  sumy+=parms.AMass[j]*vy[i][j];
	  sumz+=parms.AMass[j]*vz[i][j];
	  amsum+=parms.AMass[j];
	}
    }
  sumx=sumx/amsum;
  sumy=sumy/amsum;
  sumz=sumz/amsum;
  for(i=0;i<cons.N;i++)
    {
      vx[i][0]-=sumx;
      vy[i][0]-=sumy;
      vz[i][0]-=sumz;
      
      vx[i][1]-=sumx;
      vy[i][1]-=sumy;
      vz[i][1]-=sumz;

      vx[i][2]-=sumx;
      vy[i][2]-=sumy;
      vz[i][2]-=sumz;
    }

}//end of Maxwell


double Ensemble::KineticA(Constants cons)
{
  /*Compute Kinetic Energy Sum_i m*v[i]^2 over all atoms (O, H, H) over all N molecules*/
  /*only zeros (no Os)*/
  double KinE=0.0;
  double Kin0=0.0;
  double KinH1=0.0;
  double KinH2=0.0;
  double am0, am1, am2;
  int i;
  am0=cons.AMass[0];
  am1=cons.AMass[1];
  am2=cons.AMass[2];
  for(i=0;i<cons.NperP;i++)
    {
	  Kin0+=vx[i][0]*vx[i][0]+vy[i][0]*vy[i][0]+vz[i][0]*vz[i][0];
	  KinH1+=vx[i][1]*vx[i][1]+vy[i][1]*vy[i][1]+vz[i][1]*vz[i][1];
	  KinH2+=vx[i][2]*vx[i][2]+vy[i][2]*vy[i][2]+vz[i][2]*vz[i][2];
    }

  KinE=am0*Kin0+am1*KinH1+am2*KinH2;
  return KinE;
}

void Ensemble::ScaleVelocity(Ensemble &NNN, Constants cons, Parameters parms)
{
  /*Scale the velocities to the set temperature 
   Sets the Energy Ensemble Variable EKIN*/

  double Akin;
  double InstTemp;
  double Tf;
  int i;

  Akin=NNN.KineticA(cons);
  InstTemp=Akin/parms.ConstM;
  cout<<"Temp pre scale: "<<InstTemp<<endl;
  Tf=sqrt(cons.Temp/InstTemp);
  //scale
  for(i=0;i<cons.NperP;i++)
    {
      vx[i][0]*=Tf;
      vy[i][0]*=Tf;
      vz[i][0]*=Tf;

      vx[i][1]*=Tf;
      vy[i][1]*=Tf;
      vz[i][1]*=Tf;

      vx[i][2]*=Tf;
      vy[i][2]*=Tf;
      vz[i][2]*=Tf;

    }

  Akin=NNN.KineticA(cons);
  InstTemp=Akin/parms.ConstM;
  //Kinetic Energy
  EKIN=0.5*Akin/418.4;
  cout << "Kinetic Energy after scaling vels: "<<EKIN<<endl;
  cout <<"init Temp post scale: "<<InstTemp<<endl;


}


void Ensemble::CenterMassWater(Constants cons, WaterSend *mols)
{
  double xsum, ysum, zsum, amsum;
  int i;
  double am0=cons.AMass[0]; /* Oxygen*/
  double am1=cons.AMass[1]; /*hydrogen*/
  double am2=cons.AMass[2]; /*hydrogen*/
  amsum=am0+am1+am2;

  for(i=0;i<cons.N;i++)
    {
      xsum=ysum=zsum=0.0;
      
      xsum=mols[i].x[0]*am0+mols[i].x[1]*am1+mols[i].x[2]*am2;
      ysum=mols[i].y[0]*am0+mols[i].y[1]*am1+mols[i].y[2]*am2;
      zsum=mols[i].z[0]*am0+mols[i].z[1]*am1+mols[i].z[2]*am2;
      
      //center of mass = sum_3 (x(i)*m(i)) / sum_3 (m(i))

      XCOM[i]=xsum/amsum;
      YCOM[i]=ysum/amsum;
      ZCOM[i]=zsum/amsum;
    }
}


double Ensemble::CenterMassKinEnergy(Constants cons)
{
  double vxsum, vysum, vzsum, amsum;
  int i;
  double am0=cons.AMass[0]; /* Oxygen*/
  double am1=cons.AMass[1]; /*hydrogen*/
  double am2=cons.AMass[2]; /*hydrogen*/
  double vxCOM, vyCOM, vzCOM;
  double KEcom;
  amsum=am0+am1+am2;
  vxCOM=vyCOM=vzCOM=0.0;
  KEcom=0.0;
  for(i=0;i<cons.NperP;i++)
    {
      vxsum=vysum=vzsum=0.0;
      
      vxsum=vx[i][0]*am0+vx[i][1]*am1+vx[i][2]*am2;
      vysum=vy[i][0]*am0+vy[i][1]*am1+vy[i][2]*am2;
      vzsum=vz[i][0]*am0+vz[i][1]*am1+vz[i][2]*am2;
      
      //center of mass = sum_3 (x(i)*m(i)) / sum_3 (m(i))

      vxCOM=vxsum/amsum;
      vyCOM=vysum/amsum;
      vzCOM=vzsum/amsum;
      KEcom+=amsum*(vxCOM*vxCOM+vyCOM*vyCOM+vzCOM*vzCOM);
    }
  KEcom=KEcom/418.4;
  return KEcom;



}



void Ensemble::RattlePos(Constants cons, Parameters parms, WaterSend *local_mols)
{
  

  //updated velocities in integrateMotion call
  int i, it, j;
  int NumIt;
  double sx0, sy0, sz0;
  double sx1, sy1, sz1;
  double sx2, sy2, sz2;
  double r01s01, r01s02, r01s12, r02s01, r02s02, r02s12, r12s01, r12s02, r12s12;
  double xb, yb, zb, rb;
  double s01sq, s02sq, s12sq;
  double sig01, sig02, sig12;
  double cm11, cm21, cm31;
  double x10, y10, z10, x20, y20, z20;
  double gr01, gr02, gr12;
  double xpbc, ypbc, zpbc;
  
  double Invdt2=1.0/(2.0*cons.dt*cons.dt);
  
  //3 for x, y, z
  double * s01=new double[3];
  double * s02=new double[3];
  double * s12=new double[3];
  double * r01=new double[3];
  double * r02=new double[3];
  double * r12=new double[3];
  
  double ** Shake=new double *[3];
  double ** SInv=new double *[3];
  for(i=0;i<3;i++)
    {
      Shake[i]=new double[3];
      SInv[i]=new double[3];
    }
  
  

  //do position update step but don't move particles yet

  for(i=0;i<cons.NperP;i++)
    {
      
      sx0=local_mols[i].x[0]+cons.dt*vx[i][0];
      sy0=local_mols[i].y[0]+cons.dt*vy[i][0];
      sz0=local_mols[i].z[0]+cons.dt*vz[i][0];
      sx1=local_mols[i].x[1]+cons.dt*vx[i][1];
      sy1=local_mols[i].y[1]+cons.dt*vy[i][1];
      sz1=local_mols[i].z[1]+cons.dt*vz[i][1];
      sx2=local_mols[i].x[2]+cons.dt*vx[i][2];
      sy2=local_mols[i].y[2]+cons.dt*vy[i][2];
      sz2=local_mols[i].z[2]+cons.dt*vz[i][2];

      /*Now calculate relative vectors*/
      //O to H1 distance
      r01[0]=local_mols[i].x[0]-local_mols[i].x[1];
      r01[1]=local_mols[i].y[0]-local_mols[i].y[1];
      r01[2]=local_mols[i].z[0]-local_mols[i].z[1];
      r01[0]+=parms.BOXLENGTH*nearestInt(r01[0], parms.BOXLENGTH);
      r01[1]+=parms.BOXLENGTH*nearestInt(r01[1], parms.BOXLENGTH);
      r01[2]+=parms.BOXLENGTH*nearestInt(r01[2], parms.BOXLENGTH);
      
      /*distance O to H2*/
      r02[0]=local_mols[i].x[0]-local_mols[i].x[2];
      r02[1]=local_mols[i].y[0]-local_mols[i].y[2];
      r02[2]=local_mols[i].z[0]-local_mols[i].z[2];
      r02[0]+=parms.BOXLENGTH*nearestInt(r02[0], parms.BOXLENGTH);
      r02[1]+=parms.BOXLENGTH*nearestInt(r02[1], parms.BOXLENGTH);
      r02[2]+=parms.BOXLENGTH*nearestInt(r02[2], parms.BOXLENGTH);

      /*ditsnce H to H*/
      r12[0]=local_mols[i].x[1]-local_mols[i].x[2];
      r12[1]=local_mols[i].y[1]-local_mols[i].y[2];
      r12[2]=local_mols[i].z[1]-local_mols[i].z[2];
      r12[0]+=parms.BOXLENGTH*nearestInt(r12[0], parms.BOXLENGTH);
      r12[1]+=parms.BOXLENGTH*nearestInt(r12[1], parms.BOXLENGTH);
      r12[2]+=parms.BOXLENGTH*nearestInt(r12[2], parms.BOXLENGTH);

      /*Iterative loop begins*/
      NumIt=0;
      it=0;
      while(it<parms.iterRat)
	{
	  /*are the pbc ajdustments necessary,  only if o and h move a half a boxlength apart*/

	  s01[0]=sx0-sx1;
	  s01[1]=sy0-sy1;
	  s01[2]=sz0-sz1;
	  s01[0]+=parms.BOXLENGTH*nearestInt(s01[0], parms.BOXLENGTH);
	  s01[1]+=parms.BOXLENGTH*nearestInt(s01[1], parms.BOXLENGTH);
	  s01[2]+=parms.BOXLENGTH*nearestInt(s01[2], parms.BOXLENGTH);
	  
	  s02[0]=sx0-sx2;
	  s02[1]=sy0-sy2;
	  s02[2]=sz0-sz2;
	  s02[0]+=parms.BOXLENGTH*nearestInt(s02[0], parms.BOXLENGTH);
	  s02[1]+=parms.BOXLENGTH*nearestInt(s02[1], parms.BOXLENGTH);
	  s02[2]+=parms.BOXLENGTH*nearestInt(s02[2], parms.BOXLENGTH);
	  
	  s12[0]=sx1-sx2;
	  s12[1]=sy1-sy2;
	  s12[2]=sz1-sz2;
	  s12[0]+=parms.BOXLENGTH*nearestInt(s12[0], parms.BOXLENGTH);
	  s12[1]+=parms.BOXLENGTH*nearestInt(s12[1], parms.BOXLENGTH);
	  s12[2]+=parms.BOXLENGTH*nearestInt(s12[2], parms.BOXLENGTH);

	  s01sq=dotProd(s01, s01);
	  s02sq=dotProd(s02, s02);
	  s12sq=dotProd(s12, s12);
	  
	  //should be practically zero
	  sig01=abs(s01sq-parms.d01sq)/(2.0*parms.d01sq);
	  sig02=abs(s02sq-parms.d02sq)/(2.0*parms.d02sq);
	  sig12=abs(s12sq-parms.d12sq)/(2.0*parms.d12sq);
	  
	  if(sig01<=parms.tol&&sig02<=parms.tol&&sig12<=parms.tol)
	    it=parms.iterRat+1;
	  /*conitue out of while loop, constraints are met*/
	  
	  //calculate scalar products r is original vector, s is moved vec
	  r01s01=dotProd(r01, s01);
	  r02s01=dotProd(r02, s01);
	  r12s01=dotProd(r12, s01);
	  r01s02=dotProd(r01, s02);
	  r02s02=dotProd(r02, s02);
	  r12s02=dotProd(r12, s02);
	  r01s12=dotProd(r01, s12);
	  r02s12=dotProd(r02, s12);
	  r12s12=dotProd(r12, s12);
	  
	  //calculate matrix elements
	  
	  Shake[0][0]=r01s01*(parms.InvMass[0]+parms.InvMass[1]);
	  Shake[0][1]=r02s01*parms.InvMass[0];
	  Shake[0][2]=-r12s01*parms.InvMass[1];
	  Shake[1][0]=r01s02*parms.InvMass[0];
	  Shake[1][1]=r02s02*(parms.InvMass[0]+parms.InvMass[2]);
	  Shake[1][2]=r12s02*parms.InvMass[2];
	  Shake[2][0]=-r01s12*parms.InvMass[1];
	  Shake[2][1]=r02s12*parms.InvMass[2];
	  Shake[2][2]=r12s12*(parms.InvMass[1]+parms.InvMass[2]);
	  
	  cm11=Invdt2*(s01sq-parms.d01sq);
	  cm21=Invdt2*(s02sq-parms.d02sq);
	  cm31=Invdt2*(s12sq-parms.d12sq);
	  
	  //invert the shake matrix
	  MatrixInv(Shake, SInv);

	  //obtain solutions of linearized equation
	  gr01=SInv[0][0]*cm11+SInv[0][1]*cm21+SInv[0][2]*cm31;
	  gr02=SInv[1][0]*cm11+SInv[1][1]*cm21+SInv[1][2]*cm31;
	  gr12=SInv[2][0]*cm11+SInv[2][1]*cm21+SInv[2][2]*cm31;
	  
	  //obtain velocity adjusted
	  vx[i][0]-=cons.dt*parms.InvMass[0]*(gr01*r01[0]+gr02*r02[0]);
	  vy[i][0]-=cons.dt*parms.InvMass[0]*(gr01*r01[1]+gr02*r02[1]);
	  vz[i][0]-=cons.dt*parms.InvMass[0]*(gr01*r01[2]+gr02*r02[2]);
	  vx[i][1]-=cons.dt*parms.InvMass[1]*(-gr01*r01[0]+gr12*r12[0]);
	  vy[i][1]-=cons.dt*parms.InvMass[1]*(-gr01*r01[1]+gr12*r12[1]);
	  vz[i][1]-=cons.dt*parms.InvMass[1]*(-gr01*r01[2]+gr12*r12[2]);
	  vx[i][2]-=cons.dt*parms.InvMass[2]*(-gr02*r02[0]-gr12*r12[0]);
	  vy[i][2]-=cons.dt*parms.InvMass[2]*(-gr02*r02[1]-gr12*r12[1]);
	  vz[i][2]-=cons.dt*parms.InvMass[2]*(-gr02*r02[2]-gr12*r12[2]);
 
	  //new coords
	  sx0=local_mols[i].x[0]+cons.dt*vx[i][0];
	  sy0=local_mols[i].y[0]+cons.dt*vy[i][0];
	  sz0=local_mols[i].z[0]+cons.dt*vz[i][0];
	  sx1=local_mols[i].x[1]+cons.dt*vx[i][1];
	  sy1=local_mols[i].y[1]+cons.dt*vy[i][1];
	  sz1=local_mols[i].z[1]+cons.dt*vz[i][1];
	  sx2=local_mols[i].x[2]+cons.dt*vx[i][2];
	  sy2=local_mols[i].y[2]+cons.dt*vy[i][2];
	  sz2=local_mols[i].z[2]+cons.dt*vz[i][2];
	  
	  NumIt++;

	  it++;
	}//end while loop
      if(it==parms.iterRat)
	{
	  cerr<<"position rattle didnot converge"<<endl;
	  exit(1);
	}

      /*Now Update the position of the actual atoms, move them*/
      for(j=0;j<cons.NUMPRIM;j++)
	{
	  local_mols[i].x[j]+=cons.dt*vx[i][j];
	  local_mols[i].y[j]+=cons.dt*vy[i][j];
	  local_mols[i].z[j]+=cons.dt*vz[i][j];

	}
      //do periodic boundary conditions for Oxygen, then adjust all atoms
      xpbc=parms.BOXLENGTH*nearestInt(local_mols[i].x[0], parms.BOXLENGTH);
      ypbc=parms.BOXLENGTH*nearestInt(local_mols[i].y[0], parms.BOXLENGTH);
      zpbc=parms.BOXLENGTH*nearestInt(local_mols[i].z[0], parms.BOXLENGTH);
      
      //move O, H, H
      local_mols[i].x[0]+=xpbc;
      local_mols[i].y[0]+=ypbc;
      local_mols[i].z[0]+=zpbc;
      local_mols[i].x[1]+=xpbc;
      local_mols[i].y[1]+=ypbc;
      local_mols[i].z[1]+=zpbc;
      local_mols[i].x[2]+=xpbc;
      local_mols[i].y[2]+=ypbc;
      local_mols[i].z[2]+=zpbc;

      /*recalc the fictitious point*/
      x10=local_mols[i].x[1]-local_mols[i].x[0];
      y10=local_mols[i].y[1]-local_mols[i].y[0];
      z10=local_mols[i].z[1]-local_mols[i].z[0];
      x20=local_mols[i].x[2]-local_mols[i].x[0];
      y20=local_mols[i].y[2]-local_mols[i].y[0];
      z20=local_mols[i].z[2]-local_mols[i].z[0];
      xb=0.5*(x10+x20);
      yb=0.5*(y10+y20);
      zb=0.5*(z10+z20);
      rb=sqrt(xb*xb+yb*yb+zb*zb);
      local_mols[i].x[3]=local_mols[i].x[0]+parms.RM*xb/rb;
      local_mols[i].y[3]=local_mols[i].y[0]+parms.RM*yb/rb;
      local_mols[i].z[3]=local_mols[i].z[0]+parms.RM*zb/rb;

    }//end loop over N particles

}//end function RattlePos 

/*NON CLASS FUNCTION, NOT PRECEDED BY ENSEMBLE::*/
void MatrixInv(double ** Shake, double **Sinv)
{
  double stol=1.0E-7;
  double determ;

  //get signed cofactors, transpose and put in skinv
  
  Sinv[0][0]=Shake[1][1]*Shake[2][2]-Shake[1][2]*Shake[2][1];
  Sinv[0][1]=Shake[0][2]*Shake[2][1]-Shake[0][1]*Shake[2][2];
  Sinv[0][2]=Shake[0][1]*Shake[1][2]-Shake[0][2]*Shake[1][1];
  Sinv[1][0]=Shake[1][2]*Shake[2][0]-Shake[1][0]*Shake[2][2];
  Sinv[1][1]=Shake[0][0]*Shake[2][2]-Shake[0][2]*Shake[2][0];
  Sinv[1][2]=Shake[0][2]*Shake[1][0]-Shake[0][0]*Shake[1][2];
  Sinv[2][0]=Shake[1][0]*Shake[2][1]-Shake[1][1]*Shake[2][0];
  Sinv[2][1]=Shake[0][1]*Shake[2][0]-Shake[0][0]*Shake[2][1];
  Sinv[2][2]=Shake[0][0]*Shake[1][1]-Shake[0][1]*Shake[1][0];

  //get determinant and make sure not zero

  determ=Shake[0][0]*Sinv[0][0]+Shake[0][1]*Sinv[1][0]+Shake[0][2]*Sinv[2][0];
  
  if(abs(determ)<stol)
    {
      cerr <<"zero determinant in ShakeInv "<<endl;
      exit(1);
    }
  Sinv[0][0]/=determ;
  Sinv[0][1]/=determ;
  Sinv[0][2]/=determ;
  Sinv[1][0]/=determ;
  Sinv[1][1]/=determ;
  Sinv[1][2]/=determ;
  Sinv[2][0]/=determ;
  Sinv[2][1]/=determ;
  Sinv[2][2]/=determ;


}

void Ensemble::RattleVel(Constants cons, Parameters parms, WaterSend *local_mols)
{
  int i;
  double * r01=new double[3];
  double * r02=new double[3];
  double * r12=new double[3];
  double * v01=new double[3];
  double * v02=new double[3];
  double * v12=new double[3];

  double gv01, gv02, gv12;
  double qr11, qr21, qr31;
  double r01r01, r01r02, r01r12, r02r02, r02r12, r12r12;

  double **Rattle=new double*[3];
  double **RInv=new double*[3];
  for(i=0;i<3;i++)
    {
      Rattle[i]=new double[3];
      RInv[i]=new double[3];
    }

  for(i=0;i<cons.NperP;i++)
    {
      r01[0]=local_mols[i].x[0]-local_mols[i].x[1];
      r01[1]=local_mols[i].y[0]-local_mols[i].y[1];
      r01[2]=local_mols[i].z[0]-local_mols[i].z[1];
      r01[0]+=parms.BOXLENGTH*nearestInt(r01[0],parms.BOXLENGTH);
      r01[1]+=parms.BOXLENGTH*nearestInt(r01[1],parms.BOXLENGTH);
      r01[2]+=parms.BOXLENGTH*nearestInt(r01[2],parms.BOXLENGTH);
      
      r02[0]=local_mols[i].x[0]-local_mols[i].x[2];
      r02[1]=local_mols[i].y[0]-local_mols[i].y[2];
      r02[2]=local_mols[i].z[0]-local_mols[i].z[2];
      r02[0]+=parms.BOXLENGTH*nearestInt(r02[0],parms.BOXLENGTH);
      r02[1]+=parms.BOXLENGTH*nearestInt(r02[1],parms.BOXLENGTH);
      r02[2]+=parms.BOXLENGTH*nearestInt(r02[2],parms.BOXLENGTH);
      
      r12[0]=local_mols[i].x[1]-local_mols[i].x[2];
      r12[1]=local_mols[i].y[1]-local_mols[i].y[2];
      r12[2]=local_mols[i].z[1]-local_mols[i].z[2];
      r12[0]+=parms.BOXLENGTH*nearestInt(r12[0],parms.BOXLENGTH);
      r12[1]+=parms.BOXLENGTH*nearestInt(r12[1],parms.BOXLENGTH);
      r12[2]+=parms.BOXLENGTH*nearestInt(r12[2],parms.BOXLENGTH);

      v01[0]=vx[i][0]-vx[i][1];
      v01[1]=vy[i][0]-vy[i][1];
      v01[2]=vz[i][0]-vz[i][1];
      v02[0]=vx[i][0]-vx[i][2];
      v02[1]=vy[i][0]-vy[i][2];
      v02[2]=vz[i][0]-vz[i][2];
      v12[0]=vx[i][1]-vx[i][2];
      v12[1]=vy[i][1]-vy[i][2];
      v12[2]=vz[i][1]-vz[i][2];
      qr11=-dotProd(v01, r01);
      qr21=-dotProd(v02, r02);
      qr31=-dotProd(v12, r12);
      
      //calculate matrix elements for rattle
      r01r01=dotProd(r01, r01);
      r12r12=dotProd(r12, r12);
      r02r02=dotProd(r02, r02);
      r01r02=dotProd(r01, r02);
      r01r12=dotProd(r01, r12);
      r02r12=dotProd(r02, r12);
      
      Rattle[0][0]=r01r01*(parms.InvMass[0]+parms.InvMass[1]);
      Rattle[0][1]=r01r02*parms.InvMass[0];
      Rattle[0][2]=-r01r12*parms.InvMass[1];
      Rattle[1][0]=Rattle[0][1];
      Rattle[1][1]=r02r02*(parms.InvMass[0]+parms.InvMass[2]);
      Rattle[1][2]=r02r12*parms.InvMass[2];
      Rattle[2][0]=Rattle[0][2];
      Rattle[2][1]=Rattle[1][2];
      Rattle[2][2]=r12r12*(parms.InvMass[2]+parms.InvMass[1]);
      /*
      cout << Rattle[0][0]<<' '<<Rattle[0][1]<<' '<<Rattle[0][2]<<endl;
      cout << Rattle[1][0]<<' '<<Rattle[1][1]<<' '<<Rattle[1][2]<<endl;
      cout << Rattle[2][0]<<' '<<Rattle[2][1]<<' '<<Rattle[2][2]<<endl;
      */
      //invert the matrix
      MatrixInv(Rattle, RInv);
      /*      cout << RInv[0][0]<<' '<<RInv[0][1]<<' '<<RInv[0][2]<<endl;
      cout << RInv[1][0]<<' '<<RInv[1][1]<<' '<<RInv[1][2]<<endl;
      cout << RInv[2][0]<<' '<<RInv[2][1]<<' '<<RInv[2][2]<<endl;
      */
      //get velocity lagrange multipliers
      gv01=RInv[0][0]*qr11+RInv[0][1]*qr21+RInv[0][2]*qr31;
      gv02=RInv[1][0]*qr11+RInv[1][1]*qr21+RInv[1][2]*qr31;
      gv12=RInv[2][0]*qr11+RInv[2][1]*qr21+RInv[2][2]*qr31;

      //Evaluate the constraint velocity
      vx[i][0]+=parms.InvMass[0]*(gv01*r01[0]+gv02*r02[0]);
      vy[i][0]+=parms.InvMass[0]*(gv01*r01[1]+gv02*r02[1]);
      vz[i][0]+=parms.InvMass[0]*(gv01*r01[2]+gv02*r02[2]);
      vx[i][1]+=parms.InvMass[1]*(-gv01*r01[0]+gv12*r12[0]);
      vy[i][1]+=parms.InvMass[1]*(-gv01*r01[1]+gv12*r12[1]);
      vz[i][1]+=parms.InvMass[1]*(-gv01*r01[2]+gv12*r12[2]);
      vx[i][2]+=parms.InvMass[2]*(-gv02*r02[0]-gv12*r12[0]);
      vy[i][2]+=parms.InvMass[2]*(-gv02*r02[1]-gv12*r12[1]);
      vz[i][2]+=parms.InvMass[2]*(-gv02*r02[2]-gv12*r12[2]);
      


    }//end of loop over all molecules

}//end of RattleVel
void Ensemble::checkIdeal(char *filename, Constants cons, Parameters parms, WaterSend *mols)
{

  int i;
  double *r01=new double[3];
  double *r02=new double[3];
  double *r12=new double[3];

  double r01sq, r02sq, r12sq;
  double magr01, magr02, magr12;
  double DP;
  double cosThetaRad, ThetaRad, Thetad;
  ofstream outfile(filename);
  
  /*Now calculate relative vectors*/

  //O to H1 distance
  for(i=0;i<cons.N;i++)
    {
      r01[0]=mols[i].x[0]-mols[i].x[1];
      r01[1]=mols[i].y[0]-mols[i].y[1];
      r01[2]=mols[i].z[0]-mols[i].z[1];
      /*      r01[0]+=parms.BOXLENGTH*nearestInt(r01[0], parms.BOXLENGTH);
      r01[1]+=parms.BOXLENGTH*nearestInt(r01[1], parms.BOXLENGTH);
      r01[2]+=parms.BOXLENGTH*nearestInt(r01[2], parms.BOXLENGTH);
      */
      
      r01sq=dotProd(r01, r01);
      magr01=sqrt(r01sq);

      /*distance O to H2*/
      r02[0]=mols[i].x[0]-mols[i].x[2];
      r02[1]=mols[i].y[0]-mols[i].y[2];
      r02[2]=mols[i].z[0]-mols[i].z[2];
      /*      r02[0]+=parms.BOXLENGTH*nearestInt(r02[0], parms.BOXLENGTH);
      r02[1]+=parms.BOXLENGTH*nearestInt(r02[1], parms.BOXLENGTH);
      r02[2]+=parms.BOXLENGTH*nearestInt(r02[2], parms.BOXLENGTH);
      */
      r02sq=dotProd(r02, r02);
      magr02=sqrt(r02sq);

      /*ditsnce H to H*/
      r12[0]=mols[i].x[1]-mols[i].x[2];
      r12[1]=mols[i].y[1]-mols[i].y[2];
      r12[2]=mols[i].z[1]-mols[i].z[2];
      /*      r12[0]+=parms.BOXLENGTH*nearestInt(r12[0], parms.BOXLENGTH);
      r12[1]+=parms.BOXLENGTH*nearestInt(r12[1], parms.BOXLENGTH);
      r12[2]+=parms.BOXLENGTH*nearestInt(r12[2], parms.BOXLENGTH);
      */
      r12sq=dotProd(r12, r12);
      magr12=sqrt(r12sq);

      /*Angle between vectors OH1 and OH2*/
      
      DP=dotProd(r01, r02);
      cosThetaRad=DP/magr01/magr02;
      ThetaRad=acos(cosThetaRad);
      Thetad=ThetaRad*(cons.PiDEG/cons.Pi);
      
      outfile <<Thetad<<" "<<magr12<<" "<<magr01<<" "<<magr02<<endl;


    }

}

/*Function to remove any component of velocity along the bond, in essence
making the velocity normal to the bond in order to satisfy the constraints of 
the rigid molecule.*/

void Ensemble::VelocityNormal(Constants cons, Parameters parms, WaterSend *mols)
{
  int i;
  double * r01=new double[3];
  double * r02=new double[3];
  double * r12=new double[3];
  double * v01=new double[3];
  double * v02=new double[3];
  double * v12=new double[3];
  double m1, m2, m3, L1, L2, L3, c1, c2, c3;
 
  double r01r01, r01r02, r01r12, r02r02, r02r12, r12r12;

  double **Zmat=new double*[3];
  double **Zinv=new double*[3];
  for(i=0;i<3;i++)
    {
      Zmat[i]=new double[3];
      Zinv[i]=new double[3];
    }

  /*Subtract from velocities any components normal to the constraint surface 
allows constraint time derivatives to vanish*/
 
  

  for(i=0;i<cons.N;i++)
    {
      r01[0]=mols[i].x[0]-mols[i].x[1];
      r01[1]=mols[i].y[0]-mols[i].y[1];
      r01[2]=mols[i].z[0]-mols[i].z[1];
      r01[0]+=parms.BOXLENGTH*nearestInt(r01[0],parms.BOXLENGTH);
      r01[1]+=parms.BOXLENGTH*nearestInt(r01[1],parms.BOXLENGTH);
      r01[2]+=parms.BOXLENGTH*nearestInt(r01[2],parms.BOXLENGTH);
      
      r02[0]=mols[i].x[0]-mols[i].x[2];
      r02[1]=mols[i].y[0]-mols[i].y[2];
      r02[2]=mols[i].z[0]-mols[i].z[2];
      r02[0]+=parms.BOXLENGTH*nearestInt(r02[0],parms.BOXLENGTH);
      r02[1]+=parms.BOXLENGTH*nearestInt(r02[1],parms.BOXLENGTH);
      r02[2]+=parms.BOXLENGTH*nearestInt(r02[2],parms.BOXLENGTH);
      
      r12[0]=mols[i].x[1]-mols[i].x[2];
      r12[1]=mols[i].y[1]-mols[i].y[2];
      r12[2]=mols[i].z[1]-mols[i].z[2];
      r12[0]+=parms.BOXLENGTH*nearestInt(r12[0],parms.BOXLENGTH);
      r12[1]+=parms.BOXLENGTH*nearestInt(r12[1],parms.BOXLENGTH);
      r12[2]+=parms.BOXLENGTH*nearestInt(r12[2],parms.BOXLENGTH);

      v01[0]=vx[i][0]-vx[i][1];
      v01[1]=vy[i][0]-vy[i][1];
      v01[2]=vz[i][0]-vz[i][1];
      v02[0]=vx[i][0]-vx[i][2];
      v02[1]=vy[i][0]-vy[i][2];
      v02[2]=vz[i][0]-vz[i][2];
      v12[0]=vx[i][1]-vx[i][2];
      v12[1]=vy[i][1]-vy[i][2];
      v12[2]=vz[i][1]-vz[i][2];

      /*Calculate matrix elements*/
      
      r01r01=dotProd(r01, r01);
      r12r12=dotProd(r12, r12);
      r02r02=dotProd(r02, r02);
      r01r02=dotProd(r01, r02);
      r01r12=dotProd(r01, r12);
      r02r12=dotProd(r02, r12);

            
      Zmat[0][0]=r01r01*4.0*(parms.InvMass[0]+parms.InvMass[1]);
      Zmat[0][1]=r01r02*4.0*parms.InvMass[0];
      Zmat[0][2]=-r01r12*4.0*parms.InvMass[1];
      Zmat[1][0]=Zmat[0][1];
      Zmat[1][1]=r02r02*4.0*(parms.InvMass[0]+parms.InvMass[2]);
      Zmat[1][2]=r02r12*4.0*parms.InvMass[2];
      Zmat[2][0]=Zmat[0][2];
      Zmat[2][1]=Zmat[1][2];
      Zmat[2][2]=r12r12*4.0*(parms.InvMass[2]+parms.InvMass[1]);
      c1=2.0*dotProd(r01,v01);
      c2=2.0*dotProd(r02,v02);
      c3=2.0*dotProd(r12,v12);
      
      MatrixInv(Zmat, Zinv);

      /*Coefficients of the velocities*/
      L1=Zinv[0][0]*c1+Zinv[0][1]*c2+Zinv[0][2]*c3;
      L2=Zinv[1][0]*c1+Zinv[1][1]*c2+Zinv[1][2]*c3;
      L3=Zinv[2][0]*c1+Zinv[2][1]*c2+Zinv[2][2]*c3;

      /*Calculate constrained velocities*/
      m1=2.0*parms.InvMass[0];
      m2=2.0*parms.InvMass[1];
      m3=2.0*parms.InvMass[2];
      vx[i][0]-=m1*L1*r01[0]-m1*L2*r02[0];
      vy[i][0]-=m1*L1*r01[1]-m1*L2*r02[1];
      vz[i][0]-=m1*L1*r01[2]-m1*L2*r02[2];
      vx[i][1]+=m2*L1*r01[0]-m2*L3*r12[0];
      vy[i][1]+=m2*L1*r01[1]-m2*L3*r12[1];
      vz[i][1]+=m2*L1*r01[2]-m2*L3*r12[2];
      vx[i][2]+=m3*L2*r02[0]+m3*L3*r12[0];
      vy[i][2]+=m3*L2*r02[1]+m3*L3*r12[1];
      vz[i][2]+=m3*L2*r02[2]+m3*L3*r12[2];
      


    }
  

}/*End of VelocityNormal*/

void Ensemble::LongRangeCorrect(Constants cons, Parameters parms)
{
  /*Evaluate long range correction for the LJ interation*/
  double ELJTAIL, ELJSWF;
  double PLJSWF1, PLJSWF2, PLJSWF, PLJTAIL;
  double pf;
  
  double ru, ru2, ru3, ru4, ru5, ru6, ru7, ru8, ru9, ru10;
  double rl, rl2, rl3, rl4, rl5, rl6, rl7, rl8, rl9, rl10;
  double invRU, invRU2, invRU3, invRU4, invRU5, invRU6, invRU7, invRU8, invRU9;
  double invRL, invRL2, invRL3, invRL4, invRL5, invRL6, invRL7, invRL8, invRL9;
  double c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16, c17, c18, c19, c20;
  double TERM1, TERM2, TERM3, TERM1P, TERM2P, TERM3P, TERM4P, TERM5P, TERM6P;

  rl=cons.R_LOWER;
  ru=cons.R_UPPER;
  pf=2.0*cons.Pi*cons.N*cons.N;
  invRU=1.0/cons.R_UPPER;
  invRL=1.0/cons.R_LOWER;
  invRU2=invRU*invRU;
  invRU3=invRU2*invRU;
  invRU4=invRU3*invRU;
  invRU5=invRU4*invRU;
  invRU6=invRU5*invRU;
  invRU7=invRU6*invRU;
  invRU8=invRU7*invRU;
  invRU9=invRU8*invRU;
  invRL2=invRL*invRL;
  invRL3=invRL2*invRL;
  invRL4=invRL3*invRL;
  invRL5=invRL4*invRL;
  invRL6=invRL5*invRL;
  invRL7=invRL6*invRL;
  invRL8=invRL7*invRL;
  invRL9=invRL8*invRL;
  ru2=cons.R_UPPER2;
  ru3=ru2*ru;
  ru4=ru3*ru;
  ru5=ru4*ru;
  ru6=ru5*ru;
  ru7=ru6*ru;
  ru8=ru7*ru;
  ru9=ru8*ru;
  ru10=ru9*ru;
  rl2=cons.R_LOWER2;
  rl3=rl2*rl;
  rl4=rl3*rl;
  rl5=rl4*rl;
  rl6=rl5*rl;
  rl7=rl6*rl;
  rl8=rl7*rl;
  rl9=rl8*rl;
  rl10=rl9*rl;
  c1=1.0/3.0;
  c2=3.0/5.0;
  c3=3.0/7.0;
  c4=c1*c1;
  c5=16*c1;
  c6=4*c1;
  c7=2*c2;
  c8=4.0/7.0;
  c9=c2/3.0;
  c10=10*c1;
  c11=5.0/3.0*c3;
  c12=c3/3.0;
  c13=12*c2;
  c14=12*c3;
  c15=2*c13;
  c16=12*c8;
  c17=20*c3;
  c18=2*c3;
  c19=2.0/3.0*c2;
  c20=2.0*c19;
  TERM1=-parms.AS*(parms.ASPC*(-c1*(invRU3-invRL3)+c2*rl2*(invRU5-invRL5)-
			       c3*rl4*(invRU7-invRL7)+c4*rl6*(invRU9-invRL9))+
		   parms.CSPC*(-c1*(ru3-rl3)+3.0*rl2*(ru-rl)+
			       3.0*rl4*(invRU-invRL)-c1*rl6*(invRU3-invRL3)));
  
  TERM2=-parms.BS*(parms.ASPC*(-(invRU-invRL)+c6*rl2*(invRU3-invRL3)-
			       c7*rl4*(invRU5-invRL5)+c8*rl6*(invRU7-invRL7)-
			       c4*rl8*(invRU9-invRL9))+parms.CSPC*(-c9*(ru5-rl5)+
								   c6*rl2*(ru3-rl3)-6.0*rl4*(ru-rl)-
								   4.0*rl6*(invRU-invRL)+c1*rl8*(invRU3-invRL3)));
  
  TERM3=-parms.CS*(parms.ASPC*(ru-rl+5.0*rl2*(invRU-invRL)-
                 c10*rl4*(invRU3-invRL3)+2.0*rl6*(invRU5-invRL5)-
                 c11*rl8*(invRU7-invRL7)+c4*rl10*(invRU9-invRL9))+
                 parms.CSPC*(-c12*(ru7-rl7)+rl2*(ru5-rl5)-
                 c10*rl4*(ru3-rl3)+10.0*rl6*(ru-rl)+
		       5.0*rl8*(invRU-invRL)-c1*rl10*(invRU3-invRL3)));

  ELJSWF=pf*(TERM1+TERM2+TERM3);
  ELJTAIL=pf*(c4*parms.ASPC*invRU9-c1*parms.CSPC*invRU3);

  /*Return this value */
  ELRC=ELJSWF+ELJTAIL;

  TERM1P=-parms.AS*(parms.ASPC*(4.0*(invRU3-invRL3)-c13*rl2*(invRU5-invRL5)+
                  c14*rl4*(invRU7-invRL7)-c6*rl6*(invRU9-invRL9))+
                  parms.CSPC*(2.0*(ru3-rl3)-18.0*rl2*(ru-rl)-
			18.0*rl4*(invRU-invRL)+2.0*rl6*(invRU3-invRL3)));

  TERM2P=-parms.BS*(parms.ASPC*(12.0*(invRU-invRL)-16.0*rl2*(invRU3-invRL3)+
                  c15*rl4*(invRU5-invRL5)-c16*rl6*(invRU7-invRL7)+
                  c6*rl8*(invRU9-invRL9))+
                  parms.CSPC*(c7*(ru5-rl5)-8.0*rl2*(ru3-rl3)+
                  36.0*rl4*(ru-rl)+24.0*rl6*(invRU-invRL)-
			2.0*rl8*(invRU3-invRL3)));
  
  TERM3P=-parms.CS*(parms.ASPC*(-12.0*(ru-rl)-
                  60.0*rl2*(invRU-invRL)+40.0*rl4*(invRU3-invRL3)-
                  24.0*rl6*(invRU5-invRL5)+c17*rl8*(invRU7-invRL7)-
                  c6*rl10*(invRU9-invRL9))+
                  parms.CSPC*(c18*(ru7-rl7)-6.0*rl2*(ru5-rl5)+
                  20.0*rl4*(ru3-rl3)-60.0*rl6*(ru-rl)-
			30.0*rl8*(invRU-invRL)+2.0*rl10*(invRU3-invRL3)));

  TERM4P=-6.0*parms.AS*(parms.ASPC*(-c1*(invRU3-invRL3)+c19*rl2*(invRU5-invRL5)-
                        c12*rl4*(invRU7-invRL7))+parms.CSPC*(-c1*(ru3-rl3)+
						   2.0*rl2*(ru-rl)+rl4*(invRU-invRL)));

  TERM5P=-8.0*parms.BS*(parms.ASPC*(-(invRU-invRL)+rl2*(invRU3-invRL3)-
                        c2*rl4*(invRU5-invRL5)+c12*rl6*(invRU7-invRL7))+
                        parms.CSPC*(-c9*(ru5-rl5)+rl2*(ru3-rl3)-
			      3.0*rl4*(ru-rl)-rl6*(invRU-invRL)));
  
  TERM6P=-10.0*parms.CS*(parms.ASPC*(ru-rl+4.0*rl2*(invRU-invRL)-
                         2.0*rl4*(invRU3-invRL3)+c20*rl6*(invRU5-invRL5)-
                         c12*rl8*(invRU7-invRL7))+parms.CSPC*(-c12*(ru7-rl7)+
                         c20*rl2*(ru5-rl5)-2.0*rl4*(ru3-rl3)+
						    4.0*rl6*(ru-rl)+rl8*(invRU-invRL)));

  PLJSWF1=-pf*(TERM1P+TERM2P+TERM3P);
  PLJSWF2=-pf*(TERM4P+TERM5P+TERM6P);
  PLJSWF=PLJSWF1+PLJSWF2;
  PLJTAIL=-c6*parms.ASPC*invRU9+2.0*parms.CSPC*invRU3;
  PLJTAIL=-pf*PLJTAIL;

  /*Return this value as well*/

  PLRC=PLJSWF+PLJTAIL;


}//end of LongRangeCorrect

void Ensemble::ComputeAllEnergy(Ensemble &NNN, Constants cons, Parameters parms)
{
  /*All of these energy values of interest belong to the class Enesmble, so they
    do not need to be declared and will caryy over to main program
    EPOT set in calcEForce
    EKIN set in Kinetic energy, first set in ScaleVelocity
    RealVir_LJ and KVir_LJ set in calcEForce, RWALD and KWALD repectively
    ELRC and PLRC set in LongRangeCorrect
    kineCOM sete in CenterMassKinEnergy  

    ELRCint and PLRCint and VirTotal set here.
  */
  
  
  ELRCint=ELRC/parms.Volume;
  
  EKIN=NNN.KineticA(cons);
  EKIN=EKIN*.5/418.4;
  Etotal=EPOT+EKIN+ELRCint;
  cout << "EPOT: "<<EPOT <<" EKIN: "<<EKIN <<" ELRCint: "<<ELRCint<<endl;
  cout <<"Etot: "<<Etotal<<endl;
  
  /*Initial Center of Mass Kinteic Energy*/
  kineCOM=NNN.CenterMassKinEnergy(cons);

  VirTotal=RealVir_LJ+KVir_LJ;
  PLRCint=PLRC/parms.Volume;
  
  Pint=(kineCOM+PLRCint+VirTotal)/parms.pconstant;
  
  cout <<"RealVir_LJ: "<<RealVir_LJ <<"KVir_LJ: "<<KVir_LJ<<endl;
  cout <<"Pint: "<< Pint<<endl;
  cout <<"PLRCint: "<<PLRCint<<" kineCOM: "<<kineCOM<<" pconst: "<< parms.pconstant <<endl;


}/*End of ComputeAllEnergy*/
