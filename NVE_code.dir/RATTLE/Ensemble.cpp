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

void readParms(char* filename, Constants &cons, Parameters &parms)
{
  /*look in  parameter.h  for a description of what each value represents.*/
  
  ifstream PARMFILE(filename);
  if(!PARMFILE) cerr << " NO parm file given\n"; 

  PARMFILE >> cons.Tsteps;// 
  PARMFILE.ignore(400, '\n');
  PARMFILE >> cons.Tequil;// 
  PARMFILE.ignore(400, '\n');
  PARMFILE >> cons.Tvscale;// 
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
  PARMFILE >> parms.isave;
  PARMFILE.ignore(400, '\n');
  PARMFILE >> parms.statwrite;
  PARMFILE.ignore(400, '\n');
  PARMFILE >> parms.allsave;
  PARMFILE.ignore(400, '\n');
  PARMFILE >> parms.nperfile;
  PARMFILE.ignore(400, '\n');

  PARMFILE.close();
}

void readCoords(char *coordFileIn, WaterSend *mols)
{

  /*NATOMS=3 for water*/
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

void readVelocity(char *velFileIn, DiffSend *vels)
{
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


void InitialParms(Parameters &parms, Constants &cons, Averages &avg, Ensemble &nnn)
{
  /*Only after read in certain parameters/constants from file*/  

  double RUP_LOW;
  double RUP_LOW2;
  double RUP_LOW3;
  double RUP_LOW4;
  double RUP_LOW5;
  double RMFact=332.0637349768129609;
  double PFACT=0.0000145839756;
  double DNF=cons.N*(3.0*cons.NUMPRIM-1.0*cons.NCONST)-3.0;
  double RSPC2, RSPC6, RSPC12;
  
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
  
  /*Parameterized Values*/

  parms.iterRat=50;
  cons.ALPHA=0.35;
  cons.Alpha2=cons.ALPHA*cons.ALPHA;
  cons.SIGMA=1.0/(4.0*cons.Alpha2);
  parms.CoeffPi=2.0/sqrt(cons.Pi);
  parms.ConstM=DNF*cons.BOLTZK;

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
  parms.pconstant=1.0/(3.0*parms.Volume*PFACT); //For Pressure conversion

  cout << "MOLECULE Density: " << parms.Density <<endl;
  cout << "BOXLENGTH: "  << parms.BOXLENGTH << endl;
  cout << "VOLUME "<<parms.Volume << endl;
  
  /*Set values of parameters Q is charge 2, 3 are H  4 is Oxygen*/  
 
  parms.Q2=parms.Q2*sqrt(RMFact);
  parms.Q3=parms.Q3*sqrt(RMFact);
  parms.Q4=parms.Q4*sqrt(RMFact);

  /*Coeffieicents for LJ potential*/  
 
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

  /*Kspace Setup values*/
  cons.Kmax=11;
  cons.Kmin=-10;
  cons.Kscut=105;
  cons.Kleng=0;

  /*Initialize enery and temperature quantities*/

  avg.Eavg=0.0;
  avg.E2avg=0.0;
  avg.Epavg=0.0;
  avg.Ekavg=0.0;
  avg.Tavg=0.0;
  avg.T2avg=0.0;
  avg.Pavg=0.0;
  avg.P2avg=0.0;
  nnn.Etotal=0.0;
  nnn.InstTemp=0.0;
  nnn.Pint=0.0;
  nnn.EPOT=0.0;
  nnn.EKIN=0.0;
  nnn.KVirial=0.0;

}//End of function Initial


void RigidMol(Parameters &parms, Constants cons, WaterSend * mols, DiffSend *opos, Center *com)
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

  if(parms.BOXLOLD>2.0) cout<<" Scaling Initial config to BOXLOLD! "<<parms.BOXLOLD<<endl;
  else cout <<" no boxlength scaling of initial config. "<<endl;
  /*
    scale the old configuration IF the value in input file of
    BOXLOLD is greater than 2. Otherwise will not scale!!!!
  */
  
  //find center of mass

  double sumx=0.0;
  double sumy=0.0;
  double sumz=0.0;
  double asum=0.0;
  int i;
  double am0=parms.AMass[0];
  double am1=parms.AMass[1];
  double am2=parms.AMass[2];
  double DotProd;
  double RB, XB, YB, ZB;

  for(i=0;i<cons.N;i++)
    {
      sumx=sumy=sumz=asum=0.0;
      sumx=am0*mols[i].x[0]+am1*mols[i].x[1]+am2*mols[i].x[2];
      sumy=am0*mols[i].y[0]+am1*mols[i].y[1]+am2*mols[i].y[2];
      sumz=am0*mols[i].z[0]+am1*mols[i].z[1]+am2*mols[i].z[2];
      asum=am0+am1+am2;
      
      //COM
      com[i].x=sumx/asum;
      com[i].y=sumy/asum;
      com[i].z=sumz/asum;
      
      //new positions O, H, H
      if(parms.BOXLOLD>2.0)
	{
	  mols[i].x[0]+=(SCBOXL-1.0)*com[i].x;
	  mols[i].y[0]+=(SCBOXL-1.0)*com[i].y;
	  mols[i].z[0]+=(SCBOXL-1.0)*com[i].z;
	  
	  mols[i].x[1]+=(SCBOXL-1.0)*com[i].x;
	  mols[i].y[1]+=(SCBOXL-1.0)*com[i].y;
	  mols[i].z[1]+=(SCBOXL-1.0)*com[i].z;
	  
	  mols[i].x[2]+=(SCBOXL-1.0)*com[i].x;
	  mols[i].y[2]+=(SCBOXL-1.0)*com[i].y;
	  mols[i].z[2]+=(SCBOXL-1.0)*com[i].z;
	}
      /*Oxygen and Hydrogen positions for diffusion*/
      opos[i].x[0]=mols[i].x[0];
      opos[i].y[0]=mols[i].y[0];
      opos[i].z[0]=mols[i].z[0];
      opos[i].x[1]=mols[i].x[1];
      opos[i].y[1]=mols[i].y[1];
      opos[i].z[1]=mols[i].z[1];
      opos[i].x[2]=mols[i].x[2];
      opos[i].y[2]=mols[i].y[2];
      opos[i].z[2]=mols[i].z[2];
      /*Distance Hyd1-Oxygen*/
      XH10=mols[i].x[1]-mols[i].x[0];
      YH10=mols[i].y[1]-mols[i].y[0];
      ZH10=mols[i].z[1]-mols[i].z[0];
      
      XH10-=parms.BOXLENGTH*round(XH10*parms.InvBOXL);
      YH10-=parms.BOXLENGTH*round(YH10*parms.InvBOXL);
      ZH10-=parms.BOXLENGTH*round(ZH10*parms.InvBOXL);
      
      R10sq=XH10*XH10+YH10*YH10+ZH10*ZH10;
      R10=sqrt(R10sq);

      /*Distance Hyd2-Oxygen*/
      XH20=mols[i].x[2]-mols[i].x[0];
      YH20=mols[i].y[2]-mols[i].y[0];
      ZH20=mols[i].z[2]-mols[i].z[0];

      XH20-=parms.BOXLENGTH*round(XH20*parms.InvBOXL);
      YH20-=parms.BOXLENGTH*round(YH20*parms.InvBOXL);
      ZH20-=parms.BOXLENGTH*round(ZH20*parms.InvBOXL);
      
      R20sq=XH20*XH20+YH20*YH20+ZH20*ZH20;
      R20=sqrt(R20sq);

      /*Angle between vectors OH1 and OH2*/
      DotProd=XH20*XH10+YH20*YH10+ZH20*ZH10;
      cosThetaRad=DotProd/R20/R10;
      ThetaRad=acos(cosThetaRad);
      Thetad=ThetaRad*(cons.PiDEG/cons.Pi);
      
      /*vector of H-H (no bond there)*/
      XHH=mols[i].x[2]-mols[i].x[1];
      YHH=mols[i].y[2]-mols[i].y[1];
      ZHH=mols[i].z[2]-mols[i].z[1];
      
      XHH-=parms.BOXLENGTH*round(XHH*parms.InvBOXL);
      YHH-=parms.BOXLENGTH*round(YHH*parms.InvBOXL);
      ZHH-=parms.BOXLENGTH*round(ZHH*parms.InvBOXL);
      
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


void writeParms(Parameters parms, Constants cons)
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

void printCoords(char *filename, Constants cons, WaterSend *mols, int t)
{
  int i;
  ofstream outFile(filename);
  outFile << "Time: "<<t<<endl;
  for(i=0;i<cons.N;i++)
    {
      outFile.precision(16);
      outFile << mols[i].x[0] << "  " <<mols[i].y[0]<< "  " <<mols[i].z[0]<<endl;
      outFile << mols[i].x[1] << "  " <<mols[i].y[1]<< "  " <<mols[i].z[1]<<endl;
      outFile << mols[i].x[2] << "  " <<mols[i].y[2]<< "  " <<mols[i].z[2]<<endl;
    }
  outFile.close();
}

void printAllCoords(char *filename, Constants cons, WaterSend *mols)
{
  int i;
  ofstream outFile(filename);

  for(i=0;i<cons.N;i++)
    {
      outFile.precision(10);
      outFile << mols[i].x[0] << "  " <<mols[i].y[0]<< "  " <<mols[i].z[0]<<endl;
      outFile << mols[i].x[1] << "  " <<mols[i].y[1]<< "  " <<mols[i].z[1]<<endl;
      outFile << mols[i].x[2] << "  " <<mols[i].y[2]<< "  " <<mols[i].z[2]<<endl;
      outFile << mols[i].x[3] << "  " <<mols[i].y[3]<< "  " <<mols[i].z[3]<<endl;
    }
  outFile.close();
}

void printAllMyCoords(char *filename, Constants cons, WaterSend *mols)
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

void writeForce(char *filename, Constants cons, WaterSend *force)
{
  int i;
  ofstream outFile(filename);

  for(i=0;i<cons.NperP;i++)
    {
      outFile.precision(16);
      outFile << force[i].x[0] << "  " <<force[i].y[0]<< "  " <<force[i].z[0]<<endl;
      outFile << force[i].x[1] << "  " <<force[i].y[1]<< "  " <<force[i].z[1]<<endl;
      outFile << force[i].x[2] << "  " <<force[i].y[2]<< "  " <<force[i].z[2]<<endl;
    }
  outFile.close();
}

void writeAllForce(char *filename, Constants cons, WaterSend *force)
{
  int i;
  ofstream outFile(filename);

  for(i=0;i<cons.NperP;i++)
    {
      outFile.precision(16);
      outFile << force[i].x[0] << "  " <<force[i].y[0]<< "  " <<force[i].z[0]<<endl;
      outFile << force[i].x[1] << "  " <<force[i].y[1]<< "  " <<force[i].z[1]<<endl;
      outFile << force[i].x[2] << "  " <<force[i].y[2]<< "  " <<force[i].z[2]<<endl;
      outFile << force[i].x[3] << "  " <<force[i].y[3]<< "  " <<force[i].z[3]<<endl;
    }
  outFile.close();
}

void printVel(char *filename, Constants cons, int t, DiffSend *vel)
{
  int i;
  ofstream outFile(filename);
  
  outFile << "Veloc: "<< t<<" Rank: "<< cons.rank<<endl;
  outFile.precision(16);
  for(i=0;i<cons.NperP;i++)
    {
      
      outFile << vel[i].x[0] << "  " <<vel[i].y[0]<< "  " <<vel[i].z[0]<<endl;
      outFile << vel[i].x[1] << "  " <<vel[i].y[1]<< "  " <<vel[i].z[1]<<endl;
      outFile << vel[i].x[2] << "  " <<vel[i].y[2]<< "  " <<vel[i].z[2]<<endl;
    }
  outFile.close();
}




void RWALD(Constants cons, Parameters parms, WaterSend *mols, WaterSend *local_mols, WaterSend *force, Center *com, Ensemble &nnn, Table *terf)
{
  //Energy and derivatives in real space
  /*Want E0  EC  RealVir_LJ */
  
  //Rlower is the lower limit of the distance between the 2 molecules
 
  //do solvent-solvent interactions
  //this part must take into consideration that some of your neighbor 
  //molecules might be that way because they are close to you on the 
  //boundaries, i.e. if you are close to the edge of the box, and they are
  // opposite you then you will interact with them  
    
  //E0, EC, ALPHA, 
  //FOR WATER ASPC, CSPC, RSPC, ESPC, Q(N, 4)
  //AS, BS, CS are all set in function InitialParms

  int i;
  int j, k, n;
  double Fx0, Fy0, Fz0;
  double FXnew, FYnew, FZnew;
  double rX0, rY0, rZ0;
  double RXnew, RYnew, RZnew;
  double xCMnow, yCMnow, zCMnow;

  double r02;
  double r2, r6, r12;
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
  double rn1;
  int nindx1;
  double fact1;

  int my=cons.rank*cons.NperP;
  
  /*Return these: Energy Oxygens (LJ) and Energy Coulomb (H, H, M)*/
  nnn.E0=0.0;
  nnn.EC=0.0;

  /*Update all forces on i, loop j over all molecules, or create neighbor
   list and loop j over those.  
   do not interact with self: hence j!=(i+rank*cons.NperP) expression 
  */
  
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
	       
	       rX0=rX0-parms.BOXLENGTH*round(rX0*parms.InvBOXL);
	       rY0=rY0-parms.BOXLENGTH*round(rY0*parms.InvBOXL);
	       rZ0=rZ0-parms.BOXLENGTH*round(rZ0*parms.InvBOXL);
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
		   nnn.E0=nnn.E0+ELJ;
		   
		   force[i].x[0]+=Fx0;
		   force[i].y[0]+=Fy0;
		   force[i].z[0]+=Fz0;

		   //distance between 2 molecules Center Of Mass
		   xCMnow=com[i+my].x-com[j].x;
		   yCMnow=com[i+my].y-com[j].y;
		   zCMnow=com[i+my].z-com[j].z;
		   xCMnow-=parms.BOXLENGTH*round(xCMnow*parms.InvBOXL);
		   yCMnow-=parms.BOXLENGTH*round(yCMnow*parms.InvBOXL);
		   zCMnow-=parms.BOXLENGTH*round(zCMnow*parms.InvBOXL);

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
			   RXnew-=parms.BOXLENGTH*round(RXnew*parms.InvBOXL);
			   RYnew-=parms.BOXLENGTH*round(RYnew*parms.InvBOXL);
			   RZnew-=parms.BOXLENGTH*round(RZnew*parms.InvBOXL);
			   RijSQUARED=RXnew*RXnew+RYnew*RYnew+RZnew*RZnew;
			   Rij1=sqrt(RijSQUARED);
			   Q1=parms.Q[k]*parms.Q[n];
			   /*Add in table values*/
			   nindx1=int(Rij1*parms.irincr);
			   rn1=parms.rincr*nindx1;
			   fact1=Rij1-rn1;

			   Term1=Q1*(terf[nindx1].a+fact1*(terf[nindx1].b+0.5*fact1*terf[nindx1].c));
			   Term1A=Q1*(terf[nindx1].b+fact1*(terf[nindx1].c+0.5*fact1*terf[nindx1].d))/Rij1;
			   
			   FXnew=-Term1A*RXnew;
			   FYnew=-Term1A*RYnew;
			   FZnew=-Term1A*RZnew;
			   SumQQ+=Term1;
			   
			   force[i].x[k]+=FXnew;
			   force[i].y[k]+=FYnew;
			   force[i].z[k]+=FZnew;
			   
			   
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
	       
	       else if(r02<=cons.R_UPPER2)
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
		   nnn.E0=nnn.E0+ELJ;
		   
		   //update particle forces Oxygen
		   force[i].x[0]+=Fx0;
		   force[i].y[0]+=Fy0;
		   force[i].z[0]+=Fz0;
	   
	      //update center of mass
	      
		   xCMnow=com[i+my].x-com[j].x;
		   yCMnow=com[i+my].y-com[j].y;
		   zCMnow=com[i+my].z-com[j].z;
		   xCMnow-=parms.BOXLENGTH*round(xCMnow*parms.InvBOXL);
		   yCMnow-=parms.BOXLENGTH*round(yCMnow*parms.InvBOXL);
		   zCMnow-=parms.BOXLENGTH*round(zCMnow*parms.InvBOXL);
		   
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
			   RXnew-=parms.BOXLENGTH*round(RXnew*parms.InvBOXL);
			   RYnew-=parms.BOXLENGTH*round(RYnew*parms.InvBOXL);
			   RZnew-=parms.BOXLENGTH*round(RZnew*parms.InvBOXL);
			   RijSQUARED=RXnew*RXnew+RYnew*RYnew+RZnew*RZnew;
			   Rij1=sqrt(RijSQUARED);
			   Q1=parms.Q[k]*parms.Q[n];
			   /*Add in table values*/
			   nindx1=int(Rij1*parms.irincr);
			   rn1=parms.rincr*nindx1;
			   fact1=Rij1-rn1;
			   Term1=Q1*(terf[nindx1].a+fact1*(terf[nindx1].b+0.5*fact1*terf[nindx1].c));
			   Term1A=Q1*(terf[nindx1].b+fact1*(terf[nindx1].c+0.5*fact1*terf[nindx1].d))/Rij1;
			   Term1AA=Term1A*SZR;

			   FXnew=-Term1AA*RXnew;
			   FYnew=-Term1AA*RYnew;
			   FZnew=-Term1AA*RZnew;
			   SumQQ+=Term1*SZR;
			   force[i].x[k]+=FXnew;
			   force[i].y[k]+=FYnew;
			   force[i].z[k]+=FZnew;
			   
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
	       nnn.EC+=SumQQ;
	       Vir_LrXX+=VirijXX;
	       Vir_LrYY+=VirijYY;
	       Vir_LrZZ+=VirijZZ;
	       
	     }
	}//end of loop over second molecule in pair
    }// end of loop over pairs
  
  /*Return RealVir_LJ, E0, and EC. SUM all 3 across Processors!!!*/
  nnn.RealVir_LJ=Vir_LrXX+Vir_LrYY+Vir_LrZZ;

  nnn.Energyv[0]=nnn.E0; //Energy oxygen
  nnn.Energyv[1]=nnn.EC; //Coulomb charges
  nnn.Energyv[2]=nnn.RealVir_LJ; // Virial from real space
  
}// end of function rwald


void Primder(Constants cons, Parameters parms, WaterSend *local_mols, WaterSend *force )
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
      Fr30=X30*force[i].x[3]+Y30*force[i].y[3]+Z30*force[i].z[3];
      Ft=Fr30/Rsq30;

      //center of mass of forces
      ComFx=Gama*(force[i].x[3]-Ft*X30);
      ComFy=Gama*(force[i].y[3]-Ft*Y30);
      ComFz=Gama*(force[i].z[3]-Ft*Z30);
      
      Fx0=force[i].x[3]-ComFx;
      Fy0=force[i].y[3]-ComFy;
      Fz0=force[i].z[3]-ComFz;
      
      FxH=0.5*ComFx;
      FyH=0.5*ComFy;
      FzH=0.5*ComFz;
      
      //Now update forces of real atoms
      force[i].x[0]+=Fx0;
      force[i].y[0]+=Fy0;
      force[i].z[0]+=Fz0;
      force[i].x[1]+=FxH;
      force[i].y[1]+=FyH;
      force[i].z[1]+=FzH;
      force[i].x[2]+=FxH;
      force[i].y[2]+=FyH;
      force[i].z[2]+=FzH;

    }//end loop over NperP

}//End function primder


void Maxwell(Constants cons, Parameters parms, DiffSend *InitVel)
{
  int i, j;
  double Sigma0, Sigma1, Sigma2;
  double Rn;
  double sumx, sumy, sumz, amsum;
  Sigma0=sqrt(cons.BOLTZK*cons.Temp/parms.AMass[0]);
  Sigma1=sqrt(cons.BOLTZK*cons.Temp/parms.AMass[1]);
  Sigma2=sqrt(cons.BOLTZK*cons.Temp/parms.AMass[2]);

  for(i=0;i<cons.N;i++)
    {

      Rn=GaussV();
      InitVel[i].x[0]=Sigma0*Rn;
      Rn=GaussV();
      InitVel[i].y[0]=Sigma0*Rn;
      Rn=GaussV();
      InitVel[i].z[0]=Sigma0*Rn;


      Rn=GaussV();
      InitVel[i].x[1]=Sigma1*Rn;
      Rn=GaussV();
      InitVel[i].y[1]=Sigma1*Rn;
      Rn=GaussV();
      InitVel[i].z[1]=Sigma1*Rn;


      Rn=GaussV();
      InitVel[i].x[2]=Sigma2*Rn;
      Rn=GaussV();
      InitVel[i].y[2]=Sigma2*Rn;
      Rn=GaussV();
      InitVel[i].z[2]=Sigma2*Rn;
    }
  /*Make the Center of Mass of bulk water a Fixed point. Sum_x(everyAtom_x*hisx_velocity)
   amsum will be 18.01528*N  */
  
  sumx=0.0;
  sumy=0.0;
  sumz=0.0;
  amsum=0.0; 
  for(i=0;i<cons.N;i++)
    {
      for(j=0;j<3;j++)
	{
	  
	  sumx+=parms.AMass[j]*InitVel[i].x[j];
	  sumy+=parms.AMass[j]*InitVel[i].y[j];
	  sumz+=parms.AMass[j]*InitVel[i].z[j];
	  amsum+=parms.AMass[j];
	}
    }
  sumx=sumx/amsum;
  sumy=sumy/amsum;
  sumz=sumz/amsum;
  for(i=0;i<cons.N;i++)
    {
      InitVel[i].x[0]-=sumx;
      InitVel[i].y[0]-=sumy;
      InitVel[i].z[0]-=sumz;
      
      InitVel[i].x[1]-=sumx;
      InitVel[i].y[1]-=sumy;
      InitVel[i].z[1]-=sumz;

      InitVel[i].x[2]-=sumx;
      InitVel[i].y[2]-=sumy;
      InitVel[i].z[2]-=sumz;
    }

}//end of Maxwell


double KineticA(Constants cons, DiffSend *vel)
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
	  Kin0+=vel[i].x[0]*vel[i].x[0]+vel[i].y[0]*vel[i].y[0]+vel[i].z[0]*vel[i].z[0];
	  KinH1+=vel[i].x[1]*vel[i].x[1]+vel[i].y[1]*vel[i].y[1]+vel[i].z[1]*vel[i].z[1];
	  KinH2+=vel[i].x[2]*vel[i].x[2]+vel[i].y[2]*vel[i].y[2]+vel[i].z[2]*vel[i].z[2];
    }

  KinE=am0*Kin0+am1*KinH1+am2*KinH2;
  return KinE;
}
/*
CANNOT USE Scale Velocity unless you first do an MPIReduce to get the proper full Kinetic Energy
*/

void CenterMassWater(Constants cons, WaterSend *mols, Center *com)
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

      com[i].x=xsum/amsum;
      com[i].y=ysum/amsum;
      com[i].z=zsum/amsum;
    }
}


double CenterMassKinEnergy(Constants cons, DiffSend *vel)
{
  /*Return in units kcal/mol*/
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
      
      vxsum=vel[i].x[0]*am0+vel[i].x[1]*am1+vel[i].x[2]*am2;
      vysum=vel[i].y[0]*am0+vel[i].y[1]*am1+vel[i].y[2]*am2;
      vzsum=vel[i].z[0]*am0+vel[i].z[1]*am1+vel[i].z[2]*am2;
      
      //center of mass = sum_3 (x(i)*m(i)) / sum_3 (m(i))

      vxCOM=vxsum;// /amsum;
      vyCOM=vysum;// /amsum;
      vzCOM=vzsum;// /amsum;
      KEcom+=(vxCOM*vxCOM+vyCOM*vyCOM+vzCOM*vzCOM)/amsum;
    }
  KEcom=KEcom/418.4;
  return KEcom;

}


void RattlePos(Constants cons, Parameters parms, WaterSend *local_mols, DiffSend *myopos, DiffSend *vel, double *Mat, double *Minv)
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
  double s01[3];
  double s02[3];
  double s12[3];
  double r01[3];
  double r02[3];
  double r12[3];
  

  //do position update step but don't move particles yet

  for(i=0;i<cons.NperP;i++)
    {
      
      sx0=local_mols[i].x[0]+cons.dt*vel[i].x[0];
      sy0=local_mols[i].y[0]+cons.dt*vel[i].y[0];
      sz0=local_mols[i].z[0]+cons.dt*vel[i].z[0];
      sx1=local_mols[i].x[1]+cons.dt*vel[i].x[1];
      sy1=local_mols[i].y[1]+cons.dt*vel[i].y[1];
      sz1=local_mols[i].z[1]+cons.dt*vel[i].z[1];
      sx2=local_mols[i].x[2]+cons.dt*vel[i].x[2];
      sy2=local_mols[i].y[2]+cons.dt*vel[i].y[2];
      sz2=local_mols[i].z[2]+cons.dt*vel[i].z[2];

      /*Now calculate relative vectors*/
      //O to H1 distance
      r01[0]=local_mols[i].x[0]-local_mols[i].x[1];
      r01[1]=local_mols[i].y[0]-local_mols[i].y[1];
      r01[2]=local_mols[i].z[0]-local_mols[i].z[1];
      r01[0]-=parms.BOXLENGTH*round(r01[0]*parms.InvBOXL);
      r01[1]-=parms.BOXLENGTH*round(r01[1]*parms.InvBOXL);
      r01[2]-=parms.BOXLENGTH*round(r01[2]*parms.InvBOXL);
      
      /*distance O to H2*/
      r02[0]=local_mols[i].x[0]-local_mols[i].x[2];
      r02[1]=local_mols[i].y[0]-local_mols[i].y[2];
      r02[2]=local_mols[i].z[0]-local_mols[i].z[2];
      r02[0]-=parms.BOXLENGTH*round(r02[0]*parms.InvBOXL);
      r02[1]-=parms.BOXLENGTH*round(r02[1]*parms.InvBOXL);
      r02[2]-=parms.BOXLENGTH*round(r02[2]*parms.InvBOXL);

      /*ditsnce H to H*/
      r12[0]=local_mols[i].x[1]-local_mols[i].x[2];
      r12[1]=local_mols[i].y[1]-local_mols[i].y[2];
      r12[2]=local_mols[i].z[1]-local_mols[i].z[2];
      r12[0]-=parms.BOXLENGTH*round(r12[0]*parms.InvBOXL);
      r12[1]-=parms.BOXLENGTH*round(r12[1]*parms.InvBOXL);
      r12[2]-=parms.BOXLENGTH*round(r12[2]*parms.InvBOXL);

      /*Iterative loop begins*/
      NumIt=0;
      it=0;
      while(it<parms.iterRat)
	{
	  /*are the pbc ajdustments necessary,  only if o and h move a half a boxlength apart*/

	  s01[0]=sx0-sx1;
	  s01[1]=sy0-sy1;
	  s01[2]=sz0-sz1;
	  s01[0]-=parms.BOXLENGTH*round(s01[0]*parms.InvBOXL);
	  s01[1]-=parms.BOXLENGTH*round(s01[1]*parms.InvBOXL);
	  s01[2]-=parms.BOXLENGTH*round(s01[2]*parms.InvBOXL);
	  
	  s02[0]=sx0-sx2;
	  s02[1]=sy0-sy2;
	  s02[2]=sz0-sz2;
	  s02[0]-=parms.BOXLENGTH*round(s02[0]*parms.InvBOXL);
	  s02[1]-=parms.BOXLENGTH*round(s02[1]*parms.InvBOXL);
	  s02[2]-=parms.BOXLENGTH*round(s02[2]*parms.InvBOXL);
	  
	  s12[0]=sx1-sx2;
	  s12[1]=sy1-sy2;
	  s12[2]=sz1-sz2;
	  s12[0]-=parms.BOXLENGTH*round(s12[0]*parms.InvBOXL);
	  s12[1]-=parms.BOXLENGTH*round(s12[1]*parms.InvBOXL);
	  s12[2]-=parms.BOXLENGTH*round(s12[2]*parms.InvBOXL);

	  s01sq=s01[0]*s01[0]+s01[1]*s01[1]+s01[2]*s01[2];
	  s02sq=s02[0]*s02[0]+s02[1]*s02[1]+s02[2]*s02[2];
	  s12sq=s12[0]*s12[0]+s12[1]*s12[1]+s12[2]*s12[2];
	  
	  //should be practically zero
	  sig01=abs(s01sq-parms.d01sq)/(2.0*parms.d01sq);
	  sig02=abs(s02sq-parms.d02sq)/(2.0*parms.d02sq);
	  sig12=abs(s12sq-parms.d12sq)/(2.0*parms.d12sq);
	  
	  if(sig01<=parms.tol&&sig02<=parms.tol&&sig12<=parms.tol)
	    it=parms.iterRat+1;
	  /*conitue out of while loop, constraints are met*/
	  
	  //calculate scalar products r is original vector, s is moved vec
	  r01s01=r01[0]*s01[0]+r01[1]*s01[1]+r01[2]*s01[2];//dotProd(r01, s01);
	  r02s01=r02[0]*s01[0]+r02[1]*s01[1]+r02[2]*s01[2];//dotProd(r02, s01);
	  r12s01=r12[0]*s01[0]+r12[1]*s01[1]+r12[2]*s01[2];
	  r01s02=r01[0]*s02[0]+r01[1]*s02[1]+r01[2]*s02[2];//dotProd(r01, s02);
	  r02s02=r02[0]*s02[0]+r02[1]*s02[1]+r02[2]*s02[2];//dotProd(r02, s02);
	  r12s02=r12[0]*s02[0]+r12[1]*s02[1]+r12[2]*s02[2];//dotProd(r12, s02);
	  r01s12=r01[0]*s12[0]+r01[1]*s12[1]+r01[2]*s12[2];//dotProd(r01, s12);
	  r02s12=r02[0]*s12[0]+r02[1]*s12[1]+r02[2]*s12[2];//dotProd(r02, s12);
	  r12s12=r12[0]*s12[0]+r12[1]*s12[1]+r12[2]*s12[2];//dotProd(r12, s12);
	  
	  //calculate matrix elements
	  
	  Mat[0]=r01s01*(parms.InvMass[0]+parms.InvMass[1]);
	  Mat[1]=r02s01*parms.InvMass[0];
	  Mat[2]=-r12s01*parms.InvMass[1];
	  Mat[3]=r01s02*parms.InvMass[0];
	  Mat[4]=r02s02*(parms.InvMass[0]+parms.InvMass[2]);
	  Mat[5]=r12s02*parms.InvMass[2];
	  Mat[6]=-r01s12*parms.InvMass[1];
	  Mat[7]=r02s12*parms.InvMass[2];
	  Mat[8]=r12s12*(parms.InvMass[1]+parms.InvMass[2]);
	  
	  cm11=Invdt2*(s01sq-parms.d01sq);
	  cm21=Invdt2*(s02sq-parms.d02sq);
	  cm31=Invdt2*(s12sq-parms.d12sq);
	  
	  //invert the shake matrix
	  MatrixInv(Mat, Minv);

	  //obtain solutions of linearized equation
	  gr01=Minv[0]*cm11+Minv[1]*cm21+Minv[2]*cm31;
	  gr02=Minv[3]*cm11+Minv[4]*cm21+Minv[5]*cm31;
	  gr12=Minv[6]*cm11+Minv[7]*cm21+Minv[8]*cm31;
	  
	  //obtain velocity adjusted
	  vel[i].x[0]-=cons.dt*parms.InvMass[0]*(gr01*r01[0]+gr02*r02[0]);
	  vel[i].y[0]-=cons.dt*parms.InvMass[0]*(gr01*r01[1]+gr02*r02[1]);
	  vel[i].z[0]-=cons.dt*parms.InvMass[0]*(gr01*r01[2]+gr02*r02[2]);
	  vel[i].x[1]-=cons.dt*parms.InvMass[1]*(-gr01*r01[0]+gr12*r12[0]);
	  vel[i].y[1]-=cons.dt*parms.InvMass[1]*(-gr01*r01[1]+gr12*r12[1]);
	  vel[i].z[1]-=cons.dt*parms.InvMass[1]*(-gr01*r01[2]+gr12*r12[2]);
	  vel[i].x[2]-=cons.dt*parms.InvMass[2]*(-gr02*r02[0]-gr12*r12[0]);
	  vel[i].y[2]-=cons.dt*parms.InvMass[2]*(-gr02*r02[1]-gr12*r12[1]);
	  vel[i].z[2]-=cons.dt*parms.InvMass[2]*(-gr02*r02[2]-gr12*r12[2]);
 
	  //new coords
	  sx0=local_mols[i].x[0]+cons.dt*vel[i].x[0];
	  sy0=local_mols[i].y[0]+cons.dt*vel[i].y[0];
	  sz0=local_mols[i].z[0]+cons.dt*vel[i].z[0];
	  sx1=local_mols[i].x[1]+cons.dt*vel[i].x[1];
	  sy1=local_mols[i].y[1]+cons.dt*vel[i].y[1];
	  sz1=local_mols[i].z[1]+cons.dt*vel[i].z[1];
	  sx2=local_mols[i].x[2]+cons.dt*vel[i].x[2];
	  sy2=local_mols[i].y[2]+cons.dt*vel[i].y[2];
	  sz2=local_mols[i].z[2]+cons.dt*vel[i].z[2];
	  
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
	  local_mols[i].x[j]+=cons.dt*vel[i].x[j];
	  local_mols[i].y[j]+=cons.dt*vel[i].y[j];
	  local_mols[i].z[j]+=cons.dt*vel[i].z[j];
	  myopos[i].x[j]+=cons.dt*vel[i].x[j];
	  myopos[i].y[j]+=cons.dt*vel[i].y[j];
	  myopos[i].z[j]+=cons.dt*vel[i].z[j];
	}
      //do periodic boundary conditions for Oxygen, then adjust all atoms
      xpbc=parms.BOXLENGTH*round(local_mols[i].x[0]*parms.InvBOXL);
      ypbc=parms.BOXLENGTH*round(local_mols[i].y[0]*parms.InvBOXL);
      zpbc=parms.BOXLENGTH*round(local_mols[i].z[0]*parms.InvBOXL);
      
      //move O, H, H
      local_mols[i].x[0]-=xpbc;
      local_mols[i].y[0]-=ypbc;
      local_mols[i].z[0]-=zpbc;
      local_mols[i].x[1]-=xpbc;
      local_mols[i].y[1]-=ypbc;
      local_mols[i].z[1]-=zpbc;
      local_mols[i].x[2]-=xpbc;
      local_mols[i].y[2]-=ypbc;
      local_mols[i].z[2]-=zpbc;

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


void MatrixInv(double *Shake, double *Sinv)
{
  double stol=1.0E-7;
  double determ;
  /* Matrix  : 0 : 1 : 2 :
             : 3 : 4 : 5 :
             : 6 : 7 : 8 : */
  //get signed cofactors, transpose and put in skinv
  
  Sinv[0]=Shake[4]*Shake[8]-Shake[5]*Shake[7];
  Sinv[1]=Shake[2]*Shake[7]-Shake[1]*Shake[8];
  Sinv[2]=Shake[1]*Shake[5]-Shake[2]*Shake[4];
  Sinv[3]=Shake[5]*Shake[6]-Shake[3]*Shake[8];
  Sinv[4]=Shake[0]*Shake[8]-Shake[2]*Shake[6];
  Sinv[5]=Shake[2]*Shake[3]-Shake[0]*Shake[5];
  Sinv[6]=Shake[3]*Shake[7]-Shake[4]*Shake[6];
  Sinv[7]=Shake[1]*Shake[6]-Shake[0]*Shake[7];
  Sinv[8]=Shake[0]*Shake[4]-Shake[1]*Shake[3];

  //get determinant and make sure not zero

  determ=Shake[0]*Sinv[0]+Shake[1]*Sinv[3]+Shake[2]*Sinv[6];
  
  if(abs(determ)<stol)
    {
      cerr <<"zero determinant in ShakeInv "<<endl;
      exit(1);
    }
  Sinv[0]/=determ;
  Sinv[1]/=determ;
  Sinv[2]/=determ;
  Sinv[3]/=determ;
  Sinv[4]/=determ;
  Sinv[5]/=determ;
  Sinv[6]/=determ;
  Sinv[7]/=determ;
  Sinv[8]/=determ;


}

void RattleVel(Constants cons, Parameters parms, WaterSend *local_mols, DiffSend *vel, double *Mat, double *Minv)
{
  int i;
  double r01[3];
  double r02[3];
  double r12[3];
  double v01[3];
  double v02[3];
  double v12[3];

  double gv01, gv02, gv12;
  double qr11, qr21, qr31;
  double r01r01, r01r02, r01r12, r02r02, r02r12, r12r12;

  for(i=0;i<cons.NperP;i++)
    {
      r01[0]=local_mols[i].x[0]-local_mols[i].x[1];
      r01[1]=local_mols[i].y[0]-local_mols[i].y[1];
      r01[2]=local_mols[i].z[0]-local_mols[i].z[1];
      r01[0]-=parms.BOXLENGTH*round(r01[0]*parms.InvBOXL);
      r01[1]-=parms.BOXLENGTH*round(r01[1]*parms.InvBOXL);
      r01[2]-=parms.BOXLENGTH*round(r01[2]*parms.InvBOXL);
      
      r02[0]=local_mols[i].x[0]-local_mols[i].x[2];
      r02[1]=local_mols[i].y[0]-local_mols[i].y[2];
      r02[2]=local_mols[i].z[0]-local_mols[i].z[2];
      r02[0]-=parms.BOXLENGTH*round(r02[0]*parms.InvBOXL);
      r02[1]-=parms.BOXLENGTH*round(r02[1]*parms.InvBOXL);
      r02[2]-=parms.BOXLENGTH*round(r02[2]*parms.InvBOXL);
      
      r12[0]=local_mols[i].x[1]-local_mols[i].x[2];
      r12[1]=local_mols[i].y[1]-local_mols[i].y[2];
      r12[2]=local_mols[i].z[1]-local_mols[i].z[2];
      r12[0]-=parms.BOXLENGTH*round(r12[0]*parms.InvBOXL);
      r12[1]-=parms.BOXLENGTH*round(r12[1]*parms.InvBOXL);
      r12[2]-=parms.BOXLENGTH*round(r12[2]*parms.InvBOXL);

      v01[0]=vel[i].x[0]-vel[i].x[1];
      v01[1]=vel[i].y[0]-vel[i].y[1];
      v01[2]=vel[i].z[0]-vel[i].z[1];
      v02[0]=vel[i].x[0]-vel[i].x[2];
      v02[1]=vel[i].y[0]-vel[i].y[2];
      v02[2]=vel[i].z[0]-vel[i].z[2];
      v12[0]=vel[i].x[1]-vel[i].x[2];
      v12[1]=vel[i].y[1]-vel[i].y[2];
      v12[2]=vel[i].z[1]-vel[i].z[2];
      qr11=-(v01[0]*r01[0]+v01[1]*r01[1]+v01[2]*r01[2]);//dotProd(v01, r01);
      qr21=-(v02[0]*r02[0]+v02[1]*r02[1]+v02[2]*r02[2]);//dotProd(v02, r02);
      qr31=-(v12[0]*r12[0]+v12[1]*r12[1]+v12[2]*r12[2]);//dotProd(v12, r12);
      
      //calculate matrix elements for rattle
      r01r01=r01[0]*r01[0]+r01[1]*r01[1]+r01[2]*r01[2];//dotProd(r01, r01);
      r12r12=r12[0]*r12[0]+r12[1]*r12[1]+r12[2]*r12[2];//dotProd(r12, r12);
      r02r02=r02[0]*r02[0]+r02[1]*r02[1]+r02[2]*r02[2];//dotProd(r02, r02);
      r01r02=r01[0]*r02[0]+r01[1]*r02[1]+r01[2]*r02[2];//dotProd(r01, r02);
      r01r12=r01[0]*r12[0]+r01[1]*r12[1]+r01[2]*r12[2];//dotProd(r01, r12);
      r02r12=r02[0]*r12[0]+r02[1]*r12[1]+r02[2]*r12[2];//dotProd(r02, r12);
      
      Mat[0]=r01r01*(parms.InvMass[0]+parms.InvMass[1]);
      Mat[1]=r01r02*parms.InvMass[0];
      Mat[2]=-r01r12*parms.InvMass[1];
      Mat[3]=Mat[1];
      Mat[4]=r02r02*(parms.InvMass[0]+parms.InvMass[2]);
      Mat[5]=r02r12*parms.InvMass[2];
      Mat[6]=Mat[2];
      Mat[7]=Mat[5];
      Mat[8]=r12r12*(parms.InvMass[2]+parms.InvMass[1]);
      /*
      cout << Mat[0][0]<<' '<<Mat[0][1]<<' '<<Mat[0][2]<<endl;
      cout << Mat[1][0]<<' '<<Mat[1][1]<<' '<<Mat[1][2]<<endl;
      cout << Mat[2][0]<<' '<<Mat[2][1]<<' '<<Mat[2][2]<<endl;
      */
      //invert the matrix
      MatrixInv(Mat, Minv);
      /*      cout << Minv[0][0]<<' '<<Minv[0][1]<<' '<<Minv[0][2]<<endl;
      cout << Minv[1][0]<<' '<<Minv[1][1]<<' '<<Minv[1][2]<<endl;
      cout << Minv[2][0]<<' '<<Minv[2][1]<<' '<<Minv[2][2]<<endl;
      */
      //get velocity lagrange multipliers
      gv01=Minv[0]*qr11+Minv[1]*qr21+Minv[2]*qr31;
      gv02=Minv[3]*qr11+Minv[4]*qr21+Minv[5]*qr31;
      gv12=Minv[6]*qr11+Minv[7]*qr21+Minv[8]*qr31;

      //Evaluate the constraint velocity
      vel[i].x[0]+=parms.InvMass[0]*(gv01*r01[0]+gv02*r02[0]);
      vel[i].y[0]+=parms.InvMass[0]*(gv01*r01[1]+gv02*r02[1]);
      vel[i].z[0]+=parms.InvMass[0]*(gv01*r01[2]+gv02*r02[2]);
      vel[i].x[1]+=parms.InvMass[1]*(-gv01*r01[0]+gv12*r12[0]);
      vel[i].y[1]+=parms.InvMass[1]*(-gv01*r01[1]+gv12*r12[1]);
      vel[i].z[1]+=parms.InvMass[1]*(-gv01*r01[2]+gv12*r12[2]);
      vel[i].x[2]+=parms.InvMass[2]*(-gv02*r02[0]-gv12*r12[0]);
      vel[i].y[2]+=parms.InvMass[2]*(-gv02*r02[1]-gv12*r12[1]);
      vel[i].z[2]+=parms.InvMass[2]*(-gv02*r02[2]-gv12*r12[2]);
      


    }//end of loop over all molecules

}//end of RattleVel
void checkIdeal(char *filename, Constants cons, Parameters parms, WaterSend *mols)
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
      
      r01sq=dotProd(r01, r01);
      magr01=sqrt(r01sq);

      /*distance O to H2*/
      r02[0]=mols[i].x[0]-mols[i].x[2];
      r02[1]=mols[i].y[0]-mols[i].y[2];
      r02[2]=mols[i].z[0]-mols[i].z[2];

      r02sq=dotProd(r02, r02);
      magr02=sqrt(r02sq);

      /*ditsnce H to H*/
      r12[0]=mols[i].x[1]-mols[i].x[2];
      r12[1]=mols[i].y[1]-mols[i].y[2];
      r12[2]=mols[i].z[1]-mols[i].z[2];

      r12sq=dotProd(r12, r12);
      magr12=sqrt(r12sq);

      /*Angle between vectors OH1 and OH2*/
      
      DP=dotProd(r01, r02);
      cosThetaRad=DP/magr01/magr02;
      ThetaRad=acos(cosThetaRad);
      Thetad=ThetaRad*(cons.PiDEG/cons.Pi);
      
      outfile <<Thetad<<" "<<magr12<<" "<<magr01<<" "<<magr02<<endl;


    }

  delete[] r01;
  delete[] r02;
  delete[] r12;

}

/*Function to remove any component of velocity along the bond, in essence
making the velocity normal to the bond in order to satisfy the constraints of 
the rigid molecule.*/

void VelocityNormal(Constants cons, Parameters parms, WaterSend *mols, DiffSend *InitVel, double *Mat, double *Minv)
{
  int i;
  double r01[3];
  double r02[3];
  double r12[3];
  double v01[3];
  double v02[3];
  double v12[3];
  double m1, m2, m3, L1, L2, L3, c1, c2, c3;
  double r01r01, r01r02, r01r12, r02r02, r02r12, r12r12;

  /*Subtract from velocities any components normal to the constraint surface 
    allows constraint time derivatives to vanish*/

  for(i=0;i<cons.N;i++)
    {
      r01[0]=mols[i].x[0]-mols[i].x[1];
      r01[1]=mols[i].y[0]-mols[i].y[1];
      r01[2]=mols[i].z[0]-mols[i].z[1];
      r01[0]-=parms.BOXLENGTH*round(r01[0]*parms.InvBOXL);
      r01[1]-=parms.BOXLENGTH*round(r01[1]*parms.InvBOXL);
      r01[2]-=parms.BOXLENGTH*round(r01[2]*parms.InvBOXL);
      
      r02[0]=mols[i].x[0]-mols[i].x[2];
      r02[1]=mols[i].y[0]-mols[i].y[2];
      r02[2]=mols[i].z[0]-mols[i].z[2];
      r02[0]-=parms.BOXLENGTH*round(r02[0]*parms.InvBOXL);
      r02[1]-=parms.BOXLENGTH*round(r02[1]*parms.InvBOXL);
      r02[2]-=parms.BOXLENGTH*round(r02[2]*parms.InvBOXL);
      
      r12[0]=mols[i].x[1]-mols[i].x[2];
      r12[1]=mols[i].y[1]-mols[i].y[2];
      r12[2]=mols[i].z[1]-mols[i].z[2];
      r12[0]-=parms.BOXLENGTH*round(r12[0]*parms.InvBOXL);
      r12[1]-=parms.BOXLENGTH*round(r12[1]*parms.InvBOXL);
      r12[2]-=parms.BOXLENGTH*round(r12[2]*parms.InvBOXL);

      v01[0]=InitVel[i].x[0]-InitVel[i].x[1];
      v01[1]=InitVel[i].y[0]-InitVel[i].y[1];
      v01[2]=InitVel[i].z[0]-InitVel[i].z[1];
      v02[0]=InitVel[i].x[0]-InitVel[i].x[2];
      v02[1]=InitVel[i].y[0]-InitVel[i].y[2];
      v02[2]=InitVel[i].z[0]-InitVel[i].z[2];
      v12[0]=InitVel[i].x[1]-InitVel[i].x[2];
      v12[1]=InitVel[i].y[1]-InitVel[i].y[2];
      v12[2]=InitVel[i].z[1]-InitVel[i].z[2];

      /*Calculate matrix elements*/
      /*
      r01r01=dotProd(r01, r01);
      r12r12=dotProd(r12, r12);
      r02r02=dotProd(r02, r02);
      r01r02=dotProd(r01, r02);
      r01r12=dotProd(r01, r12);
      r02r12=dotProd(r02, r12);
      */
      r01r01=r01[0]*r01[0]+r01[1]*r01[1]+r01[2]*r01[2];
      r12r12=r12[0]*r12[0]+r12[1]*r12[1]+r12[2]*r12[2];
      r02r02=r02[0]*r02[0]+r02[1]*r02[1]+r02[2]*r02[2];
      r01r02=r01[0]*r02[0]+r01[1]*r02[1]+r01[2]*r02[2];
      r01r12=r01[0]*r12[0]+r01[1]*r12[1]+r01[2]*r12[2];
      r02r12=r02[0]*r12[0]+r02[1]*r12[1]+r02[2]*r12[2];
      
            
      Mat[0]=r01r01*4.0*(parms.InvMass[0]+parms.InvMass[1]);
      Mat[1]=r01r02*4.0*parms.InvMass[0];
      Mat[2]=-r01r12*4.0*parms.InvMass[1];
      Mat[3]=Mat[1];
      Mat[4]=r02r02*4.0*(parms.InvMass[0]+parms.InvMass[2]);
      Mat[5]=r02r12*4.0*parms.InvMass[2];
      Mat[6]=Mat[2];
      Mat[7]=Mat[5];
      Mat[8]=r12r12*4.0*(parms.InvMass[2]+parms.InvMass[1]);
      c1=2.0*(r01[0]*v01[0]+r01[1]*v01[1]+r01[2]*v01[2]);
      c2=2.0*(r02[0]*v02[0]+r02[1]*v02[1]+r02[2]*v02[2]);
      c3=2.0*(r12[0]*v12[0]+r12[1]*v12[1]+r12[2]*v12[2]);
      
      MatrixInv(Mat, Minv);

      /*Coefficients of the velocities*/
      L1=Minv[0]*c1+Minv[1]*c2+Minv[2]*c3;
      L2=Minv[3]*c1+Minv[4]*c2+Minv[5]*c3;
      L3=Minv[6]*c1+Minv[7]*c2+Minv[8]*c3;

      /*Calculate constrained velocities*/
      m1=2.0*parms.InvMass[0];
      m2=2.0*parms.InvMass[1];
      m3=2.0*parms.InvMass[2];
      InitVel[i].x[0]-=m1*L1*r01[0]+m1*L2*r02[0];
      InitVel[i].y[0]-=m1*L1*r01[1]+m1*L2*r02[1];
      InitVel[i].z[0]-=m1*L1*r01[2]+m1*L2*r02[2];
      InitVel[i].x[1]+=m2*L1*r01[0]-m2*L3*r12[0];
      InitVel[i].y[1]+=m2*L1*r01[1]-m2*L3*r12[1];
      InitVel[i].z[1]+=m2*L1*r01[2]-m2*L3*r12[2];
      InitVel[i].x[2]+=m3*L2*r02[0]+m3*L3*r12[0];
      InitVel[i].y[2]+=m3*L2*r02[1]+m3*L3*r12[1];
      InitVel[i].z[2]+=m3*L2*r02[2]+m3*L3*r12[2];
      
    }

}/*End of VelocityNormal*/

void LongRangeCorrect(Constants cons, Parameters parms, Ensemble &nnn)
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
  nnn.ELRC=ELJSWF+ELJTAIL;

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

  nnn.PLRC=PLJSWF+PLJTAIL;


}//end of LongRangeCorrect


void TableSet(Table *terf, Constants cons, Parameters &parms)
{
  int i;
  double fct=2.0/sqrt(cons.Pi);
  double A2=cons.ALPHA*cons.ALPHA;
  double r, r2, r3;
  double aerc, axp;
  parms.rincr=parms.BOXLENGTH/500000.0;
  parms.irincr=1.0/parms.rincr;

  for(i=1;i<500001;i++)
    {
      r=parms.rincr*i*1.0;
      r2=r*r;
      r3=r2*r;
      aerc=erfc(cons.ALPHA*r);
      axp=exp(-A2*r2);
      terf[i].a=aerc/r;
      terf[i].b=-(fct*axp*cons.ALPHA+terf[i].a)/r;
      terf[i].c=fct*cons.ALPHA*axp*2.0*(1.0/r2 +A2)+2.0*terf[i].a/r2;
      terf[i].d=-fct*cons.ALPHA*axp*2.0*(3.0/r3+(1.0/r2+A2)*(A2*r))-6.0*aerc/(r2*r2);
    }

}//end of TableSet
