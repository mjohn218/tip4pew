/*Program to compute radial distribution file for Oxygen Oxygen
  in a water simulation
  Maggie Johnson
  6/2005
*/

#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#define Np 1728

using namespace std;

class Coord
{
public:
  double x[Np];
  double y[Np];
  double z[Np];

};

void calcRDF(Coord c, double *gr, double boxl, int N, double delr);
void readConfig(ifstream &cfile, Coord &c, int N);

int main(int argc, char *argv[])
{
  /*Read in coordinates of oxygens for each time point
   Average over ~10,000 configurations
  */
  int n, NumFiles, k;
  int i; //num particles
  //N=1728;
  Coord c;
  int numperfile;
  ifstream cfile(argv[1]);
  double boxl, rho;
  double  Temp; //for naming output file
  double r, vb;
  int tag;
  int start;
  cout << "enter rho (g/cc) : "<<endl;
  cin >> rho;
  cout << "enter temperature : "<<endl;
  cin >> Temp;
  cout << "enter any integer to tag output file : "<<endl;
  cin >> tag;
  
  cout << "NumFiles : "<<endl;
  cin >> NumFiles;

  cout <<"Num configs per file:"<<endl;
  cin >> numperfile;

  cout <<"start with file save_pbc_pos.? enter int: "<<endl;
  cin >> start;

  int NumConfigs;
  NumConfigs=numperfile*NumFiles;
  cout << "Averaging over: "<<NumConfigs <<" configurations"<<endl;
  
  double ndensity;
  ndensity=rho*.602214199/18.01528;
  boxl=pow(Np*1.0/ndensity, 1.0/3.0);
  
  cout<< "ndensity "<< ndensity<<endl;
  cout <<" boxlength "<< boxl << endl;
  double delr=.025;
  int nbins;

  nbins=int(boxl/(2.0*delr));
  
  double * gr=new double[nbins];
  char str[20];
  double dub;
  int num;
  char fname[80];
  int pt=start;
  ofstream file1("pos_1.out");
  for(k=0;k<1;k++)
    {
      
      for(n=0;n<1;n++)
	{
	  
	  //cfile >> str >>str>> num >> str >>dub;
	  for(i=0;i<Np;i++)
	    {
	      cfile >> c.x[i] >>c.y[i]>>c.z[i];
	      cfile >> dub >>dub >>dub;
	      cfile >> dub >>dub >>dub;
	      file1<<c.x[i]<<' '<<c.y[i]<<' '<<c.z[i]<<endl;
	    }
	  
	  calcRDF(c, gr, boxl, Np, delr);
	  //cout << c.x[0] <<endl;
	}
      pt+=1;
      cfile.close();
    }
  /*  cout << "x[0]: "<<c.x[0]<<endl;
  cout << "z[511]: "<<c.z[511]<<endl;
  cout << "x[0]: "<<c.x[0]<<endl;
  cout << "z[511]: "<<c.z[511]<<endl;
  */  


  /*Average rdf over all Numconfig configurations*/



}/*End main*/

void readConfig(ifstream &cfile, Coord &c, int N)
{
  int i;
  cfile.ignore(400,'\n');
  cfile.ignore(400,'\n');
  cfile.ignore(400,'\n');
  cfile.ignore(400,'\n');
  cfile.ignore(400,'\n');
  
  for(i=0;i<N;i++)
    {
      cfile >> c.x[i] >>c.y[i]>>c.z[i];
    }
  cout << "In read: "<<c.x[0]<<endl;
}

void calcRDF(Coord c, double *gr, double boxl, int N, double delr )
{
  int i, j;
  double rx, ry, rz;
  double r2, r;
  int ig;
  double iboxl=1.0/boxl;
  for(i=0;i<N-1;i++)
    {
      for(j=i+1;j<N;j++)
	{
	  rx=c.x[i]-c.x[j];
	  ry=c.y[i]-c.y[j];
	  rz=c.z[i]-c.z[j];
	  rx-=boxl*round(rx*iboxl);
	  ry-=boxl*round(ry*iboxl);
	  rz-=boxl*round(rz*iboxl);
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
