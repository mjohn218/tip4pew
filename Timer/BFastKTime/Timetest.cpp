#include <iostream.h>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <math.h>
#include "md_timer.h"

struct MD_Timer sintime;
struct MD_Timer erftime;
struct MD_Timer multtime;
struct MD_Timer exptime;


int main()
{
  double t;

  double start, end;
  cout.precision(16);

  int i;
  int num;
  cout << " num: "<<endl;
  cin >> num;
  double * y=new double[num];
  double * ye=new double[num];

  //  start=clock();
  initialize_timer(&sintime);
  initialize_timer(&erftime);
  initialize_timer(&multtime);
  initialize_timer(&exptime);
  
  start_timer(&sintime);
  for(i=0;i<num;i++)
    {
      t=i*.01;
      y[i]=sin(t);
    }
  stop_timer(&sintime);

  cout << timer_duration(sintime) <<" Time sin" <<endl;
  
  start_timer(&erftime);
  for(i=0;i<num;i++)
    {
      t=i*.01;
      ye[i]=erf(t);
    }
  stop_timer(&erftime);
  cout << timer_duration(erftime) <<" Time erf" <<endl;
  
  start_timer(&multtime);
  for(i=0;i<num;i++)
    {
      t=i*.01;
      ye[i]=t*t+t*t+2*t*t;
    }
  stop_timer(&multtime);
  cout << timer_duration(multtime) <<" Time multiply" <<endl;

  start_timer(&exptime);
  for(i=0;i<num;i++)
    {
      t=i*.01;
      ye[i]=exp(t);
    }
  stop_timer(&exptime);
  cout << timer_duration(exptime) <<" Time exp" <<endl;

}
