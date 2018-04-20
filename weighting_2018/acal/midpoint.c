/* calculate the mid-point of the minor arc between two points */

#include <math.h>

#define TINY 0.0001
#define D double

void midpoint(t1,p1,t2,p2,t,p)
float t1,p1,t2,p2,*t,*p;
{

double x,y,z,r,arg,sin_theta;

  x=sin((D)t1)*cos((D)p1)+sin((D)t2)*cos((D)p2);
  y=sin((D)t1)*sin((D)p1)+sin((D)t2)*sin((D)p2);
  z=cos((D)t1)+cos((D)t2);
  
  r=sqrt(x*x+y*y+z*z);
  x/=r;
  y/=r;
  z/=r;
  
  *t=arg=acos(z);
  sin_theta=sin(arg);
  
  if(sin_theta<TINY)
  { *p=0.5*(p1+p2);
     return;
  }
   
  *p=acos(x/sin_theta);
  if(y<0.0) *p=2.*3.14159265-(*p);
   
}
  
