#include <math.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "acal/const.h"

extern void stop();

float sphdist();
float ctbtaper();

float pathprojectn1(dep1,ths1,phs1,thr1,phr1,
                   dep2,ths2,phs2,thr2,phr2,depres,disresr,disress)
float dep1,ths1,phs1,thr1,phr1,dep2,ths2,phs2,thr2,phr2,depres,disresr,disress;
{
float prj,ddep,contrib,dist;

  prj=0.0;
  
  ddep=fabs((D)(dep1-dep2));
  if (ddep<depres) {
    contrib=ctbtaper(ddep,depres);
    dist=sphdist(thr1,phr1,thr2,phr2);
    if (dist<disresr) {
      contrib*=ctbtaper(dist,disresr);
      dist=sphdist(ths1,phs1,ths2,phs2);
      if (dist<disress) prj=contrib*ctbtaper(dist,disress);
    }
  }
  return(prj);
}

float ctbtaper(x,res)
float x,res;
{ 
float y;
double arg;

  if(x<0.0 || x>res || res<=0.0) stop("ctbtaper: error wrong argument");
  arg=PI*x/res;
  y=0.5+0.5*cos(arg);
  return(y);
}
