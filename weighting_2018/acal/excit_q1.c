/* calculates Green's functions of source-receiver functions R,S,D*/
/* S-T could not be considered */

#include <math.h>
#include "acal/structprem.h"
#include "acal/structwavket_A3.h"

extern void stop();

#define FPI 12.566371
#define CONST 9.52278e-26 /* rn^{-4}*wn^{-2}*rhobar^{-1} */
#define SQRTCONST 3.08590e-13  /* sqrt(CONST) */

void excit_q(r,eig,trace1,R,S,dSdr,D,flag)
tracehdrH_st *trace1;
sprem_st *eig;
float r; /* radius of the source location -- normalized by a */
float *R;  /* R[3]  */
float *S; /* S[4]: M_RR, M_+, (M_RL,M_RT), (M_-,M_LT) */
float *dSdr;
float *D; /* D[0]: D0--M_RL
           * D[1]: D1--M_- 
           * D[2]: D1--M_RR
           * D[3]: D1--M_+  
           * D[4]: D2--M_RL 
           * D[5]: D3--M_-  */
int flag; /* 0--> no dSdr and D calculated */
{
static int first=1;
static float h,hn,hn2,hn3;
float a,b,v,vp,u,up,x,L,L1,k0,k1,k2,k3,vpp,upp,xp,f,fp,w;
int l;

  if(first)
  { if(r<0.0 || r>=1.0 ) stop("excit_q: depth of source outoff range");
    h=r-eig->r1;
    hn=eig->r2-eig->r1;
    hn2=1./hn/hn;
    hn3=hn2/hn;
/*    first =0;*/
  }
  
  l=eig->l;
  L=l*(l+1);

  a=hn3*( hn*(eig->v1p+eig->v2p)+2.*(eig->v1-eig->v2) );
  b=hn2*(-hn*(eig->v2p+2.*eig->v1p)+3.*(eig->v2-eig->v1) );
  v=eig->v1 + h*(eig->v1p+ h*(b+h*a) );
  vp=eig->v1p+ h*(2.*b+3.*h*a);
  vpp=2.*b+6.*h*a;
  x=vp-v/r;
  xp=vpp-vp/r+v/r/r;
  
  if(eig->typ=='S')
  { a=hn3*( hn*(eig->u1p+eig->u2p)+2.*(eig->u1-eig->u2) );
    b=hn2*(-hn*(eig->u2p+2.*eig->u1p)+3.*(eig->u2-eig->u1) );
    u=eig->u1 + h*(eig->u1p+ h*(b+h*a) );
    up=eig->u1p + h*(2.*b+3.*h*a);
    upp=2.*b+6.*h*a;
    x+=u/r;
    xp+=up/r-u/r/r;
    f=(2.*u-L*v)/r;
    fp=(2.*up-L*vp)/r-(2.*u-L*v)/r/r;
  }
  
  k0=SQRTCONST*(float)sqrt( (double)(2*l+1)/FPI );
  k1=0.5*k0*(float)sqrt( (double)L );
  k2=0; if(l>0) k2=0.5*k1*(float)sqrt( (double)((l-1)*(l+2)) );
  
  w=eig->w;

 /* modification BR 10/15/96  */

    R[0]=0.0;
    R[1]=0.0;R[2]=0.0;

  if(eig->typ=='S'&& trace1->comp=='Z')
   R[0]=-k0*eig->va/w;      /*va=-U*w*w; <w*U,w*U>=1 */

 if(trace1->comp=='L') {  R[1]=-2.*k1*eig->ha/w; /* different sign
			from transverse because vt=-vtheta 10/15/97 br*/
			  
			 }
  
 if(trace1->comp=='T') { R[1]=2.*k1*eig->ha/w;
			 }


/* end modification */

  S[2]=2.*k1*x*w;              /* <w*Z,w*Z>=1 */

  S[3]=4.*k2*v*w/r;

 
  if(eig->typ=='S')
  { S[0]=k0*up*w;
    S[1]=k0*f*w;
  }
  else S[0]=S[1]=0.0;
  
  if(!flag) return;
  
  k3=0; if(l>1) k3=0.5*k2*(float)sqrt( (double)((l-2)*(l+3)) );
  
  dSdr[2]=2.*k1*xp*w;
  dSdr[3]=4.*k2*w*(vp-v/r)/r;
  if(eig->typ=='S')
  { dSdr[0]=k0*upp*w;
    dSdr[1]=k0*fp*w;
  }
  else dSdr[0]=dSdr[1]=0.0;
  
  L1=0.0; if(l>1) L1=0.5*(float)((l+2)*(l-1));
  if(eig->typ=='S')
  { D[0]=k0*(-0.5*L*x-f+2.*up)*w/r;
    D[1]=2.*k1*(-L1*v/r+x)*w/r;
    D[2]=2.*k1*(-x+up)*w/r;
    D[3]=2.*k1*(x+f)*w/r;
    D[4]=2.*k2*(x-2.*v/r)/r;
    D[5]=4.*k3*v/r/r;
  }
  else
  { D[0]=k0*0.5*L*x*w/r;       /* - */
    D[1]=2.*k1*(L1*v/r-x)*w/r; /* - */
    D[2]=2.*k1*(-x)*w/r;       /* + */
    D[3]=2.*k1*(x)*w/r;        /* + */
    D[4]=2.*k2*(x-2.*v/r)/r;   /* + */
    D[5]=4.*k3*v/r/r;          /* + */
  }
  
}           
