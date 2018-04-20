/* subs for streamlined acal program 

   excit_q(r,&prem,trace,R[i],S[i],dSdr,Di,0) is called only
   in the first run, in which comp info in trace is used, 
    --- > mixed comp data in the same wp files is not allowed

  ADD 

  #define EIG  X for XDeigenfuncion; Y for yannos eigenfunction
  in EIG_XY.h to switch for different eigenfunction 

  use the QC80 kernel, setmode is modified    Dec 1/2003 gung


  from subs_q36fbry_auto_test702.c 
  remove the skip in 1S0, since there is no 1S0 in yannos
  sum_modes in setmodes_A3 to be replaced by TMODE and SMODE 

  version with Q1  BR 11/09/96 
  
  patch for autopicking bug added (Federica Oct 26, 2006)*/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include "endian.h"
#include "acal/dimensiony.h"
#include "acal/EIG_XY.h"

#include "acal/const.h"
#include "acal/structwavket_A3.h"
#include "acal/structbodinv_A3_aniso.h"
#include "acal/structknl.h"

#include "acal/structsynth.h"
#include "acal/structprem.h"

#define SRH 0.7071068
#define NWP 360
#define CONST 9.52278e-26 /* rn^{-4}*wn^{-2}*rhobar^{-1} */
#define VERBOSE 1
#define EPS 1.0e-5

#ifdef __hpux
#define LEGENDRE(x1,x2,x3,x4,x5,x6) legndr(x1,x2,x3,x4,x5,x6)
#else
#define LEGENDRE(x1,x2,x3,x4,x5,x6) legndr_(x1,x2,x3,x4,x5,x6)
#endif

FILE *prop_file;	/* has proportionality coefficient : dlnVs/dlnVp */

extern void legndr_();
extern int fgetrecord();
extern void excit_q();
extern void Euler();
extern void rotmatrix();
extern void getprem();
extern int geteigys_ellrot();
extern float ***farray3();
extern float **farray2();
extern float *free_farray2();
extern float *farray1();
extern float *free_farray1();
extern char *carray1();
extern int *iarray1();
extern int **iarray2();
extern void stop();
void loadmdlpar_A3();
void loadmdlpar_A3_anis();
datparbr_st *getdatpar_auto();
void seta();
int dataselector_auto();
int dataselector_auto_reg();
int dataselector_auto_reg_sem();

wpinfoved_st *loadtrace_auto();
wpinfoved_st *loadtrace_auto_reg();
wpinfoved_st *loadtrace_auto_all();
wpinfoved_st *loadtrace_auto_syn();
wpinfoved_st *loadtrace_auto_sem_reg();
int **setmodes_A3_anis();
source_st *transource();
void locaknl();
window_st *setwindow();
int slkrcd();
int whichorbit();
void collapsX20_q();
void collapsX20();
void rotcoef_q();
void Gelement_anis();
void a_term_anis();
float **garray();
void rotcoef();


/**********************************************/
void loadmdlpar_A3(fpname,par)
mdlparA3_st *par;
char *fpname;
{
FILE *fp;
char string[80];
int d1,d2;
float f;

  if((fp=fopen(fpname,"r"))==0)
    stop("loadmdlpar_A3: cannot open model file");
  fgetrecord(fp,string); sscanf(string,"%d%d%f",&d1,&d2,&f);
    par->nknotA2=d1; par->nknotA1=d2; par->damp=f;

   fgetrecord(fp,string); sscanf(string,"%d%d%f",&d1,&d2,&f);
    par->nknots=d1; par->seafl=d2; par->damps=f;
  fgetrecord(fp,string); sscanf(string,"%d%d%f",&d1,&d2,&f);
    par->nknotm=d1; par->moho=d2;  par->dampm=f;
  fgetrecord(fp,string); sscanf(string,"%d%d%f",&d1,&d2,&f);
    par->nknotd=d1; par->d670=d2;  par->dampd=f;
  fgetrecord(fp,string); sscanf(string,"%d%d%f",&d1,&d2,&f);
    par->nknotc=d1; par->cmb=d2;   par->dampc=f;
    


    if( (par->seafl!=0&& (par->nknots!=par->nknotA2)) ||
	(par->moho !=0&& (par->nknotm!=par->nknotA2)) ||
	(par->d670 !=0&& (par->nknotd!=par->nknotA2)) ||
	(par->cmb  !=0&& (par->nknotc!=par->nknotA2)) )
     stop("loadmdlpar_A3: error -  nonconformal model not available yet");
							 

  par->dim=(par->nknotA2)*(par->nknotA1)
         +(par->nknots)*par->seafl
         +(par->nknotm)*par->moho
         +(par->nknotd)*par->d670
         +(par->nknotc)*par->cmb;
}

/**********************************************/
void loadmdlpar_A3_anis(fpname,par)
gmdlparA3_st *par;
char *fpname;
{
FILE *fp;
char string[80];
int d1,d2,ii;
float f;

  if ((fp=fopen(fpname,"r"))==0)
    stop("loadmdlpar_A3_anis: cannot open model file");
  fgetrecord(fp,string); sscanf(string,"%d",&d1);
    par->nquant=d1;
    
/* Assumes radial parametrization is the same for all */

  for (ii=0;ii<par->nquant;ii++) {
    fgetrecord(fp,string);
    sscanf(string,"%s",&par->descript[ii]);
    fgetrecord(fp,string); sscanf(string,"%d%d%f",&d1,&d2,&f);
    par->nknotA2[ii]=d1; par->nknotA1[ii]=d2; par->damp[ii]=f;
  }
  for (ii=par->nquant;ii<PARAM;ii++) {
    par->nknotA2[ii]=0; par->nknotA1[ii]=0; par->damp[ii]=0.;
  }

  if (par->descript[0]!='S')
    stop ("loadmdlpar_A3_anis: First parameter must be the isotropic velocity structure");
  for (ii=0;ii<par->nquant;ii++) {
/*    if (par->nknotA1[0]!=par->nknotA1[ii] || par->nknotA2[0]!=par->nknotA2[ii])*/
    if (par->nknotA1[0]!=par->nknotA1[ii])
       stop ("loadmdlpar_A3_anis: Vertical parametrization must be the same for all parameters");
  }
  fgetrecord(fp,string); sscanf(string,"%d%d%f",&d1,&d2,&f);
    par->nknots=d1; par->seafl=d2; par->damps=f;
  fgetrecord(fp,string); sscanf(string,"%d%d%f",&d1,&d2,&f);
    par->nknotm=d1; par->moho=d2; par->dampm=f;
  fgetrecord(fp,string); sscanf(string,"%d%d%f",&d1,&d2,&f);
    par->nknotd=d1; par->d670=d2; par->dampd=f;
  fgetrecord(fp,string); sscanf(string,"%d%d%f",&d1,&d2,&f);
    par->nknotc=d1; par->cmb=d2; par->dampc=f;

  if((par->seafl!=0 && (par->nknots!=par->nknotA2[0])) ||
     (par->moho !=0 && (par->nknotm!=par->nknotA2[0])) ||
     (par->d670 !=0 && (par->nknotd!=par->nknotA2[0])) ||
     (par->cmb  !=0 && (par->nknotc!=par->nknotA2[0])))
     stop("loadmdlpar_A3_anis: Discontinuities must have the same parametrization as the isotropic structure");
  par->ndisc=par->seafl+par->moho+par->d670+par->cmb;
							 
  par->dim=0;
  for (ii=0;ii<par->nquant;ii++) {
    par->dim+=(par->nknotA2[ii])*(par->nknotA1[ii]);
  }
  par->dim+=(par->nknots)*par->seafl
         +(par->nknotm)*par->moho
         +(par->nknotd)*par->d670
         +(par->nknotc)*par->cmb;
}

/**********************************************/
datparbr_st *getdatpar_auto()
{
static datparbr_st datpar;
int i,n;
char string[80];
float delta1,delta2;
FILE *file;

  if((file=fopen("lu_data.par","r"))==0)
    stop("getdatpar: cannot open lu_data.par");

  n=-5; /* n indicates the number of lines in lu_data.par before the phase names */
  while(fgetrecord(file,string)>0)  n++;
  datpar.n=n;
  
  datpar.code=iarray1(0,n-1);
  datpar.delta1=farray1(0,n-1);
  datpar.delta2=farray1(0,n-1);
  
  rewind(file);
  fgetrecord(file,string); sscanf(string,"%f",&datpar.vcutoffr);
  fgetrecord(file,string); sscanf(string,"%f",&datpar.vcutoffs);
  fgetrecord(file,string); sscanf(string,"%s",datpar.netwk);
  fgetrecord(file,string); sscanf(string,"%s",datpar.chnnl);
  fgetrecord(file,string); sscanf(string,"%s",datpar.comp);
	
  for(i=0;i<n;i++)
  { if(fgetrecord(file,string)<=0) stop("getdatpar: error 1");
    if(string[0]=='*') 
    { datpar.n=-1;
      fclose(file); 
      return(&datpar);
    }
    sscanf(string,"%d %f %f", &(datpar.code[i]),&delta1,&delta2);
 datpar.delta1[i]=delta1*CONV;
    datpar.delta2[i]=delta2*CONV;
  }
  fclose(file);
  return(&datpar);
}

/*************/

void seta(smax,coef,ptl,A)
float *A,*ptl;
int smax;
float **coef[];
{
int s,t0,t1,t2;
float add;

  for(s=0;s<=smax;s++)
  { if(s%2) t0=1; else t0=3;
    for(t1=0;t1<=2*s;t1++)
    { if(s%2) add=0.0;
      else add=coef[s][0][t1]*ptl[0];
      for(t2=t0;t2<2*s;t2+=4)
        add+=coef[s][t2][t1]*ptl[t2]+coef[s][t2+1][t1]*ptl[t2+1];
      A[s*s+t1]=add;
    }
  }
}
/************************************/
int dataselector_auto(packhdr,d,netwk,chnnl,comp)
/* compare data to lu_data.par 
   if within limits return 1, else return 0 */

packhdrbr_st *packhdr;
float d; /*delta*/
char *netwk,*chnnl;
char comp;
{
int i;
static int first=1;
static datparbr_st *dp;
  
 if(first)
   {
     first=0; 
     dp=getdatpar_auto();
   }
 
 if(packhdr->rmsr>dp->vcutoffr||packhdr->rmss>dp->vcutoffs) {    
   if(packhdr->rmsr==999)  printf("RMSR = 999 Wavepacket Rejected \n"); 
    return(0); 
 }
  
  if(dp->n==-1) return(1);
  if(netwk[0]!='*' && dp->netwk[0]!='*')
    for(i=0;i<4;i++) 
      if(netwk[i]!=dp->netwk[i] && netwk[i])   
        return(0);    
  
  if(chnnl[0]!='*' && dp->chnnl[0]!='*')
    for(i=0;i<4;i++) if(chnnl[i]!=dp->chnnl[i]) return(0);
   
  if(dp->comp[0]!='*')
    if(comp!=dp->comp[0]) {
      return(0);
    }

  if((dp->code)[0]=='*') return(1);
  for(i=0;i<dp->n;i++){
    if(packhdr->phase==(short)((dp->code)[i]) && d>=(dp->delta1)[i] && 
       d<=(dp->delta2)[i]) return(1);
  }
  
  return(0);
}
/************************************/
int dataselector_auto_reg(packhdr,d,netwk,chnnl,comp)
/* compare data to lu_data.par and check the weight
   if within limits return 1, else return 0 */

packhdrbr_st *packhdr;
float d; /*delta*/
char *netwk,*chnnl;
char comp;
{
int i;
static int first=1;
static datparbr_st *dp;
  
 if(first)
   {
     first=0; 
     dp=getdatpar_auto();
   }
 
  if(packhdr->rmsr>dp->vcutoffr||packhdr->rmss>dp->vcutoffs||packhdr->weight==0)     
    return(0); 
  
  if(dp->n==-1) return(1);
  if(netwk[0]!='*' && dp->netwk[0]!='*')
    for(i=0;i<4;i++) 
      if(netwk[i]!=dp->netwk[i] && netwk[i])   
        return(0);    
  
  if(chnnl[0]!='*' && dp->chnnl[0]!='*')
    for(i=0;i<4;i++) if(chnnl[i]!=dp->chnnl[i]) return(0);
   
  if(dp->comp[0]!='*')
    if(comp!=dp->comp[0]) {
      return(0);
    }

  if((dp->code)[0]=='*') return(1);
  for(i=0;i<dp->n;i++){
    if(packhdr->phase==(short)((dp->code)[i]) && d>=(dp->delta1)[i] && 
       d<=(dp->delta2)[i]) return(1);
  }
  
  return(0);
}
/************************************/
/************************************/
int dataselector_auto_reg_sem(packhdr,d,netwk,chnnl,comp)
/* compare data to lu_data.par and check the weight
   if within limits return 1, else return 0 */

packhdrbr_st *packhdr;
float d; /*delta*/
char *netwk,*chnnl;
char comp;
{
int i;
static int first=1;
static datparbr_st *dp;
  
 if(first)
   {
     first=0; 
     dp=getdatpar_auto();
   }
 
  if(packhdr->rmsr>dp->vcutoffr||packhdr->rmss>dp->vcutoffs)     
    return(0); 
  
  if(dp->n==-1) return(1);
  if(netwk[0]!='*' && dp->netwk[0]!='*')
    for(i=0;i<4;i++) 
      if(netwk[i]!=dp->netwk[i] && netwk[i])   
        return(0);    
  
  if(chnnl[0]!='*' && dp->chnnl[0]!='*')
    for(i=0;i<4;i++) if(chnnl[i]!=dp->chnnl[i]) return(0);
   
  if(dp->comp[0]!='*')
    if(comp!=dp->comp[0]) {
      return(0);
    }

  if((dp->code)[0]=='*') return(1);
  for(i=0;i<dp->n;i++){
    if(packhdr->phase==(short)((dp->code)[i]) && d>=(dp->delta1)[i] && 
       d<=(dp->delta2)[i]) return(1);
  }
  
  return(0);
}
/************************************/
wpinfoved_st *loadtrace_auto(wpd,trace,title,nwp,maxorbit)
int wpd,*nwp,*maxorbit;
tracehdrH_st *trace;
title_st *title;
{
  int nread,size,i,orbit;
  packhdrbr_st packhdr;
  static wpinfoved_st wpinfo[NWP];
  static int prnwp=0;
  float gv1,gv2,dist; 
  double dum,diff;

  for(i=0;i<prnwp;i++) free_farray1(wpinfo[i].data,0);
  
  *nwp=0;
  *maxorbit=1;
  do
  { 
    if(!(nread=read(wpd,&packhdr,sizeof(packhdr)))) break;
    if(nread!=sizeof(packhdr)) stop("loadtrace: error 1");
    if(packhdr.id!=trace->id) break;
    size=packhdr.ndata*sizeof(float);
    
    diff=modf((packhdr.t0+title->dt)/trace->smplintv,&dum);
    if(dataselector_auto(&packhdr,trace->delta,trace->netwk,
			 trace->chnnl,trace->comp)
       && (fabs(diff)<EPS||fabs(1.0-diff)<EPS))
      { if(*nwp>NWP) stop("loadtrace:not enough working space");
      wpinfo[*nwp].ndata=packhdr.ndata;
      wpinfo[*nwp].pv1=packhdr.pv1;
      wpinfo[*nwp].pv2=packhdr.pv2;
      wpinfo[*nwp].gv1=packhdr.gv1;
      wpinfo[*nwp].gv2=packhdr.gv2;
      wpinfo[*nwp].t0=packhdr.t0;
      wpinfo[*nwp].rmsd=packhdr.rmsd;
      wpinfo[*nwp].rmsr=packhdr.rmsr;
      wpinfo[*nwp].rmss=packhdr.rmss;
      wpinfo[*nwp].phase=packhdr.phase;
      wpinfo[*nwp].data=farray1(0,packhdr.ndata);
      wpinfo[*nwp].weight=packhdr.weight;
      wpinfo[*nwp].orbit=orbit=whichorbit(packhdr.phase);
      
      /* readjust gv window for synthetic data */
      if(wpinfo[*nwp].gv1==0.0&&wpinfo[*nwp].gv2==1.0) {	
	if(orbit==1||orbit==3||orbit==5) 
	  dist=((orbit-1)*PI+trace->delta);
	else if (orbit==2||orbit==4||orbit==6) 
	  dist=(orbit*PI-trace->delta);	      
	
	wpinfo[*nwp].gv1=dist/(packhdr.t0+trace->smplintv*(packhdr.ndata-1));
	wpinfo[*nwp].gv2=dist/packhdr.t0;	
      }
      
      if(orbit>*maxorbit) *maxorbit=orbit;
      
      if(read(wpd,wpinfo[*nwp].data,size)!=size)
	stop("loadtrace: error 2");      
      (*nwp)++;
      }
    else {
      lseek(wpd,(long)size,1);
    }
    
  } while(packhdr.id==trace->id);
  if(*maxorbit>MAXORBIT) 
    stop("loadtrace: error -- not enough space");
  prnwp=*nwp;
  
  return(wpinfo);
}

/************************************/
wpinfoved_st *loadtrace_auto_reg(wpd,trace,title,nwp,maxorbit)
int wpd,*nwp,*maxorbit;
tracehdrH_st *trace;
title_st *title;
{
  int nread,size,i,orbit;
  packhdrbr_st packhdr;
  static wpinfoved_st wpinfo[NWP];
  static int prnwp=0;
  float gv1,gv2,dist; 
  double dum,diff;

  for(i=0;i<prnwp;i++) free_farray1(wpinfo[i].data,0);
  
  *nwp=0;
  *maxorbit=1;
  do
  { 
    if(!(nread=read(wpd,&packhdr,sizeof(packhdr)))) break;
    if(nread!=sizeof(packhdr)) stop("loadtrace: error 1");
    if(packhdr.id!=trace->id) break;
    size=packhdr.ndata*sizeof(float);
    
    diff=modf((packhdr.t0+title->dt)/trace->smplintv,&dum);
    if(dataselector_auto_reg(&packhdr,trace->delta,trace->netwk,
			 trace->chnnl,trace->comp)
       && (fabs(diff)<EPS||fabs(1.0-diff)<EPS))
      { if(*nwp>NWP) stop("loadtrace:not enough working space");
      wpinfo[*nwp].ndata=packhdr.ndata;
      wpinfo[*nwp].pv1=packhdr.pv1;
      wpinfo[*nwp].pv2=packhdr.pv2;
      wpinfo[*nwp].gv1=packhdr.gv1;
      wpinfo[*nwp].gv2=packhdr.gv2;
      wpinfo[*nwp].t0=packhdr.t0;
      wpinfo[*nwp].rmsd=packhdr.rmsd;
      wpinfo[*nwp].rmsr=packhdr.rmsr;
      wpinfo[*nwp].rmss=packhdr.rmss;
      wpinfo[*nwp].phase=packhdr.phase;
      wpinfo[*nwp].data=farray1(0,packhdr.ndata);
      wpinfo[*nwp].weight=packhdr.weight;
      wpinfo[*nwp].orbit=orbit=whichorbit(packhdr.phase);
      
      /* readjust gv window for synthetic data */
      if(wpinfo[*nwp].gv1==0.0&&wpinfo[*nwp].gv2==1.0) {	
	if(orbit==1||orbit==3||orbit==5) 
	  dist=((orbit-1)*PI+trace->delta);
	else if (orbit==2||orbit==4||orbit==6) 
	  dist=(orbit*PI-trace->delta);	      
	
	wpinfo[*nwp].gv1=dist/(packhdr.t0+trace->smplintv*(packhdr.ndata-1));
	wpinfo[*nwp].gv2=dist/packhdr.t0;	
      }
      
  
      if(read(wpd,wpinfo[*nwp].data,size)!=size)
	stop("loadtrace: error 2");      
      (*nwp)++;
      }
    else {
      lseek(wpd,(long)size,1);
    }
    
  } while(packhdr.id==trace->id);
  if(*maxorbit>MAXORBIT) 
    stop("loadtrace: error -- not enough space");
  prnwp=*nwp;
  
  return(wpinfo);
}

/************************************/
/************************************/
wpinfoved_st *loadtrace_auto_sem_reg(wpd,smd,trace,title,nwp,maxorbit)
int wpd,smd,*nwp,*maxorbit;
tracehdrH_st *trace;
title_st *title;

{  int nread,size,i,orbit;
  packhdrbr_st packhdrwp;
  packhdrbr_st packhdrsm;
  static wpinfoved_st wpinfo[NWP];
  static int prnwp=0;
  float gv1,gv2,dist; 
  double dum,diff1,diff2;

  for(i=0;i<prnwp;i++) free_farray1(wpinfo[i].data,0);
  for(i=0;i<prnwp;i++) free_farray1(wpinfo[i].semdata,0);

  *nwp=0;
  *maxorbit=1;
  do
  { 
    if(!(nread=read(wpd,&packhdrwp,sizeof(packhdrwp)))) break;
    if(!(nread=read(smd,&packhdrsm,sizeof(packhdrsm)))) break;

    if(nread!=sizeof(packhdrwp)) stop("loadtrace: error 1");
    
    if(packhdrwp.id!=trace->id) break;
    if(packhdrsm.id!=trace->id) break;

    size=packhdrwp.ndata*sizeof(float);
       
    diff1=modf((packhdrwp.t0+title->dt)/trace->smplintv,&dum);
    diff2=modf((packhdrsm.t0+title->dt)/trace->smplintv,&dum);

    if(dataselector_auto_reg_sem(&packhdrsm,trace->delta,trace->netwk,
			 trace->chnnl,trace->comp)
       && dataselector_auto_reg(&packhdrwp,trace->delta,trace->netwk,
			 trace->chnnl,trace->comp)
       && (fabs(diff1)<EPS||fabs(1.0-diff1)<EPS)
       && (fabs(diff2)<EPS||fabs(1.0-diff2)<EPS)  )
      { if(*nwp>NWP) stop("loadtrace:not enough working space");

      wpinfo[*nwp].ndata=packhdrwp.ndata;
      wpinfo[*nwp].pv1=packhdrwp.pv1;
      wpinfo[*nwp].pv2=packhdrwp.pv2;
      wpinfo[*nwp].gv1=packhdrwp.gv1;
      wpinfo[*nwp].gv2=packhdrwp.gv2;
      wpinfo[*nwp].t0=packhdrwp.t0;
      wpinfo[*nwp].rmsd=packhdrwp.rmsd;
      wpinfo[*nwp].rmsr=packhdrwp.rmsr;
      wpinfo[*nwp].rmss=packhdrwp.rmss;
      wpinfo[*nwp].phase=packhdrwp.phase;
      wpinfo[*nwp].data=farray1(0,packhdrwp.ndata);
      wpinfo[*nwp].semdata=farray1(0,packhdrsm.ndata);
      wpinfo[*nwp].weight=packhdrwp.weight;
      wpinfo[*nwp].orbit=orbit=whichorbit(packhdrwp.phase);      

      /* readjust gv window for synthetic data */
      
      if(wpinfo[*nwp].gv1==0.0&&wpinfo[*nwp].gv2==1.0) {	
	if(orbit==1||orbit==3||orbit==5) 
	  dist=((orbit-1)*PI+trace->delta);
	else if (orbit==2||orbit==4||orbit==6) 
	  dist=(orbit*PI-trace->delta);	      
	
	wpinfo[*nwp].gv1=dist/(packhdrwp.t0+trace->smplintv*(packhdrwp.ndata-1));
	wpinfo[*nwp].gv2=dist/packhdrwp.t0;	
      }
 
      if(orbit>*maxorbit) *maxorbit=orbit;
      /* printf("Before reading sem file nwp = %d\n",*nwp); */
      if(read(smd,wpinfo[*nwp].semdata,size)!=size) {	stop("loadtrace: error semsyn 2"); }
      
      /* printf("Before reading wpd file nwp = %d\n",*nwp); */
      if(read(wpd,wpinfo[*nwp].data,size)!=size) {	stop("loadtrace: error data 2"); }
 
      (*nwp)++;
      }
    else {
      lseek(wpd,(long)size,1);
      lseek(smd,(long)size,1);
    }
    
  } while(packhdrwp.id==trace->id && packhdrsm.id==trace->id);
  if(*maxorbit>MAXORBIT) 
    stop("loadtrace: error -- not enough space");
  prnwp=*nwp;

 return(wpinfo);
  
}

/************************************/


wpinfoved_st *loadtrace_auto_all(wpd,trace,title,nwp,maxorbit)
int wpd,*nwp,*maxorbit;
tracehdrH_st *trace;
title_st *title;
{
  int nread,size,i,orbit;
  packhdrbr_st packhdr;
  static wpinfoved_st wpinfo[NWP];
  static int prnwp=0;
  float gv1,gv2,dist; 

  for (i=0;i<prnwp;i++) free_farray1(wpinfo[i].data,0);
  *nwp=0;
  *maxorbit=1;
  do
  { 
    if (!(nread=read(wpd,&packhdr,sizeof(packhdr)))) break;
    if (nread!=sizeof(packhdr)) stop("loadtrace: error 1");
    if (packhdr.id!=trace->id) break;
    size=packhdr.ndata*sizeof(float);
    
    if (*nwp>NWP) stop("loadtrace:not enough working space");
    wpinfo[*nwp].ndata=packhdr.ndata;
    wpinfo[*nwp].pv1=packhdr.pv1;
    wpinfo[*nwp].pv2=packhdr.pv2;
    wpinfo[*nwp].gv1=packhdr.gv1;
    wpinfo[*nwp].gv2=packhdr.gv2;
    wpinfo[*nwp].t0=packhdr.t0;
    wpinfo[*nwp].rmsd=packhdr.rmsd;
    wpinfo[*nwp].rmsr=packhdr.rmsr;
    wpinfo[*nwp].rmss=packhdr.rmss;
    wpinfo[*nwp].phase=packhdr.phase;
    wpinfo[*nwp].data=farray1(0,packhdr.ndata);
    wpinfo[*nwp].weight=packhdr.weight;
    wpinfo[*nwp].orbit=orbit=whichorbit(packhdr.phase);
      
    /* readjust gv window for synthetic data */
    if (wpinfo[*nwp].gv1==0.0&&wpinfo[*nwp].gv2==1.0) {	
      if (orbit==1||orbit==3||orbit==5) 
        dist=((orbit-1)*PI+trace->delta);
      else if (orbit==2||orbit==4||orbit==6) 
        dist=(orbit*PI-trace->delta);	      
	
      wpinfo[*nwp].gv1=dist/(packhdr.t0+trace->smplintv*(packhdr.ndata-1));
      wpinfo[*nwp].gv2=dist/packhdr.t0;	
    }
      
    if (orbit>*maxorbit) *maxorbit=orbit;
      
    if (read(wpd,wpinfo[*nwp].data,size)!=size)
      stop("loadtrace: error 2");      
    (*nwp)++;
    
  } while(packhdr.id==trace->id);
  if (*maxorbit>MAXORBIT) 
    stop("loadtrace: error -- not enough space");
  prnwp=*nwp;
  
  return(wpinfo);
}

/************************************/
wpinfoved_st *loadtrace_auto_syn(wpd,trace,title,nwp,maxorbit)
     /* for synthetic dummp wpd file only , data reading is removed*/
     int wpd,*nwp,*maxorbit;
     tracehdrH_st *trace;
     title_st *title;
{
  int nread,size,i,orbit,ii;
  packhdrbr_st packhdr;
  static wpinfoved_st wpinfo[NWP];
  static int prnwp=0;
  float gv1,gv2,dist; 
  
  for(i=0;i<prnwp;i++) free_farray1(wpinfo[i].data,0);
  *nwp=0;  *maxorbit=1;
  do
    { 
      if(!(nread=read(wpd,&packhdr,sizeof(packhdr)))) break;
      if(nread!=sizeof(packhdr)) stop("loadtrace: error 1");
      if(packhdr.id!=trace->id) break;
      size=packhdr.ndata*sizeof(float);
      if(dataselector_auto(&packhdr,trace->delta,trace->netwk,
			   trace->chnnl,trace->comp))
	{ if(*nwp>NWP) stop("loadtrace: NWP too small");
	wpinfo[*nwp].ndata=packhdr.ndata;
	wpinfo[*nwp].pv1=packhdr.pv1;
	wpinfo[*nwp].pv2=packhdr.pv2;
	wpinfo[*nwp].gv1=packhdr.gv1;
	wpinfo[*nwp].gv2=packhdr.gv2;
	wpinfo[*nwp].t0=packhdr.t0;
	wpinfo[*nwp].rmsd=packhdr.rmsd;
	wpinfo[*nwp].rmsr=packhdr.rmsr;
	wpinfo[*nwp].rmss=packhdr.rmss;
	wpinfo[*nwp].phase=packhdr.phase;
	wpinfo[*nwp].data=farray1(0,packhdr.ndata);
	wpinfo[*nwp].weight=packhdr.weight;
	wpinfo[*nwp].orbit=orbit=whichorbit(packhdr.phase);

	/* readjust gv window for synthetic data */
	if(wpinfo[*nwp].gv1==0.0&&wpinfo[*nwp].gv2==1.0) {	
	  if(orbit==1||orbit==3||orbit==5) 
	    dist=((orbit-1)*PI+trace->delta);
	  else if (orbit==2||orbit==4||orbit==6) 
	    dist=(orbit*PI-trace->delta);	
	  
	  wpinfo[*nwp].gv1=dist/(packhdr.t0+trace->smplintv*(packhdr.ndata-1));
	  wpinfo[*nwp].gv2=dist/packhdr.t0;	
	}
	
	if(orbit>*maxorbit) *maxorbit=orbit;  
	for(ii=0;ii<packhdr.ndata;ii++)wpinfo[*nwp].data[ii]=1.e-7;
	(*nwp)++; 
	}
    } while(packhdr.id==trace->id);
  
  if(*maxorbit>MAXORBIT) stop("loadtrace: error -- not enough space");
  prnwp=*nwp;
  return(wpinfo);
}



/*******************************************************/
int **setmodes_A3_anis(mp,mdl,title,trace,neffk,Aknot,mode,swtch,ellswtch,unconf,nstart_l1)
                       
gmdlparA3_st *mp;
float *mdl;
title_st *title;
tracehdrH_st *trace;
sphA1reg_st *Aknot;
mode_anis_st *mode;
int swtch,ellswtch,*neffk,*nstart_l1,unconf;  /* swtch not used */
{
static int firstT=1,firstS=1,exist;
static int **start;
static float *knl;
static float ***ptl;
static float prop,prop1,prop2,prop3;
static int nmode[2],pmax;
static int ll[2][NMODE],nn[2][NMODE];
float **R;  
float **S;	/* S[4]: M_RR, M_+, (M_RL,M_RT), (M_-,M_LT) */

sort_knlrechdr_st rcdhdr;
source_st *src;         /* source parameters */
sprem_st prem[2][NMODE];

int i,ii,j,p,ob,kid,lmax,ik;
int end,hjump,vjump,begin;
float **f,f1,f2,delta,arg,r;
float *dSdr;
float *Di;		/* dummy fed to excit() */

float **cX20,**pavaX20;
void collapsX20_q();
char mode_type;

  if (trace->comp=='L' || trace->comp=='Z') {           /* if spheroidal */
    mode_type='S';
    ik=0;
  } else if (trace->comp=='T') {                        /* if toroidal */
    mode_type='T';
    ik=1;
  }	
  
  if ((firstT) && (firstS)) {
    pmax=0;
    for (i=0;i<mp->nquant;i++) {
      pmax+=mp->nknotA1[i];
    }
    pmax+=mp->seafl+mp->moho+mp->d670+mp->cmb;
    ptl=farray3(0,1,0,NMODE,0,pmax);
    start=iarray2(0,1,0,999);
  }

  if (((firstT) && mode_type=='T') || ((firstS) && mode_type=='S')) {

    int fd,pp,rcdlen,pdim,sum_modes;
    short kmax,ndisc;
    int firstl1;
    float fac,maxw;
    
    /* prop=s4, prop1=s3, prop2=s2, prop3=s1 */
    if ((prop_file=fopen("lu_prop.par","r"))==0) 
      stop("setmodes_A3_anis: cannot open prop file");
    fscanf(prop_file,"%f %f %f %f",&prop,&prop1,&prop2,&prop3);
    fclose(prop_file);
    prop=1./prop;
    prop1=1./prop1;

    if (trace->comp=='L' || trace->comp=='Z') {   /* if spheroidal */
      sum_modes=SMODE;
      firstS=0;
      if ((fd=open("bu_selfcplS",0))<0) 
        stop("setmodes_A3_anis: cannot open selfcplS.dat");

    } else if (trace->comp=='T') {		 /* if toroidal */
      firstT=0;
      sum_modes=TMODE;
      if ((fd=open("bu_selfcplT",0))<0) 
        stop("setmodes_A3_anis: cannot open selfcplT.dat");

    } else stop("setmodes_A3_anis: non-existent component");
        
    read(fd,&rcdlen,4);
    read(fd,&maxw,4);
    read(fd,&kmax,2);  

    if (kmax!=mp->nknotA1[0]) {
      printf("kmax %d nknotA1 %d\n",kmax,mp->nknotA1[0]); 
      stop("setmodes_A3_anis: dimension error in bu_selfcpl 1");
    }

    read(fd,&ndisc,2);
    read(fd,&lmax,4);
    
    read(fd,start[ik],sizeof(int)*(lmax+1));

    /* limited to ACFLNrho partial files (6 parameters) */
    /* pdim needs to stay like this even for the case of only viso
    because of reading anisotropic kernels */
    
    pdim=6*kmax+ndisc;
    
    if (pdim*sizeof(float)+sizeof(rcdhdr)!=rcdlen) 
      stop("setmodes_A3_anis: dimension error in bu_selfcpl 2");
    knl=farray1(0,pdim-1);
        
    i=0;               /*  In original anisotropic code i=i0[ik]; */
    firstl1=1;    

    /* loading the kernels into mode structure */
    while (read(fd,&rcdhdr,sizeof(rcdhdr))) {
      if (read(fd,knl,pdim*sizeof(float))!=pdim*sizeof(float))
        stop("setmodes_A3_anis: reading error");
    
      nn[ik][i]=rcdhdr.n; 
      ll[ik][i]=rcdhdr.l;
      
      if (rcdhdr.l==1 && firstl1==1) {
	firstl1=0;
	*nstart_l1=rcdhdr.n;
      }  
      
      if (rcdhdr.l==1) {
	if (i!=start[ik][(int)rcdhdr.l]+rcdhdr.n-*nstart_l1) {
	  printf("l=%d n=%d i=%d , start[%d]=%d"
		 ,rcdhdr.l,rcdhdr.n,i,rcdhdr.l,start[ik][(int)rcdhdr.l]);
	  stop("setmodes_A3_anis: error 21");
	}
      }           	
      else if (i!=start[ik][(int)rcdhdr.l]+rcdhdr.n) {	
	printf("l=%d n=%d i=%d , start[%d]=%d"
	       ,rcdhdr.l,rcdhdr.n,i,rcdhdr.l,start[ik][(int)rcdhdr.l]);
	stop("setmodes_A3: error 22");
      } 
            
      if (ellswtch)
        ptl[ik][i][pmax]=4.*sqrt(8.*PI/5.)*rcdhdr.ell*rcdhdr.w; 
      else  
        ptl[ik][i][pmax]=0.;
        
      pp=0;
      fac=0.5/rcdhdr.w;
    
/* create the combined kernels for elastic structure in N,L and F  */

/* mp->nquant==1 inversion for Vsiso */

      if (mp->nquant==1) { 

/* VSiso kernel */ 
	if (mp->descript[0]!='S')
	  stop("setmodes_A3_anis wrong order of variables (S)");

        for (p=0;p<kmax;p++) {
           ptl[ik][i][pp++]=fac*(knl[p+3*kmax]      /* Kvs */
	                    +prop*knl[p]            /* Kvp */
                            +prop1*knl[p+5*kmax]);  /* Krho */
        }
/* discontinuities */
        if (mp->seafl) { ptl[ik][i][pp++]=fac*knl[6*kmax];} 	  
        if (mp->moho)  { ptl[ik][i][pp++]=fac*knl[6*kmax+1];}
        /* printf("n = %d  l = %d  pp = %d  seafloor  = %g  fac = %g\n",nn[ik][i],ll[ik][i],pp,ptl[ik][i][16],fac);
         printf("n = %d  l = %d  pp = %d  moho  = %g  fac = %g\n",nn[ik][i],ll[ik][i],pp,ptl[ik][i][17],fac); } */
	if (mp->d670)  ptl[ik][i][pp++]=fac*knl[6*kmax+2];	
        if (mp->cmb)   ptl[ik][i][pp++]=fac*knl[6*kmax+3];
	if (pp!=pmax)
	  stop("setmodes_A3_anis pp!=pmax");


/* mp->nquant==2 inversion for Vs and xi */

      } else if (mp->nquant==2) { 

/* VSiso kernel */ 
	if (mp->descript[0]!='S')
	  stop("setmodes_A3_anis wrong order of variables (S)");

        for (p=0;p<kmax;p++) {
           ptl[ik][i][pp++]=fac*(knl[p+3*kmax]      /* Kvs */
	                    +prop*knl[p]            /* Kvp */
                            +prop1*knl[p+5*kmax]);  /* Krho */
	
        }

/* xi kernel */
	if (mp->descript[1]!='X')
	  stop("setmodes_A3_anis wrong order of variables (X)");   

        for (p=0;p<kmax;p++) {
           ptl[ik][i][pp++]=fac*(knl[p+4*kmax]      /* Kxi */
	                   +prop2*knl[p+2*kmax]     /* Keta */
	                   +prop3*knl[p+kmax]);     /* Kphi */
	   /*	   if(pp==24) printf("ik = %d  i = %d  pp = %d Xi kernel = %g\n",ik,i,pp,ptl[ik][i][pp]);*/
        }

/* discontinuities */
        if (mp->seafl) ptl[ik][i][pp++]=fac*knl[6*kmax];
        if (mp->moho)  ptl[ik][i][pp++]=fac*knl[6*kmax+1];
        if (mp->d670)  ptl[ik][i][pp++]=fac*knl[6*kmax+2];
        if (mp->cmb)   ptl[ik][i][pp++]=fac*knl[6*kmax+3];
	if (pp!=pmax)
	  stop("setmodes_A3_anis pp!=pmax");

      } else 
        stop ("setmodes_A3_anis: this option has not been implemented yet");

      mode[i].ptl=farray1(0,pmax);  /* ptl[pmax] used for ell */
	
      i++;
    } /* end  while read(fd) */
  
    close(fd);
    nmode[ik]=i;
      
    if (nmode[ik]>=sum_modes || nmode[ik]>=NMODE) stop("setmodes_A3_anis: parameter sum_modes too small");
  }	/* end first time around */
  
/* load eigenfunctions appropriate for modetype ik */
 
  for (i=0;i<nmode[ik];i++)
    for (p=0;p<=pmax;p++) /* Loop goes to pmax instead of one less than pmax because ptl[ik][i][pmax] is where the ellipticity is stored */ 
       mode[i].ptl[p]=ptl[ik][i][p];

  r=(6371000.-title->depth)/6371000.0;
  R=farray2(0,nmode[ik]-1,0,2);
  S=farray2(0,nmode[ik]-1,0,3);
	
  for (i=0;i<nmode[ik];i++) {
  /* mode 0 is 0S0 : prem_st returns junk <-- non e' vero!*/
  if (EIG=='X')
    getprem(nn[ik][i],mode_type,ll[ik][i],title->depth,0,&exist,&prem[ik][i]);
  else if (EIG=='Y')
   exist=geteigys_ellrot(nn[ik][i],mode_type,ll[ik][i],title->depth,&prem[ik][i]);
  else stop("setmodes_A3: error EIG setup"); 
     
   /* if (nn[ik][i]==0 && ll[ik][i]==0 && mode_type=='T') {
      if (i!=0) stop("setmodes_A3: index error");
      R[i][0]=R[i][1]=R[i][2]=0.0;
      S[i][0]=S[i][1]=S[i][2]=S[i][3]=0.0;
    } else excit_q(r,&prem[ik][i],trace,R[i],S[i],dSdr,Di,0);*/
    excit_q(r,&prem[ik][i],trace,R[i],S[i],dSdr,Di,0);
  }
   
  pavaX20=farray2(0,3,1,2);
  cX20=farray2(0,3,0,4);

  /* Transform source parameters */
  src=transource(title,trace);

  /* ellipticity correction : modified for focusing */ 
  collapsX20_q(cX20,title->theta,title->phi,trace->theta,trace->phi);

  /* ob=1 minor arc averages
     ob=2 great circle averages */

  for (ob=1;ob<=2;ob++) {

    if (ob==1)  /*minor arc */
      delta=trace->delta;
    else       /* 2pi circle */
      delta=PI*2;      
    
    pavaX20[0][ob]=cX20[0][0];     
    arg=delta*2; /*h=2 for ell correction */
    f1=sin(arg)/(arg);
    f2=(1.0-cos(arg))/(arg);
    pavaX20[0][ob]+=cX20[0][3]*f1+cX20[0][4]*f2;    
  }		/* end ob loop */

  for (i=0;i<nmode[ik];i++) {  /*0*/
    double add1,add2; 
      
    add1=add2=0;
    if (unconf==0) {
      for (p=0;p<pmax;p++) {
        for (j=0;j<neffk[0];j++) {
	  kid = Aknot[j].index+mp->nknotA2[0]*p;
	  add1+=mode[i].ptl[p]*Aknot[j].pava[1]*mdl[kid];
	  add2+=mode[i].ptl[p]*Aknot[j].pava[0]*mdl[kid];
       }
      }
    } else {
      begin=0;
      end=0;
      hjump=0;    
      vjump=0;
      if (mp->nquant!=2) 
        stop("setmodes_A3: This option has not been implemented yet");
      for (ii=0;ii<mp->nquant;ii++) {
        end+=neffk[ii];
	if (ii!=0) {
	  hjump+=mp->nknotA2[ii-1]*mp->nknotA1[ii-1];
	  vjump+=mp->nknotA1[ii-1];
	}
        for (p=0;p<mp->nknotA1[ii];p++) {
          for (j=begin;j<end;j++) {
	    kid = Aknot[j].index+hjump+mp->nknotA2[ii]*p;
	    add1+=mode[i].ptl[vjump+p]*Aknot[j].pava[1]*mdl[kid];
	    add2+=mode[i].ptl[vjump+p]*Aknot[j].pava[0]*mdl[kid];
	  }
        }
	begin+=neffk[ii];
      }
      
      hjump+=mp->nknotA2[1]*mp->nknotA1[1];
      vjump+=mp->nknotA1[1];

      if (mp->seafl) {
        for (j=0;j<neffk[0];j++) {
	  kid = Aknot[j].index+hjump;
	  add1+=mode[i].ptl[vjump]*Aknot[j].pava[1]*mdl[kid];
	  add2+=mode[i].ptl[vjump]*Aknot[j].pava[0]*mdl[kid];
	}
	hjump+=mp->nknotA2[0];
        vjump+=1;
      }

      if (mp->moho) {
        for (j=0;j<neffk[0];j++) {
	  kid = Aknot[j].index+hjump;
	  add1+=mode[i].ptl[vjump]*Aknot[j].pava[1]*mdl[kid];
	  add2+=mode[i].ptl[vjump]*Aknot[j].pava[0]*mdl[kid];
	}
	hjump+=mp->nknotA2[0];
        vjump+=1;
      }

      if (mp->d670) {
        for (j=0;j<neffk[0];j++) {
	  kid = Aknot[j].index+hjump;
	  add1+=mode[i].ptl[vjump]*Aknot[j].pava[1]*mdl[kid];
	  add2+=mode[i].ptl[vjump]*Aknot[j].pava[0]*mdl[kid];
	}
	hjump+=mp->nknotA2[0];
        vjump+=1;
      }
       
      if (mp->cmb) {
        for (j=0;j<neffk[0];j++) {
	  kid = Aknot[j].index+hjump;
	  add1+=mode[i].ptl[vjump]*Aknot[j].pava[1]*mdl[kid];
	  add2+=mode[i].ptl[vjump]*Aknot[j].pava[0]*mdl[kid];
	}
	hjump+=mp->nknotA2[0];
        vjump+=1;
      }
    }
       
    mode[i].dw[1]=add1/trace->delta+mode[i].ptl[pmax]*pavaX20[0][1];
    mode[i].dw[2]=add2/PI/2.+mode[i].ptl[pmax]*pavaX20[0][2];
        
    /* mode[i].rs[0] = R0~
       mode[i].rs[1] = -R1~
       mode[i].rs[2] = -S0~
       mode[i].rs[3] = -S1~
       mode[i].rs[4] = -S2~
       mode[i].rs[5] = R-1~
       mode[i].rs[6] = S-1~ 
       mode[i].rs[7] = S-2~   modification BR 10/16/96 */
    
    mode[i].rs[0]= R[i][0];
    mode[i].rs[1]= R[i][1];
    mode[i].rs[2]= S[i][0]*src->moment[0]+
    S[i][1]*.5*(src->moment[1]+src->moment[2]);
    
    if (mode_type=='S') {
      mode[i].rs[4]=-S[i][3]*.5*(src->moment[2]-src->moment[1]);
      mode[i].rs[3]= S[i][2]*src->moment[4];
      mode[i].rs[6]=-S[i][2]*src->moment[3];
      mode[i].rs[7]= S[i][3]*src->moment[5];
    } else if(mode_type=='T') { 
      mode[i].rs[4]=-S[i][3]*src->moment[5];
      mode[i].rs[3]= S[i][2]*src->moment[3];
      mode[i].rs[6]= S[i][2]*src->moment[4];
      mode[i].rs[7]= S[i][3]*.5*(src->moment[1]-src->moment[2]);
    } else stop("setmodes_A3_anis: non-existent modetype");
    
    mode[i].rs[5]=0.; /* for self-coupling */
  } 	/* 0, end sum on modes */

free_farray2(cX20,0,3,0);
free_farray1(pavaX20,0,3,1);  
free_farray2(R,0,nmode[ik]-1,0);
free_farray2(S,0,nmode[ik]-1,0);

return(start);
}

/*********************************************************************/

/* transform source parameters into equatorial coordinate frame
   the azimuth (arg) is counted anti-clockwise from south */

source_st *transource(title,trace)
title_st *title;
tracehdrH_st *trace;
{
float ss,cc,cs;
double c,s,arg;
static source_st new;

  new.theta=HPI;
  new.phi=0.0;
  new.dep=title->depth;
  
  arg=(double)trace->az;
  c=cos(arg);
  s=sin(arg);
  cc=c*c;
  ss=s*s;
  cs=c*s;
  
  new.moment[0]=title->moment.rr;
  new.moment[1]=ss*title->moment.tt+cc*title->moment.pp
                -2.*cs*title->moment.tp;
  new.moment[2]=cc*title->moment.tt+ss*title->moment.pp
                +2.*cs*title->moment.tp;
  new.moment[3]=s*title->moment.rt-c*title->moment.rp;
  new.moment[4]=c*title->moment.rt+s*title->moment.rp;
  new.moment[5]=cs*(title->moment.tt-title->moment.pp)
                +(ss-cc)*title->moment.tp;
                
  return(&new);

}
/*******************************************************************/
void locaknl(file,v,rclen)
int file;
float v;
long rclen;
{
int i;
long jump;
/* struct bdrcd knlrcd;  */
knlrechdr_st knlrcd;
int hdrlen;
float grv[100];
long lct[100];

  lseek(file,200L,0);
  read(file,grv,400);
  read(file,lct,400);
  
  for(i=0;i<100;i++) if(grv[i]>v) break;
  if(!i) i++;
  
  jump=rclen*lct[i-1];
  if(lseek(file,jump,1)!=jump+1000) stop("locaknl: reading error");
  
  hdrlen=sizeof(knlrcd);
  do
  { if(read(file,&knlrcd,hdrlen)<=0) 
         stop("locaknl: specified group velocity too large");
    lseek(file,(long)(rclen-hdrlen),1);
  } while(knlrcd.grv1<v||knlrcd.grv2<v);
  
  lseek(file,-rclen,1);

}
/****************************************************************************/
window_st *setwindow(nwp,wpinfo,gvmin,gvmax,knlhdr)
wpinfoved_st *wpinfo;
int nwp;
float *gvmin,*gvmax;
/* struct hdr1_2 *knlhdr; */
knlhdrA1_st *knlhdr;
{
static window_st window[10];
float dgrvel;
int i;
/* float grvlobe(),lobewidth; */
float lobewidth;
       
  if (nwp>10) stop("setwindow: not enough working space");
  *gvmin=1.0;
  *gvmax=0.0;

  for (i=0;i<nwp;i++) {  
     if (wpinfo[i].gv1<knlhdr->minvel) wpinfo[i].gv1=knlhdr->minvel;
     if (wpinfo[i].gv2>knlhdr->maxvel) wpinfo[i].gv2=knlhdr->maxvel;
     dgrvel=wpinfo[i].gv2-wpinfo[i].gv1;
/* lobewidth=GRVLOBE*grvlobe(wpinfo[i].t0,(int)wpinfo[i].phase); */
     
     lobewidth=GRVLOBE;
     window[i].pv1=wpinfo[i].pv1;
     window[i].pv2=wpinfo[i].pv2;	  	  
     window[i].gv1=wpinfo[i].gv1-dgrvel*lobewidth;
	  
     if (window[i].gv1<knlhdr->minvel) {
	window[i].gv1=knlhdr->minvel;
	window[i].gv2=wpinfo[i].gv1+dgrvel*lobewidth;
     } else
	window[i].gv2=wpinfo[i].gv1;

     window[i].gv3=wpinfo[i].gv2;
     window[i].gv4=wpinfo[i].gv2+dgrvel*lobewidth;
	  
     if (window[i].gv4>knlhdr->maxvel) {
	window[i].gv4=knlhdr->maxvel;
	window[i].gv3=knlhdr->maxvel-dgrvel*lobewidth;
     }

     if (window[i].gv4>*gvmax) *gvmax=window[i].gv4;
     if (window[i].gv1<*gvmin) *gvmin=window[i].gv1;
  }
        
  return(window);  
}
/*************************************************************************/
int slkrcd(n,knlrcd,window,w1,w2,gvmax)
int n;
/* struct bdrcd *knlrcd; */
knlrechdr_st *knlrcd;
window_st *window;
float gvmax,w1,w2;
{
int status,i;
float pv1,pv2,gv1,gv2;
  status=-1;
  pv1=knlrcd->w1/(0.5+(float)knlrcd->l1);
  pv2=knlrcd->w2/(0.5+(float)knlrcd->l2);
  gv1=knlrcd->grv1;
  gv2=knlrcd->grv2;
  for(i=0;i<n;i++)
  { if(gv1>=window[i].gv1 && gv1<=window[i].gv4 &&
             gv2>=window[i].gv1 && gv2<=window[i].gv4 &&
             pv1>=window[i].pv1 && pv1<=window[i].pv2 &&
             pv2>=window[i].pv1 && pv2<=window[i].pv2 &&
             knlrcd->w1>=w1 && knlrcd->w1<=w2 &&
             knlrcd->w2>=w1 && knlrcd->w2<=w2) return(1);
    if(gv1<=gvmax||gv2<=gvmax) status=0;
  }
  return(status);
}


/***********************************************************************/
int whichorbit(phase)
short phase;
{
int orbit;
  if(phase>0) orbit=1;
  else if(phase<0)
    { 
      if(phase==-11|| phase==-21||phase==-51||phase==-52) orbit=1;
      else if(phase==-31 || phase==-41||phase==-61|phase==-81) orbit=1;
      else if(phase==-12|| phase==-22||phase==-53||phase==-54||phase==-71) orbit=2;
      else if(phase==-32|| phase==-42|| phase==-62|| phase==-63) orbit=2;
      else if(phase==-13|| phase==-23||phase==-33|| phase== -43||phase==-64) orbit=3;  
      else if(phase==-14|| phase==-24||phase==-34|| phase==-44) orbit=4;
      else if(phase==-15|| phase==-25|| phase==-35|| phase==-45) orbit=5;
      else if(phase==-16|| phase==-26|| phase==-36|| phase==-46) orbit=6;      
	
      /* else stop("whichorbit: error 1 -- undefined phase"); */

      else orbit=1;

  }
  else stop("whichorbit: error 2 -- undefined phase");
  
  return(orbit);
}


/***************************************************************************/

void collapsX20_q(c,ths,phs,thr,phr)
float **c,ths,phs,thr,phr;
{
float **coef,**coef1,**coef2;

/* c[0] = model ; c[1]=dm/dtheta; c[2]=d2m/dtheta2; c[3]=dm/dfi  */
    /* CHM 10/31/96 */
    coef =farray2(0,4,0,4);
    coef1=farray2(0,4,0,4);
    coef2=farray2(0,4,0,4);

    rotcoef_q(1,2,2,ths,phs,thr,phr,coef,coef1,coef2);
    
    c[0][0]=coef[0][0];
    c[0][1]=0.0;
    c[0][2]=0.0;
    c[0][3]=coef[3][0];
    c[0][4]=coef[4][0];
    free_farray2(coef,0,4,0);

    c[1][0]=coef1[0][0];
    c[1][1]=0.0;
    c[1][2]=0.0;
    c[1][3]=coef1[3][0];
    c[1][4]=coef1[4][0];
    free_farray2(coef1,0,4,0);

    c[2][0]=coef2[0][0];
    c[2][1]=0.0;
    c[2][2]=0.0;
    c[2][3]=coef2[3][0];
    c[2][4]=coef2[4][0];

    c[3][0]=coef2[0][0];
    c[3][1]=0.0;
    c[3][2]=0.0;
    c[3][3]=-2.*c[0][4];
    c[3][4]=2.*c[0][3];
    free_farray2(coef2,0,4,0);


}
/***************************************************************/

void collapsX20(c,ths,phs,thr,phr)
float *c,ths,phs,thr,phr;
{
float **coef;
void rotcoef();

    coef=farray2(0,4,0,4);
    rotcoef(1,2,2,ths,phs,thr,phr,coef);
    
    c[0]=coef[0][0];
    c[1]=0.0;
    c[2]=0.0;
    c[3]=coef[3][0];
    c[4]=coef[4][0];
    free_farray2(coef,0,4,0);
}
/***************************************************************/
void rotcoef_q(newpath,smax,s,ths,phs,thr,phr,m,m1,m2)
int newpath;
int smax;
int s;
float ths,phs,thr,phr;
float **m,**m1,**m2;
{
static int first=0;
static float *cosg,*sing,*cosa,*sina;
static float alpha,beta,gamma;
float **d1,**d2;
float *x,*xp,*xcosec,*x2p;
float fac,f,f1,f2,cc,ss,cs,sc;
float f11,f12,f21,f22,hpi;
int t1,t2,t,t0;
  
  hpi=(float)HPI;
  if (s>smax) stop("rotcoef_q: illegal argument");
  
  if (newpath) {
    float calpha,salpha,cgamma,sgamma;
    Euler(ths,phs,thr,phr,&alpha,&beta,&gamma);
    if (first) {
      free_farray1(cosg,0);
      free_farray1(sing,0);
      free_farray1(cosa,0);
      free_farray1(sina,0);
    }
    first=1;
    cosg=farray1(0,smax);
    sing=farray1(0,smax);
    cosa=farray1(0,smax);
    sina=farray1(0,smax);
    cosg[0]=cosa[0]=1.0;
    sing[0]=sina[0]=0.0;
    calpha=cos((D)alpha);
    salpha=sin((D)alpha);
    cgamma=cos((D)gamma);
    sgamma=sin((D)gamma);
    for (t=1;t<=smax;t++) {
      cosa[t]=cosa[t-1]*calpha-sina[t-1]*salpha;
      sina[t]=cosa[t-1]*salpha+sina[t-1]*calpha;
      cosg[t]=cosg[t-1]*cgamma-sing[t-1]*sgamma;
      sing[t]=cosg[t-1]*sgamma+sing[t-1]*cgamma;
    }
  }
   
  d1=farray2(-s,s,-s,s);
  rotmatrix(s,s,beta,d1);
  d2=farray2(0,0,-s,s);
  rotmatrix(0,s,HPI,d2);

  fac=sqrt((D)(2*s+1)/(4.*PI));

/* modif. 10/02/96 to compute derivatives of Plm */
/* xp is first derivative, x2p is second derivative */


  x=farray1(0,s);
  xp=farray1(0,s);
  x2p=farray1(0,s);
  xcosec=farray1(0,s);
  LEGENDRE(&hpi,&s,&s,x,xp,xcosec); 

  for (t1=0;t1<=s;t1++) {
    x[t1]=x[t1]/fac;
    xp[t1]=xp[t1]/fac;
    x2p[t1]=-(s*(s+1)-t1*t1)*x[t1];
  }

/* end modif 10/02/96 */  

/* selection rule t2+s=even otherwise Plm(pi/2)=0 */
  
  if (!(s%2)) {
    m[0][0]=SRH*fac*d1[0][0]*d2[0][0];
    m1[0][0]=SRH*fac*d1[0][0]*xp[0];
    m2[0][0]=SRH*fac*d1[0][0]*x2p[0];

    for (t1=1;t1<=s;t1++) {
      f=fac*d1[0][t1]*d2[0][0];
      f1=fac*d1[0][t1]*xp[0];
      f2=fac*d1[0][t1]*x2p[0];
      m[0][2*t1-1]=f*cosa[t1];
      m1[0][2*t1-1]=f1*cosa[t1];
      m2[0][2*t1-1]=f2*cosa[t1];
      m[0][2*t1  ]=f*sina[t1];
      m1[0][2*t1  ]=f1*sina[t1];
      m2[0][2*t1  ]=f2*sina[t1];
    }
  }
  
  if (s%2) t0=1;
  else t0=2;

  for (t2=t0;t2<=s;t2+=2) {
    f1=fac*d1[t2 ][0]*d2[0][ t2];
    f2=fac*d1[-t2][0]*d2[0][-t2];
    m[2*t2-1][0]= (f1+f2)*cosg[t2]*SRH;
    m[2*t2  ][0]=-(f1+f2)*sing[t2]*SRH;

    f11=fac*d1[t2][0]*xp[t2];

/* sign changes only for s odd in f12 and f22 when t2 goes to -t2: 
BR 08/20/99 */

    if (!(s%2)) f12=fac*d1[-t2][0]*xp[t2];
    else f12=-fac*d1[-t2][0]*xp[t2];
    f21=fac*d1[t2][0]*x2p[t2];
    if (!(s%2)) f22=fac*d1[-t2][0]*x2p[t2];
    else f22=-fac*d1[-t2][0]*x2p[t2];
    m1[2*t2-1][0]= (f11+f12)*cosg[t2]*SRH;
    m1[2*t2  ][0]=-(f11+f12)*sing[t2]*SRH;
    m2[2*t2-1][0]= (f21+f22)*cosg[t2]*SRH;
    m2[2*t2  ][0]=-(f21+f22)*sing[t2]*SRH;

    for (t1=1;t1<=s;t1++) {
      f1=fac*d1[ t2][t1]*d2[0][ t2];
      f2=fac*d1[-t2][t1]*d2[0][-t2];
      f11=fac*d1[ t2][t1]*xp[ t2];
      if (!(s%2)) f12=fac*d1[-t2][t1]*xp[t2];
      else f12=-fac*d1[-t2][t1]*xp[t2];
      f21=fac*d1[ t2][t1]*x2p[ t2];
      if (!(s%2)) f22=fac*d1[-t2][t1]*x2p[t2];
      else f22=-fac*d1[-t2][t1]*x2p[t2];

      cc=cosa[t1]*cosg[t2];
      ss=sina[t1]*sing[t2];
      cs=cosa[t1]*sing[t2];
      sc=sina[t1]*cosg[t2];
      m[2*t2-1][2*t1-1]= f1*(cc-ss)+f2*(cc+ss);
      m[2*t2-1][2*t1  ]= f1*(sc+cs)+f2*(sc-cs);
      m[2*t2  ][2*t1-1]=-f1*(sc+cs)+f2*(sc-cs);
      m[2*t2  ][2*t1  ]= f1*(cc-ss)-f2*(cc+ss);

      m1[2*t2-1][2*t1-1]= f11*(cc-ss)+f12*(cc+ss);
      m1[2*t2-1][2*t1  ]= f11*(sc+cs)+f12*(sc-cs);
      m1[2*t2  ][2*t1-1]=-f11*(sc+cs)+f12*(sc-cs);
      m1[2*t2  ][2*t1  ]= f11*(cc-ss)-f12*(cc+ss);

      m2[2*t2-1][2*t1-1]= f21*(cc-ss)+f22*(cc+ss);
      m2[2*t2-1][2*t1  ]= f21*(sc+cs)+f22*(sc-cs);
      m2[2*t2  ][2*t1-1]=-f21*(sc+cs)+f22*(sc-cs);
      m2[2*t2  ][2*t1  ]= f21*(cc-ss)-f22*(cc+ss);

    }
  }
  
  free_farray2(d1,-s,s,-s);
  free_farray2(d2,0,0,-s);
  free_farray1(x,0);
  free_farray1(xp,0);
  free_farray1(x2p,0);
  free_farray1(xcosec,0);
  
}
/************************************/

void Gelement_anis(g,a,mode,l,fac,A,RS)
/* struct bdrcd *mds; */
/*knlrechdr_st *mds; */
mode_anis_st mode;
float **g,**a,**A,*RS;      /* a[3][3] should be defined in selfterm_q */
float fac;        
int l;
{
/*static float P[4]; *//* P_{kk'}^1, P_{kk'}^2, P_{k'k}^1, P_{k'k}^2 */
/* static float **A; */
int l1,l2,j,j1;
float rs0,rs1,rs2,rs3,rs4,rs5,rs6,rs7,k,diff;  
  
  l1=l2=l;
  /* l1=mds->l1; structure mds not passsed anynore 11/29/96
  l2=mds->l2;*/


/* modified 10/11/96 by BR to include rs1[0] and rs2[0] terms for
spheroidal modes and new definition of rs2 */

  rs1=mode.rs[1];
  rs2=mode.rs[2]; 
  rs3=mode.rs[3]; 
  rs4=mode.rs[4]; 
  rs0=mode.rs[0];
  rs5=mode.rs[5]; 
  rs6=mode.rs[6];
  rs7=mode.rs[7];
  
/*  printf("\nrs in Gelement: rs0 %g rs1 %g rs2 %g rs3 %g rs4 %g rs5 %g rs6 %g rs7 %g\n",rs0,rs1,rs2,rs3,rs4,rs5,rs6,rs7);*/

/* j=0  zeroth order terms, defined as previously T[0],T[1] not needed in this
	version 
   j=1-3: F1 -> F6 
   j=5-6: F7 -> F10 */
	RS[0]=rs1*rs3+rs0*rs2+rs0*rs4; /* T0 */
	RS[1]=rs1*rs4+rs1*rs2-rs0*rs3; /* T1 */
	RS[2]=2.*rs0*rs7+rs1*rs6; /* T2 */
	RS[3]=-rs0*rs6+2.*rs1*rs7; /* T3 */
	RS[4]=rs5*rs3; /*T4 */
	RS[5]=-rs5*(rs2+rs4); /*T5*/
	RS[6]= 2.*(rs1*rs3+2.*rs4*rs0);/* T6 */

	RS[7]= -rs0*rs3+rs1*rs2+5.*rs1*rs4;/* T7 */
	RS[8]= rs5*rs6; /* T8 */
	RS[9]= 2.*rs5*rs7; /* T9 */

/*	for(i=0;i<6;i++)RS[i]=-RS[i];  modif 10/15/97 */

/* expressions for A[i][j] last modified 08/26/99 br */

/* zeroth order terms: */

	A[0][0]=RS[0];
	A[1][0]=RS[1];
/*	printf(" A00 %g A10 %g\n",A[0][0],A[1][0]);fflush(stdout); */

/* additional 1/l order terms: */

for(j=1;j<=3;j++)
{  A[0][j]= a[0][j]*RS[0] ;
	A[0][j]+=a[1][j]*RS[8]-a[2][j]/2.*RS[6];
   A[1][j]= a[0][j]*RS[1] ;
	A[1][j]+=a[1][j]*RS[9]-a[2][j]/2.*RS[7];
/*	printf("a0 %g a1 %g i %d\n",A[0][j],A[1][j],j); */ }

for(j=4;j<=5;j++)
{  j1=j-3;
   A[0][j]=a[3][j1]*RS[2]-a[4][j1]*RS[4];
   A[1][j]=a[3][j1]*RS[3]-a[4][j1]*RS[5];
/*	printf("a0 %g a1 %g i %d\n",A[0][j],A[1][j],j); */ }

}
/******************************************************************/

void a_term_anis(delta,mode,a)
mode_anis_st mode;
float delta,**a;
{
float sd,cd;

sd=sin(delta);
cd=cos(delta);

/* zeroth-order terms */

a[0][0]=1.;
a[1][0]=0.;
a[2][0]=0.;

/* X-terms */

a[0][1]=(mode.ddw[2]-mode.ddw[1])*0.5;
a[0][1]+=(mode.dw[2]-mode.dw[1])/sd*cd/8.;
a[1][1]=(mode.dw[2]-mode.dw[1])/sd;
a[2][1]=a[1][1]*cd;
a[3][1]=mode.e1w[2]-mode.e1w[1];
a[4][1]=mode.e2w[2]-mode.e2w[1];
/* printf("a01 %g a31 %g a41 %g\n",a[0][1],a[3][1],a[4][1]); */

/* A-terms */

a[0][2]=mode.ddw[2]*0.5;
a[0][2]+=mode.dw[2]/sd*cd/8.;
a[1][2]=mode.dw[2]/sd;
a[2][2]=a[1][2]*cd;
a[3][2]=mode.e1w[2];
a[4][2]=mode.e2w[2];

/* alpha-terms */
/* last modified 08/.26/99 br */

a[0][3]=-cd/sd/8;
/*	printf ("a03 %g\n",a[0][3]); */
a[1][3]=-1./sd;
a[2][3]=a[1][3]*cd;

}
/***************************************************************************/
float **garray(delta)
float delta;
{
static int first=1;
static float **g;
float lambda,fac,arg;
int l;
  
  if(first)
  { first=0;
    g=farray2(0,LMAX,0,1);
  }
  
  for(l=0;l<=LMAX;l++)
  { lambda=0.5+(F)l;
  
    fac=lambda*HPI*(F)sin((D)delta);

/*    fac=CONST*(F)sqrt(1./fac); */
/* modification 10/11/96 BR */
    fac= (F)sqrt(1./fac);
    
    arg=lambda*delta-0.5*HPI;
    
    g[l][0]=fac*(F)cos((D)arg);
    g[l][1]=fac*(F)sin((D)arg);
  }
  
  return(g);
}





/****************************************************/

void rotcoef(newpath,smax,s,ths,phs,thr,phr,m)
int newpath;
int smax;
int s;
float ths,phs,thr,phr;
float **m;
{
static int first=0;
static float *cosg,*sing,*cosa,*sina;
static float alpha,beta,gamma;
float **d1,**d2;
float fac,f,f1,f2,cc,ss,cs,sc;
int t1,t2,t,t0;
  
  if(s>smax) stop("rotcoef: illegal argument");
  
  if(newpath)
  { float calpha,salpha,cgamma,sgamma;
    Euler(ths,phs,thr,phr,&alpha,&beta,&gamma);
    if (first) {
      free_farray1(cosg,0);
      free_farray1(sing,0);
      free_farray1(cosa,0);
      free_farray1(sina,0);
    }
    first=1;
    cosg=farray1(0,smax);
    sing=farray1(0,smax);
    cosa=farray1(0,smax);
    sina=farray1(0,smax);
    cosg[0]=cosa[0]=1.0;
    sing[0]=sina[0]=0.0;
    calpha=cos((D)alpha);
    salpha=sin((D)alpha);
    cgamma=cos((D)gamma);
    sgamma=sin((D)gamma);
    for(t=1;t<=smax;t++)
    { cosa[t]=cosa[t-1]*calpha-sina[t-1]*salpha;
      sina[t]=cosa[t-1]*salpha+sina[t-1]*calpha;
      cosg[t]=cosg[t-1]*cgamma-sing[t-1]*sgamma;
      sing[t]=cosg[t-1]*sgamma+sing[t-1]*cgamma;
    }
  }

/*  m=farray2(0,2*s,0,2*s);*/
  
  d1=farray2(-s,s,-s,s);
  rotmatrix(s,s,beta,d1);
  d2=farray2(0,0,-s,s);
  rotmatrix(0,s,HPI,d2);
  
  fac=sqrt((D)(2*s+1)/(4.*PI));

/* selection rule t2+s=even (otherwise Plm(pi/2)=0) */
  
  if(!(s%2))
  { m[0][0]=SRH*fac*d1[0][0]*d2[0][0];
    for(t1=1;t1<=s;t1++)
    { f=fac*d1[0][t1]*d2[0][0];
      m[0][2*t1-1]=f*cosa[t1];
      m[0][2*t1  ]=f*sina[t1];
    }
  }
  
  if(s%2) t0=1;
  else t0=2;
 for(t2=t0;t2<=s;t2+=2)
  { f1=fac*d1[t2 ][0]*d2[0][ t2];
    f2=fac*d1[-t2][0]*d2[0][-t2];
    m[2*t2-1][0]= (f1+f2)*cosg[t2]*SRH;
    m[2*t2  ][0]=-(f1+f2)*sing[t2]*SRH;
    for(t1=1;t1<=s;t1++)
    { f1=fac*d1[ t2][t1]*d2[0][ t2];
      f2=fac*d1[-t2][t1]*d2[0][-t2];
      cc=cosa[t1]*cosg[t2];
      ss=sina[t1]*sing[t2];
      cs=cosa[t1]*sing[t2];
      sc=sina[t1]*cosg[t2];
      m[2*t2-1][2*t1-1]= f1*(cc-ss)+f2*(cc+ss);
      m[2*t2-1][2*t1  ]= f1*(sc+cs)+f2*(sc-cs);
      m[2*t2  ][2*t1-1]=-f1*(sc+cs)+f2*(sc-cs);
      m[2*t2  ][2*t1  ]= f1*(cc-ss)-f2*(cc+ss);
    }
  }
  
  free_farray2(d1,-s,s,-s);
  free_farray2(d2,0,0,-s);
  
}
