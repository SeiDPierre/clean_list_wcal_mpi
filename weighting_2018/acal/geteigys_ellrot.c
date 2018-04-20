/* This routine reads reformatted eigenfunction files used for the
   computation of the source terms written by yannos_ST_format800km_3_ellrot.f
   and fills a variable of structure sprem_st. 
   The ellipticity and rotational splitting parameters are stored in the read
   eigenfunction files. The terms a6,a7,a8,a9,a10 and aker[i] are set to 0.
   This routine makes use of the modified structure eigenh_ellrot_st.
   
   Modified after geteigys.c (Federica 2/25/2005)

   bu_premy.bin is an fortran output file, in which	     
   (ACCESS=sequential and FORM=unformatted ). Each record is     
   preceded and terminated with an integer*4 count, making each 
   record 8 characters longer than normal 	
   that's why a 4 byte dummy read is used before and after 	
   the data string 						    
*/

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>

#include "acal/dimensiony.h"
#include "acal/const.h"
#include "acal/structprem.h" 
#include "acal/eigeny_st.h"

extern void stop();

int geteigys_ellrot(n,modetype,l,depth,prem)

int n,l;
char modetype;
float depth;
sprem_st *prem;

{
  eigenh_ellrot_st eigh;
  float *farray1();

  static int nots,nsurf;
  static int first=1;
  static float prvdep=-1.;
  static long jump;
  static int iq1,iq2;

  /* the default values of static variable is 0 */
  static int nrec[2][MAXL],maxl[2],minl[2];
  static int minn[2][MAXL],maxn[2][MAXL];
  static int total_l[2],irlast[2],eigsize[2];
  static int eigfile,ksize;
  static float *r;  

  int type,nrecord,position,i;
  int kjump;
  float rs,fdumy;

  if (first) /* read header for eigenfunctions */
  { 
    int i;
    int dummy,sz1,sz2;

    first=0;
    if ((eigfile=open("bu_premy.bin",0))<0) 
      stop("geteigys_ellrot: cannot open file bu_premy.bin");

    /* The 1st record in yannos_ST_format writing */
    read(eigfile,&dummy,4); 
    read(eigfile,total_l,8);
    read(eigfile,irlast,8);
    read(eigfile,&nots,4);    

    r=farray1(0,nots-1);
  
    read(eigfile,&nsurf,4);

    read(eigfile,eigsize,8);
    read(eigfile,r,nots*4);
    read(eigfile,minl,8);
    read(eigfile,maxl,8);
 
    sz1=(maxl[0]+1)*4; /* the record of maxl starts from 0 */
    sz2=(maxl[1]+1)*4;
    read(eigfile,minn[0],sz1);
    read(eigfile,minn[1],sz2);
    read(eigfile,maxn[0],sz1);
    read(eigfile,maxn[1],sz2);
    read(eigfile,nrec[0],sz1);
    read(eigfile,nrec[1],sz2);
    read(eigfile,&dummy,4); 
    /* End of the 1st record in yannos_ST_format writing */

    jump=lseek(eigfile,0L,1);
    jump+=4;    
  
  }

  if (depth!=prvdep) {
    prvdep=depth;
    rs=1.-depth/RADIUS;
    for (i=0;i<nots;i++) if (r[i]>rs) break;
    if (!i)stop("geteigys_ellrot: error in search of source depth");
    iq1=i-1;
    iq2=i;
  }
  
  if (n==0 && l==1 && modetype=='T') return(0);
    
  prem->r1=r[iq1];
  prem->r2=r[iq2];
  
  prem->n=n;
  prem->l=l;
  prem->typ=modetype;
 
  prem->a6=0.;
  prem->a7=0.;
  prem->a8=0.;
  prem->a9=0.;
  prem->a10=0.;
  for (i=0;i<=9;i++)prem->aker[i]=0.;

  /* eigsize[0] = eigsize of S mode
     eigsize[1] = eigsize of T mode */

  ksize=nots*4;

  type=0;
  if (modetype=='S') type=1;
  
  if (l>maxl[type] || l<minl[type]) return(0);
  if (n>maxn[type][l] || n<minn[type][l]) return(0);
  
  nrecord=(n-minn[type][l])+nrec[type][l];
  if (type==0) position=jump+(nrecord-1)*(eigsize[0]+8);
  else {    
    position=jump+(nrecord-1)*(eigsize[1]+8)
      -irlast[0]*(eigsize[1]-eigsize[0]);   
  }
  
  lseek(eigfile,(long)(position),0);
       
  read(eigfile,&eigh,sizeof(eigh));

  if (n!=eigh.n||l!=eigh.l) stop("geteigys_ellrot: ERROR 1");
     
  prem->w=eigh.w;
  prem->q=eigh.q;
  prem->va=eigh.a_vert;
  prem->ha=eigh.a_hori;
  prem->grv=eigh.cgp;
  prem->ell=eigh.ell;
  prem->rot=eigh.rot;
  
  kjump=iq1*sizeof(float);

  if (type==0) {
    lseek(eigfile,(long)(kjump),1);
    read(eigfile,&prem->v1,sizeof(float));    
    read(eigfile,&prem->v2,sizeof(float));   
    lseek(eigfile,(long)(ksize-8),1);
    read(eigfile,&prem->v1p,sizeof(float));    
    read(eigfile,&prem->v2p,sizeof(float));   
    prem->u1=0;prem->u1p=0;
    prem->u2=0;prem->u2p=0;
  
  } else {
    lseek(eigfile,(long)(kjump),1);
    read(eigfile,&prem->u1,sizeof(float));    
    read(eigfile,&prem->u2,sizeof(float));    
    lseek(eigfile,(long)(ksize-8),1);
    read(eigfile,&prem->u1p,sizeof(float));    
    read(eigfile,&prem->u2p,sizeof(float));   
    lseek(eigfile,(long)(ksize-8),1);
    read(eigfile,&prem->v1,sizeof(float));    
    read(eigfile,&prem->v2,sizeof(float));    
    lseek(eigfile,(long)(ksize-8),1);
    read(eigfile,&prem->v1p,sizeof(float));    
    read(eigfile,&prem->v2p,sizeof(float));   
  }

  if (modetype=='T' && eigh.type!=2) 
    stop("geteigys_ellrot: got wrong mode error 2");
  if (modetype=='S') {   
    if (l && eigh.type!=3) stop("geteigys_ellrot: got wrong mode error 3");
    if (!l && eigh.type!=1) stop("geteigys_ellrot: got wrong mode error 4");
  }

  return(1);
}
