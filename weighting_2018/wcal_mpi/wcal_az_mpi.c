#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <unistd.h>
#include <string.h>
#include <fcntl.h>
#include <sys/types.h>
#include <mpi.h>

/* local libs from acal */
#include "acal/const.h"
#include "acal/structwavket_A3z.h"

/* Max number of phases allowed */
#define MAXPHS 500000

extern void   wcal_az(char *event);
extern float  pathprojectn1();
extern void   stop();
extern void   fgetrecord();
extern int    dataselector_auto();     

int	main(int argc, char **argv)
{
        int rank, size;
	/* Initialize env */
        MPI_Init(&argc, &argv);
	/* Get the rank of procs */
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	/* Get the number of procs */
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        printf("[%i/%i] Threads initialized\n", rank, size - 1);
        if (rank == 0)
	{
		int msg;
                MPI_Status status;
                char event[9];

                /* Opening lu_events determining how many threads */
		FILE *f_evt = fopen("lu_events", "r");
                while (fscanf(f_evt, "%8s", event) == 1)
		{
                        MPI_Recv(&msg,  1, MPI_INT, MPI_ANY_SOURCE,     0, MPI_COMM_WORLD, &status);
                        MPI_Send(event, 9, MPI_CHAR, status.MPI_SOURCE, 1, MPI_COMM_WORLD);
                }

                int ndone = 0;
                while (ndone < size - 1)
		{
			MPI_Recv(&msg,  1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
                        MPI_Send(event, 9, MPI_CHAR, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
                        ndone += 1;
                }

                fclose(f_evt);
        } 
	else 
	{
		int msg;
                char event[9];
                MPI_Status status;

                /* check in with root and receive first event */
                MPI_Send(&msg, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                MPI_Recv(event, 9, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

                /* loop over events */
                while ( status.MPI_TAG )
		{
			double wt0, wt1;
                        /* calculate weights */
                        wt0 = MPI_Wtime();
			wcal_az(event);
			wt1 = MPI_Wtime();
                        printf("Rank [%2i]: Weighting of event %8s took %f s\n", rank, event, wt1 - wt0);

                        /* get new event */
                        MPI_Send(&msg, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                        MPI_Recv(event, 9, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                }

                printf("Rank [%2i] is done\n", rank);

        }

        MPI_Finalize();
}

/* Weight calculation according to Li & Romanowicz, 1996 */
/* Weighting with 2 different horizontal resolutions for receiver and source */
void
wcal_az(event)
  char *event;
{
  int   wpd,wph,wpN;
  FILE  *list,*par;

  title_st        title;
  tracehdrHH_st   trace;
  tracehdrAreg_st traceA;
  packhdrbr_st    packheader;

  float dep,thr,phr,ths,phs;
  float Depth,depth,stheta,sphi,rtheta,rphi;
  float disress,disresr,depres;    /* distance resolution (source, receiver) and depth resolution */
  char  string[120];
  int   nread,thsize,phsize,ndata,phase,p,ifphase;
  float rwgt,prj,contrib,ddep,dist;
  float vdep[MAXPHS],vths[MAXPHS],vphs[MAXPHS],vthr[MAXPHS],vphr[MAXPHS];
  int   vp[MAXPHS];
  int   i,nphase;

  /* Reading wp[dHN] files */
  strcpy(string,event);
  strcpy(&string[8],".wpH\0");

  wph = open(string,0);
  string[11]='d';
  wpd = open(string,2);
  string[11]='N';
  wpN = open(string,0);
  
  read( wph, &title, sizeof(title) );
  printf("Event = %s\n",string);
  
  Depth   = title.depth;
  stheta  = title.theta;
  sphi    = title.phi;

  if (sphi<0.0) sphi+=2.*PI;

  /* Reading distance resolution */
  if ( ( par=fopen("lu_wcal.par","r") ) == NULL ) 
  {
    stop("wcal_az : cannot open lu_wcal.par");
  }

  fgetrecord(par,string); sscanf(string,"%f",&disress);
  fgetrecord(par,string); sscanf(string,"%f",&disresr);
  fgetrecord(par,string); sscanf(string,"%f",&depres);
  
  /* Reading listed phases */
  if ( (list=fopen("lu_phase.lst","r")) == NULL ) 
  {
    stop("wcal_az : cannot open lu_phase.lst");
  }

  i = 0;
  /* Need to read source [depth,theta,phi], receiver[theta,phi], phase code */
  while (fscanf(list,"%f %f %f %f %f %d",&vdep[i],&vths[i],&vphs[i],&vthr[i],&vphr[i],&vp[i])!=EOF) 
  {
     if (vphs[i]<0.0) vphs[i]+=2.*PI; 
     if (vphr[i]<0.0) vphr[i]+=2.*PI; 
     i++;
     if (!i%1000) printf("i %d ",i);
  }
  printf("List of phases read correctly\n");
  fclose(list);

  nphase = i; 
  if (nphase>=MAXPHS) 
  {
    stop("wcal_az: hard coded max number of phases too small, see line 13");
  }

  printf("Number of listed phases is %d\n",nphase); 

  /* Loading traces in wp files */
  thsize = sizeof(trace);
  phsize = sizeof(packheader);
  while (nread=read(wph,&trace,thsize)) 
  {
    if (nread!=thsize) 
    {
      stop("wcal_az: error while reading event %s\n",string);
    }

    printf("Station = %.4s\n",trace.stn);
  
    /* Loading effective knots along paths */
    lseek(wpN,trace.locatnA,0);  
    read(wpN,&traceA,sizeof(traceA));     
    if ( traceA.neffkreg != 0 )
    {
      rtheta    = trace.theta;
      rphi      = trace.phi;

      if (rphi<0.0) rphi+=2.*PI;

      if (lseek(wpd,trace.locatn,0)<0) 
      {
        stop("wcal_az : Error while lseeking trace in %s.wpd\n",string);
      }

      while (1) 
      {
        if (!(nread=read(wpd,&packheader,phsize))) break;  
        if (nread!=phsize) 
        {
          stop("wcal_az : error while reading packheader in %s.wpd file\n",string);
        }

        ndata   = packheader.ndata;
        if (packheader.id!=trace.id) break;
	      if (dataselector_auto(&packheader,trace.delta,trace.netwk,trace.chnnl,trace.comp))
        { 
          phase=packheader.phase;
          if (phase==121) phase=20;
          if (phase==124) phase=40;
    
          /* Interaction with other path */
          rwgt    = 0.0;
          ifphase = 0;

          for ( i = 0; i < nphase; i++) 
          {
            dep = vdep[i]; ths = vths[i];
            phs = vphs[i]; thr = vthr[i];
            phr = vphr[i]; p   = vp[i];

            depth = Depth;

            if (p<0)    {dep=0.0; depth=0.0;} /* surface wave case */
            if (p==121) p=20;
            if (p==124) p=40;

            prj = 0.0;        
            
            if (p==phase) 
            {
              prj=pathprojectn1(dep*1000.,ths,phs,thr,phr,
                depth,stheta,sphi,rtheta,rphi,depres,disresr,disress);

              if ((fabs((D)(rtheta-thr))<0.002) && (fabs((D)(rphi-phr))<0.002) && 
                  (fabs((D)(stheta-ths))<0.002) && (fabs((D)(sphi-phs))<0.002)) 
              {
		            ifphase=1;
              }
                rwgt+=prj;
            } /* end of if (p==phase) */ 
          } /* for ( i = 0; i < nphase; i++) */
      
          if (ifphase==1) 
          {
            printf( " Phase = %4d | eff number of traces = %5.2f\n", packheader.phase, rwgt );
  	    
            if (rwgt==0.0) 
            {
              printf("stheta %f sphi %f rtheta %f rphi %f\n",stheta,sphi,rtheta,rphi);
              stop(" wcal_az : an error occured when lseeking ");            
            }

            rwgt  = sqrt(rwgt);
            rwgt  *= packheader.rmsd*ndata; /* one packet one datum */
            rwgt  /= sqrt((D)(ndata*trace.smplintv));

            packheader.weight=1./rwgt;
  	         
            if (packheader.weight>1e10) 
            {
              /* To avoid the model being influenced by few phases with a weight much higher than the average */
              packheader.weight=1e10;
              printf("Warning: weight for %d is large!\n",p);
            }
          } /* End of if (iphase==1) */ 
	  
          if (lseek(wpd,(off_t)(-phsize),1)<0)        stop("wcal_az : lseek error ");
          if (write(wpd,&packheader,phsize)!=phsize)  stop("wcal_az : writing error ");
        } /* End of if (dataselector_auto(&packheader,trace.delta,trace.netwk,trace.chnnl,trace.comp)) */

        if (lseek(wpd,(off_t)(4*ndata),1)<0) stop("wcal_az: lseek error n2");
      } /* End of while(1) for (lseek(wpd,trace.locatn,0)<0) condition */
    } /* End of if ( traceA.neffkreg != 0 ) */
    else 
    {
      if (lseek(wpd,trace.locatn,0)<0) 
      {
        stop("wcal_az : lseek error, could not search for trace");
      }
      while (1) 
      {

        if (!(nread=read(wpd,&packheader,phsize)))  break;
        if (nread!=phsize)                          stop("wcal_az : error reading wpd file");
        if (packheader.id!=trace.id)                break;
        packheader.weight = 0.;
        ndata             = packheader.ndata;
        if (lseek(wpd,(off_t)(-phsize),1)<0)        stop("wcal_az : lseek error 1");
        if (write(wpd,&packheader,phsize)!=phsize)  stop("wcal_az : writing error");
        if (lseek(wpd,(off_t)(4*ndata),1)<0)        stop("wcal_az : lseek error 2");
      } /* End while(1) for condition (lseek(wpd,trace.locatn,0)<0) */
    } /* End for Else in ( traceA.neffkreg != 0 ) condition */
  } /* End of while (nread=read(wph,&trace,thsize)) */
} /* End of main */

