/* Make an ASCII file to list the relevant phases crossing the target region.
It distinguishes between phases for the minor and major arc 
It reads an input file with a station list:
if flag=0 it will list only phases recorded at stations in the list
if flag=1 it will list only phases recorded at stations not in the list 
patch for autopicking bug added (Federica Nov 10, 2006) */
#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <assert.h>
#include <mpi.h>
#include "acal/const.h"
#include "acal/structwavket_A3z.h"
#include "acal/dimensiony.h"
#define EPS 1.0e-5

extern void   listing_hdr_az(char *event, int flag);
extern void   stop();
extern void   midpoint();
extern int    fgetrecord();
extern int    dataselector_auto();
extern float  ***farray3();
extern char   **carray2();
extern void   free_farray3();
int           check_station();
extern void 	exit_error();

int     main(int argc, char **argv)
{

        int rank, size, flag;

	if (argc < 1) 
	{
		printf( "Usage: %s flag for check_station\n"
			"      flag = 0 -> from the list considered\n"
			"      flag = 1 -> not from the list consiered\n", argv[0]);
		exit_error(stdout, "Error [%s]: check usage\n", __func__);
  	} else
	{
		sscanf(argv[1], "%i", &flag);
	}
        
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
                FILE *f_evt;
		if ((f_evt = fopen("lu_events", "r")) == NULL)
		{
			exit_error(stdout, "Error: cannot open lu_events - %s\n", strerror(errno));
		}
		else 
		{
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
			listing_hdr_az( event, flag );
                        wt1 = MPI_Wtime();
                        printf("Rank [%2i]: Listing of event %8s took %f s\n", rank, event, wt1 - wt0);

                        /* get new event */
                        MPI_Send(&msg, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                        MPI_Recv(event, 9, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                }

                printf("Rank [%2i] is done\n", rank);

        }

        MPI_Finalize();
}

void
listing_hdr_az(event,flag)
     char *event;
     int  flag;
{
     char    string[120],**stat, evt[120], *file_out_name, fl_for_out[120];
     int     wph,wpd,wpN,nwp,nread,nstat,ifstat;
     int     size;
     FILE    *file, *file_out; //file_out is for writing out
     float   ths1,phs1,dep1,thr1,phr1,dis1,smp1,t01,rmsd1,rmsr1,wgt1;
     float   ths2,phs2,dep2,thr2,phr2,dis2,smp2,t02,rmsd2,rmsr2,wgt2;
     float   rmss1,rmss2;
     float   thm1,thm2,thm,phm1,phm2,phm;
     float   ***B;
     int     i,JJ,ndata1,ndata2,iwp1,iwp2,minor,major,first_orbit,ksmax;
     int     ifcode,ifmoment,ifdepth,ifsloctn,ifrloctn,ifdist,ifsmplint,ifcomp;
     int     ifnetwk,ifchnnl,ifphaseid,ifndata,ifwwindow,ifgwindow,ifpwindow;
     int     ifrmsd,ifrmsr,ifweight,ifrmss,onlysrc,ifstn;
     double  dum,diff;

     title_st        title;
     tracehdrHH_st    trace;
     tracehdrAreg_st traceA;
     packhdrbr_st    packhdr;
     sphA1zreg_st    Aknot[MAXKNOT];

     /* Naming output file from flag and event name */
     sprintf( fl_for_out, "%i", flag);
     strncpy( evt, event, 8 );
     file_out_name  = strcat(evt,"_flag_");
     file_out_name  = strcat(file_out_name,fl_for_out);
     file_out_name  = strcat(file_out_name,"_list");
     file_out 	    = fopen(file_out_name, "w");
    
     /* opening wp files */ 
     strncpy( string, event, 8 );
     strncpy( &string[8], ".wpH\0", 5 );

     wph            = open( string, 0 );
     string[11]     = 'd';
     wpd            = open(string,0); 
     string[11]     = 'N';
     wpN            = open(string,0); 
       
     file           = fopen("lu_filter.par","r");
     fgetrecord(file,string);sscanf(string,"%g %g",&ths2,&ths1);
     fgetrecord(file,string);sscanf(string,"%g %g",&phs1,&phs2);
     fgetrecord(file,string);sscanf(string,"%g %g",&dep1,&dep2);
     fgetrecord(file,string);sscanf(string,"%g %g",&thr2,&thr1);
     fgetrecord(file,string);sscanf(string,"%g %g",&phr1,&phr2);
     fgetrecord(file,string);sscanf(string,"%g %g",&dis1,&dis2);
     fgetrecord(file,string);sscanf(string,"%g %g",&smp1,&smp2);
     fgetrecord(file,string);sscanf(string,"%g %g",&t01,&t02);
     fgetrecord(file,string);sscanf(string,"%g %g",&rmsd1,&rmsd2);
     fgetrecord(file,string);sscanf(string,"%g %g",&rmsr1,&rmsr2);
     fgetrecord(file,string);sscanf(string,"%g %g",&rmss1,&rmss2);
     fgetrecord(file,string);sscanf(string,"%g %g",&wgt1,&wgt2);
     fgetrecord(file,string);sscanf(string,"%d %d",&ndata1,&ndata2);
     fgetrecord(file,string);sscanf(string,"%d %d",&iwp1,&iwp2);
     fgetrecord(file,string);sscanf(string,"%g %g",&thm2,&thm1);
     fgetrecord(file,string);sscanf(string,"%g %g",&phm1,&phm2);
     fclose(file);
       
     ths1=(90.-ths1)*CONV; ths2=(90.-ths2)*CONV;
     thr1=(90.-thr1)*CONV; thr2=(90.-thr2)*CONV;
     thm1=(90.-thm1)*CONV; thm2=(90.-thm2)*CONV;
     phs1*=CONV; phs2*=CONV; phr1*=CONV; phr2*=CONV; 
     phm1*=CONV; phm2*=CONV; 
     dis1*=CONV; dis2*=CONV; 
         
     onlysrc=0;

     file=fopen("lu_list.par","r");
     fgetrecord(file,string);sscanf(string,"%d",&ifcode);
     fgetrecord(file,string);sscanf(string,"%d",&ifmoment);
     fgetrecord(file,string);sscanf(string,"%d",&ifdepth);
     fgetrecord(file,string);sscanf(string,"%d",&ifsloctn);
     fgetrecord(file,string);sscanf(string,"%d",&ifrloctn);
     fgetrecord(file,string);sscanf(string,"%d",&ifdist);
     fgetrecord(file,string);sscanf(string,"%d",&ifsmplint);
     fgetrecord(file,string);sscanf(string,"%d",&ifcomp);
     fgetrecord(file,string);sscanf(string,"%d",&ifnetwk);
     fgetrecord(file,string);sscanf(string,"%d",&ifchnnl);
     fgetrecord(file,string);sscanf(string,"%d",&ifphaseid);
     fgetrecord(file,string);sscanf(string,"%d",&ifndata);
     fgetrecord(file,string);sscanf(string,"%d",&ifwwindow);
     fgetrecord(file,string);sscanf(string,"%d",&ifgwindow);
     fgetrecord(file,string);sscanf(string,"%d",&ifpwindow);
     fgetrecord(file,string);sscanf(string,"%d",&ifrmsd);
     fgetrecord(file,string);sscanf(string,"%d",&ifrmsr);
     fgetrecord(file,string);sscanf(string,"%d",&ifrmss);
     fgetrecord(file,string);sscanf(string,"%d",&ifweight);
     fgetrecord(file,string);sscanf(string,"%d",&ifstn);
     fclose(file);

     /* Read station list */
     if( (file=fopen("lu_station.par","r"))==0 ) {
          stop("listing_hdr_az: cannot open lu_station.par");
     }

     nstat     = 0;
     while(fgetrecord(file,string)>0) nstat++;
     rewind(file);

     stat      = carray2(0,nstat-1,0,3);

     for ( i = 0; i < nstat; i++ ) {
          fgetrecord(file,string);
          sscanf(string,"%s",&stat[i][0]);
     }
     
     fclose(file);
       
     read( wph, &title, sizeof(title) );
     if (title.phi < 0.0) title.phi += 2.*PI;

     if (title.theta<ths1 || title.theta>ths2 || title.phi<phs1 
          || title.phi>phs2 || title.depth<dep1 || title.depth>dep2) {

          stop("listing_hdr_az: error");
     }

     if ( ifrloctn==0 && ifdist==0 && ifsmplint==0
          && ifcomp==0 && ifnetwk==0 && ifchnnl==0
          && ifphaseid==0 && ifndata==0 && ifwwindow==0
          && ifgwindow==0 && ifpwindow==0 && ifrmsd==0
          && ifrmsr==0 && ifrmss==0 && ifweight==0 ) onlysrc=1;
       
     if (onlysrc) {
          if (ifcode)    printf("%.8s ",title.event);
          if (ifmoment)  printf("%9.3e %9.3e %9.3e %9.3e %9.3e %9.3e ",
               title.moment.rr, title.moment.tt, title.moment.pp,
               title.moment.rt, title.moment.rp, title.moment.tp);

          if (ifdepth) printf("%5.1f ",title.depth/1000.);
          if (ifsloctn) printf("%6.3f %6.3f ",title.theta,title.phi);

          printf(" %6.3f ",title.dt);
          printf("\n");

          stop("listing_hdr_az: error");
     }
      
     while ( read(wph, &trace, sizeof(trace))==sizeof(trace) ) {
          lseek( wpN, trace.locatnA, 0 );  
          read( wpN, &traceA, sizeof(traceA) );
          minor = 0;
          major = 0;    
          
          if (traceA.neffkreg!=0) {

               B = farray3( 0, traceA.neffk-1, 0, SMAX*2, 0, 4 );

               for ( i = 0; i < traceA.neffk; i++ ) {

                    if ( read( wpN, &Aknot[i], sizeof(Aknot[i]))!=sizeof(Aknot[i]) ) {
                         printf("%.8s\n",title.event);
                         stop("listing_hdr_az: Error in loading wpN n1");
                    }

                    ksmax     = Aknot[i].smax;
                    if ( ksmax > SMAX ) stop("SMAX in include/structwavket_A3z is declared too small");     
                    if ( ksmax > 0) {
                         for( JJ = 1; JJ <= ksmax*2; JJ++ ) {
                              if ( read( wpN, &B[i][JJ][0], 5*sizeof(float)) != 5*sizeof(float) ) stop("listing_hdr_az: error in loading wpN 2");  
                         }   
                    }
          	
                    if (Aknot[i].reg == 1) {
                         if ( Aknot[i].rphi < 0 ) Aknot[i].rphi = Aknot[i].rphi + 2*PI;
                         if ( Aknot[i].rphi >= 0 && Aknot[i].rphi <= trace.delta) {
                              minor++;
                         } 
                         else if ( Aknot[i].rphi > trace.delta && Aknot[i].rphi < 2*PI ) {
                              major++;
                         } 
                         else {
                              printf("rphi %f, delta %f\n",Aknot[i].rphi,trace.delta);
                              stop(" listing_hdr_az: Error invalid values for knot position");
                         }
                    } /* End if (Aknot[i].reg == 1) */
               } /* for ( i = 0; i < traceA.neffk; i++ ) */

               //sscanf(argv[2],"%d",&flag);

               ifstat = check_station( nstat, trace, stat, flag );
               midpoint( title.theta, title.phi, trace.theta, trace.phi, &thm, &phm );
               if ( trace.theta   >=thr1  && trace.theta   <=thr2 &&
                    trace.phi     >=phr1 && trace.phi     <=phr2 &&
                    thm           >=thm1 && thm           <=thm2 &&
                    phm           >=phm1 && phm           <=phm2 &&
                    trace.delta   >=dis1 && trace.delta   <=dis2 &&
                    trace.smplintv>=smp1 && trace.smplintv<=smp2 &&
                    ifstat == 1 ) {

                    lseek(wpd,trace.locatn,0);
                    nwp = 0;
                    do {
                         if ( !( nread = read( wpd, &packhdr, sizeof(packhdr) ) ) ) break;
                         if ( nread != sizeof(packhdr) ) stop("listing_hdr_az: Error 1");
                         if ( packhdr.id != trace.id ) break;
                         size = (off_t)packhdr.ndata*sizeof(float);  
                         lseek( wpd, (off_t)size, 1 );
                              if (packhdr.t0   >=t01    && packhdr.t0    <=t02   &&
                                 packhdr.rmsd  >=rmsd1  && packhdr.rmsd  <=rmsd2 &&
                                 packhdr.rmsr  >=rmsr1  && packhdr.rmsr  <=rmsr2 &&
                                 packhdr.rmss  >=rmss1  && packhdr.rmss  <=rmss2 &&
                                 packhdr.weight>=wgt1   && packhdr.weight<=wgt2  &&
                                 packhdr.ndata >=ndata1 && packhdr.ndata <=ndata2)
     	     
                                   if (packhdr.phase>0 || packhdr.phase==-11 || packhdr.phase==-21 ||
     	                             packhdr.phase==-31 || packhdr.phase==-41 || packhdr.phase==-51 ||
     	                             packhdr.phase==-61) {
                                        
                                        first_orbit = 1;
                                   } 
                              else first_orbit = 0;

                              diff = modf( (packhdr.t0+title.dt)/trace.smplintv, &dum );

                              if ( dataselector_auto( &packhdr, trace.delta, trace.netwk, trace.chnnl, trace.comp ) && 
                                   (fabs(diff)<EPS || fabs(1.0-diff)<EPS) && 
                                   ( (first_orbit==1 && minor!=0) || (first_orbit==0 && major!=0) ) ) { 
                                   nwp++;
                                   if ( nwp > iwp2 ) break;
                                   if ( nwp >= iwp1 )
				   {
                                        if (ifcode)    fprintf(file_out, "%.8s ",title.event);
                                        if (ifmoment)  fprintf(file_out, "%9.3e %9.3e %9.3e %9.3e %9.3e %9.3e ",
                                             title.moment.rr, title.moment.tt, title.moment.pp, 
                                             title.moment.rt, title.moment.rp, title.moment.tp);

                                        if (ifdepth)   fprintf(file_out, "%5.1f ",         title.depth/1000.);
                                        if (ifsloctn)  fprintf(file_out, "%6.3f %6.3f ",   title.theta,title.phi);
                                        if (ifrloctn)  fprintf(file_out, "%6.3f %6.3f ",   trace.theta,trace.phi);
                                        if (ifdist)    fprintf(file_out, "%6.3f ",         trace.delta);
                                        if (ifsmplint) fprintf(file_out, "%4.0f ",         trace.smplintv);
                                        if (ifcomp)    fprintf(file_out, "%c ",            trace.comp);
                                        if (ifnetwk)   fprintf(file_out, "%.4s ",          trace.netwk);
                                        if (ifchnnl)   fprintf(file_out, "%.4s ",          trace.chnnl);
                                        if (ifphaseid) fprintf(file_out, "%3d ",           packhdr.phase);
                                        if (ifndata)   fprintf(file_out, "%3d ",           packhdr.ndata);
                                        if (ifwwindow) fprintf(file_out, "%9.3e %9.3e %9.3e %9.3e ",
                                             trace.w1,trace.w2,trace.w3,trace.w4);
                                        if (ifgwindow) fprintf(file_out, "%9.3e %9.3e ",   packhdr.gv1,packhdr.gv2);
                                        if (ifpwindow) fprintf(file_out, "%9.3e %9.3e ",   packhdr.pv1,packhdr.pv2);
                                        if (ifrmsd)    fprintf(file_out, "%9.3e ",         packhdr.rmsd);
                                        if (ifrmsr)    fprintf(file_out, "%5.2f ",         packhdr.rmsr);
                                        if (ifrmss)    fprintf(file_out, "%5.2f ",         packhdr.rmss);
                                        if (ifweight)  fprintf(file_out, "%9.3e ",         packhdr.weight);
                                        if (ifstn)     fprintf(file_out, "%.4s ",          trace.stn);
                                        fprintf(file_out,"\n");
                                   } 
                              }
                         } 
                    while (packhdr.id==trace.id);
               }

          free_farray3(B,0,traceA.neffk-1,0,2*SMAX,0);
          }/* End if (traceA.neffkreg!=0) */

     } /* end while ( read(wph, &trace, sizeof(trace))==sizeof(trace) ) */

     fclose(file_out);
} /* End main */

/******************************************************************************/
int check_station(nstat,trace,stat,flag)
/* flag=0 looking for stations in the list,
   flag=1 looking for stations not in the list 
   return=0 station should not be considered
   return=1 station should be considered */

	int		nstat,flag;
	char		**stat;
	tracehdrHH_st 	trace;
{
	int i;
  	char station[4];
	
	sscanf(trace.stn,"%s",station);
	for (i=0;i<nstat;i++) {
		if (station[0]==stat[i][0] && station[1]==stat[i][1]
		&& station[2]==stat[i][2] && station[3]==stat[i][3]) {

		if (flag==0) return(1);
  		else return(0);
		}
	}
	if (flag==0) return(0);
	else return(1);
}
