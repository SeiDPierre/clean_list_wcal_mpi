/* Put 0 weighting to all wavepackets */

#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <unistd.h>
#include <string.h>
#include <fcntl.h>
#include <sys/types.h>
#include "acal/const.h"
#include "acal/structwavket_A3z.h"

#define MAXPHS 80000

extern void 	stop();
extern void 	fgetrecord();

int	main(int argc, char **argv)
{

	int	wpd,wph;
	FILE	*list,*par;
	char 	string[120];
	int 	nread,thsize,phsize,ndata;

	title_st 	title;
	tracehdrHH_st 	trace;
	packhdrbr_st 	packheader;

	strcpy(string,argv[1]);
	strcpy(&string[8],".wpH\0");
	wph = open(string,0);
 	string[11]='d';
	wpd = open(string,2);
  
	read(wph,&title,sizeof(title));
  
	thsize = sizeof(trace);
	phsize = sizeof(packheader);
  
	while (nread=read(wph,&trace,thsize)) {
		if ( nread != thsize) stop( " Clean : Error Reading wph file " );
  		if ( lseek( wpd, trace.locatn, 0) <0 ) stop(" Clean : lseek error 0 ");
		do {
			if ( !( nread = read( wpd, &packheader, phsize )) ) break;
			if ( nread != phsize) stop(" Clean : error reading wpd file");
		
			ndata = packheader.ndata;
			packheader.weight = 0.;
			if ( lseek( wpd, (off_t)(-phsize), 1 ) < 0 ) stop(" Clean : lseek error 1 " );
			if ( write( wpd, &packheader, phsize) != phsize ) stop( " Clean : writing error " );
			if ( lseek( wpd, (off_t)(4*ndata), 1 ) <0 ) stop( " Clean : lseek error 2");
		} while (packheader.id==trace.id);
	}
}
