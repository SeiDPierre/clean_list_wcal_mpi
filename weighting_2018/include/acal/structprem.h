typedef struct { 
float w,q,va,ha,grv,ell,a6,a7,a8,a9,a10,rot;
/* w   : frequency of the mode
   q   : 1/Q 
   va  : U(r=seafloor)
   ha  : V(r=seafloor) spheroidal
         W(r=seafloor) toroidal
   grv : group velocity (rad/s)
   ell : ellipticity
   rot : rotation parameter (splitting due to rot)
*/
float aker[10],u1,u1p,v1,v1p,u2,u2p,v2,v2p;
/*  u1  : U(r1) 
    u1p : u'(r1)
    vi  : V or W 
*/
float r1,r2;	/* r1 : r under source; r2 : r above source */
short n,l;	/* branch, angular order */
char  typ;  	/* 'S' or 'T' */
} sprem_st;

/* Old prem file format */
typedef struct { 
int nmode;			/* # of modes */
int lnts,lntt;			/* max angular order for S & T */
float r[48];	  		/* normalized radius 48  */
int indsfr[330],kntsfr[330];		
int indtor[300],knttor[300];	
} spremhdr_st;

/* New prem file format */
typedef struct { 
int nmode;			/* # of modes */
int lnts,lntt;			/* max angular order for S & T */
int nlayers;			/* number of layers in truncated model */
int nrecheadsph,nrecheadtor;
int total;			/* number of layers in non-truncated model */
} header_head_st;

typedef struct { 
header_head_st head;
float rad[166];	  		/* normalized radius 166  */
int indsfr[2048],kntsfr[2048];	/* index & knot? spheroidal modes */	
int indtor[2048],knttor[2048];	/* idem toroidal */
} spremhdr_750_st;
           
