#ifndef STRUCTBODINV_A3_ANISO
#define STRUCTBODINV_A3_ANISO

#define NDEP 48
#define MAXORBIT 6
#define PARAM 5		/* Max number of blocks */


typedef struct {float vcutoffr,vcutoffs;  /* cutoff rmsr and rmss*/
                char netwk[4],chnnl[4],comp[1];
                int n;       /* number of phases to select; -1 -> all */
                int *code;
                float *delta1,*delta2;} datparbr_st;

typedef struct {short nknotA1,nknotA2,damp; /* vertical & horizontal knots */
                short nknots,seafl; float damps;	/* seafloor */
                short nknotm,moho;  float dampm;	/* moho */
                short nknotd,d670;  float dampd;	/* 670 disc */
                short nknotc,cmb;   float dampc;	/* cmb */
                int dim; 
                char extr[8];	
               } mdlparA3_st; /*damp can be used to adjust amplitude*/

typedef struct {int nquant; /* number of physical parameters */
	char descript[PARAM]; /* A C F L N or Vp Vs.. */
	int ndisc;
	int nknotA2[PARAM],nknotA1[PARAM];float damp[PARAM]; /* vertical & horizontal knots */
	int nknots, seafl; float damps;	/* seafloor + flag */
	int nknotm, moho;  float dampm;	/* moho */
	int nknotd, d670;  float dampd;	/* 670 disc */
	int nknotc, cmb;   float dampc;	/* cmb */
	int dim;			/* dim : total # of parameters */
	} gmdlparA3_st;

typedef struct {float vcutoff;  /* cut-off value for rmsr */
                char netwk[4],chnnl[4];
                int n;     /* number of phases to select; -1 -> all */
                int *code;
                float *delta1,*delta2;} datpar_st;

typedef struct {short n,l;float R,S[2][2][NDEP-4];}
                  RS_st; /* ocean & first depth nodal skipped */

typedef struct {short n,l;float ha,W[2][NDEP];} temp_st;


typedef struct {int ndata;
                float pv1,pv2,gv1,gv2;
                float t0;
                float rmsd; /* rms of original data*/
                float rmsr; /* rms of residual by PREM + ell */
 		float rmss; /* rms of residual by PREM + ell */
                short phase,orbit;
                float *data;
                float weight;} wpinfobr_st;

typedef struct {int ndata;
                float pv1,pv2,gv1,gv2;
                float t0;
                float rmsd; /* rms of original data*/
                float rmsr; /* rms of residual by PREM + ell */
 		float rmss; /* rms of residual by PREM + ell */
                short phase,orbit;
                float *data;
                float *semdata;
                float weight;} wpinfoved_st;

typedef struct {float pv1,pv2,gv1,gv2,gv3,gv4;} window_st;

/* ptl : partial d (extracted from kernel file) */
/* ptlq : partial d for Q */

typedef struct {float rs[8],dw[MAXORBIT+2],da[MAXORBIT+2],ddw[MAXORBIT+2],e1w[MAXORBIT+2],e2w[MAXORBIT+2],*ptl;} mode_anis_st; 

typedef struct {float rs[8],dw[MAXORBIT+2],da[MAXORBIT+2],ddw[MAXORBIT+2],
e1w[MAXORBIT+2],e2w[MAXORBIT+2],*ptl,*ptlq;} modeq1_st; 

typedef struct {int nwp;        /* number of wave packets */
                float smplintv; /* sampling interval of data */
                float ths,phs,thr,phr;
                float delta;
                char  pad[4];
               } trchdr_st;

typedef struct {int phase;   /* phase code */
                int ndata;
                float rms;  /* rms of original data */
                float weight;
               } wphdr_st;

#endif
