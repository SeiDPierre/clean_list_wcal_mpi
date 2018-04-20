#ifndef STRUCTSYNTH_H
#define STRUCTSYNTH_H

typedef struct {float w1,w2,w3,w4;
               }              band_st;

typedef struct {float moment[6];       /* Newton.meters */
		float theta, phi;      /* in radians */
		float dep;             /* in meters */
	       }source_st;
	       
typedef struct {float theta,phi;       /* in radians */
               }              recvlctn_st;
               
typedef struct {float t0;       /* relative to event time in seconds */
		float wlowcut,wlowcorner;
		float whghcut,whghcorner;
		char  modetype; /* 'S', 'T', or 'B' (both) for N and/or
		                    E components*/
		char  whichbrnch;       /* 'f'=fundament only; 'a'=all;
		                           'o'=overtone only */
		char  datatype; /* 'a'=accel.; 'c'=counts; etc. */
		char  com;      /* 'Z','L','T','N','E' */
		float pv1,pv2;  /* phase velocity range in rad/s */
		float gv1,gv2;  /* group velocity range in rad/s */
		float dt;       /* sampling interval in seconds */
		int  ndata;    /* number of data points */
		int orbnum;	/* orbit number */
	       }              syntrace_st;

/* XD's synth/struct.h does not have syntrace_st.orbnum */
#endif
