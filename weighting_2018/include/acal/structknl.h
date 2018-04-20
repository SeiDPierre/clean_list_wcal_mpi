                       
typedef struct { short n1,q1,l1,pad1; float grv1,w1,alpha1;
                 short n2,q2,l2,pad2; float grv2,w2,alpha2;
                 float dw;        /* 2*(w2-w1)/(w1+w2) */
                 float dw2_ell;   /* ell effect, dw^2 ; added 7/21/95*/
               } knlrechdr_st;
               
typedef struct { short n,q,l,pad; float gv,w,alp,ell;
               } sort_knlrechdr_st;
               
typedef struct { int rclength,hdrlen,nknl,seclength,maxs;
                 float maxvel,minvel; /* group vel */
                 float maxw;   /* in rad/sec */
                 float dwcut;  /* in percentage */
                 short pmaxu,pmaxl,ndisc,pad;
                 float ru1,ru2; /* radius range for upper layer, r1<r2*/
                 float rl1,rl2;
                 char junk[140];
               } knlhdr_st;

typedef struct { int rclength,hdrlen,nknl,seclength,maxs; /* */
                 float maxvel,minvel; /* group vel */
                 float maxw;   /* in rad/sec       */
                 float dwcut;  /* in percentage    */
                 short pmax,ndisc,pad;		  
                 float r1,r2; /* radius range for the kernels, r1<r2 */
               } knlhdrA1_st;
