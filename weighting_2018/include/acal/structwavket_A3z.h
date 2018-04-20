#define MAXKNOT    13000

typedef struct {int level,smax;  /*smax used in nact */
		int index;
                float theta,phi;   /* colat(0,pi), phi (0,2pi) */
                float rtheta,rphi; /*after rotation (theta == pi/2 ) */
                float pava[2][5];     /* 0 for great circle, 1 for minor arc
					0->4 for const,2psi,4psi terms */
               } sphA1z_st;

typedef struct {int level,smax;  /*smax used in nact */
		int reg;         /* 1 when knot in the region, otherwise 0 */
		int index;
		int indexreg;      /* index of all the knots in the region */
                float theta,phi;   /* colat(0,pi), phi (0,2pi) */
                float rtheta,rphi; /*after rotation (theta == pi/2 ) */
                float pava[2][5];     /* 0 for great circle, 1 for minor arc
					0->4 for const,2psi,4psi terms */
               } sphA1zreg_st;

typedef struct {char stn[4];
		float theta;       /* in radian */
		float phi;         /* in radian */
		float delta;       /* in radian */
		float az;	   /* in epicentral coordinates measured
				      from south counterclockwisely */
		short id;	   /* trace id */
		int neffk;	   /* number of effective knots */
		} tracehdrA_st;

typedef struct {char stn[4];
		float theta;       /* in radian */
		float phi;         /* in radian */
		float delta;       /* in radian */
		float az;	   /* in epicentral coordinates measured
				      from south counterclockwisely */
		short id;	   /* trace id */
		int neffk;	   /* number of effective knots */
		int neffkreg;	   /* number of effective knots inside the region of interest */
		} tracehdrAreg_st;


typedef struct {float rr,tt,pp,rt,rp,tp;} tensor_st;
    
typedef struct {char event[8];
                float theta;            /* in radian */
                float phi;              /* in radian */
                float depth;		/* in meters */
                tensor_st moment;       /* N.m */
                float dt;               /* refined cmt t0 relative to 
                                           Harvard t0 in sec */
                }    title_st;      /*48bytes*/ /* wph file header */
	        
typedef struct {short phase;	   /* code for P sS etc.*/
		short ndata;       /* number of data points in packet*/
		short id;          /* belong to which trace */
		char pad[2];	   /* dummy */
		float t0;          /* starting time from reftime */
		float gv1,gv2;     /* group velocity range in rad/s */
		float pv1,pv2;     /* phase velocity range in rad/s */
		float rmsd;        /* rms of data*/
		float rmsr;        /* rms residual/rmsd  prem+ell */
		float rmss;	   /* rms of synthetic */
		float weight;}       /* weights */     
	        packhdrbr_st;        /*44bytes*/

		
typedef struct {char stn[4];
                int locatn;       /* location in *.wpd file (bytes) */
		float theta;       /* in radian */
		float phi;         /* in radian */
		float delta;       /* in radian */
		float az;	   /* in epicentral coordinates measured
				      from south counterclockwisely */
		float dip;	   /* -pi/2 is vertical*/
		float smplintv;	   /* in second */
		float w1;	   /* cut-off omega at low end (rad/s)*/
		float w2;          /* corner  omega at low end (rad/s)*/
		float w3;          /* corner  omega at hi  end (rad/s)*/
		float w4;          /* cut-off omega at hi  end (rad/s)*/
		short id;	   /* trace id */
		char  reftime;	   /* c=cmt;b=bulletin;r=refined cmt*/
		char  comp;        /* 'Z','L','T','E','N' */
		char  netwk[4];    /* GEO, GSN, etc -- added 9/10/92 */
		char  chnnl[4];    /* VLP, LP, etc -- added 9/10/92 */
		char  extr[4];}   tracehdr_st;       /*64bytes*/
		/* wph file for each station */


typedef struct {char stn[4];
                int locatn;        /* location in *.wpd file (bytes) */
		float theta;       /* in radian */
		float phi;         /* in radian */
		float delta;       /* in radian */
		float az;	   /* in epicentral coordinates measured
				      from south counterclockwisely */
		float dip;	   /* -pi/2 is vertical*/
		float smplintv;	   /* in second */
		float w1;	   /* cut-off omega at low end (rad/s)*/
		float w2;          /* corner  omega at low end (rad/s)*/
		float w3;          /* corner  omega at hi  end (rad/s)*/
		float w4;          /* cut-off omega at hi  end (rad/s)*/
		short id;	   /* trace id */
		char  reftime;	   /* c=cmt;b=bulletin;r=refined cmt*/
		char  comp;        /* 'Z','L','T','E','N' */
		char  netwk[4];    /* GEO, GSN, etc -- added 9/10/92 */
		char  chnnl[4];    /* VLP, LP, etc -- added 9/10/92 */
		char  extr[4];    
		int locatnA;  /* location in *.wpA file (bytes) */
		}   tracehdrH_st;  

		/* wpH file for each station */

typedef struct {char stn[4];
                int locatn;        /* location in *.wpd file (bytes) */
                float theta;       /* in radian */
                float phi;         /* in radian */
                float delta;       /* in radian */
                float az;          /* in epicentral coordinates measured
                                      from south counterclockwisely */
                float dip;         /* -pi/2 is vertical*/
                float smplintv;    /* in second */
                float w1;          /* cut-off omega at low end (rad/s)*/
                float w2;          /* corner  omega at low end (rad/s)*/
                float w3;          /* corner  omega at hi  end (rad/s)*/
                float w4;          /* cut-off omega at hi  end (rad/s)*/
                short id;          /* trace id */
                char  reftime;     /* c=cmt;b=bulletin;r=refined cmt*/
                char  comp;        /* 'Z','L','T','E','N' */
                char  netwk[4];    /* GEO, GSN, etc -- added 9/10/92 */
                char  chnnl[4];    /* VLP, LP, etc -- added 9/10/92 */
                char  extr[4];
                long locatnA;  /* location in *.wpA file (bytes) */
                }   tracehdrHH_st;




/* seismic phase code: (Toroidal)
 * 11=P;
 * 20=Sdif; 21=S; 22=SS; 23=SSS; 24=S4; 25=ScS; 26=(ScS)2; 27=(ScS)3;
 * 40=sSdif;41=s; 42=sS; 43=sSS; 44=s(S3);45=sScS;46=s(ScS)2;47=s(ScS)3;
 *-11= Rayleigh Wave 1st orbit; -12= 2nd orbit; and so on; 
 *-21= Love Wave first orbit; -22=second orbit ; an so on;
 *-31= Rayleigh X phase 1st orbit, -32= XR2; -33=XR3
 *-41= Love X phase 1st orbit, -42=XG2; -43=XG3 ; -44=XG4 
 *-51=R1+X1; -52=R1+X2; -53=R2+X2; -54=R2+X3; -55=R3+X4
 *-61=L1+GX1; -62=L1+GX2; -63=L2+GX2; -64=G2+GX3; -65=G3+GX4
 *-71=RX1+RX2; 
 *-81=GX1+GX2;
 * 121=S+ScS; 122=SS+(ScS)2; 123=SSS+(ScS)2; 124=sS+sScS; 125=SS+sS;
 * 126=SSS+sSS; 127=S4+sS3; 128=SSS+sScS2; 129=SS+sScS; 
 * 140=sSdif+sScS; 141=ScS2+sSS; 142=ScS+sS; 143=S4+ScS3;144=SS+ScS;
 * 145=sS3+(ScS)2; 146=s(ScS)2+sS3; 147=s(ScS)2+sSS; 148=S4+(ScS)2;
 * 149=S4+sScS2; 160=Sdif+sS;
 * 1000=unknown phase;
 *
 * Spheroidal :	2--:Shallow phase; 3--:Deep phase; 4--:2 phases; 
 *		6--:3 phases; 9-- quadruple; 10-- quintuple;
 *
 * Shallow single phases
 * 200=Sdiff; 201=S; 202=SS; 203=; 204=ScS; 205=PKiKS; 206=PKiKP; 207=PcP; 
 * 208=SP; 209=P; 210=PcS; 211=Pdiff; 212=PP; 213=ScP; 214=PS; 215=S4;
 * 216=PKKPab; 217=SKKPbc; 218=P'P'ab; 219=ScS2; 220=ScS3; 221=SKKSac;
 * 222=PKSab; 223=PKKSbc; 224=P'P'bc; 225=SKiKP; 226=PKKPdf; 227=SKKPdf;
 * 228=PKKSdf; 229=SKKSdf; 230=P'P'df; 231=S'S'df; 232=SKPab;
 * 233=SKPdf; 234=PKSab; 235=PKSbc; 236=PKSdf; 237=SKSac; 238=SKSdf;
 * 239=PKPab; 240=PKPdf; 241=PKPbc; 242=S'S'ac; 243=SKPbc; 244=PKKPbc;
 * 245=X1; 246=SSS; 
 *
 * Deep single phases
 * 302=sSdiff; 303=s; 304=sS; 305=sSS; 306=sS3; 307=sScS; 308=sScS2;
 * 309=pPKiKP; 310=sPKiKP; 311=pP; 312=sP; 313=pS; 314=pPdiff; 
 * 315=sPdiff; 316=pPKPdf; 317=sPKPdf; 318=; 319=pSdiff;
 * 320=pSKSac; 321=sSKSac; 322=sScS3;
 *
 * Paired phases
 * 401=PcP+P; 402=S+ScS; 403=SP+PS; 404=pP+sP; 405=sS+pS; 406=ScS+SKSac;
 * 407=sPKiKP+pPKiKP; 408=pPdiff+sPdiff; 409=sSdiff+pSdiff; 410=SKSac+SKSdf;
 * 411=PKPab+PKPdf; 412=PcS+ScP; 413=PP+pP; 414=ScP+sS; 415=ScS+PKiKP;
 * 416=PP+sP; 417=SKKPdf+PKKSdf; 418=PKiKP+PP; 419=Sdiff+sSKSac;
 * 420=pSKSac+sSKSac; 421=PcP+PP; 422=sS+ScS; 423=Sdiff+SKSac;
 * 424=ScS+pPKiKP; 425=ScS+sPKiKP; 426=PcP+pP; 427=PcP+sP; 428=ScS+SKiKP;
 * 429=PKiKS+SKiKP; 430=SS+sS; 431=P+pP; 432=PcS+sS; 433=PS+PKiKP;
 * 434=SP+PKiKP; 435=Sdiff+PKKPdf; 436=S+ScP; 437=S+PcS; 438=PKiKP+ScS;
 * 439=S+PKiKP; 440=S+SKiKP; 441=SS+SKiKP; 442=Pdiff+pPdiff; 443=sS+pPKiKP;
 * 444=sS+sPKiKP; 445=PS+PKKPdf; 446=SP+PKKPdf; 447=PKKPab+PKKPdf;
 * 448=SKiKP+sPKiKP; 449=sS+pSKSac; 450=SS+sSKSac; 451=Sdiff+SKKSac;
 * 452=SKKSac+SKSdf; 453=SP+SKSac; 454=SP+sSKSac; 455=Sdiff+PKKPab;
 * 456=PKiKP+PKPdf; 457=SS+ScS; 458=P+sP; 459=S+PcP; 460=P+PP; 461=PKiKP+pPKiKP;
 * 462=SKKSac+SKSac; 463=SS+PcS; 464=SS+PKiKP; 465=S+sS; 466=SS+pPKiKP;
 * 467=SP+pSKSac; 468=PS+sS; 469=PS+sSdiff; 470=PKKPab+PKKPbc; 
 * 471=SKKPbc+PKKSbc; 472=SP+pSdiff; 473=PP+SKiKP; 474=PP+pPKiKP; 
 * 475=SKKSac+pSKSac; 476=Sdiff+sSdiff; 477=Sdiff+PKiKP; 478=PP+sPKiKP;
 * 479=PS+SKKPbc; 480=SKSac+pSKSac; 481=SKiKP+pPKiKP; 482=PKiKP+sS;
 * 483=PKiKP+SKiKP; 484=SKiKP+sS; 485=S+pPKiKP; 486=S+sPKiKP; 487=sS+sSKSac;
 * 488=SS+PKKPbc; 489=SS+SKKPbc; 490=pP+sPdiff; 491=Sdiff+pSKSac; 
 * 492=PKiKP+pPdiff; 493=PKiKP+PKPab; 494=PcS+PP; 495=S+X1; 496=PKKPbc+SSS
 * 497=SS+SSS; 498=SSS+sS3; 499=ScS2+sScS2; 501=ScS3+sScS3; 502=S4+SSS; 
 * 503=S4+PKKPbc; 504=S4+ScS2; 505=PKiKP+PcP; 506=S4+sScS2; 507=SSS+sScS2;
 * 
 * Triple phases
 * 601=S+SKSac+ScS; 602=sP+pP+PP; 603=pPKiKP+sPKiKP+ScS; 604=sP+pP+PcP
 * 605=sS+ScS+SKiKP; 606=SP+PS+PKiKP; 607=S+PcS+ScP; 608=PcS+ScP+sS;
 * 609=Pdiff*sPdiff+sSdiff; 610=sS+pPKiKP+sPKiKP; 611=SP+PS+PKKPdf;
 * 612=S+ScS+SKKSac; 613=P+pP+sP; 614=P+PP+pP; 615=PcP+P+PP; 616=SS+ScS+PKiKP;
 * 617=PcP+P+pP; 618=PKiKP+pPKiKP+sPKiKP; 619=PP+PKiKP+pPKiKP; 
 * 620=Sdiff+SKKSac+SKSac; 621=SP+PS+PKKPab; 622=SKKSac+SKSac+pSKSac;
 * 623=SS+PcS+ScP; 624=PcP+PP+pP; 625=ScS+PKiKP+pPKiKP; 626=ScS+PKiKP+sS;
 * 627=SS+pPKiKP+sPKiKP; 628=S+PKiKP+pPKiKP; 629=S+ScS+PKiKP; 630=sS+pS+pSKSac;
 * 631=PS+sS+pS; 632=PS+sSdiff+pSdiff; 633=Sdiff+SKKSac+pSKSac;
 * 634=SP+PS+sSdiff; 635=SP+pSKSac+sSKSac; 636=SP+pSdiff+pSKSac; 
 * 637=PP+SKiKP+pPKiKP; 638=PP+pPKiKP+sPKiKP; 639=PKiKP+PP+sPKiKP;
 * 640=SKKS+pSKSac+sSKSac; 641=Sdiff+pSdiff+sSdiff; 642=Sdiff+PKiKP+SKKSac;
 * 643=Sdiff+PKKPab+SKKSac; 644=Sdiff+PKKPab+PKKPbc; 645=PS+SKKPbc+sSdiff;
 * 646=SP+PS+SKKPbc; 647=SKSac+pSKSac+sSKSac; 648=SS+PcS+sS;
 * 649=SKiKP+pPKiKP+sPKiKP; 650=SS+ScS+sS; 651=SS+SKiKP+sPKiKP; 652=SS+SKiKP+sS;
 * 653=SS+ScS+SKiKP; 654=S+pPKiKP+sPKiKP; 655=SS+sS+pSKSac; 
 * 656=sS+pSKSac+sSKSac; 657=PKKPab+PKKPbc+sSdiff; 658=SS+SKKPbc+PKKSbc;
 * 659=SP+PS+SKSac; 660=Sdiff+pSKSac+sSKSac; 661=PKiKP+pPdiff+sPdiff;
 * 662=PKiKP+PKPab+PKPbc; 663=S+PcP+sS; 664=PcS+PP+ScP; 665=S+PKiKP+sS; 
 * 666=S+sS+sPKiKP; 667=SP+PS+PKKPbc; 668=SS+SKiKP+SSS; 669=SKKPbc+PKKSbc+SSS;
 * 670=S4+SSS+sS3 ; 671=S4+PKKPbc+SSS; 672=S4+PKKPbc+sS3; 673=PKiKP+PcP+P;
 *
 * Quadruple phases
 * 901=S+ScS+SKSac+SKKSac; 902=P+PP+pP+sP; 903=PcP+P+pP+sP; 904=PcP+P+PP+pP;
 * 905=PP+PKiKP+pPKiKP+sPKiKP; 906=SKKSac+SKSac+pSKSac+sSKSac; 
 * 907=SS+PcS+ScP+sS; 908=PcP+PP+pP+sP; 909=SS+ScS+PKiKP+sS; 	
 * 910=ScS+PKiKP+pPKiKP+sPKiKP; 911=SS+ScS+PKiKP+pPKiKP; 
 * 912=S+PKiKP+pPKiKP+sPKiKP; 913=sS+pS+pSKSac+sSKSac; 914=PS+sS+pS+sSKSac;
 * 915=SP+PS+sSdiff+pSdiff; 916=SP+PS+PKKPab+PKKPbc; 
 * 917=SP+pSdiff+pSKSac+sSKSac; 918=Sdiff+SKKSac+pSKSac+sSKSac;
 * 919=PP+SKiKP+pPKiKP+sPKiKP; 920=Sdiff+PKKPab+SKKSac+PKKPbc;
 * 921=Sdiff+PKKPab+PKKPbc+pSKSac; 922=SP+PS+SKKPbc+PKKSbc; 
 * 923=SS+ScS+sS+pPKiKP; 924=ScS+PKiKP+sS+pPKiKP; 925=SS+SKiKP+sS+pPKiKP;
 * 926=SS+SKiKP+sS+sPKiKP; 927=SS+sS+pSKSac+sSKSac; 928=S+ScS+SKKSac+SKiKP;
 * 929=SP+PS+SKSac+pSdiff; 930=S+PcS+ScP+sS; 931=S+PKiKP+sS+pPKiKP;
 * 932=S+sS+sPKiKP+pS; 933=S+PcS+ScP+X1; 934=S4+PKKPbc+SSS+sS3; 
 * 935=PKiKP+PcP+P+PP;
 *
 * Quintuple phases 
 * 1001=PcP+P+PP+pP+sP; 1002=SS+ScS+PKiKP+pPKiKP+sPKiKP;
 * 1003=Sdiff+PKKPab+PKKPbc+pSKSac+sSKSac; 1004=SS+SKiKP+sS+pPKiKP+sPKiKP;
 * 1005=S+ScS+SKKSac+SKiKP+SKSac; 1006=S+PKiKP+sS+pPKiKP+sPKiKP; 
 * 1007=PKiKP+PcP+P+PP+SKiKP;
 *
 * 0    = no phases in wavepacket
 * 400  = 1 unknown phase in wavepacket
 * 500  = 2 unknown phases in wavepacket
 * 800  = 3 unknown phases in wavepacket
 * 999  = 4 unknown phases in wavepacket
 * 1000 = 5 unknown phases in wavepacket
 * 1100 = more than 5 phases in wavepacket
 * 
 * 
 */
