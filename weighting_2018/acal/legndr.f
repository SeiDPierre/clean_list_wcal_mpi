      SUBROUTINE LEGNDR(THETA,L,M,X,XP,XCOSEC)
      DIMENSION X(2),XP(2),XCOSEC(2)
      REAL*8 CT,ST,FCT,COT,FPI,BOEING,X1,X2,X3,F1,F2,XM,TH
      REAL*8 SMALL,SUM,COMPAR,DFLOAT
      DATA FPI/12.56637062D0/,BOEING/.70710678D0/
      DFLOAT(I)=FLOAT(I)
      SUM=0.D0
      TH=THETA
      CT=DCOS(TH)
      ST=DSIN(TH)
      LP1=L+1
      FCT=DSQRT(FLOAT(2*L+1)/FPI)
      SFL3=SQRT(FLOAT(L*(L+1)))
      COMPAR=DFLOAT(2*L+1)/FPI
      DSFL3=SFL3
      SMALL=1.D-16*COMPAR
      COT=CT/ST
      MP1=M+1
      DO 1 I=1,MP1
      X(I)=0.
      XCOSEC(I)=0.
1      XP(I)=0.
      IF(L.GT.1.AND.ABS(THETA).GT.1.E-5)GO TO 3
      IF(L.GT.1)GO TO 3
      X(1)=FCT
      IF(L.EQ.0)RETURN
      X(1)=CT*FCT
      X(2)=-ST*FCT/DSFL3
      XP(1)=-ST*FCT
      XP(2)=-CT*FCT*DSFL3/2
      IF(ABS(THETA).LT.1.E-5)XCOSEC(2)=XP(2)
      IF(ABS(THETA).GE.1.E-5)XCOSEC(2)=X(2)/ST
      RETURN
3      X1=1.D0
      X2=CT
      DO 4 I=2,L
      X3=(DFLOAT(2*I-1)*CT*X2-DFLOAT(I-1)*X1)/DFLOAT(I)
      X1=X2
4      X2=X3
      COT=CT/ST
      COSEC=1./ST
      X3=X2*FCT
      X2=DFLOAT(L)*(X1-CT*X2)*FCT/ST
      X(1)=X3
      X(2)=X2
      XP(1)=-X2
      SUM=X3*X3
      XP(2)=FLOAT(L*(L+1))*X3-COT*X2
      X(2)=-X(2)/SFL3
      XP(2)=-XP(2)/SFL3
      XCOSEC(2)=X(2)*COSEC
      SUM=SUM+2.D0*X(2)*X(2)
      IF(SUM-COMPAR.GT.SMALL)RETURN
      X1=X3
      X2=-X2/SQRT(FLOAT(L*(L+1)))
      DO 5 I=3,MP1
      K=I-1
      F1=DSQRT(DFLOAT((L+I-1)*(L-I+2)))
      F2=DSQRT(DFLOAT((L+I-2)*(L-I+3)))
      XM=K
      X3=-(2.D0*COT*(XM-1.D0)*X2+F2*X1)/F1
      SUM=SUM+2.D0*X3*X3
      IF(SUM-COMPAR.GT.SMALL.AND.I.NE.LP1)RETURN
      X(I)=X3
      X1=X2
      XCOSEC(I)=X(I)*COSEC
      XP(I)=-(F1*X2+XM*COT*X3)
5      X2=X3
      RETURN
      END
