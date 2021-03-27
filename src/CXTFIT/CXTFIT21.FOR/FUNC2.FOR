C    ------------------------------------------------------------------
      REAL*8 FUNCTION DBEXP(A)
C
C     PURPOSE:  TO CALCULATE DEXP(A)
C
      IMPLICIT REAL*8(A-H,O-Z)
      DBEXP=0.D00
      IF(A.LT.-100) RETURN
      DBEXP=DEXP(A)
      RETURN
      END
C
C    ------------------------------------------------------------------
      REAL*8 FUNCTION EXF(A,B)
C
C     PURPOSE:  TO CALCULATE EXP(A) ERFC(B)
C
      IMPLICIT REAL*8(A-H,O-Z)
      EXF=0.D00
      IF((DABS(A).GT.170.).AND.(B.LE.0.)) RETURN
      C=A-B*B
      IF((DABS(C).GT.170.).AND.(B.GE.0.)) RETURN
      IF(C.LT.-170.) GO TO 3
      X=DABS(B)
      IF(X.GT.3.0) GO TO 1
      T=1./(1.+.3275911*X)
      Y=T*(.2548296-T*(.2844967-T*(1.421414-T*(1.453152-1.061405*T))))
      GO TO 2
    1 Y=.5641896/(X+.5/(X+1./(X+1.5/(X+2./(X+2.5/X+1.)))))
    2 EXF=Y*DEXP(C)
    3 IF(B.LT.0.0) EXF=2.*DEXP(A)-EXF
      RETURN
      END
C
C     ------------------------------------------------------------------
      REAL*8 FUNCTION EXPBI0(X,Z)
      IMPLICIT REAL*8(A-Z)
C
C     returns  EXP(Z)*Io(X)  for any real X and Z.
C
      DATA P1,P2,P3,P4,P5,P6,P7/1.0D0,3.5156229D0,3.0899424D0,
     $     1.2067492D0,0.2659732D0,0.360768D-1,0.45813D-2/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1,
     $     0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1,
     $     0.2635537D-1,-0.1647633D-1,0.392377D-2/
C
      IF (DABS(X).LT.3.75) THEN
        Y=(X/3.75)**2
        EXPBI0=DBEXP(Z)*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
        AX=DABS(X)
        Y=3.75/AX
       EXPBI0=DBEXP(AX+Z)/DSQRT(AX)*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*
     $         (Q6+Y*(Q7+Y*(Q8+Y*Q9))))))))
      ENDIF
      RETURN
      END
C     ----------------------------------------------------------
      REAL*8 FUNCTION EXPBI1(X,Z)
      IMPLICIT REAL*8(A-Z)
C
C     returns  EXP(Z)*I1(X)  for any real X and Z.
C
      DATA P1,P2,P3,P4,P5,P6,P7/0.5D0,0.87890594D0,0.51498869D0,
     $     0.15084984D0,0.2658733D-1,0.301532D-2,0.32411D-3/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,-0.3988024D-1,
     $     -0.3662018D-2,0.163801D-2,-0.1031555D-1,0.2282967D-1,
     $     -0.2895312D-1,0.1787654D-1,-0.420059D-2/
C
      IF (DABS(X).LT.3.75) THEN
        Y=(X/3.75)**2
        EXPBI1=DBEXP(Z)*X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
        AX=DABS(X)
        Y=3.75/AX
        EXPBI1=DBEXP(AX+Z)/DSQRT(AX)*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*
     $         (Q6+Y*(Q7+Y*(Q8+Y*Q9))))))))
        IF(X.LT.0.) EXPBI1=-EXPBI1
      ENDIF
      RETURN
      END
C
C--------------------------------------------------------------
      REAL*8 FUNCTION PHI1(TAU)
C
C     PURPOSE: TO CALCULATE SERIES IN EXPONENTIAL BVP
C               FOR EQUILLIBRIUM PHASE
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON MODC,MM,A,B,DA,CX,P,BETA,BETR,OMEGA,R,DMU1,DMU2,MODE,
     & MCON,T,Z,STOPER,LEVEL,ICHEB
      COMMON/BOUN/MODB,NPULSE,TPULSE(50),PULSE(50),MASSST
      PARAMETER(ZERO=0.0D+00,ONE=1.0D+00)
C
      OMOM=A*B
      C=TPULSE(1)-B
C
      PHI1=ZERO
      FN=ONE
      DO 12 N=1,25
       SUMK=ZERO
       K=1
       CK=1.0
       FNK=FN
       FN=FN*N
    9  CK=-CK
       SUMK=SUMK+CK*(T-TAU)**(N-K)/(C**K)/FNK
       IF (N.EQ.K) GOTO 12
       FNK=FNK/(N-K)
       K=K+1
       GOTO 9
   12  PHI1=PHI1+SUMK*((OMOM*TAU)**N)/FN
C
       RETURN
       END
C------------------------------------------------------------------
      REAL*8 FUNCTION PHI2(TAU)
C
C     PURPOSE: TO CALCULATE SERIES IN EXPONENTIAL BVP
C               FOR NONEQUILLIBRIUM PHASE
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON MODC,MM,A,B,DA,CX,P,BETA,BETR,OMEGA,R,DMU1,DMU2,MODE,
     & MCON,T,Z,STOPER,LEVEL,ICHEB
      COMMON/BOUN/MODB,NPULSE,TPULSE(50),PULSE(50),MASSST
      PARAMETER(ZERO=0.0D+00,ONE=1.0D+00)
C
C
      OMOM=A*B
      C=TPULSE(1)-B
C
      PHI2=ZERO
      FN=ONE
      DO 12 N=1,25
       SUMK=ZERO
       K=1
       CK=1.0
       FNK=FN
       FN=FN*N
    9  CK=-CK
       SUMK=SUMK+CK*(T-TAU)**(N-K+1)/(C**K)/FNK/(N-K+1)
       IF (N.EQ.K) GOTO 12
       FNK=FNK/(N-K)
       K=K+1
       GOTO 9
   12  PHI2=PHI2+SUMK*((OMOM*TAU)**N)/FN
C
       RETURN
       END
C     ------------------------------------------------------------------
      REAL*8 FUNCTION GOLD(X,Y)
C
C     PURPOSE: TO CALCULATE GOLDSTEIN'S J-FUNCTION J(X,Y)
C      introduce RT  (8/24/94  Nobuo Toride)
      IMPLICIT REAL*8(A-H,O-Z)
      GOLD=0.0
      BF=0.0
      E=2.*DSQRT(DMAX1(1.D-35,X*Y))
      Z=X+Y-E
      IF(Z.GT.17.) GO TO 8
      IF(E.GT.1.D-15) GO TO 2
      GOLD=DEXP(-X)
      RETURN
    2 A=DMAX1(X,Y)
      B=DMIN1(X,Y)
      RT=11.+2.*B+0.3*A
      IF(RT.GT.25) GO TO 6
      NT=INT(RT)
      I=0
      IF(X.LT.Y) I=1
      GXY=1.+I*(B-1.)
      GXYO=GXY
      GX=1.0
      GY=GXY
      GZ=1.0
      DO 4 K=1,NT
      GX=GX*A/K
      GY=GY*B/(K+I)
      GZ=GZ+GX
      GXY=GXY+GY*GZ
      IF((GXY-GXYO)/GXY.LT.1.D-08) GO TO 5
    4 GXYO=GXY
    5 GOLD=GXY*EXF(-X-Y,0.0D00)
      GO TO 8
    6 DA=DSQRT(A)
      DB=DSQRT(B)
      P=3.75/E
      B0=(.3989423+P*(.01328592+P*(.00225319-P*(.00157565-P*(.00916281-P
     1*(.02057706-P*(.02635537-P*(.01647633-.00392377*P))))))))/DSQRT(E)
      BF=B0*EXF(-Z,0.0D00)
      P=1./(1.+.3275911*(DA-DB))
      ERF=P*(.2548296-P*(.2844967-P*(1.421414-P*(1.453152-P*1.061405))))
      P=0.25/E
      C0=1.-1.772454*(DA-DB)*ERF
      C1=0.5-Z*C0
      C2=0.75-Z*C1
      C3=1.875-Z*C2
      C4=6.5625-Z*C3
      SUM=.1994711*(A-B)*P*(C0+1.5*P*(C1+1.666667*P*(C2+1.75*P*(C3+P*(C4
     1*(1.8-3.3*P*Z)+97.45313*P)))))
      GOLD=0.5*BF+(.3535534*(DA+DB)*ERF+SUM)*BF/(B0*DSQRT(E))
    8 IF(X.LT.Y) GOLD=1.+BF-GOLD
      RETURN
      END
C

