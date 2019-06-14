C
C     ******************************************************************
C     *                                                                *
C     *       NON-LINEAR LEAST-SQUARES ANALYSIS             CFITIM     *
C     *                                                                *
C     ******************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION Y(90),X(90),F(90),R(90),DELZ(90,5),B(10),E(5),TH(10),
     1P(5),PHI(5),Q(5),LSORT(90),TB(10),A(5,5),D(5,5),INDEX(5)
      CHARACTER TITLE*60,BI(10)*3
      DATA STOPCR/0.0005/
C
C     ----- OPEN I/O FILES -----
      OPEN(6,FILE='CFITIM.OUT')
      CLOSE(6,STATUS='DELETE')
      OPEN(5,FILE='CFITIM.IN')
      OPEN(6,FILE='CFITIM.OUT',STATUS='NEW')
C
C     ----- READ NUMBER OF CASES -----
      READ(5,1006) NC,KP
      IF(KP.NE.6) KP=0
      DO 120 NCASE=1,NC
      WRITE(KP,1000)
C
C     ----- READ INPUT PARAMETERS -----
      READ(5,1006) MODE,NDATA,MIT,NOB
      M=(MODE-1)/2
      IF(M.EQ.0) WRITE(KP,1021)
      IF(M.EQ.1) WRITE(KP,1022)
      IF(M.EQ.2) WRITE(KP,1023)
      N=MODE-2*M
      IF(N.EQ.1) WRITE(KP,1024)
      IF(N.EQ.2) WRITE(KP,1025)
      NU=3
      IF(MODE.GT.2) NU=5
      NU1=NU+1
      NU2=2*NU
      READ(5,1001) TITLE
      WRITE(KP,1002) TITLE
C
C     ----- READ COEFFICIENTS NAMES -----
      READ(5,1004) (BI(I),I=1,NU2)
C
C     ----- READ INITIAL ESTIMATES -----
      READ(5,1005) (B(I),I=NU1,NU2)
      IF(M.EQ.2) B(8)=1./B(7)
C
C     ----- READ INDICES -----
      READ(5,1006) (INDEX(I),I=1,NU)
      IF(M.EQ.2) INDEX(3)=0
      WRITE(KP,1007)
      DO 4  I=1,NU
      J=2*I-1
    4 WRITE(KP,1008) I,BI(J),BI(J+1),B(I+NU)
      IF(NDATA.EQ.0) GO TO 10
C
C     ----- READ AND WRITE EXPERIMENTAL DATA -----
      DO 6 I=1,NOB
    6 READ(5,1005) X(I),Y(I)
   10 WRITE(KP,1009)
      DO 12 I=1,NOB
   12 WRITE(KP,1010) I,X(I),Y(I)
C
C     ----------
      NP=0
      DO 14 I=NU1,NU2
      TB(I)=B(I)
      IF(INDEX(I-NU).EQ.0) GO TO 14
      NP=NP+1
      K=2*NP-1
      J=2*(I-NU)-1
      BI(K)=BI(J)
      BI(K+1)=BI(J+1)
      B(NP)=B(I)
      TB(NP)=B(I)
      TH(NP)=B(NP)
   14 TH(I)=B(I)
C
C     ----------
      GA=0.02
      NIT=0
      NP2=2*NP
      CALL MODEL(TH,F,NOB,X,INDEX,MODE)
      SSQ=0.
      DO 32 I=1,NOB
      R(I)=Y(I)-F(I)
   32 SSQ=SSQ+R(I)*R(I)
      WRITE(KP,1011) (BI(J),BI(J+1),J=1,NP2,2)
      WRITE(KP,1012) NIT,SSQ,(B(I),I=1,NP)
C
C     ----- BEGIN OF ITERATION -----
   34 NIT=NIT+1
      GA=0.1*GA
      DO 38 J=1,NP
      TEMP=TH(J)
      TH(J)=1.01*TH(J)

      Q(J)=0
      CALL MODEL(TH,DELZ(1,J),NOB,X,INDEX,MODE)
      DO 36 I=1,NOB
      DELZ(I,J)=DELZ(I,J)-F(I)
   36 Q(J)=Q(J)+DELZ(I,J)*R(I)
      Q(J)=100.*Q(J)/TH(J)
C
C     ----- Q=XT*R (STEEPEST DESCENT) -----
   38 TH(J)=TEMP
      DO 44 I=1,NP
      DO 42 J=1,I
      SUM=0
      DO 40 K=1,NOB
   40 SUM=SUM+DELZ(K,I)*DELZ(K,J)
      D(I,J)=10000.*SUM/(TH(I)*TH(J))
   42 D(J,I)=D(I,J)
   44 E(I)=DSQRT(D(I,I))
   50 DO 52 I=1,NP
      DO 52 J=1,NP
   52 A(I,J)=D(I,J)/(E(I)*E(J))
C
C     ----- A IS THE SCALED MOMENT MATRIX -----
      DO 54 I=1,NP
      P(I)=Q(I)/E(I)
      PHI(I)=P(I)
   54 A(I,I)=A(I,I)+GA
      CALL MATINV(A,NP,P)
C
C     ----- P/E IS THE CORRECTION VECTOR -----
      STEP=1.0
   56 DO 58 I=1,NP
   58 TB(I)=P(I)*STEP/E(I)+TH(I)
      DO 62 I=1,NP
      IF(TH(I)*TB(I))66,66,62
   62 CONTINUE
      SUMB=0
      CALL MODEL(TB,F,NOB,X,INDEX,MODE)
      DO 64 I=1,NOB
      R(I)=Y(I)-F(I)
   64 SUMB=SUMB+R(I)*R(I)
   66 SUM1=0.0
      SUM2=0.0
      SUM3=0.0
      DO 68 I=1,NP
      SUM1=SUM1+P(I)*PHI(I)
      SUM2=SUM2+P(I)*P(I)
   68 SUM3=SUM3+PHI(I)*PHI(I)
      ARG=SUM1/DSQRT(SUM2*SUM3)
      ARG1=0.
      IF(NP.GT.1) ARG1=DSQRT(1.-ARG*ARG)
      ANGLE=57.29578*DATAN2(ARG1,ARG)
C
C     ----------
      DO 72 I=1,NP
      IF(TH(I)*TB(I))74,74,72
   72 CONTINUE
      IF((SUMB-SSQ).LT.1.D-08)GO TO 80
   74 IF(ANGLE-30.0)76,76,78
   76 STEP=0.5*STEP
      GO TO 56
   78 GA=10.*GA
      GO TO 50
C
C     ----- PRINT COEFFICIENTS AFTER EACH ITERATION -----
   80 CONTINUE
      DO 82 I=1,NP
   82 TH(I)=TB(I)
      WRITE(KP,1012) NIT,SUMB,(TH(I),I=1,NP)
      DO 86 I=1,NP
      IF(DABS(P(I)*STEP/E(I))/(1.0D-20+DABS(TH(I)))-STOPCR) 86,86,94
   86 CONTINUE
      GO TO 96
   94 SSQ=SUMB
      IF(NIT.LE.MIT) GO TO 34
C
C     ----- END OF ITERATION LOOP -----
   96 CONTINUE
      CALL MATINV(D,NP,P)
C
C     ----- WRITE CORRELATION MATRIX -----
      DO 98 I=1,NP
   98 E(I)=DSQRT(DMAX1(D(I,I),1.D-20))
      WRITE(KP,1013) (I,I=1,NP)
      DO 102 I=1,NP
      DO 100 J=1,I
  100 A(J,I)=D(J,I)/(E(I)*E(J))
  102 WRITE(KP,1014) I,(A(J,I),J=1,I)
C
C     ----- CALCULATE 95% CONFIDENCE INTERVAL -----
      Z=1./FLOAT(NOB-NP)
      SDEV=DSQRT(Z*SUMB)
      TVAR=1.96+Z*(2.3779+Z*(2.7135+Z*(3.187936+2.466666*Z**2)))
      WRITE(KP,1015)
      DO 108 I=1,NP
      SECOEF=E(I)*SDEV
      TVALUE=TH(I)/SECOEF
      TSEC=TVAR*SECOEF
      TMCOE=TH(I)-TSEC
      TPCOE=TH(I)+TSEC
      J=2*I-1
  108 WRITE(KP,1016) I,BI(J),BI(J+1),TH(I),SECOEF,TVALUE,TMCOE,TPCOE
C
C     ----- PREPARE FINAL OUTPUT -----
      LSORT(1)=1
      DO 116 J=2,NOB
      TEMP=R(J)
      K=J-1
      DO 111 L=1,K
      LL=LSORT(L)
      IF(TEMP-R(LL)) 112,112,111
  111 CONTINUE
      LSORT(J)=J
      GO TO 116
  112 KK=J
  113 KK=KK-1
      LSORT(KK+1)=LSORT(KK)
      IF(KK-L) 115,115,113
  115 LSORT(L)=J
  116 CONTINUE
      WRITE(KP,1017)
      DO 118 I=1,NOB
      J=LSORT(NOB+1-I)
  118 WRITE(KP,1018) I,X(I),Y(I),F(I),R(I),J,X(J),Y(J),F(J),R(J)
      WRITE(KP,1020)
  120 CONTINUE
C
C     ----- END OF PROBLEM -----
 1000 FORMAT(//5X,67(1H*)/5X,1H*,65X,1H*/5X,1H*,5X,'NONLINEAR LEAST SQUA
     1RES ANALYSIS',28X,1H*/5X,1H*65X,1H*)
 1001 FORMAT(A60)
 1002 FORMAT(5X,1H*,A60,5X,1H*/5X,1H*,65X,1H*/5X,67(1H*))
 1004 FORMAT(5(2A3,4X))
 1005 FORMAT(5F10.0)
 1006 FORMAT(5I5)
 1007 FORMAT(//5X,'INITIAL VALUES OF COEFFICIENTS'/5X,30(1H=)/5X,'NO',
     16X,'NAME',7X,'INITIAL VALUE')
 1008 FORMAT(4X,I3,5X,2A3,4X,F12.3)
 1009 FORMAT(//5X,'OBSERVED DATA'/5X,13(1H=)/5X,'OBS. NO.',5X,'PORE VOLU
     1ME',5X,'CONCENTRATION')
 1010 FORMAT(5X,I5,5X,F12.4,4X,F12.4)
 1011 FORMAT(//5X,'ITERATION',5X,'SSQ',2X,5(4X,2A3))
 1012 FORMAT(5X,I5,6X,F9.7,2X,3F9.5,2F10.5)
 1013 FORMAT(//5X,'CORRELATION MATRIX'/5X,18(1H=)/8X,10(2X,I2,7X))
 1014 FORMAT(3X,I3,10(2X,F7.4,2X))
 1015 FORMAT(///5X,'NON-LINEAR LEAST SQUARES ANALYSIS, FINAL RESULTS'
     1/5X,48(1H=)//54X,'95% CONFIDENCE LIMITS'/5X,'VAR',2X,'NAME'
     2,6X,'VALUE',6X,'S.E.COEFF.',3X,'T-VALUE',5X,'LOWER',8X,'UPPER')
 1016 FORMAT(5X,I2,2X,2A3,F12.5,3X,F9.4,4X,F8.2,2X,F9.4,4X,F9.4)
 1017 FORMAT(//5X,4(1H-),'ORDERED BY COMPUTER INPUT',5(1H-), 3X,8(1H-
     1),'ORDERED BY RESIDUALS',6(1H-)/11X,'PORE',3X,'CONCENTRATION',
     23X,'RESI-',9X,'PORE',3X,'CONCENTRATION',3X,'RESI-'/5X,'NO',2X,
     3'VOLUME',4X,'OBS.',2X,'FITTED',4X,'DUAL',3X,'NO',2X,'VOLUME',4X,
     4'OBS.',2X,'FITTED',4X,'DUAL')
 1018 FORMAT(5X,I2,4F8.3,3X,I2,4F8.3)
 1020 FORMAT(///5X,'END OF PROBLEM'/5X,14(1H=))
 1021 FORMAT(5X,1H*,5X,'EQUILIBRIUM TRANSPORT (MODEL A)',29X,1H*)
 1022 FORMAT(5X,1H*,5X,'NON-EQUILIBRIUM TRANSPORT (MODEL B)',25X,1H*)
 1023 FORMAT(5X,1H*,5X,'ONE-SITE KINETIC ADSORPTION (MODEL D)',23X,1H*)
 1024 FORMAT(5X,1H*,5X,'FIRST-TYPE BOUNDARY CONDITION',31X,1H*)
 1025 FORMAT(5X,1H*,5X,'THIRD-TYPE BOUNDARY CONDITION',31X,1H*)
C
C     ----- CLOSE FILES -----
      CLOSE(5)
      CLOSE(6)
      STOP
      END
C
C     ------------------------------------------------------------------
      SUBROUTINE MATINV(A,NP,B)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(5,5),B(10),INDEX(5,2)
      DO 2 J=1,5
    2 INDEX(J,1)=0
      I=0
    4 AMAX=-1.0
      DO 11 J=1,NP
      IF(INDEX(J,1)) 11,6,11
    6 DO 10 K=1,NP
      IF(INDEX(K,1)) 10,8,10
    8 P=DABS(A(J,K))
      IF(P.LE.AMAX) GO TO 10
      IR=J
      IC=K
      AMAX=P
   10 CONTINUE
   11 CONTINUE
      IF(AMAX) 30,30,14
   14 INDEX(IC,1)=IR
      IF(IR.EQ.IC) GO TO 18
      DO 16 L=1,NP
      P=A(IR,L)
      A(IR,L)=A(IC,L)
   16 A(IC,L)=P
      P=B(IR)
      B(IR)=B(IC)
      B(IC)=P
      I=I+1
      INDEX(I,2)=IC
   18 P=1./A(IC,IC)
      A(IC,IC)=1.0
      DO 20 L=1,NP
   20 A(IC,L)=A(IC,L)*P
      B(IC)=B(IC)*P
      DO 24 K=1,NP
      IF(K.EQ.IC) GO TO 24
      P=A(K,IC)
      A(K,IC)=0.0
      DO 22 L=1,NP
   22 A(K,L)=A(K,L)-A(IC,L)*P
      B(K)=B(K)-B(IC)*P
   24 CONTINUE
      GO TO 4
   26 IC=INDEX(I,2)
      IR=INDEX(IC,1)
      DO 28 K=1,NP
      P=A(K,IR)
      A(K,IR)=A(K,IC)
   28 A(K,IC)=P
      I=I-1
   30 IF(I) 26,32,26
   32 RETURN
      END
C
C     ------------------------------------------------------------------
      SUBROUTINE MODEL(B,Y,NOB,X,INDEX,MODE)
C
C     PURPOSE: TO CALCULATE CONCENTRATIONS FOR GIVEN PORE VOLUME
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION B(10),Y(90),X(90),INDEX(5),XG(20),W(20)
      DATA XG/.03877242,.1160841,.1926976,.2681522,.3419941,.4137792,
     1.4830758,.5494671,.6125539,.6719567,.7273183,.7783057,.8246122,
     2.8659595,.9020988,.9328128,.9579168,.9772599,.9907262,.9982377/
      DATA W/.07750595,.07703982,.07611037,.07472317,.07288658,.07061165
     1,.06791204,.06480401,.06130624,.05743977,.05322785,.04869581,.0438
     27091,.03878217,.03346019,.02793701,.02224585,.01642106,.01049828,.
     300452128/
      K=0
      IF(MODE.LE.2) GO TO 12
C
C     ----- SOLUTION FOR NON-EQUILIBRIUM TRANSPORT (MODEL B) -----
      DO 2 I=6,10
      IF(INDEX(I-5).EQ.0) GO TO 2
      K=K+1
      B(I)=B(K)
    2 CONTINUE
      P=B(6)
      R=B(7)
      IF(MODE.GE.5) B(8)=1./R
      BETA=DMIN1(B(8),.9999D00)
      OMEGA=B(9)
      DO 10 J=1,NOB
      DO 8 M=1,2
      C=0.0
      T=X(J)+(1-M)*B(10)
      IF(T.LE.0.) GO TO 6
      A=DSQRT(1.+.05*P)
      T2=DMIN1(T,BETA*R*(1.+40.*(1.+A)/P))
      T1=DMAX1(0.D00,BETA*R*(1.+40.*(1.-A)/P))
      IF(T2.LE.T1) GO TO 6
      DO 4 I=1,20
      TAU=0.5*(T1+T2+(T2-T1)*XG(I))
      C=C+W(I)*CCO(P,R,BETA,OMEGA,T,TAU,MODE)
      TAU=0.5*(T1+T2+(T1-T2)*XG(I))
    4 C=C+W(I)*CCO(P,R,BETA,OMEGA,T,TAU,MODE)
      C=0.5*(T2-T1)*C
    6 IF(M.EQ.2) GO TO 8
      Y(J)=C
    8 CONTINUE
   10 Y(J)=Y(J)-C
      RETURN
C
C     ----- SOLUTION FOR EQUILIBRIUM TRANSPORT (MODEL A) -----
   12 DO 14 I=4,6
      IF(INDEX(I-3).EQ.0) GO TO 14
      K=K+1
      B(I)=B(K)
   14 CONTINUE
      E=0.0
      P=B(4)
      R=B(5)
      DO 18 J=1,NOB
      DO 16 M=1,2
      C=0.0
      T=X(J)+(1-M)*B(6)
      IF(T.LE.0.) GO TO 18
      CM=0.5*(R-T)*DSQRT(P/(R*T))
      CP=0.5*(R+T)*DSQRT(P/(R*T))
      C=0.5*EXF(E,CM)+0.5*EXF(P,CP)
      IF(MODE.EQ.2) C=C+DSQRT(.3183099*P*T/R)*EXF(-CM*CM,E)-0.5*(2.+P+P*
     1T/R)*EXF(P,CP)
      IF(M.EQ.2) GO TO 18
      Y(J)=C
   16 CONTINUE
   18 Y(J)=Y(J)-C
      RETURN
      END
C
C     ------------------------------------------------------------------
      FUNCTION CCO(P,R,BETA,OMEGA,T,TAU,MODE)
C
C     PURPOSE: TO CALCULATE THE ARGUMENT UNDER THE INTEGRAL SIGN -----
C
      IMPLICIT REAL*8 (A-H,O-Z)
      CCO=0.0
      BER=BETA*R
      CM=P*(BER-TAU)**2/(4.*BER*TAU)
      C=.2820948*DSQRT(P*BER/TAU**3)*DEXP(-CM)
      IF((MODE.EQ.3).OR.(MODE.EQ.5)) GO TO 2
      CP=(BER+TAU)*DSQRT(P/(4.*BER*TAU))
      C=2.*C*TAU/BER-0.5*P*EXF(P,CP)/BER
    2 IF(C.LT.1.D-07) RETURN
      EPSI=OMEGA*TAU/BER
      ETHA=OMEGA*(T-TAU)/(R-BER)
      CCO=C*GOLD(EPSI,ETHA)
      RETURN
      END
C
C     ------------------------------------------------------------------
      FUNCTION GOLD(X,Y)
C
C     PURPOSE: TO CALCULATE J(X,Y)
C
      IMPLICIT REAL*8(A-H,O-Z)
      GOLD=0.0
      BF=0.0
      E=2.*DSQRT(X*Y)
      Z=X+Y-E
      IF(Z.GT.17.) GO TO 8
      IF(E.NE.0.) GO TO 2
      GOLD=DEXP(-X)
      RETURN
    2 A=DMAX1(X,Y)
      B=DMIN1(X,Y)
      NT=11.+2.*B+0.3*A
      IF(NT.GT.25) GO TO 6
      I=0
      IF(X.LT.Y) I=1
      GXY=1.+I*(B-1.)
      GX=1.0
      GY=GXY
      GZ=1.0
      DO 4 K=1,NT
      GX=GX*A/K
      GY=GY*B/(K+I)
      GZ=GZ+GX
    4 GXY=GXY+GY*GZ
      GOLD=GXY*DEXP(-X-Y)
      GO TO 8
    6 DA=DSQRT(A)
      DB=DSQRT(B)
      P=3.75/E
      B0=(.3989423+P*(.01328592+P*(.00225319-P*(.00157565-P*(.00916281-P
     1*(.02057706-P*(.02635537-P*(.01647633-.00392377*P))))))))/DSQRT(E)
      BF=B0*DEXP(-Z)
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
C     ------------------------------------------------------------------
      FUNCTION EXF(A,B)
C
C     PURPOSE: TO CALCULATE EXP(A) ERFC(B)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      EXF=0.0
      IF((DABS(A).GT.170.).AND.(B.LE.0.)) RETURN
      IF(B.NE.0.0) GO TO 1
      EXF=DEXP(A)
      RETURN
    1 C=A-B*B
      IF((DABS(C).GT.170.).AND.(B.GT.0.)) RETURN
      IF(C.LT.-170.) GO TO 4
      X=DABS(B)
      IF(X.GT.3.0) GO TO 2
      T=1./(1.+.3275911*X)
      Y=T*(.2548296-T*(.2844967-T*(1.421414-T*(1.453152-1.061405*T))))
      GO TO 3
    2 Y=.5641896/(X+.5/(X+1./(X+1.5/(X+2./(X+2.5/(X+1.))))))
    3 EXF=Y*DEXP(C)
    4 IF(B.LT.0.0) EXF=2.*DEXP(A)-EXF
      RETURN
      END
