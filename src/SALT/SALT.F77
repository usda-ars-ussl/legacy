C     ******************************************************************
C     *                                                                *
C     *       NON-LINEAR LEAST-SQUARES ANALYSIS          SALT          *
C     *       OF SALINITY-RESPONSE CURVES                              *
C     *                                                                *
C     *       OCTOBER 15, 1985                                         *
C     *                                                                *
C     ******************************************************************
C
      DIMENSION C(50),Y(50),F(50),R(50),DELZ(50,8),B(8),E(8),TH(8),P(8)
     1,PHI(8),Q(8),TB(8),A(8,8),D(8,8),IND(50),NPA(13)
      CHARACTER BI(16)*6, TITLE*60
      DATA MIT/20/,STOPCR/0.00001/,NPA/2,1,2,2,3,4,5,6,7,8,2,3,2/
C
C     ----- OPEN I/O FILES -----
      OPEN(6,FILE='SALT.OUT',STATUS='OLD')
      CLOSE(6,STATUS='DELETE')
      OPEN(5,FILE='SALT.IN',STATUS='OLD')
      OPEN(6,FILE='SALT.OUT',STATUS='NEW')
C
C     ----- READ NUMBER OF CASES -----
      READ(5,1002) NC,KP
      IF(KP.NE.6) KP=0
      DO 120 NCASE=1,NC
C
C     ----- READ TITLE AND INPUT PARAMETERS -----
      READ(5,1000) TITLE
      READ(5,1002) KNOB,NOPT,KBI,KB
      IF(NOPT.GE.18) KBI=1
      IF(NOPT.GE.11) KB=1
      I=NOPT
      IF(I.GT.10) I=I-7
      NP=NPA(I)
      IF(KNOB.GT.0) WRITE(KP,1003) TITLE,NOPT,NP
      NOB=IABS(KNOB)
C
C     ----- READ INITIAL ESTIMATES -----
      IF(NCASE.EQ.1.OR.KBI.EQ.1) READ(5,1006) (BI(I),I=9,16)
      IF(KB.EQ.1) READ(5,1008) (B(I),I=1,8)
      CMAX=0.0
      CMIN=10000.0
      YMAX=0.0
      YMIN=10000.0
C
C     ----- READ INPUT DATA -----
      DO 2 I=1,NOB
      READ(5,1010) C(I),Y(I),IND(I)
      IF(NP.LE.3) IND(I)=1
      CMIN=AMIN1(CMIN,C(I))
      YMIN=AMIN1(YMIN,Y(I))
      IF(Y(I).LT.0.001) GO TO 2
      CMAX=AMAX1(CMAX,C(I))
    2 YMAX=AMAX1(YMAX,Y(I))
      IF(KNOB.LT.0) GO TO 120
      IF(NOPT-2) 4,10,22
C
C     ----- NOPT = 1: LINEAR REGRESSION -----
    4 SX=0.0
      SY=0.0
      SXY=0.0
      SXX=0.0
      SYY=0.0
      ROB=FLOAT(NOB)
      DO 6 I=1,NOB
      SX=SX+C(I)
      SY=SY+Y(I)
      SXY=SXY+C(I)*Y(I)
      SXX=SXX+C(I)*C(I)
    6 SYY=SYY+Y(I)*Y(I)
      SD=ROB*SXX-SX*SX
      YZERO=(SXX*SY-SX*SXY)/SD
      SLOPE=(SX*SY-ROB*SXY)/SD
      SR=ABS(SLOPE)*SQRT(ABS(SD/(ROB*SYY-SY*SY)))
      SIGMA=(SYY-YZERO*SY+SLOPE*SXY)/(ROB-2.)
      STDA=SQRT(SIGMA*SXX/SD)
      STDB=SQRT(ROB*SIGMA/SD)
      WRITE(KP,1012) YZERO,STDA,SLOPE,STDB,SR
      YM=Y(1)
      IF(KB.EQ.1) YM=B(3)
      CT=(YZERO-YM)/SLOPE
      WRITE(KP,1014) YM,CT
      SREL=SLOPE/YM
      CZERO=YZERO/SLOPE
      WRITE(KP,1034) SLOPE,SREL,YZERO,CZERO
      TH(3)=YM
      DO 8 I=1,NOB
      F(I)=YZERO-SLOPE*C(I)
    8 R(I)=Y(I)-F(I)
      GO TO 109
C
C     ----- NOPT = 2: SLOPE IS FITTED ONLY -----
   10 YM=Y(1)
      CT=C(1)
      IF(NOPT.EQ.2) GO TO 12
      NOPT=2
      WRITE(KP,1003) TITLE,NOPT,NP
   12 IF(KB.EQ.0) GO TO 14
      YM=B(3)
      CT=B(2)
   14 SA=0.0
      SB=0.0
      DO 16 I=1,NOB
      SA=SA+Y(I)-YM
   16 SB=SB+C(I)-CT
      SLOPE=-SA/SB
      SREL=SLOPE/YM
      YZERO=YM+SLOPE*CT
      CZERO=CT+YM/SLOPE
      WRITE(KP,1014) YM,CT
      WRITE(KP,1034) SLOPE,SREL,YZERO,CZERO
      TH(3)=YM
      DO 18 I=1,NOB
      F(I)=YM-SLOPE*(C(I)-CT)
   18 R(I)=Y(I)-F(I)
      IR=0
      DO 20 I=1,NOB
      IF(F(I).GE.0.) GO TO 20
      IR=1
      C(I)=CT
      Y(I)=YM
   20 CONTINUE
      IF(IR.GT.0) GO TO 14
      GO TO 109
C
C     ----- NONLINEAR LEAST-SQUARES ANALYSIS -----
   22 IF(KB.EQ.1) GO TO 26
      B(1)=(1.-YMIN/YMAX)/(CMAX-CMIN)
      B(2)=AMAX1(1.1*CMIN,0.25*CMAX)
      DO 24 I=3,8
   24 B(I)=0.98*YMAX
      IF(NOPT.EQ.4.OR.NOPT.EQ.11) B(3)=Y(1)
   26 IP=MAX0(NP,3)
      IF(NOPT.EQ.22) IP=2
      NP2=2*NP
      IP2=2*IP
      DO 28 I=1,8
   28 BI(I)=BI(I+8)
      NT=0
      B1=B(1)
      B2=B(2)
      B3=B(3)
      IF(NOPT.NE.3) GO TO 30
      B2=B(3)
      B(3)=B(2)
      B(2)=B2
      BI(2)=BI(11)
      BI(3)=BI(10)
   30 CONTINUE
      DO 32 I=1,8
      TB(I)=B(I)
   32 TH(I)=B(I)
C
C     ----------
      GA=0.02
      NIT=0
      CALL MODEL(TH,F,NOB,C,IND,NOPT)
      SSQ=0.
      DO 34 I=1,NOB
      R(I)=Y(I)-F(I)
   34 SSQ=SSQ+R(I)*R(I)
      WRITE(KP,1016) (BI(I),I=1,IP)
      WRITE(KP,1018) NIT,SSQ,(B(I),I=1,IP)
C
C     ----- BEGIN OF ITERATION -----
   36 NIT=NIT+1
      GA=0.1*GA
      DO 40 J=1,NP
      TEMP=TH(J)
      TH(J)=1.01*TH(J)
      Q(J)=0.0
      CALL MODEL(TH,DELZ(1,J),NOB,C,IND,NOPT)
      DO 38 I=1,NOB
      DELZ(I,J)=DELZ(I,J)-F(I)
   38 Q(J)=Q(J)+DELZ(I,J)*R(I)
      Q(J)=100.*Q(J)/TH(J)
   40 TH(J)=TEMP
      DO 46 I=1,NP
      DO 44 J=1,I
      SUM=0.0
      DO 42 K=1,NOB
   42 SUM=SUM+DELZ(K,I)*DELZ(K,J)
      D(I,J)=10000.*SUM/(TH(I)*TH(J))
   44 D(J,I)=D(I,J)
   46 E(I)=SQRT(D(I,I))
   50 DO 52 I=1,NP
      DO 52 J=1,NP
   52 A(I,J)=D(I,J)/(E(I)*E(J))
      DO 54 I=1,NP
      P(I)=Q(I)/E(I)
      PHI(I)=P(I)
   54 A(I,I)=A(I,I)+GA
      CALL MATINV(A,NP,P)
      STEP=1.0
   56 DO 58 I=1,NP
   58 TB(I)=P(I)*STEP/E(I)+TH(I)
      DO 62 I=1,NP
      IF(TH(I)*TB(I)) 66,66,62
   62 CONTINUE
      SUMB=0
      CALL MODEL(TB,F,NOB,C,IND,NOPT)
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
      ARG=SUM1/SQRT(SUM2*SUM3)
      ANGLE=57.29578*ATAN2(SQRT(1.-ARG*ARG),ARG)
      DO 72 I=1,NP
      IF(TH(I)*TB(I)) 74,74,72
   72 CONTINUE
      IF(SUMB/SSQ-1.0) 80,80,74
   74 IF(ANGLE-30.0) 76,76,78
   76 STEP=0.5*STEP
      GO TO 56
   78 GA=10.*GA
      GO TO 50
C
C     ----- PRINT COEFFICIENTS AFTER EACH ITERATION -----
   80 CONTINUE
      DO 82 I=1,NP
   82 TH(I)=TB(I)
      WRITE(KP,1018) NIT,SUMB,(TH(I),I=1,IP)
      IF(NOPT.GT.10.AND.NOPT.LT.18) GO TO 88
      IF(NOPT.EQ.3.OR.NOPT.EQ.20) GO TO 88
      IF(NOPT.LT.11) GO TO 85
      IF(ABS(TH(1)).GT.1.0E-06) GO TO 88
      IF(NT.EQ.1) GO TO 83
      NT=1
      B(1)=-B1
      B(2)=B2
      B(3)=B3
      GO TO 30
   83 B(1)=B(2)
      B(2)=B(3)
      DO 84 I=1,3
   84 BI(I)=BI(I+1)
      NOPT=20
      NP=2
      WRITE(KP,1020)
      WRITE(KP,1003) TITLE,NOPT,NP
      GO TO 30
   85 IF(TH(2).GT.0.001) GO TO 88
      IF(NOPT.EQ.5) GO TO 86
      WRITE(KP,1022)
      GO TO 119
   86 NP=2
      WRITE(KP,1024)
      GO TO 10
   88 DO 90 I=1,NP
      IF(ABS(P(I)*STEP/E(I))/(1.0E-20+ABS(TH(I)))-STOPCR) 90,90,92
   90 CONTINUE
      GO TO 94
   92 SSQ=SUMB
      IF(NIT.LE.MIT) GO TO 36
   94 CONTINUE
      IF(NOPT.LT.5.OR.NOPT.GT.10) GO TO 96
      IF(TH(2).GT.CMIN) GO TO 96
      WRITE(KP,1025)
      GO TO 110
   96 CALL MATINV(D,NP,P)
C
C     ----- WRITE CORRELATION MATRIX -----
      DO 98 I=1,NP
   98 E(I)=SQRT(AMAX1(D(I,I),1.E-20))
      WRITE(KP,1026) (I,I=1,NP)
      DO 102 I=1,NP
      DO 100 J=1,I
  100 A(J,I)=D(J,I)/(E(I)*E(J))
  102 WRITE(KP,1028) I,(A(J,I),J=1,I)
C
C     ----- CALCULATE 95% CONFIDENCE INTERVAL -----
      Z=1./FLOAT(NOB-NP)
      SDEV=SQRT(Z*SUMB)
      TVAR=1.96+Z*(2.3779+Z*(2.7135+Z*(3.187936+2.466666*Z**2)))
      WRITE(KP,1030)
      DO 108 I=1,NP
      SEC=E(I)*SDEV
      TVALUE=TH(I)/SEC
      TSEC=TVAR*SEC
      TMC=TH(I)-TSEC
      TPC=TH(I)+TSEC
  108 WRITE(KP,1032) BI(I),TH(I),SEC,TVALUE,TMC,TPC
C
C     ----- CALCULATE VARIOUS PARAMETERS -----
      IF(NOPT.GT.10) GO TO 109
  110 SLOPE=TH(1)*TH(3)
      IF(NOPT.EQ.3) SLOPE=TH(1)*TH(2)
      YZERO=TH(3)*(1.+TH(1)*TH(2))
      IF(NOPT.EQ.3) YZERO=B(2)*(1.+TH(1)*TH(3))
      CZERO=TH(2)+1./TH(1)
      IF(NOPT.EQ.3) CZERO=TH(3)+1./TH(1)
      WRITE(KP,1034) SLOPE,TH(1),YZERO,CZERO
C
C     ----- PREPARE FINAL OUTPUT -----
  109 WRITE(KP,1036)
      DO 118 I=1,NOB
      K=2+IND(I)
      IF(NOPT.EQ.3.OR.NOPT.EQ.20) K=2
      RY=Y(I)/TH(K)
      WRITE(KP,1038) C(I),Y(I),F(I),R(I),RY,IND(I)
  118 CONTINUE
  119 IF(NOPT.NE.5) GO TO 120
      IF(TH(2).GT.CMIN) GO TO 120
      WRITE(KP,1024)
      NP=1
      GO TO 10
  120 WRITE(KP,'(A)') CHAR(12)
C
C     ----- END OF PROBLEM -----
 1000 FORMAT(A60)
 1002 FORMAT(6I5)
 1003 FORMAT(5X,67(1H*)/5X,1H*,65X,1H*/5X,1H*,5X,'LEAST SQUARES ANALYSIS
     1 OF SALINITY RESPONSE CURVE',2X,'SALT',5X,1H*/5X,1H*,65X,1H*/5X,1H
     2*,A60,5X,1H*/5X,1H*,5X,'NOPT =',I3,': NP =',I3,42X,1H*/5X,1H*,65X,
     31H*/5X,67(1H*))
 1006 FORMAT(8(A6,4X))
 1008 FORMAT(8F10.0)
 1010 FORMAT(2F10.0,I10)
 1012 FORMAT(//5X,'LINEAR REGRESSION RESULTS FOR Y=YZERO-SLOPE*C'/5X,45(
     11H=)//5X,'YZERO =',F10.4,' WITH STANDARD ERROR OF',F10.4/5X,'SLOPE
     2 =',F10.4,' WITH STANDARD ERROR OF',F10.4/5X,'CORRELATION COEFFICI
     3ENT =',F10.4)
 1014 FORMAT(//5X,'CONTROL YIELD (YM) =',F10.4/5X,'THRESHOLD (CT) =',F10
     1.4)
 1016 FORMAT(///5X,'NIT',4X,'SSQ',5X,9(1X,A6))
 1018 FORMAT(4X,I3,F10.4,F9.4,F8.4,9F7.2)
 1020 FORMAT(//5X,'ALPHA IS TOO SMALL, NEW NOPT =20'/5X,32(1H=))
 1022 FORMAT(//5X,'THRESHOLD IS TOO SMALL, USE ANOTHER NOPT'/5X,40(1H=))
 1024 FORMAT(//5X,'CHANGED TO OPTION NUMBER 2')
 1025 FORMAT(//5X,'WARNING: THRESHOLD IS LESS THAN CMIN, USE ANOTHER NOP
     1T'/5X,55(1H=))
 1026 FORMAT(//7X,'CORRELATION MATRIX'/7X,18(1H=)/8X,10(3X,I2,3X))
 1028 FORMAT(4X,I4,10(1X,F6.3,1X))
 1030 FORMAT(//50X,'95% CONFIDENCE LIMITS'/6X,'VARIABLE',4X,'VALUE',5X,
     1'S.E.COEFF',5X,'T-VALUE',4X,'LOWER',8X,'UPPER')
 1032 FORMAT(6X,A6,F12.6,2X,F10.5,2X,F10.3,2X,F9.4,3X,F9.4)
 1034 FORMAT(/5X,'SLOPE FOR ORIGINAL DATA (S*YM) =',F10.6/5X,'SLOPE FOR
     1RELATIVE YIELD DATA (S) =',F10.6/5X,'INTERSECTION AT ZERO SALINITY
     2 (YZERO) =',F10.5/5X,'SALINITY EXTRAPOLATED TO ZERO YIELD (CZERO)=
     3',F10.5)
 1036 FORMAT(/7X,'CONC',7X,'Y-OBS',6X,'Y-FITTED',3X,'DEVIATION',3X,'REL
     1YIELD',2X,'INDEX')
 1038 FORMAT(5F12.4,5X,I2)
      CLOSE(5)
      CLOSE(6)
      STOP
      END
C
C     __________________________________________________________________
      SUBROUTINE MODEL(B,Y,NOB,C,IND,NOPT)
C
C     PURPOSE: TO CALCULATE Y(C)
C
      DIMENSION B(8),Y(50),C(50),IND(50)
C
C     ---------------
      IF(NOPT.GT.10) GO TO 10
    2 B2=B(2)
      IF(NOPT.EQ.3) B2=B(3)
      IF(NOPT.EQ.3) B3=B(2)
      DO 8 I=1,NOB
      IF(NOPT.NE.3) B3=B(2+IND(I))
      IF(C(I)-B2) 4,4,6
    4 Y(I)=B3
      GO TO 8
    6 Y(I)=B3-B(1)*B3*(C(I)-B2)
      Y(I)=AMAX1(Y(I),0.)
    8 CONTINUE
      RETURN
   10 IF(NOPT.GT.17) GO TO 14
      DO 12 I=1,NOB
   12 Y(I)=B(2+IND(I))/(1.+(C(I)/B(1))**B(2))
      RETURN
   14 IF(NOPT.EQ.20) GO TO 18
      DO 16 I=1,NOB
   16 Y(I)=B(3)*EXP(B(1)*C(I)-B(2)*C(I)**2)
      RETURN
   18 DO 20 I=1,NOB
   20 Y(I)=B(2)*EXP(-B(1)*C(I)**2)
      RETURN
      END
C
C     __________________________________________________________________
      SUBROUTINE MATINV(A,NP,B)
      DIMENSION A(8,8),B(8),INDEX(8,2)
      DO 2 J=1,NP
    2 INDEX(J,1)=0
      I=0
    4 AMAX=-1.0
      DO 12 J=1,NP
      IF(INDEX(J,1)) 12,6,12
    6 DO 10 K=1,NP
      IF(INDEX(K,1)) 10,8,10
    8 P=ABS(A(J,K))
      IF(P.LE.AMAX) GO TO 10
      IR=J
      IC=K
      AMAX=P
   10 CONTINUE
   12 CONTINUE
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
