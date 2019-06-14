C CXTFIT2.FOR 3/3/95
C
C     ******************************************************************
C     *                                                                *
C     *   CXTFIT VERSION 2.1                                           *
C     *                        LAST MODIFIED:       Feb 12, 2004       *
C     *                                                                *
C     *   NON-LINEAR LEAST-SQUARES ANALYSIS OF C(X,T) DATA             *
C     *   FOR  ONE-DIMENSIONAL  DETERMINISTIC OR STOCHASTIC            *
C     *   EQUILIBRIUM AND NONEQUILIBRIUM CONVECTIVE DISPERSIVE         *
C     *   EQUATION                                                     *
C     *                                                                *
C     *                                                                *
C     *                       NOBUO TORIDE                             *
C     *                       U.S. SALINITY LABORATORY                 *
C     *                       450 West Big Spring Road                 *
C     *                       RIVERSIDE CA 92507-4617                  *
C     *                                                                *
C     *           CURRENT ADDRESS                                      *
C     *               DEPT. OF AGRICULTURAL SCIENCES                   *
C     *               SAGA UNIV., SAGA 840, JAPAN                      *
C     *               TEL. +81-952-28-8756                             *
C     *               FAX. +81-952-28-8709                             *
C     *               E-MAIL NOBUO@CC.SAGA-U.AC.JP                     *
C     ******************************************************************
C   
      IMPLICIT REAL*8(A-H,O-Z)
C    MAXOB; No. OF OBSERVATION. 
      PARAMETER (MAXOB=405)
C
C
      DIMENSION T(MAXOB),Z(MAXOB),C(MAXOB),C1(MAXOB),C2(MAXOB),
     & VAR1(MAXOB),VAR2(MAXOB),F(MAXOB),R(MAXOB),DELZ(MAXOB,15),B(30),
     & E(15),TH(30),P(15),PHI(15),Q(15),LSORT(MAXOB),TB(30),A(15,15),
     & D(15,15),BMAX(15),BMIN(15),DERL(30)
      CHARACTER BI(15)*6
      CHARACTER FILEIN*15,FILEOUT*15
      COMMON MODC,MM,AA,BB,DA,CX,PEC,BETA,BETR,OMEGA,RE,DMU1,DMU2,
     &   MODE,MCON,TT,ZZ,STOPPER,LEVEL,ICHEB
      COMMON/MODAT/INDEX(15),INVERSE,NREDU,NVAR,ZL,MIT,MDEG,MASS,
     &     DUMTP(50),DUMGA1(50),MNEQ,ISKIP,PHIM,PHIIM
      COMMON/STOCH/MODD,MODK,MODS,MSTOCH,CORR,MCORR,MSD,SDLNK,
     & V,SDLNV,VMAX,VMIN,DIS,AVEY,SDLNY,YMAX,YMIN,ALPHA,SDLND
     & ,DK,RHOTH,MD56,MK34,MAL8
      COMMON/BOUN/MODB,NPULSE,TPULSE(50),PULSE(50),MASSST
      COMMON/INITI/MODI,NINI,CINI(50),ZINI(50)
      COMMON/PROD/MODP,NPRO1,NPRO2,GAMMA1(50),ZPRO1(50),GAMMA2(50),
     &    ZPRO2(50)
C
C     ----- OPEN I/O FILES ----------- 
       KP=7
       WRITE(*, 1000)
       READ(*, 1001) FILEIN
       IF(FILEIN.EQ.'  ')FILEIN='CXTFIT2.IN'
       OPEN(5, FILE = FILEIN, STATUS = 'OLD')
C
       WRITE(*, 1002)
       READ(*, 1001) FILEOUT
       IF(FILEOUT.EQ.'  ')FILEOUT='CXTFIT2.OUT'
       OPEN(KP, FILE = FILEOUT, STATUS = 'UNKNOWN')
C
C     ----- READ NUMBER OF CASES -----
      READ(5,*) NC
      DO 155 NCASE=1,NC
       ISKIP=0
C     ----  READ, CHECK, AND WRITE DATA -------
C      
      CALL DATAIN(NOB,NU1,NU2,T,Z,C,B,BMAX,BMIN,BI,ILMT,
     &   KP,MAXOB,MPRINT,DT,DZ,NT,NZ,FILEIN,INRHO)
C
      CALL CONST1(MAXTRY,STOPCR,GA,GD,DERL,STSQ,MM,ICHEB,NU1,NU2,OMMAX,
     &           MIT)
C
C    ------ DIRECT PROBLEM ----------------
      IF(INVERSE.LE.0) THEN
         CALL DIRECT(T,Z,B,C1,C2,VAR1,VAR2,MPRINT,DT,DZ,NT,NZ,MAXOB,KP)
        GOTO 154
      END IF  
C
C    ------ PREPARE INVERSE PROBLEM -------
C
C
C   WRITE EXPERIMENTAL DATA
C      WRITE(KP,1009)
C      DO 14 I=1,NOB
C        WRITE(KP,1010) I,C(I),Z(I),T(I)
C   14 CONTINUE
C
C
C  TWO-SITE MODEL; BETA SHOULD BE LARGER THAN 1/R
C   (MODE=2, MNEQ=2, R=B(NVAR+3),BETA=B(NVAR+4))
C
C  TWO-REGION MODEL; PHIM/R < BETA < (PHIM+R-1)/R
C   (MODE=2, MNEQ=3, R=B(NVAR+3),BETA=B(NVAR+4))
C
C  BMINTS; INITIAL VALUE OF A MINIMUM CONSTRAINT FOR BETA
C  BMAXTS; INITIAL VALUE OF A MINIMUM CONSTRAINT FOR BETA
         IF(MNEQ.LE.1.OR.ILMT.NE.1) GOTO 210
           BMINTS=BMIN(4)
           BMAXTS=BMAX(4)
C
C     BETA SHOULD BE LESS THAN 1 FOR THE NONEQUILIBRIUM CDE
C     OMEGA SHOULD BE LESS THAN OMMAX
C         (DEFAULT VALUE IS 100, SEE CONST1)
  210    NB=0
         NR=0
         NOMEG=0
         IF(MODE.EQ.2) THEN
           NB=INDEX(1)+INDEX(2)+INDEX(3)+INDEX(4)
           NR=INDEX(1)+INDEX(2)+INDEX(3)
           NOMEG=INDEX(1)+INDEX(2)+INDEX(3)+INDEX(4)+INDEX(5)
         END IF
         IF(MODE.LE.2) GOTO 212
C
         IF(MOD(MODE,2).EQ.1) GOTO 212
           NOMEG=INDEX(1)+INDEX(2)+INDEX(3)+INDEX(4)
C
C     ----- REARRANGE VARIABLE ARRAYS -----------------
  212 NP=0
      DO 20 I=NU1,NU2
      TB(I)=B(I)
      IF(INDEX(I-NVAR).EQ.0) GO TO 20
      NP=NP+1
      DERL(NP)=DERL(I)
      BI(NP)=BI(I-NVAR)
      B(NP)=B(I)
      TB(NP)=B(I)
      TH(NP)=B(NP)
      IF(ILMT.NE.0) BMIN(NP)=BMIN(I-NVAR)
      IF(ILMT.NE.0) BMAX(NP)=BMAX(I-NVAR)
   20 TH(I)=B(I)
C
C      ALL PARAMETERS SHOULD BE POSITIVE (EXCEPT RHO FOR MODE>=3).
C       [MODE=1]>1.D-10  [MODE>=2] > 1.D-7
C       [MODE>=3]  RHO < -1.D5, RHO > 1.D5 
          DO 22 I=1,NP
          IF(I.EQ.NP.AND.INRHO.EQ.1) THEN
             IF(TH(I).LT.-1.0) TH(I)=-1.D+00
             IF(TH(I).GT.1.0) TH(I)=1.D+00
             IF(DABS(TH(I)).LT.1.D-05) THEN
               IF((TH(I)-1.D-30).GT.0.D+00) THEN
                 TH(I)=1.D-05
               ELSE
                 TH(I)=-1.D-05
               END IF
             END IF
             GOTO 22
          END IF 
          IF(MODE.EQ.1.AND.DABS(TH(I)).LT.1.D-10) TH(I)=1.D-10
          IF(MODE.GE.2.AND.DABS(TH(I)).LT.1.D-07) TH(I)=1.D-07
   22     CONTINUE
C
C    ------ START INVERSE PROBLEM -------
      NIT=0
      NP2=2*NP
      ISKIP=0
      DO 200 I=1,NOB
      TT=T(I)
      ZZ=Z(I)
      CALL MODEL(TH,TT,ZZ,CXT1,CXT2,DUM1,DUM2)
      ISKIP=1
      F(I)=CXT1
C  TOTAL RESIDANT CONCENTRATIONS FOR THE NONEQUILIBRIUM TRANSPORT
         IF(MODE.EQ.2.AND.(MODC.EQ.4.OR.MODC.EQ.6)) THEN
          F(I)=F(I)+CXT2
         END IF
  200 CONTINUE
      IF (MIT.EQ.0) GO TO 140
      SSQ=0.
      DO 32 I=1,NOB
      R(I)=C(I)-F(I)
   32 SSQ=SSQ+R(I)*R(I)
      WRITE(KP,1011) (BI(I),I=1,NP)
      WRITE(KP,1012) NIT,SSQ,(B(I),I=1,NP)
      IF(MODE.GE.2) THEN
      WRITE(*,1011) (BI(I),I=1,NP)
      WRITE(*,1012) NIT,SSQ,(B(I),I=1,NP)
      END IF
C
C     ----- BEGIN ITERATION ----FOR NIT=1 GA=GA/GD (=0.001)
      MFINAL=0
   34 NIT=NIT+1
      ISIG=1
      NTRIAL=0
      GA=GA/GD
      DO 38 J=1,NP
         IF(J.EQ.NP.AND.INRHO.EQ.1) THEN
           IF(DABS(TH(J)*(1+DERL(J))).GT.1.0D+00) ISIG=-1
         END IF
         IF(MODE.EQ.2.AND.INDEX(4).EQ.1) THEN
           IF((TH(J)*(1+DERL(J))).GE.1.0D+00) ISIG=-1
         END IF
      TEMP=TH(J)
C
           TH(J)=(1.0+ISIG*DERL(J))*TH(J)
C  
      Q(J)=0.
C
      ISKIP=0
      DO 201 I=1,NOB
      TT=T(I)
      ZZ=Z(I)
      CALL MODEL(TH,TT,ZZ,CXT1,CXT2,DUM1,DUM2)
      ISKIP=1
      DELZ(I,J)=CXT1
         IF(MODE.EQ.2.AND.(MODC.EQ.4.OR.MODC.EQ.6)) THEN
           DELZ(I,J)=DELZ(I,J)+CXT2
         END IF
  201 CONTINUE
      DO 36 I=1,NOB
      DELZ(I,J)=ISIG*(DELZ(I,J)-F(I))
   36 Q(J)=Q(J)+DELZ(I,J)*R(I)
      Q(J)=Q(J)/TH(J)/DERL(J)
C
C     ----- Q=XT*R (STEEPEST DESCENT) -----
C            (DERL(I) INCREMENT FOR PARAMETER I)
   38 TH(J)=TEMP
      DO 44 I=1,NP
      DO 42 J=1,I
      SUM=0.
      DO 40 K=1,NOB
   40 SUM=SUM+DELZ(K,I)*DELZ(K,J)
      D(I,J)=SUM/(TH(I)*TH(J))/DERL(I)/DERL(J)
   42 D(J,I)=D(I,J)
      E(I)=DSQRT(D(I,I))
   44 E(I)=DMAX1(E(I),1.D-30)
   50 DO 52 I=1,NP
      DO 52 J=1,NP
   52 A(I,J)=D(I,J)/(E(I)*E(J))
C
C     ----- A IS THE SCALED MOMENT MATRIX -----
C           PHI IS THE SCALED VECOTR 
      DO 54 I=1,NP
      P(I)=Q(I)/E(I)
      PHI(I)=P(I)
   54   A(I,I)=A(I,I)+GA
C
        IF(MFINAL.EQ.0) THEN
           CALL MATINV(A,NP,P)
	  ELSE
           CALL MATINV(D,NP,P)
           GOTO 97
	  END IF
C
C     ----- P/E IS THE CORRECTION VECTOR -----
C           A = THE INVERSE OF A (COVARIANCE MATRIX)
      STEP=1.0
   56 CONTINUE
C
C  TWO-SITE MODEL; BETA SHOULD BE LARGER THAN 1/R
C   (MODE=2, MNEQ=2)
C
C  TWO-REGION MODEL; PHIM/R < BETA < (PHIM+R-1)/R
C   (MODE=2, MNEQ=3)
         IF(MNEQ.GE.2.AND.INDEX(4).EQ.1) THEN
           IF(INDEX(3).EQ.0)  THEN
             CR=B(NVAR+3)
           ELSE
             TB(NR)=P(NR)*STEP/E(NR)+TH(NR)
             CR=TB(NR)
           END IF
         IF(ILMT.EQ.0) GOTO 255
            IF(MNEQ.EQ.2) THEN
              BMIN(NB)=1/CR
              IF(BMIN(NB).LT.BMINTS) BMIN(NB)=BMINTS
            ELSE IF(MNEQ.EQ.3) THEN
              BMIN(NB)=PHIM/CR
              IF(BMIN(NB).LT.BMINTS) BMIN(NB)=BMINTS
              BMAX(NB)=(PHIM+CR-1.0)/CR
              IF(BMAX(NB).GT.BMAXTS) BMAX(NB)=BMAXTS
            END IF
         END IF
C
  255 DO 59 I=1,NP
      TB(I)=P(I)*STEP/E(I)+TH(I)
C     STOCHASTIC MODEL RHO CANNOT BE ZERO TO EVALUATE DERIVATIVES
         IF(I.EQ.NP.AND.INRHO.EQ.1) THEN
           IF(DABS(TB(I)).LT.1.D-05) THEN
             IF((TB(I)-1.D-30).GT.0.D+00) THEN
               TB(I)=1.D-05
             ELSE
               TB(I)=-1.D-05
             END IF
               GOTO 58
           END IF
         END IF
       IF(ILMT.EQ.0) THEN
C      STOCHASTIC MODEL   -1=<RHO<=1
           IF(I.EQ.NP.AND.INRHO.EQ.1) THEN
             IF(TB(I).LE.-1.0) TB(I)=-1.D+00
             IF(TB(I).GE.1.0) TB(I)=1.D+00 
             GOTO 58
           END IF
C  PARAMETERS SHOULD BE POSITIVE
       IF(MODE.EQ.1) THEN
         IF(TB(I).LT.1.D-10) THEN
             TB(I)=1.D-10
           GOTO 58
         ELSE
          GO TO 59
         END IF
       ELSE
         IF(TB(I).LT.1.D-07) THEN
             TB(I)=1.D-07
             
           GOTO 58
         ELSE
          GO TO 59
         END IF
       END IF
      END IF
      IF(DABS(BMAX(I)-BMIN(I)).LT.1.D-10) GO TO 59
      IF((TB(I)-BMAX(I)).GT.1.D-20) GO TO 57
      IF((TB(I)-BMIN(I)).GT.1.D-20) GO TO 59
      TB(I)=BMIN(I)
      GO TO 58
   57 TB(I)=BMAX(I)
   58 P(I)=(TB(I)-TH(I))*E(I)/STEP
   59 CONTINUE
C
C  BETA SHOULD BE LESS THAN ONE
         IF(MODE.EQ.2) THEN 
          IF(INDEX(4).EQ.1.AND.TB(NB).GT.0.9999) THEN
            TB(NB)=0.9999
             P(NB)=(TB(NB)-TH(NB))*E(NB)/STEP
          END IF
          IF(INDEX(4).EQ.1.AND.TB(NB).LT.0.0001) THEN
           TB(NB)=0.0001
           P(NB)=(TB(NB)-TH(NB))*E(NB)/STEP
          END IF
          IF(ILMT.NE.0) GO TO 60 
C
         IF(MNEQ.GE.2.AND.INDEX(4).EQ.1) THEN
           IF(MNEQ.EQ.2) THEN
             IF(TB(NB).LE.1/CR) THEN 
               TB(NB)=1/CR
               P(NB)=(TB(NB)-TH(NB))*E(NB)/STEP
             END IF
           ELSE IF(MNEQ.EQ.3) THEN
             IF(TB(NB).LE.PHIM/CR)  THEN 
              TB(NB)=PHIM/CR
               P(NB)=(TB(NB)-TH(NB))*E(NB)/STEP
             ELSE IF(TB(NB).GE.(PHIM+CR-1.0)/CR) THEN
              TB(NB)=(PHIM+CR-1.0)/CR
               P(NB)=(TB(NB)-TH(NB))*E(NB)/STEP
             END IF            
           END IF
         END IF
C  OMEGA IS LESS THAN OMMAX (=100)
         IF(INDEX(5).EQ.1.AND.TB(NOMEG).GT.OMMAX) THEN
             TB(NOMEG)=OMMAX
             P(NOMEG)=(TB(NOMEG)-TH(NOMEG))*E(NOMEG)/STEP
           END IF
         END IF
          IF(MODE.EQ.4.OR.MODE.EQ.6) THEN
            IF(INDEX(4).EQ.1.AND.TB(NOMEG).GT.1.E04) THEN
              TB(NOMEG)=OMMAX
              P(NOMEG)=(TB(NOMEG)-TH(NOMEG))*E(NOMEG)/STEP
            END IF
          END IF   
         IF(MODE.EQ.8.AND.INDEX(4).EQ.1) THEN 
              DOMEGA=TB(NOMEG)*RHOTH*TB(NVAR+4)*ZL/TB(NVAR+1) 
           IF(DOMEGA.GT.OMMAX) THEN
              TB(NOMEG)=OMMAX/RHOTH/TB(NVAR+4)/ZL*TB(NVAR+1) 
              P(NOMEG)=(TB(NOMEG)-TH(NOMEG))*E(NOMEG)/STEP
           END IF
         END IF
   60 DO 62 I=1,NP
       IF(I.EQ.NP.AND.INRHO.EQ.1) GOTO 62
      IF(TH(I)*TB(I))66,66,62
   62 CONTINUE
      SUMB=0.0
C
      ISKIP=0
      DO 205 I=1,NOB
      TT=T(I)
      ZZ=Z(I)
      CALL MODEL(TB,TT,ZZ,CXT1,CXT2,DUM1,DUM2)
      ISKIP=1
      F(I)=CXT1
         IF(MODE.EQ.2.AND.(MODC.EQ.4.OR.MODC.EQ.6)) THEN
           F(I)=F(I)+CXT2
         END IF
  205 CONTINUE
C ----- DEVIATION FROM THE STEEPEST DESCENT DIRECTION --------
C         ANGLE BETWEEN PHI AND P (CORRECTION VECTOR) 
C            AND PHI (STEEPEST DECENT DIRECTION)  
      DO 64 I=1,NOB
      R(I)=C(I)-F(I)
   64 SUMB=SUMB+R(I)*R(I)
   66 SUM1=0.0
      SUM2=0.0
      SUM3=0.0
      DO 68 I=1,NP
      SUM1=SUM1+P(I)*PHI(I)
      SUM2=SUM2+P(I)*P(I)
   68 SUM3=SUM3+PHI(I)*PHI(I)
      ARG=SUM1/DSQRT(DMAX1(SUM2*SUM3,1D-20))
      ARG1=0.0D+00
      IF(NP.GT.1) ARG1=DSQRT(1.-ARG*ARG)
C   (57.29578=180/pi)
      ANGLE=57.29578*DATAN2(ARG1,ARG)
C
C     ----------
      DO 72 I=1,NP
         IF(I.EQ.NP.AND.INRHO.EQ.1) GOTO 72
      IF(TH(I)*TB(I))74,74,72
   72 CONTINUE
      NTRIAL=NTRIAL+1
C      write(7,*) 'ntrial step ga ',ntrial,step,ga
C      write(7,*) 'ssq new/ssq angle',sumb/ssq,angle
C
C IF THE TRIAL SUCCEED, DECREASE GA. 
C FAILURE, IF ANGLE < 30, DECREASE THE SCALE OF THE CORRECTION VECTOR.
C           OTHERWISE, INCREASE GA.
       IF(NTRIAL.GT.MAXTRY) GO TO 95
      IF((SUMB/SSQ-1.0).LE.1.D-10.OR.STEP.LT.1.D-20) GOTO 80
   74 IF(ANGLE-30.0)76,76,78
C  IF 
   76 STEP=0.5*STEP
      GO TO 56
   78 GA=GA*GD
      GO TO 50
C   ----- PRINT COEFFICIENTS AFTER EACH ITERATION -----
   80 CONTINUE
      DO 82 I=1,NP
   82 TH(I)=TB(I)
      WRITE(KP,1012) NIT,SUMB,(TH(I),I=1,NP)
      IF(MODE.GE.2) THEN
      WRITE(*,1012) NIT,SUMB,(TH(I),I=1,NP)
      END IF
C  IF CHANGES OF ALL PARAMETE VALUES < STOPCR, STOP THE ITERATION.       
      DO 86 I=1,NP
       DUMX=DABS(P(I)*STEP/E(I))/DABS(1.0D-20+TH(I))
C      WRITE(*,*) 'I STOP',I,DUMX
      IF(DUMX-STOPCR) 86,86,94
   86 CONTINUE
      GO TO 96
   94 SSQ=SUMB
C STSQ=STOP CRITERIA TO EVALUATE THE IMPROVEMENT OF SSQ
         IF(NIT.EQ.1) SSQ3=SSQ
         IF(NIT.EQ.2) SSQ2=SSQ
         IF(NIT.EQ.3) SSQ1=SSQ
         IF(NIT.GE.4) THEN
           DSSQ1=1.0-SSQ/SSQ1
           DSSQ2=1.0-SSQ1/SSQ2
           DSSQ3=1.0-SSQ2/SSQ3 
           IF(DSSQ3.LT.0.0) GOTO 91 
           IF(DSSQ2.LT.0.0) GOTO 91 
           IF(DSSQ1.LT.0.0) GOTO 91 
           IF(DSSQ1.LT.STSQ.AND.DSSQ2.LT.STSQ.AND.DSSQ3.LT.STSQ) THEN
             WRITE(KP,1035) NIT-3,NIT
             GO TO 96
           END IF
   91      SSQ3=SSQ2
           SSQ2=SSQ1
           SSQ1=SSQ
         END IF
         IF(NIT.LT.MIT) GO TO 34
         IF(NIT.EQ.MIT) WRITE(KP,1034) MIT
         GO TO 96
   95 WRITE(KP,1038) MAXTRY
C
C     ----- END OF ITERATION LOOP -----
C     EVALUATE COVARIANCE MATRIX USING THE ESTIMATED PARAMETERS
   96 MFINAL=1
      GOTO 34
C
C     ----- WRITE COVARIANCE MATRIX -----
   97 DO 98 I=1,NP 
      IF(D(I,I).LT.0.0+00) THEN
        WRITE(KP,1020)
        GOTO 104
      END IF
      E(I)=DSQRT(D(I,I))
   98 E(I)=DMAX1(E(I),1.D-30)
      IF(NP.EQ.1) GO TO 104
      WRITE(KP,1013) (BI(I),I=1,NP)
      DO 102 I=1,NP
      DO 100 J=1,I
  100 A(J,I)=D(J,I)/(E(I)*E(J))
  102 WRITE(KP,1014) BI(I),(A(J,I),J=1,I)
C   ------- COEEFICIENT OF DETERMINATION (r2) ------
  104 SUMC=0.0
      DO 105 I=1,NOB
      SUMC=SUMC+C(I)
  105 CONTINUE
      SUMC=SUMC/NOB
      SUMC2=0.0
      SUMCF=0.0
      DO 106 I=1,NOB
      SUMC2=SUMC2+(C(I)-SUMC)**2
      SUMCF=SUMCF+(C(I)-F(I))**2
  106 CONTINUE
      RSQ=1.0-SUMCF/SUMC2
      WRITE(KP,1041) RSQ
C  TMSE: THE MEAN SQURE FOR ERROR (MSE)
	TMSE=SUMB/FLOAT(NOB-NP)
	WRITE(KP,1050) TMSE
C
C     ----- CALCULATE 95% CONFIDENCE INTERVAL -----
      XZ=1./FLOAT(NOB-NP)
	SDEV=DSQRT(XZ*SUMB)
      TVAR=1.96+XZ*(2.3779+XZ*(2.7135+XZ*(3.187936+2.466666*XZ**2)))
      IF(NP.EQ.1)WRITE(KP,1042)
      IF(NP.GT.1)WRITE(KP,1015)
      DO 108 I=1,NP
      SECOEF=E(I)*SDEV
      TVALUE=TH(I)/SECOEF
      TSEC=TVAR*SECOEF
      TMCOE=TH(I)-TSEC
      TPCOE=TH(I)+TSEC
      IF(NP.EQ.1) WRITE(KP,1043) BI(I),TH(I),SECOEF,TMCOE,TPCOE
      IF(NP.GT.1) WRITE(KP,1016) BI(I),TH(I),SECOEF,TVALUE,TMCOE,TPCOE
  108 CONTINUE
C  
C   ------ COMMENT FOR THE NONEQUILIBRIUM CDE ---------
C     FOR BETA = 0.9999  
           IF(MODE.EQ.2) THEN 
            IF(INDEX(4).EQ.1.AND.ABS(TB(NB)-0.9999).LT.1.E-5) THEN
             WRITE(KP,1045)
            END IF
           END IF
C  FOR OMEGA = OMMAX (=100)
         IF(ILMT.NE.0) GO TO 250
         IF(MODE.EQ.2) THEN 
           IF(INDEX(5).EQ.1.AND.ABS(TB(NOMEG)-OMMAX).LT.
     &          OMMAX*1.E-04) THEN
             WRITE(KP,1047)
           END IF
         END IF
         IF(MODE.EQ.4.OR.MODE.EQ.6) THEN
           IF(INDEX(4).EQ.1.AND.(TB(NOMEG)-OMMAX).LT.OMMAX*1.E-04) THEN
            WRITE(KP,1047)
           END IF
         END IF   
         IF(MODE.EQ.8.AND.INDEX(4).EQ.1) THEN 
            DOMEGA=TB(NOMEG)*RHOTH*TB(NVAR+4)*ZL/TB(NVAR+1) 
         IF((DOMEGA-OMMAX).LT.OMMAX*1.E-04) THEN
            WRITE(KP,1049)
         END IF
         END IF  
C     ----- PREPARE FINAL OUTPUT -----
  250 LSORT(1)=1
      DO 115 J=2,NOB
      TEMP=R(J)
      K=J-1
      DO 111 L=1,K
      LL=LSORT(L)
      IF(TEMP-R(LL)) 112,112,111
  111 CONTINUE
      LSORT(J)=J
      GO TO 115
  112 KK=J
  113 KK=KK-1
      LSORT(KK+1)=LSORT(KK)
      IF(KK-L) 114,114,113
  114 LSORT(L)=J
  115 CONTINUE
      WRITE(KP,1017)
      DO 116 I=1,NOB
  116 WRITE(KP,1018) I,Z(I),T(I),C(I),F(I),R(I)
C      WRITE(KP,1019)
C      DO 117 I=1,NOB
C      J=LSORT(NOB+1-I)
C  117 WRITE(KP,1018) J,Z(J),T(J),C(J),F(J),R(J)
      GO TO 150
  140 WRITE(KP,1030)
      DO 145 I=1,NOB
  145 WRITE(KP,1033)I,Z(I),T(I),F(I)
  150 CONTINUE
C      WRITE(KP,'(A)') CHAR(12)
  154 CONTINUE
C
  155 CONTINUE
C
C     ----- END OF PROBLEM -----
 1000 FORMAT(' Enter input file name (default = CXTFIT2.IN)')
 1001 FORMAT(A15)
 1002 FORMAT(' Enter output file name (default = CXTFIT2.OUT)')
 1011 FORMAT(//5X,'ITER',5X,'SSQ',3X,5(4X,A6))
 1012 FORMAT(4X,I3,1X,E12.4,5(E10.3))
 1013 FORMAT(///,5X,'COVARIANCE MATRIX FOR FITTED PARAMETERS'/5X,39(1H=)
     &/15X,10(A6,1X))
 1014 FORMAT(8X,A6,10(F7.3))
 1015 FORMAT(5X,'NON-LINEAR LEAST SQUARES ANALYSIS, FINAL RESULTS'
     1/5X,48(1H=)//48X,'95% CONFIDENCE LIMITS'/5X,1X,'NAME',6X,'VALUE',
     25X,'S.E.COEFF.',1X,'T-VALUE',6X,'LOWER',8X,'UPPER')
 1016 FORMAT(5X,A6,E11.4,E12.4,1X,E9.4,E12.4,E13.4)
 1017 FORMAT(//5X,18(1H-),'ORDERED BY COMPUTER INPUT',19(1H-)/40X,
     1'CONCENTRATION',9X,'RESI-'/'$',4X,'NO',4X,'DISTANCE',7X,'TIME',
     28X,'OBS',8X,'FITTED',7X,'DUAL')
 1018 FORMAT(4X,I3,F12.4,4F12.4)
 1019 FORMAT('$',//5X,18(1H-),'ORDERED BY RESIDUAL',25(1H-)/40X,
     1'CONCENTRATION',9X,'RESI-'/5X,'NO',4X,'DISTANCE',7X,'TIME',
     28X,'OBS',8X,'FITTED',7X,'DUAL')
 1020 FORMAT(///,5X,'COVARIANCE MATRIX FOR FITTED PARAMETERS'/5X,39(1H=)
     &/15X,' OUT OF RANGE !!')
 1030 FORMAT(//5X,6(1H-),'RESULTS FOR INITIAL COEFFICIENT VALUES',6(1H-)
     1/5X,'NO',8X,'DISTANCE',9X,'TIME',6X,'CONCENTRATION')
 1033 FORMAT(4X,I3,3(3X,F12.4))
 1034 FORMAT(/5X,'CONVERGENCE CRITERIA NOT MET IN',I4,' ITERATIONS')
 1035 FORMAT(/5X,'NO FURTHER DECREASE IN SSQ OBTAINED FROM ',I3,' TO ',
     &I3,' ITERATIONS')
 1038 FORMAT(/5X,'NO FURTHER DECREASE IN SSQ OBTAINED AFTER ',I3,' TRIAL
     1S')
 1041 FORMAT(/7X,'RSQUARE FOR REGRESSION OF OBSERVED VS PREDICTED =',
     1F10.8,/, 10X,'(COEFFICIENT OF DETERMINATION)'/)
 1042 FORMAT(5X,'NON-LINEAR LEAST SQUARES ANALYSIS, FINAL RESULTS'
     1/5X,48(1H=)//47X,'95% CONFIDENCE LIMITS'/5X,'NAME',8X,'VALUE',
     28X,'S.E.COEFF.',8X,'LOWER',10X,'UPPER')
 1043 FORMAT(5X,A6,1X,E13.4,3X,E13.4,1X,E13.4,2X,E13.4)
 1045 FORMAT(/7X,'BETA = 0.9999  THE EQUILIBRIUM CDE SHOULD BE USED ! ')
 1047 FORMAT(/7X,'OMEGA = 100    THE EQUILIBRIUM CDE SHOULD BE USED ! ')
 1049 FORMAT(/7X,'<OMEGA> = 100  THE EQUILIBRIUM CDE SHOULD BE USED ! ')
 1050 FORMAT(7X,'MEAN SQUARE FOR ERROR (MSE) =',1E12.4//)

C
C     ----- CLOSE FILES -----
      CLOSE(5)
      STOP
      END
C     ------------------------------------------------------------------
      SUBROUTINE MATINV(A,NP,B)
C
C     PURPOSE: PERFORM MATRIX INVERSION FOR PARAMETER ESTIMATION
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(15,15),B(15),INDEX(15,2)
      DO 2 J=1,15
    2 INDEX(J,1)=0
      I=0
    4 AMAX=-1.0
      DO 12 J=1,NP
      IF(INDEX(J,1)) 12,6,12
    6 DO 10 K=1,NP
      IF(INDEX(K,1)) 10,8,10
    8 P=DABS(A(J,K))
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
C 