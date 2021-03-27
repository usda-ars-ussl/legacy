C    -----------------------------------------------------------
      SUBROUTINE DETCDE(C1,C2)
C
C   PURPOSE: CALCULATE DETERMINISTIC EQUILIBRIUM & NONEQUILIBRIUM CDE
C
      IMPLICIT REAL*8(A-H,O-Z)
C
C
      COMMON MODC,MM,A,B,DA,CX,P,BETA,BETR,OMEGA,R,DMU1,DMU2,MODE,
     & MCON,TT,ZZ,STOPER,LEVEL,ICHEB
      COMMON/MODAT/INDEX(15),INVERSE,NREDU,NVAR,ZL,MIT,MDEG,MASS,
     & DUMTP(50),DUMGA1(50),MNEQ,ISKIP,PHIM,PHIIM
      COMMON/STOCH/MODD,MODK,MODS,MSTOCH,CORR,MCORR,MSD,SDLNK,
     & V,SDLNV,VMAX,VMIN,DIS,AVEY,SDLNY,YMAX,YMIN,ALPHA,SDLND
     & ,DK,RHOTH,MD56,MK34,MAL8
      COMMON/BOUN/MODB,NPULSE,TPULSE(50),PULSE(50),MASSST
      COMMON/INITI/MODI,NINI,CINI(50),ZINI(50)
      COMMON/PROD/MODP,NPRO1,NPRO2,GAMMA1(50),ZPRO1(50),GAMMA2(50),
     & ZPRO2(50)
C
C
       C1=0.0
       C2=0.0
       BETR=BETA*R
       MMODE=MOD(MODE,2)
       IF(MMODE.EQ.1.OR.(OMEGA+DMU2).LT.1.D-07) THEN
          A=0.0
          DA=0.0
          CX=0.0
        ELSE
          A=OMEGA*OMEGA/(OMEGA+DMU2)/BETR
          DA=OMEGA*DMU2/(OMEGA+DMU2)
          CX=OMEGA/(OMEGA+DMU2)
       END IF
      IF(MMODE.EQ.1.OR.
     &     (BETA.GE.0.9999999.AND.OMEGA.LT.1.D-07)) THEN
        B=0.0
      ELSE
        IF((OMEGA+DMU2).LT.1.D-07) THEN
         B=0.0
        ELSE
         B=(OMEGA+DMU2)/(R-BETR)
       END IF
      END IF
C    --------- BOUNDARY VALUE PROBLEM ---------------
      IF (MODB.EQ.0) GO TO 10
      CALL BOUND(CBOU1,CBOU2)
      GO TO 15
   10 CBOU1=0.0
      CBOU2=0.0
C    --------  INITIAL VALUE PROBLEM ----------
   15 IF (MODI.EQ.0) GO TO 20
      CALL INITIAL(CINT1,CINT2)
      GO TO 25
   20 CINT1=0.0
      CINT2=0.0
C    --------  PRODUCTION VALUE PROBLEM ----------
   25  IF (MODP.EQ.0) GO TO 30
      CALL PRODUC(CPRO1,CPRO2)
      GO TO 35
   30 CPRO1=0.0
      CPRO2=0.0
C
   35 C1=CBOU1+CINT1+CPRO1
      IF(MOD(MODE,2).EQ.1) GO TO 99
      C2=CBOU2+CINT2+CPRO2
C
   99 RETURN
      END
C    -----------------------------------------------------------
      SUBROUTINE BOUND(C1,C2)
C
C       PURPOSE: BOUNDARY VALUE PROBLEM FOR EQUILIBRIUM AND
C                   NONEQUILIBRIUM CDE
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON MODC,MM,A,B,DA,CX,P,BETA,BETR,OMEGA,R,DMU1,DMU2,MODE,
     & MCON,TT,Z,STOPER,LEVEL,ICHEB
      COMMON/MODAT/INDEX(15),INVERSE,NREDU,NVAR,ZL,MIT,MDEG,MASS,
     &        DUMTP(50),DUMGA1(50),MNEQ,ISKIP,PHIM,PHIIM
      COMMON/BOUN/MODB,NPULSE,TPULSE(50),PULSE(50),MASSST
      COMMON/STOCH/MODD,MODK,MODS,MSTOCH,CORR,MCORR,MSD,SDLNK,
     & V,SDLNV,VMAX,VMIN,DIS,AVEY,SDLNY,YMAX,YMIN,ALPHA,SDLND
     & ,DK,RHOTH,MD56,MK34,MAL8
      EXTERNAL CBJ,CBAL,CTRAN,CBEXP,CBIN1,CBIN2
C
C     MODB=0  ZERO INPUT (SOLUTE FREE WATER)
C     MODB=1  DIRAC DELTA INPUT
C     MODB=2  STEP INPUT
C     MODB=3  SINGLE PULSE INPUT
C     MODB=4  MULTIPLE PULSE INPUT
C     MODB=5  EXPONENTIAL INPUT
C     MODB=6  ARBITRARY INPUT (DEFINE FUNCTION CINPUT(T))
C
C       MCON=0,1 CALCULATE PHASE 1 CONC.
C       MCON=0,2 CALCULATE PHASE 2 CONC.
C
C MODJH: Calculation control code for step input
C     MODJH=0; evaluate Eq.(3.23) or (3.24);
C          =1 Eq.(3.21) or (3.22) based on Goldstein's J-function
      MODJH=1
C
      C1=0.0
      C2=0.0
      A1=0.0
      A2=0.0
      T=TT
      IF (TT.GT.1.D-07) GO TO 1
      RETURN
C
    1 GO TO (10,20,20,20,30,40) MODB
C ---- DELTA INPUT -----------------------------------
C  >>>  ADJUSTMENT FOR  DIMENSIONAL DELTA INPUT  <<<<
   10       IF(NREDU.LE.1) THEN
              BMASS=PULSE(1)*V/ZL
            ELSE
              BMASS=PULSE(1)
            END IF
C
        C1=CC0(TT)*DBEXP(-OMEGA*TT/BETR)
        IF(MODE.EQ.1) THEN
           C1=C1*BMASS
           GOTO 99
        END IF
          AP=DSQRT(1.0+P*Z/30.0)
          TMAX=DMIN1(TT,BETR*(Z+60.0*(1.0+AP)/P))
          TMIN=DMAX1(0.D00,BETR*(Z+60.0*(1.0-AP)/P))
        IF(MCON.EQ.0.OR.MCON.EQ.1) THEN
           MC=1
           CALL CHEBYCON(CTRAN,A1,TMIN,TMAX,MC)
           C1=(C1+A1)*BMASS
        END IF
        IF(MCON.EQ.0.OR.MCON.EQ.2) THEN
           MC=2
           CALL CHEBYCON(CTRAN,C2,TMIN,TMAX,MC)
           C2=C2*BMASS
        END IF
        GO TO 99
C
C ---- STEPWISE INPUT -----------------------------
   20 DO 25 I=1,NPULSE
          TT=TT-TPULSE(I)
        IF (TT.LE.0.D00) GO TO 99
C   EQUILIBRIUM CDE
        IF(MODE.EQ.1) THEN
           A1=CC4(TT,DA)
           GOTO 27
        END IF
C   NONEQUILIBRIUM CDE
C     Solution based on J-function Eq.(3.21) and (3.22)
        IF(MODJH.EQ.1) THEN
          AP=DSQRT(1.0+P*Z/30.0)
          TMAX=DMIN1(TT,BETR*(Z+60.0*(1.0+AP)/P))
          TMIN=DMIN1(0.D00,BETR*(Z+60.0*(1.0-AP)/P))
           IF(MCON.EQ.0.OR.MCON.EQ.1) THEN
              MC=1
              CALL CHEBYCON(CBJ,A1,TMIN,TMAX,MC)
             IF(MODC.NE.3.AND.MODC.NE.4) THEN
                IF(Z.LT.1.D-07.AND.TT.GT.0.D00) A1=1.0
             END IF
           END IF
           IF(MCON.EQ.0.OR.MCON.EQ.2) THEN
              MC=2
              CALL CHEBYCON(CBJ,A2,TMIN,TMAX,MC)
           END IF
C   Alternative solution Eq.(3.23) and (3.24)
        ELSE
           IF(MCON.EQ.0.OR.MCON.EQ.1) THEN
              MC=1
              CALL CHEBYCON(CBAL,A11,0.D00,TT,MC)
              A1=CC4(TT,DA)*DBEXP(-A*TT)+A11
           END IF
           IF(MCON.EQ.0.OR.MCON.EQ.2) THEN
              MC=2
              CALL CHEBYCON(CBAL,A2,0.D00,TT,MC)
           END IF
        END IF
C
   27 IF (I.EQ.1) THEN
        C1=PULSE(I)*A1
        C2=PULSE(I)*A2
      ELSE
        C1=C1+(PULSE(I)-PULSE(I-1))*A1
        C2=C2+(PULSE(I)-PULSE(I-1))*A2
      END IF
      TT=T
   25 CONTINUE
      GO TO 99
C
C ----  EXPONENTIAL INPUT  ------------------------
C  CONSTANT TERM
   30 A2=0.0
      IF(MODE.EQ.1) THEN
           A1=CC4(TT,DA)
           GOTO 32
        END IF
C
       IF(MCON.EQ.0.OR.MCON.EQ.1) THEN
           MC=1
           CALL CHEBYCON(CBAL,A11,0.D00,TT,MC)
           A1=CC4(TT,DA)*DBEXP(-A*TT)+A11
        END IF
        IF(MCON.EQ.0.OR.MCON.EQ.2) THEN
            MC=2
            CALL CHEBYCON(CBAL,A2,0.D00,TT,MC)
        END IF
C
   32  C1=PULSE(1)*A1
       C2=PULSE(1)*A2
C
C EXPONENTIAL TERM  TPULSE(1)= lamda^b
      IF(ABS(TPULSE(1)).LT.1.D-07) THEN
         C1=C1+PULSE(2)*A1
         C2=C2+PULSE(2)*A2
         GO TO 99
       END IF
      IF(MODE.EQ.1) THEN
           DUM=-R*TPULSE(1)
           DUM1=1.0+4.0*(DMU1+DUM)/P
           IF(DUM1.LT.1.D-07) THEN
              CALL CHEBYCON(CBEXP,A1,0.D00,TT,MC)
             ELSE
              A1=DBEXP(-TPULSE(1)*TT)*CC4(TT,DUM)
            END IF
              GOTO 35
      END IF
      IF(MCON.EQ.0.OR.MCON.EQ.1) THEN
              MC=1
              CALL CHEBYCON(CBEXP,A1,0.D00,TT,MC)
      END IF
      IF(MCON.EQ.0.OR.MCON.EQ.2) THEN
              MC=2
              CALL CHEBYCON(CBEXP,A2,0.D00,TT,MC)
              A2=OMEGA/(R-BETR)*A2
       END IF
   35    C1=C1+PULSE(2)*A1
         C2=C2+PULSE(2)*A2
C
         GO TO 99
C
C ---- ARBITARY FUNCTION   (SPECIFY CINPUT(T)) ---------------
   40  IF(MCON.EQ.0.OR.MCON.EQ.1) THEN
           CALL CHEBY(CBIN1,C1,0.D00,TT)
             IF(MODC.NE.3.AND.MODC.NE.4) THEN
                IF(Z.LT.1.D-07.AND.TT.GT.0.D00) C1=CINPUT(T)
             END IF
        END IF
C IF YOU DON'T NEED C2 WHEN MODE=2 & MODB=6, USE NEXT LINE TO SKIP C2
C   IT IS NECESARRY TO EVALUATE TRIPLE INTEGRALS IN (3.14)
C   NOTE THAT SIMILAR MODIFICATION IS NEEDED FOR IVP AND PVP
C
C       IF(MODE.EQ.1.OR.MODB.EQ.6) GOTO 99
       IF(MODE.EQ.1) GOTO 99
        IF(MCON.EQ.0.OR.MCON.EQ.2) THEN
           CALL CHEBY(CBIN2,C2,0.D00,TT)
           C2=C2*OMEGA/(R-BETR)
        END IF
C
   99 IF(C1.LT.1.0D-50) C1=0.0
      IF(C2.LT.1.0D-50) C2=0.0
      TT=T
      RETURN
      END
C   -----------------------------------------------
      SUBROUTINE INITIAL(C1,C2)
C
C       PURPOSE: INITIAL VALUE PROBLEM FOR EQUILIBRIUM AND
C                   NONEQUILIBRIUM CDE
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/MODAT/INDEX(15),INVERSE,NREDU,NVAR,ZL,MIT,MDEG,MASS,
     & DUMTP(50),DUMGA1(50),MNEQ,ISKIP,PHIM,PHIIM
      COMMON/STOCH/MODD,MODK,MODS,MSTOCH,CORR,MCORR,MSD,SDLNK,
     & V,SDLNV,VMAX,VMIN,D,AVEY,SDLNY,YMAX,YMIN,ALPHA,SDLND
     & ,DK,RHOTH,MD56,MK34,MAL8
      COMMON MODC,MM,A,B,DA,CX,P,BETA,BETR,OMEGA,R,DMU1,DMU2,MODE,
     & MCON,T,Z,STOPER,LEVEL,ICHEB
      COMMON/INITI/MODI,NINI,CINI(50),ZINI(50)
      EXTERNAL CIVP
C
C      MODI=0   ZERO INITIAL CONCENTRATION
C      MODI=1   CONSTANT INITIAL CONCENTRATION
C      MODI=2   STEPWISE INITIAL CONCENTRATION
C      MODI=3   EXPONENTIAL INITIAL CONCENTRATION
C      MODI=4   DELTA(Z-Z0) & CONSTANT BACKGROUND
C
        C1=0.0
        C2=0.0
        A1=0.0
        A2=0.0
C
       MMODI=MODI
       IF(MODI.EQ.4) THEN
         IF(CINI(1).GT.1.E-7) THEN
           MODI=1
         ELSE
           IF (T.LT.1.D-07) RETURN
           DG=DBEXP(-(DMU1+DA)/BETR*T)
           GOTO 40
         END IF
       END IF
C --- Initial concentration for small T --------
      IF (T.GT.1.D-07) GO TO 11
      GO TO (1,1,9) MODI
    1 IF (NINI.EQ.1) THEN
        C1=CINI(1)
        C2=CINI(1)
      ELSE IF(Z-ZINI(NINI).GE.1.D-08
     &        .OR.DABS(Z-ZINI(NINI)).LE.1.D-07) THEN
           C1=CINI(NINI)
           C2=CINI(NINI)
         ELSE
       DO 2 I=2,NINI
        IF((Z-ZINI(I-1)).GE.1.D-08.AND.
     &         ZINI(I)-Z.GT.1.D-08) THEN
            C1=CINI(I-1)
            C2=CINI(I-1)
           RETURN
           END IF
         IF (DABS(ZINI(I-1)-Z).LT.1.D-07) THEN
            C1=CINI(I-1)
            C2=CINI(I-1)
           RETURN
         END IF
    2 CONTINUE
        END IF
       RETURN
C
    9   C1=CINI(1)+CINI(2)*DBEXP(-ZINI(1)*Z)
        C2=C1
       RETURN
C -------------------------------------------------
   11  DG=DBEXP(-(DMU1+DA)/BETR*T)
      IF(MODE.EQ.1) THEN
        A1=0.0
        GOTO 7
      END IF
      IF(MCON.EQ.0.OR.MCON.EQ.1) THEN
              MC=1
              CALL CHEBYCON(CIVP,A1,0.D00,T,MC)
      END IF
      IF(MCON.EQ.0.OR.MCON.EQ.2) THEN
              MC=2
              CALL CHEBYCON(CIVP,A2,0.D00,T,MC)
      END IF
    7   GO TO (10,10,30) MODI
C ---- UNIFORM OR STEPWISE INITIAL CONDITION ---------
   10 IF(MCON.EQ.2) GO TO 17
      DO 15 I=1,NINI
      IF (I.EQ.1) THEN
        G=-CINI(I)*DG*(CC1(T)-1.0)
C
       ELSE
        G=G+(CINI(I-1)-CINI(I))*DG*(CC2(T,ZINI(I))-1.0)
       END IF
   15 CONTINUE
      C1=G*DBEXP(-A*T)+A1
      IF(MODE.EQ.1) GOTO 40
      IF(MCON.EQ.1) GOTO 40
C
   17 IF (NINI.EQ.1) THEN
        C2=CINI(1)*DBEXP(-B*T)+A2
      ELSE IF(Z.GE.ZINI(NINI)) THEN
        C2=CINI(NINI)*DBEXP(-B*T)+A2
      ELSE
       DO 20 I=2,NINI
          IF(Z.GE.ZINI(I-1).AND.Z.LT.ZINI(I) ) THEN
              C2=CINI(I-1)*DBEXP(-B*T)+A2
              GO TO 40
          END IF
   20  CONTINUE
      END IF
      GO TO 40
C ---- EXPONENTIAL INITIAL CONDITION ------------
   30 IF(MCON.EQ.2) GOTO 32
       C1=(CINI(1)*DG*(1.0-CC1(T))+CINI(2)*DG*CC3(T,ZINI(1)))
     &      *DBEXP(-A*T)+A1
      IF(MODE.EQ.1) GOTO 40
      IF(MCON.EQ.1) GOTO 40
   32   C2=(CINI(1)+CINI(2)*DBEXP(-ZINI(1)*Z))*DBEXP(-B*T)+A2
C
   40 IF(MMODI.NE.4) GOTO 100
C --- DELTA INITIAL CONDITION ---------------------
      MODI=4
      DMASS=CINI(2)
C  >> NOTICE DELTA(z)=DELTA(Z)/ZL
      IF(NREDU.NE.2) DMASS=DMASS/ZL
      IF(DMASS.LT.1.E-10) GOTO 100
      IF(MODE.EQ.1) THEN
        A1=0.0
        GOTO 47
      END IF
      IF(MCON.EQ.0.OR.MCON.EQ.1) THEN
              MC=1
              CALL CHEBYCON(CIVP,A1,0.D00,T,MC)
      END IF
      IF(MCON.EQ.0.OR.MCON.EQ.2) THEN
              MC=2
              CALL CHEBYCON(CIVP,A2,0.D00,T,MC)
      END IF
   47 CONTINUE
C
      IF(MCON.EQ.2) GO TO 49
      IF(MODC.LE.4.AND.ABS(ZINI(2)).LT.1.E-10) THEN
       G=DMASS*DG*CC0(T)*BETA*R
      ELSE
       G=DMASS*DG*CC5(T,ZINI(2))
      END IF
        C1=C1+G*DBEXP(-A*T)+DMASS*A1
      IF(MODE.EQ.1) GOTO 100
      IF(MCON.EQ.1) GOTO 100
C
   49 C2=C2+DMASS*A2
C  Note that IVP and PVP allow negative values,e.g.,negative flux conc.
C           (See Toride et al. 1993 SSSAJ, p1406-)
  100 IF(DABS(C1).LT.1.0D-40) C1=0.0
      IF(DABS(C2).LT.1.0D-40) C2=0.0
      RETURN
      END
C   ------------------------------------------------
      SUBROUTINE PRODUC(C1,C2)
C
C       PURPOSE: PRODUCTION VALUE PROBLEM FOR EQUILIBRIUM AND
C                   NONEQUILIBRIUM CDE
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/STOCH/MODD,MODK,MODS,MSTOCH,CORR,MCORR,MSD,SDLNK,
     & V,SDLNV,VMAX,VMIN,D,AVEY,SDLNY,YMAX,YMIN,ALPHA,SDLND
     & ,DK,RHOTH,MD56,MK34,MAL8
      COMMON MODC,MM,A,B,DA,CX,P,BETA,BETR,OMEGA,R,DMU1,DMU2,MODE,
     &       MCON,T,Z,STOPER,LEVEL,ICHEB
      COMMON/PROD/MODP,NPRO1,NPRO2,GAMMA1(50),ZPRO1(50),GAMMA2(50),
     &   ZPRO2(50)
      EXTERNAL C1PRO,C2PRO
C
C      MODP=0   ZERO PRODUCTION
C      MODP=1   CONSTANT PRODUCTION
C      MODP=2   STEPWISE PRODUCTION
C      MODP=3   EXPONENTIAL PRODUCTION
C
        C1=0.0
        C2=0.0
      IF (T.GT.1.D-07) GO TO 5
      RETURN
C
    5  IF(MODE.EQ.1.AND.(NPRO1.EQ.1.OR.MODP.EQ.1)) THEN
C
C  MODP1:  Calculation control code for constant production
C                     for the equilibrium CDE
C     MODP1=0; evaluate the integral in Eq.(2.32)
C          =1   Eq.(2.33) or (2.34) (current setting)
          MODP1=1
C
          GO TO (6,7) MODP1
    6      IF((OMEGA+DMU1).LT.1.D-07) THEN
             C1=PROD0(T)
             RETURN
          ELSE
           C1=GAMMA1(1)/DMU1*(1-DBEXP(-(DMU1*T)/BETR)*
     &                               (1-CC1(T))-CC4(T,DA))
          RETURN
          END IF
    7     CALL CHEBY(C1PRO,A1,0.D00,T)
          C1=A1
            RETURN
C
        END IF
C  ---- Equilibrium (phase 1) concnetration ------
      IF(MCON.EQ.0.OR.MCON.EQ.1) THEN
        CALL CHEBY(C1PRO,A1,0.D00,T)
           C1=A1
      END IF
        IF(MODE.EQ.1) RETURN
C
C  ---- Nonquilibrium (phase 2) concnetration ------
      IF((OMEGA+DMU2).LT.1.D-07) THEN
c
C
       GO TO (10,10,20) MODP
C
   10 IF (NPRO2.EQ.1) THEN
        C2=GAMMA2(1)*T/(R-BETR)
        GO TO 40
      ELSE IF(Z.GE.ZPRO2(NPRO2)) THEN
        C2=GAMMA2(NPRO2)*T/(R-BETR)
        GO TO 40
         ELSE
           DO 15 I=2,NPRO2
           IF(Z.GE.ZPRO2(I-1).AND.Z.LT.ZPRO2(I) ) THEN
             C2=GAMMA2(I-1)*T/(R-BETR)
             GO TO 40
          END IF
   15  CONTINUE
        END IF
       GO TO 40
C
   20   C2=(GAMMA2(1)+GAMMA2(2)*DBEXP(-ZPRO2(1)*Z))*T/(R-BETR)
   40      RETURN
C
      ELSE
      IF(MCON.EQ.0.OR.MCON.EQ.2) THEN
         CALL CHEBY(C2PRO,A2,0.D00,T)
      END IF
C
      GO TO (50,50,60) MODP
   50  IF (NPRO2.EQ.1) THEN
        C2=GAMMA2(1)/(OMEGA+DMU2)*
     &           (1.0-DBEXP(-(OMEGA+DMU2)*T/(R-BETR)))+A2
      ELSE IF(Z.GE.ZPRO2(NPRO2)) THEN
        C2=GAMMA2(NPRO2)/(OMEGA+DMU2)*
     &           (1.0-DBEXP(-(OMEGA+DMU2)*T/(R-BETR)))+A2
         ELSE
           DO 55 I=2,NPRO2
           IF(Z.GE.ZPRO2(I-1).AND.Z.LT.ZPRO2(I) ) THEN
        C2=GAMMA2(I-1)/(OMEGA+DMU2)*
     &         (1.0-DBEXP(-(OMEGA+DMU2)*T/(R-BETR)))+A2
              RETURN
          END IF
   55  CONTINUE
        END IF
       RETURN
C
   60   GAMMA=GAMMA2(1)+GAMMA2(2)*DBEXP(-ZPRO2(1)*Z)
        C2=GAMMA/(OMEGA+DMU2)*(1.0-DBEXP(-(OMEGA+DMU2)*T/(R-BETR)))+A2
C
      END IF
      RETURN
      END
