C     ------------------------------------------------------------------
      SUBROUTINE DATAIN(NOB,NU1,NU2,T,Z,C,B,BMAX,BMIN,BI,ILMT,KP,MAXOB,
     & MPRINT,DT,DZ,NT,NZ,FILEIN,INRHO)
C
C      PURPOSE; READ INPUT PARAMETERS
C               CALL CHECK TO CHECK INPUT PARAMETERS
C               CALL DATAOUT TO WRITE PARAMETERS
C               CALL DIRECT FOR DIRECT PROBLEMS
C
      IMPLICIT REAL*8(A-H,O-Z)
C
C
      DIMENSION T(MAXOB),Z(MAXOB),C(MAXOB),B(30),BMAX(15),BMIN(15)
      CHARACTER TITLE1*60,TITLE2*60,BI(15)*6,FILEIN*15
      COMMON MODC,MM,A,BB,DA,CX,PEC,BETA,BETR,OMEGA,R,DMU1,DMU2,
     &   MODE,MCON,TT,ZZ,STOPER,LEVEL,ICHEB
      COMMON/MODAT/INDEX(15),INVERSE,NREDU,NVAR,ZL,MIT,MDEG,MASS,
     &        DUMTP(50),DUMGA1(50),MNEQ,ISKIP,PHIM,PHIIM
      COMMON/STOCH/MODD,MODK,MODS,MSTOCH,CORR,MCORR,MSD,SDLNK,
     & V,SDLNV,VMAX,VMIN,DIS,AVEY,SDLNY,YMAX,YMIN,ALPHA,SDLND
     & ,DK,RHOTH,MD56,MK34,MAL8
      COMMON/BOUN/MODB,NPULSE,TPULSE(50),PULSE(50),MASSST
      COMMON/INITI/MODI,NINI,CINI(50),ZINI(50)
      COMMON/PROD/MODP,NPRO1,NPRO2,GAMMA1(50),ZPRO1(50),GAMMA2(50),
     &    ZPRO2(50)
C  ----------------------------------------------------------------
C  INITIALIZE PARAMETERS
       MCON=0
       DO 1 I=1,30
       B(I)=0.0
    1  CONTINUE
       DO 2 I=1,15
       INDEX(I)=0
       BMAX(I)=0.0
    2  BMIN(I)=0.0
       DO 3 I=1,10
       DUMTP(I)=0.0
       DUMGA1(I)=0.0
       TPULSE(I)=0.0
       PULSE(I)=0.0
       CINI(I)=0.0
       ZINI(I)=0.0
       DUMGA1(I)=0.0
       GAMMA1(I)=0.0
       ZPRO1(I)=0.0
       GAMMA2(I)=0.0
       ZPRO2(I)=0.0
    3  CONTINUE
C----- READ MODEL DESCRIPTION (BLOCK A)--------------------------
      READ(5,*)
      READ(5,1001) TITLE1
      READ(5,1001) TITLE2
C
C   MODE=1  EQUILIBRIUM  MODEL, MODE=2 NONEQUILIBRIUM  MODEL
C   MODE>=3 STOCHASTIC   MODEL
C     3; Kd&V EQUILIBRIUM      4; Kd&V NONEQUILIBRIUM
C     5; D&V EQUILIBRIUM      6; D&V NONEQUILIBRIUM
C     8; ALPHA & V NONEQUILIBRIUM
C
C   MODC=1; flux conc. or area-averaged flux conc.
C   MODC=2; field-scale flux conc.
C   MODC=3; third-type resident conc.
C   MODC=4; third-type total resident concentration.
C   MODC=5; first-type resident conc.
C   MODC=6; first-type total resident conc.
C
C   MODC= 1,2,3,5  PHASE 1 CONC. IS USED FOR PARAMETER ESTIMATION.
C
C   NREDU=0  REAL t ,x ,<c2>  FOR MODE=4,6,8
C   NREDU=1; REAL t & x  C2=<c2/k>
C   NREDU=2; DIMENSIONLESS T & Z
C   NREDU=3; DIMENSIONLESS T & REAL x
C
C
      READ(5,*)
      READ(5,*) INVERSE,MODE,NREDU
      IF(MODE.EQ.7) MODE=8
C   [EQUILIBRIUM CDE WITH NREDU=1, ZL=MAX OF Z]
       IF((MOD(MODE,2).EQ.1).AND.NREDU.LE.1) THEN
        READ(5,*)
        READ(5,*) MODC
       ELSE
        READ(5,*)
        READ(5,*) MODC,ZL
       END IF
C  ----- READ PARAMETERS FOR AN INVERSE PROBLEM (BLOCK B) ------------
      MIT=-1
      MNEQ=0
      MDEG=0
      ILMT=0
      IF(INVERSE.EQ.1) THEN
       READ(5,*)
       READ(5,*)
       READ(5,*) MIT,ILMT,MASS
        IF(MOD(MODE,2).EQ.0) THEN
          READ(5,*)
          READ(5,*) MNEQ,MDEG
              IF(ILMT.EQ.0.AND.MNEQ.EQ.3) THEN
                 WRITE(*,1020)
                 WRITE(KP,1020)
                 STOP
              END IF
            IF(MNEQ.EQ.3.OR.(MNEQ.EQ.0.AND.MDEG.GE.2)) THEN
              READ(5,*)
              READ(5,*) PHIM
                PHIIM=1.0-PHIM
            END IF
        END IF
      END IF
C
      IF(INVERSE.LE.0) MASS=0
C   [ MODE=4,6,8; ONE-SITE MODEL]
      IF(MODE.EQ.4) MNEQ=1
      IF(MODE.EQ.6) MNEQ=1
      IF(MODE.EQ.8) MNEQ=1
C
C
C  ----  READ TRANSPORT PARAMETERS  (BLOCK C)--------------
      READ(5,*)
C   NTP: NUMBER OF TRANSPORT PARAMETERS
C   NVAR:TOTAL NUMBER OF PARAMETERS INCLUDING MASS ESTIMATION
      IF(MODE.EQ.1) NTP=4
      IF(MODE.EQ.2) NTP=7
      IF(MODE.EQ.3) NTP=8
      IF(MODE.EQ.4) NTP=10
      IF(MODE.EQ.5) NTP=8
      IF(MODE.EQ.6) NTP=10
      IF(MODE.EQ.8) NTP=11
C
      NVAR=NTP
C
C     ----- COEFFICIENTS NAMES -----
      READ(5,*)
      GOTO (70,71,72,73,72,73,74,73) MODE
   70 BI(1)=' V....'
      BI(2)=' D....'
      BI(3)=' R....'
      BI(4)=' mu...'
      GO TO 20
   71 BI(1)=' V....'
      BI(2)=' D....'
      BI(3)=' R....'
      BI(4)=' beta.'
      BI(5)=' omega'
      BI(6)=' mu1..'
      BI(7)=' mu2..'
      GO TO 20
   72 BI(1)=' <V>..'
      BI(2)=' <D>..'
      BI(3)=' <Kd>.'
      BI(4)=' mu1..'
      BI(5)=' SD.v.'
      BI(6)=' SD.Kd'
      BI(7)=' SD.D.'
       IF(MODE.EQ.3) THEN
         BI(8)='RhovKd'
       ELSE IF(MODE.EQ.5) THEN
         BI(8)=' RhovD'
       END IF
      GO TO 20
   73 BI(1)=' <V>..'
      BI(2)=' <D>..'
      BI(3)=' <Kd>.'
      IF(MODE.NE.8) THEN
         BI(4)=' omega'
      ELSE
         BI(4)=' alpha'
      END IF
      BI(5)=' mu1..'
      BI(6)=' mu2..'
      BI(7)=' SD.v.'
      BI(8)=' SD.Kd'
      BI(9)=' SD.D.'
       IF(MODE.EQ.4) THEN
         BI(10)='RhovKd'
       ELSE IF(MODE.EQ.6) THEN
         BI(10)=' RhovD'
       ELSE IF(MODE.EQ.8) THEN
         BI(10)='SD.alp'
         BI(11)='RhovAl'
       END IF
      GOTO 20
   74 WRITE(*,*) 'WARNING ! PLEASE CHECK MODE.'
      STOP
C
C     ----- READ INITIAL ESTIMATES -----
   20 IF(MODE.LE.2) THEN
        READ(5,*) (B(I),I=1,NTP)
      ELSE
        READ(5,*) (B(I),I=1,NTP),RHOTH
      END IF
C
C     ----- READ INDICES -----
      IF(INVERSE.LE.0) GO TO 21
      READ(5,*) (INDEX(I),I=1,NTP)
C
C     ----- READ CONSTRAINTS ON PARAMETER VALUES -----
      IF(ILMT.EQ.0) GO TO 21
      READ(5,*) (BMIN(I),I=1,NTP)
      READ(5,*) (BMAX(I),I=1,NTP)
C
C     ------- PARAMETERS FOR BOUNDARY VALUE PROBLEMS (BLOCK D)-----------
C     MODB=0  ZERO INPUT (SOLUTE FREE WATER)
C     MODB=1  DIRAC DELTA INPUT
C     MODB=2  STEP INPUT
C     MODB=3  SINGLE PULSE INPUT
C     MODB=4  MULTIPLE PULSE INPUT
C     MODB=5  EXPONENTIAL INPUT
C     MODB=6  ARBITRARY INPUT (DEFINE FUNCTION CINPUT(T))
   21 READ(5,*)
      READ(5,*)
      MASSST=0
      READ(5,*) MODB
      DO 81 I=1,10
      PULSE(I)=0.0
      TPULSE(I)=0.0
   81 CONTINUE
      IF(MODB.EQ.0.OR.MODB.GE.4) MASS=0
      IF(MODB.EQ.6) GO TO 28
      IF(MODB.EQ.0) THEN
           NPULSE=1
           GO TO 28
      END IF
        IF (MODB.EQ.1) THEN
           IF(MODE.LE.2) THEN
              READ(5,*) PULSE(1)
           ELSE
C  MASSST=1 AMOUNT OF SOLUTE IS CONSTANT REGARDELESS IF V
C        =0                  DEPENDS V (pulse*<v>/v)
              READ(5,*) PULSE(1),MASSST
           END IF
            IF(MASS.EQ.1) THEN
              NVAR=NTP+1
              B(NVAR)=PULSE(1)
              BI(NVAR)=' MASS.'
              READ(5,*) INDEX(NVAR)
                IF(ILMT.EQ.1) THEN
                 READ(5,*) BMIN(NVAR)
                 READ(5,*) BMAX(NVAR)
                END IF
             END IF
        ELSE IF (MODB.EQ.2) THEN
           READ(5,*) PULSE(1)
           TPULSE(1)=0.0
           NPULSE=1
           DUMTP(1)=0.0
            IF(MASS.EQ.1) THEN
              NVAR=NTP+1
              B(NVAR)=PULSE(1)
              BI(NVAR)=' Cin..'
              READ(5,*) INDEX(NVAR)
                IF(ILMT.EQ.1) THEN
                 READ(5,*) BMIN(NVAR)
                 READ(5,*) BMAX(NVAR)
                END IF
             END IF
        ELSE IF (MODB.EQ.3) THEN
           IF(MODE.LE.2) THEN
              READ(5,*) PULSE(1),TPULSE(2)
           ELSE
C  MASSST=1 AMOUNT OF SOLUTE IS CONSTANT REGARDLESS OF V
C              DURATION TIME IS <T>*<v>/v
C              [MODE=7 AND 8 IN CXTFIT]
C        =0  DURATION TIME IS CONSTANT FOR ALL STREAM TUNES.
C           AMOUNT OF SOLUTE IS PROPORTIONAL TO V.
C              [MODE=5 AND 6 IN CXTFIT]
C
              READ(5,*) PULSE(1),TPULSE(2),MASSST
           END IF
           TPULSE(1)=0.0
           PULSE(2)=0.0
           NPULSE=2
           DUMTP(1)=0.0
           DUMTP(2)=TPULSE(2)
            IF(MASS.EQ.1) THEN
              NVAR=NTP+2
              B(NVAR-1)=PULSE(1)
              B(NVAR)=DUMTP(2)
              BI(NVAR-1)=' Cin..'
              IF(MASSST.EQ.0) THEN
                 BI(NVAR)=' T2...'
              ELSE IF(MASSST.EQ.1) THEN
                 BI(NVAR)=' Mo.'
              END IF
              READ(5,*) INDEX(NVAR-1),INDEX(NVAR)
                IF(ILMT.EQ.1) THEN
                 READ(5,*) BMIN(NVAR-1),BMIN(NVAR)
                 READ(5,*) BMAX(NVAR-1),BMAX(NVAR)
                END IF
             END IF
        ELSE IF (MODB.EQ.4) THEN
           READ(5,*) NPULSE
           DO 25 I=1,NPULSE
            READ(5,*) PULSE(I),TPULSE(I)
            DUMTP(I)=TPULSE(I)
   25      CONTINUE
        ELSE IF (MODB.EQ.5) THEN
           READ(5,*) PULSE(1),PULSE(2),TPULSE(1)
           DUMTP(1)=TPULSE(1)
      END IF
C
C --- REARRANGE B(I)----------
   28  NU1=NVAR+1
       NU2=NVAR*2
       DO 29 I=1,NVAR
       B(I+NVAR)=B(I)
       B(I)=0.0
   29  CONTINUE
C
C     ------- PARAMETERS FOR INITIAL VALUE PROBLEMS (BLOCK E)---------
C      MODI=0   ZERO INITIAL CONCENTRATION
C      MODI=1   CONSTANT INITIAL CONCENTRATION
C      MODI=2   STEPWISE INITIAL CONCENTRATION
C      MODI=3   EXPONENTIAL INITIAL CONCENTRATION
C      MODI=4   DELTA(Z-Z1) & CONSTANT INITIAL CONCENTRATION
      READ(5,*)
      READ(5,*)
      READ(5,*) MODI
      DO 82 I=1,10
      CINI(I)=0.0
      ZINI(I)=0.0
   82 CONTINUE
      IF(MODI.EQ.0) THEN
        NINI=1
        GO TO 31
      END IF
      IF(MODI.EQ.1) THEN
         READ(5,*) CINI(1)
        ZINI(1)=0.0
        NINI=1
      ELSE IF (MODI.EQ.2) THEN
         READ(5,*) NINI
         DO 26 I=1,NINI
         READ(5,*) CINI(I),ZINI(I)
   26 CONTINUE
       ELSE IF (MODI.EQ.3) THEN
          READ(5,*) CINI(1),CINI(2),ZINI(1)
       ELSE IF (MODI.EQ.4) THEN
          NINI=1
          ZINI(1)=0
C DELTA CONDITION;  MASS:CINI(2)
C      Z1 (position of a spike):ZINI(2)  UNIFORM CONC: CINI(1)
          READ(5,*) CINI(2),ZINI(2),CINI(1)
      END IF
C
C     ------- PARAMETERS FOR PRODUCTION VALUE PROBLEMS (BLOCK F)
C      MODP=0   ZERO PRODUCTION
C      MODP=1   CONSTANT PRODUCTION
C      MODP=2   STEPWISE PRODUCTION
C      MODP=3   EXPONENTIAL PRODUCTION
C
C      MPRO=0; SAME CONDITION FOR PHASE 1 & 2
C      MPRO=1; DIFFERENT CONDITION FOR PHASE 1 & 2
   31 READ(5,*)
      READ(5,*)
      MPRO=0
      READ(5,*) MODP
      IF(MODP.EQ.0) THEN
         NPRO1=1
         NPRO2=1
         GO TO 32
      ELSE IF (MODP.EQ.1) THEN
         READ(5,*) GAMMA1(1)
         NPRO1=1
         NPRO2=1
         DUMGA1(1)=GAMMA1(1)
          IF(MOD(MODE,2).EQ.0) THEN
            READ(5,*) GAMMA2(1)
          END IF
       ELSE IF (MODP.EQ.2) THEN
           READ(5,*) MPRO
           READ(5,*) NPRO1
            DO 27 I=1,NPRO1
            READ(5,*) GAMMA1(I),ZPRO1(I)
            DUMGA1(I)=GAMMA1(I)
   27    CONTINUE
           IF(MPRO.EQ.0) THEN
              NPRO2=NPRO1
            DO 38 I=1,NPRO2
             GAMMA2(I)=GAMMA1(I)
             ZPRO2(I)=ZPRO1(I)
   38       CONTINUE
           END IF
           IF(MOD(MODE,2).EQ.0.AND.MPRO.EQ.1) THEN
            READ(5,*) NPRO2
            DO 39 I=1,NPRO2
            READ(5,*) GAMMA2(I),ZPRO2(I)
   39    CONTINUE
            END IF
        ELSE IF (MODP.EQ.3) THEN
            READ(5,*) GAMMA1(1),GAMMA1(2),ZPRO1(1)
            DUMGA1(1)=GAMMA1(1)
            DUMGA1(2)=GAMMA1(2)
           IF(MPRO.EQ.0) THEN
             GAMMA2(1)=GAMMA1(1)
             GAMMA2(2)=GAMMA1(2)
             ZPRO2(1)=ZPRO1(1)
             DUMPLA2=ZPRO1(1)
           END IF
           IF(MOD(MODE,2).EQ.0.AND.MPRO.EQ.1) THEN
              READ(5,*) GAMMA2(1),GAMMA2(2),ZPRO2(1)
           END IF
      END IF
C
C
   32 IF(INVERSE.EQ.1) THEN
C ----- READ DATA FOR INVERSE PROBLEM (BLOCK G) ----------
C     INPUTM = 0  Z(I), T(I), C(I)
C            = 1  T(I), C(I) FOR SAME Z (Breakthrough Curve)
C            = 2  Z(I), C(I) FOR SAME T
C            = 3  C(I), Z(I), T(I) FIXED FORMAT, 3F10.0, FOR CXTFIT DATA
C                    AT THE END OF DATA SET, GIVE "0,0,0"
C                    TO INDICATE THE LAST LINE
        READ(5,*)
        READ(5,*)
        NOB=0
        READ(5,*) INPUTM
         IF(INPUTM.EQ.1.OR.INPUTM.EQ.2) READ(5,*)  DUMTZ
        READ(5,*)
        DO 10 I=1,MAXOB
         IF(INPUTM.EQ.1) THEN
          IF(MIT.GE.1) THEN
            READ(5,*) T(I),C(I)
            Z(I)=DUMTZ
            IF(ABS(T(I)).LT.1.E-7.AND.ABS(C(I)).LT.1.E-7) THEN
             IF(NOB.GE.1) GOTO 11
             NOB=NOB+1
            ELSE
              NOB=NOB+1
            END IF
          ELSE IF(MIT.EQ.0) THEN
            READ(5,*) T(I)
            Z(I)=DUMTZ
            IF(ABS(T(I)).LT.1.E-10.AND.NOB.GE.1) THEN
              GOTO 11
            ELSE
              NOB=NOB+1
            END IF
          END IF
         ELSE IF (INPUTM.EQ.2) THEN
          IF(MIT.GE.1) THEN
            READ(5,*) Z(I),C(I)
            T(I)=DUMTZ
            IF(ABS(Z(I)).LT.1.E-7.AND.ABS(C(I)).LT.1.E-7) THEN
             IF(NOB.GE.1) GOTO 11
             NOB=NOB+1
            ELSE
              NOB=NOB+1
            END IF
          ELSE IF(MIT.EQ.0) THEN
            READ(5,*) Z(I)
            T(I)=DUMTZ
            IF(ABS(Z(I)).LT.1.E-10.AND.NOB.GE.1) THEN
              GOTO 11
            ELSE
              NOB=NOB+1
            END IF
          END IF
         ELSE IF (INPUTM.EQ.0) THEN
          IF(MIT.GE.1) THEN
            READ(5,*) Z(I),T(I), C(I)
            IF(ABS(Z(I)).LT.1.E-7.AND.ABS(T(I)).LT.1.E-7.AND.
     &        ABS(C(I)).LT.1.E-7) THEN
             IF(NOB.GE.1) GOTO 11
             NOB=NOB+1
            ELSE
              NOB=NOB+1
            END IF
          ELSE IF(MIT.EQ.0) THEN
            READ(5,*) Z(I),T(I)
            IF(ABS(Z(I)).LT.1.E-7.AND.ABS(T(I)).LT.1.E-7.AND.
     &        NOB.GE.1) THEN
              GOTO 11
            ELSE
              NOB=NOB+1
            END IF
          END IF
         ELSE
            READ(5,1010) C(I),Z(I),T(I)
            IF(ABS(Z(I)).LT.1.E-7.AND.ABS(T(I)).LT.1.E-7.AND.
     &        ABS(C(I)).LT.1.E-7) THEN
             IF(NOB.GE.1) GOTO 11
             NOB=NOB+1
            ELSE
              NOB=NOB+1
            END IF
         END IF
   10   CONTINUE
   11   END IF
      IF(INVERSE.LE.0) THEN
C-------T & Z FOR DIRECT PROBLEM (BLOCK H) ---------------
C   MPRINT=1 PRINT CONC VS. TIME
C   MPRINT=2       CONC VS.DEPTH
        READ(5,*)
        READ(5,*)
       READ(5,*) NZ,DZ,ZI,NT,DT,TI,MPRINT
C
         Z(1)=ZI
         T(1)=TI
      DO 33 I=1,NZ-1
   33  Z(I+1)=Z(I)+DZ
      DO 34 I=1,NT-1
   34  T(I+1)=T(I)+DT
C
      END IF
C   [EQUILIBRIUM CDE WITH NREDU=1, ZL=MAX OF Z]
       IF(MOD(MODE,2).EQ.1.AND.NREDU.LE.1) THEN
          IF(INVERSE.LE.0) THEN
             ZL=Z(NZ)
          ELSE
             ZL=0.0
             DO 37 I=1,NOB
   37        ZL=DMAX1(ZL,Z(I))
          END IF
        END IF
        IF(DABS(ZL).LT.1.D-20) ZL=1.0
C
C  [RETARDATION FACTOR  FOR STOCASTIC MODEL]
      IF(MODE.GE.3) R=1.0+RHOTH*B(NVAR+3)
C ------  CHECK INPUT DATA  -------------------------------
C
      CALL CHECK(ILMT,B,BMAX,BMIN,BI,KP,INRHO)
C
C ------  WRITE INPUT DATA TO FILE ------------------------
C
      CALL DATAOUT(TITLE1,TITLE2,ILMT,B,BMAX,BMIN,BI,KP,FILEIN)
C
C--------CHANGE FROM ACTUAL TO REDUCED DISTSNCE OF INITIAL
C                 & PRODUCTION DISTRIBUTION -------------------
      IF(NREDU.LE.1.OR.NREDU.EQ.3) THEN
       V=B(NVAR+1)
C
        IF(MODI.EQ.1.OR.MODI.EQ.2) THEN
          DO 52 I=1,NINI
           ZINI(I)=ZINI(I)/ZL
   52     CONTINUE
        END IF
        IF(MODI.EQ.3) ZINI(1)=ZINI(1)*ZL
        IF(MODI.EQ.4) ZINI(2)=ZINI(2)/ZL
C
        IF(MODP.EQ.1.OR.MODP.EQ.2) THEN
          DO 53 I=1,NPRO1
           ZPRO1(I)=ZPRO1(I)/ZL
   53     CONTINUE
          DO 55 I=1,NPRO2
           ZPRO2(I)=ZPRO2(I)/ZL
   55     CONTINUE
       END IF
       IF(MODP.EQ.3) THEN
           ZPRO1(1)= ZPRO1(1)*ZL
           ZPRO2(1)= ZPRO2(1)*ZL
       END IF
      END IF
C
C
 1001 FORMAT(A60)
 1010 FORMAT(7F10.0)
 1020 FORMAT(5X' ERROR ! ILMT SHOULD BE 1 FOR MNEQ =3')
C
      RETURN
      END
C
C    -----------------------------------------------------------
      SUBROUTINE CHECK(ILMT,B,BMAX,BMIN,BI,KP,INRHO)
C
C       PURPOSE: CHECK INPUT DATA
C
      IMPLICIT REAL*8(A-H,O-Z)
C
C
      CHARACTER BI(15)*6
      DIMENSION B(30),BMAX(15),BMIN(15)
      COMMON MODC,MM,A,BB,DA,CX,P,BETA,BETR,OMEGA,R,DMU1,DMU2,MODE,
     &   MCON,TT,ZZ,STOPER,LEVEL,ICHEB
      COMMON/MODAT/INDEX(15),INVERSE,NREDU,NVAR,ZL,MIT,MDEG,MASS,
     &     DUMTP(50),DUMGA1(50),MNEQ,ISKIP,PHIM,PHIIM
      COMMON/STOCH/MODD,MODK,MODS,MSTOCH,CORR,MCORR,MSD,SDLNK,
     & V,SDLNV,VMAX,VMIN,D,AVEY,SDLNY,YMAX,YMIN,ALPHA,SDLND
     & ,DK,RHOTH,MD56,MK34,MAL8
      COMMON/BOUN/MODB,NPULSE,TPULSE(50),PULSE(50),MASSST
      COMMON/INITI/MODI,NINI,CINI(50),ZINI(50)
      COMMON/PROD/MODP,NPRO1,NPRO2,GAMMA1(50),ZPRO1(50),GAMMA2(50),
     &        ZPRO2(50)
C
C --- CHECK TRANSPORT PARAMETER VALUES -----------------------
      NERR=0
C ---- DELTA IVP FOR FIRST-TYPE --------------------------------
C  Corrected on 2/12/05
       IF(MODC.GE.5.AND.MODI.EQ.4) THEN
         IF(ABS(ZINI(2)).LE.1.E-10) THEN
           WRITE(*,*) ' WARNING ! Z0 > 0 FOR DIRAC IC WITH FIRST TYPE IN
     &LET'
           WRITE(KP,*) ' WARNING ! Z0 > 0 FOR DIRAC IC WITH FIRST TYPE I
     &NLET'
           NERR=NERR+1
          END IF
	  END IF
      GOTO (10,20,30,40,30,40,50,40) MODE
C
C  ----- MODE =1 -------------------------
   10 DO 11 I=2,3
       IF(B(NVAR+I).LE.1.D-09) THEN
           WRITE(*,1000) BI(I),B(NVAR+I)
           WRITE(KP,1000) BI(I),B(NVAR+I)
           NERR=NERR+1
       END IF
   11 CONTINUE
      IF(B(NVAR+4).LT.0.0D+00) THEN
           WRITE(*,1002) BI(4),B(NVAR+4)
           WRITE(KP,1002) BI(4),B(NVAR+4)
           NERR=NERR+1
       END IF
       GOTO 100
C
C  -----  MODE 2 --------------------------
   20 IF(B(NVAR+4).GE.1.0D+00.OR.B(NVAR+4).LE.1.D-09) THEN
           WRITE(*,1004) BI(4),B(NVAR+4)
           WRITE(KP,1004) BI(4),B(NVAR+4)
           NERR=NERR+1
      END IF
      DO 21 I=2,5
         IF(I.EQ.4) GO TO 21
         IF(B(NVAR+I).LE.1.D-09) THEN
             WRITE(*,1000) BI(I),B(NVAR+I)
             WRITE(KP,1000) BI(I),B(NVAR+I)
             NERR=NERR+1
         END IF
   21 CONTINUE
      DO 22 I=6,7
         IF(B(NVAR+I).LT.0.0D+00) THEN
           WRITE(*,1002) BI(I),B(NVAR+I)
           WRITE(KP,1002) BI(I),B(NVAR+I)
           NERR=NERR+1
         END IF
   22 CONTINUE
      IF(B(NVAR+5).GT.100.0) THEN
           WRITE(*,1042)
           WRITE(KP,1042)
      END IF
       GOTO 100
C
C  -------- MODE 3,5 (STOCHASTIC EQUILIBRIUM CDE) ---------------------
   30 IF(B(NVAR+2).LE.1.D-09) THEN
             WRITE(*,1000) BI(2),B(NVAR+2)
             WRITE(KP,1000) BI(2),B(NVAR+2)
             NERR=NERR+1
         END IF
      DO 32 I=3,7
      IF(B(NVAR+I).LT.0.0D+00) THEN
         WRITE(*,1002) BI(I),B(NVAR+I)
         WRITE(KP,1002) BI(I),B(NVAR+I)
         NERR=NERR+1
      END IF
   32 CONTINUE
      IF(ABS(B(NVAR+3)).LE.1.D-05.AND.B(NVAR+6).GT.1.D-05) THEN
          WRITE(*,1036)
          WRITE(KP,1036)
          NERR=NERR+1
      END IF
      IF(ABS(B(NVAR+8)).GT.1.0D+00) THEN
         WRITE(*,1006) BI(8),B(NVAR+8)
         WRITE(KP,1006) BI(8),B(NVAR+8)
         NERR=NERR+1
      END IF
      IF(RHOTH.LT.0.0D+00) THEN
         WRITE(*,1001) RHOTH
         WRITE(KP,1001) RHOTH
         NERR=NERR+1
      END IF
      IF(MODE.EQ.3) THEN
        IF(B(NVAR+5).LT.1.D-07.AND.B(NVAR+6).LT.1.D-07) THEN
         WRITE(*,1007)
         WRITE(KP,1007)
         NERR=NERR+1
        END IF
        IF(B(NVAR+5).LT.1.D-07.AND.B(NVAR+7).GE.1.D-07) THEN
         WRITE(*,1020)
         WRITE(KP,1020)
         B(NVAR+7)=0.0
        END IF
      ELSE IF(MODE.EQ.5) THEN
        IF(B(NVAR+5).LT.1.D-07.AND.B(NVAR+7).LT.1.D-07) THEN
         WRITE(*,1008)
         WRITE(KP,1008)
         NERR=NERR+1
        END IF
        IF(B(NVAR+5).LT.1.D-07.AND.B(NVAR+6).GE.1.D-07) THEN
         WRITE(*,1022)
         WRITE(KP,1022)
         B(NVAR+6)=0.0
        END IF
      END IF
      GOTO 100
C
C   ------MODE 4,6,8 (STOCHASTIC ONE-SITE MODEL)-----------------
   40 IF (MODB.EQ.6) THEN
        WRITE(*,1024)
        WRITE(KP,1024)
        NERR=NERR+1
      END IF
      DO 41 I=2,4
         IF(I.EQ.3) GOTO 41
         IF(B(NVAR+I).LE.1.D-09) THEN
             WRITE(*,1000) BI(I),B(NVAR+I)
             WRITE(KP,1000) BI(I),B(NVAR+I)
             NERR=NERR+1
         END IF
   41 CONTINUE
      IF(B(NVAR+3).LT.0.0D+00) THEN
         WRITE(*,1002) BI(3),B(NVAR+3)
         WRITE(KP,1002) BI(3),B(NVAR+3)
         NERR=NERR+1
      END IF
      DO 42 I=5,9
      IF(B(NVAR+I).LT.0.0D+00) THEN
         WRITE(*,1002) BI(I),B(NVAR+I)
         WRITE(KP,1002) BI(I),B(NVAR+I)
         NERR=NERR+1
      END IF
   42 CONTINUE
      IF(ABS(B(NVAR+3)).LE.1.D-05.AND.B(NVAR+8).GT.1.D-05) THEN
          WRITE(*,1036)
          WRITE(KP,1036)
          NERR=NERR+1
      END IF
      IF(RHOTH.LT.0.0D+00) THEN
         WRITE(*,1001) RHOTH
         WRITE(KP,1001) RHOTH
         NERR=NERR+1
      END IF
      IF(MODE.LE.6) THEN
        IF(ABS(B(NVAR+10)).GT.1.0D+00) THEN
           WRITE(*,1006) BI(10),B(NVAR+10)
           WRITE(KP,1006) BI(10),B(NVAR+10)
           NERR=NERR+1
        END IF
      ELSE IF(MODE.EQ.8) THEN
        IF(B(NVAR+10).LT.0.0D+00) THEN
           WRITE(*,1002) BI(10),B(NVAR+10)
           WRITE(KP,1002) BI(10),B(NVAR+10)
           NERR=NERR+1
        END IF
        IF(ABS(B(NVAR+11)).GT.1.0D+00) THEN
           WRITE(*,1006) BI(11),B(NVAR+11)
           WRITE(KP,1006) BI(11),B(NVAR+11)
           NERR=NERR+1
        END IF
      END IF
C
      IF(MODE.EQ.4) THEN
      IF (B(NVAR+7).LT.1.D-07.AND.B(NVAR+8).LT.1.D-07) THEN
         WRITE(*,1007)
         WRITE(KP,1007)
         NERR=NERR+1
      END IF
        IF(B(NVAR+7).LT.1.D-07.AND.B(NVAR+9).GE.1.D-07) THEN
         WRITE(*,1020)
         WRITE(KP,1020)
         B(NVAR+9)=0.0
        END IF
      ELSE IF(MODE.EQ.6) THEN
      IF (B(NVAR+7).LT.1.D-07.AND.B(NVAR+9).LT.1.D-07) THEN
         WRITE(*,1008)
         WRITE(KP,1008)
         NERR=NERR+1
      END IF
        IF(B(NVAR+7).LT.1.D-07.AND.B(NVAR+8).GE.1.D-07) THEN
         WRITE(*,1022)
         WRITE(KP,1022)
         B(NVAR+8)=0.0
        END IF
      ELSE IF(MODE.EQ.8) THEN
      IF (B(NVAR+7).LT.1.D-07.AND.B(NVAR+10).LT.1.D-07) THEN
         WRITE(*,1009)
         WRITE(KP,1009)
         NERR=NERR+1
      END IF
        IF(B(NVAR+7).LT.1.D-07.AND.B(NVAR+9).GE.1.D-07) THEN
         WRITE(*,1020)
         WRITE(KP,1020)
         B(NVAR+9)=0.0
        END IF
        IF(B(NVAR+7).LT.1.D-07.AND.B(NVAR+8).GE.1.D-07) THEN
         WRITE(*,1022)
         WRITE(KP,1022)
         B(NVAR+8)=0.0
        END IF
      END IF
      IF(MODE.EQ.4.OR.MODE.EQ.6) THEN
        IF(B(NVAR+4).GT.100.0) THEN
           WRITE(*,1042)
           WRITE(KP,1042)
        END IF
      ELSE IF(MODE.EQ.8) THEN
         DOMEGA=B(NVAR+4)*RHOTH*B(NVAR+3)*ZL/B(NVAR+1)
        IF(DOMEGA.GT.100.0) THEN
           WRITE(*,1044)
           WRITE(KP,1044)
        END IF
      END IF
      GOTO 100
   50 NERR=NERR+1
C ------CHECK MAX & MIN COMSTRAINTS  ----------------------------
  100   INRHO=0
        IF(MODE.EQ.3.OR.MODE.EQ.5) THEN
           IF(INDEX(8).EQ.1) INRHO=1
        ELSE IF (MODE.EQ.4.OR.MODE.EQ.6) THEN
           IF(INDEX(10).EQ.1) INRHO=1
        ELSE IF (MODE.EQ.8) THEN
           IF(INDEX(11).EQ.1) INRHO=1
        END IF
C
        IF(MIT.LE.0) GO TO 102
         DO 52 I=1,NVAR
          IF(ILMT.NE.1) GO TO 52
          IF(INDEX(I).EQ.0) GO TO 52
          IF(I.EQ.NVAR.AND.INRHO.EQ.1) THEN
            IF(DABS(BMAX(I)-BMIN(I)).LT.1.E-10) GOTO 52
            IF(BMIN(I).LT.-1) THEN
               WRITE(*,1037) BI(I)
               WRITE(KP,1037) BI(I)
               BMIN(I)=-1.D+00
            END IF
            IF(BMAX(I).GT.1) THEN
               WRITE(*,1038) BI(I)
               WRITE(KP,1038) BI(I)
               BMAX(I)= 1.D+00
            END IF
            GOTO 56
          END IF
            IF(BMIN(I).LT.0.0) THEN
               WRITE(*,1026) BI(I)
               WRITE(KP,1026) BI(I)
               BMIN(I)=1.D-25
            END IF
   56       IF(BMAX(I).LT.BMIN(I)) THEN
               WRITE(*,1025) BI(I)
               WRITE(KP,1025) BI(I)
               NERR=NERR+1
            END IF
   52  CONTINUE
C  CONSTRAINT FOR BETA
C     BETA SHOULD BE LESS THAN 1 FOR THE NONEQUILIBRIUM CDE
      IF(MODE.NE.2) GOTO 92
      IF(ILMT.NE.1) GOTO 92
      IF(INDEX(4).NE.1) GOTO 92
C  CHECK BMAX
      IF(MNEQ.LE.2) THEN
         IF(BMAX(4).LT.0.9999) GOTO 90
         IF(BMAX(4).GT.0.9999) THEN
           WRITE(*,1027)
           WRITE(KP,1027)
           BMAX(4)=0.9999
           GOTO 90
         END IF
      ELSE IF(MNEQ.EQ.3) THEN
         IF(INDEX(3).EQ.1) THEN
               WRITE(*,1040)
               WRITE(KP,1040)
         END IF
         BBMAX=(PHIM+B(NVAR+3)-1.0)/B(NVAR+3)
         IF(BMAX(4).LE.BBMAX) THEN
           GOTO 90
         ELSE
           WRITE(*,1031)
           WRITE(KP,1031)
           IF(INDEX(3).EQ.0) BMAX(4)=BBMAX
              IF(BMAX(4).GT.0.9999) THEN
               WRITE(*,1027)
               WRITE(KP,1027)
               BMAX(4)=0.9999
             END IF
         END IF
      END IF
C   CHECK BMIN
   90 IF(MNEQ.LE.1) GOTO 92
      IF(MNEQ.EQ.2) THEN
       IF(BMIN(4).GE.1.0/B(NVAR+3)) THEN
          GOTO 92
       ELSE
          WRITE(*,1029)
          WRITE(KP,1029)
       END IF
      ELSE IF(MNEQ.EQ.3) THEN
       IF(BMIN(4).GE.PHIM/B(NVAR+3)) THEN
          GOTO 92
       ELSE
          WRITE(*,1028)
          WRITE(KP,1028)
          IF(INDEX(3).EQ.0) BMIN(4)=PHIM/B(NVAR+3)
       END IF
      END IF
C  CONSTRAINT FOR OMEGA
C     OMEGA SHOULD BE LESS THAN 1 FOR THE NONEQUILIBRIUM CDE
   92 IF(MODE.EQ.2.AND.ILMT.EQ.1) THEN
      IF(INDEX(5).NE.1) GOTO 102
         IF(BMAX(5).GT.100.0) THEN
           WRITE(*,1046)
           WRITE(KP,1046)
           GOTO 102
         END IF
      ELSE IF(MODE.EQ.4.OR.MODE.EQ.6) THEN
         IF(INDEX(4).NE.1.OR.ILMT.EQ.0) GOTO 102
         IF(BMAX(4).GT.100.0) THEN
           WRITE(*,1046)
           WRITE(KP,1046)
           GOTO 102
         END IF
      ELSE IF(MODE.EQ.8) THEN
         IF(INDEX(4).NE.1.OR.ILMT.EQ.0) GOTO 102
          DUMAX=BMAX(4)/RHOTH/B(NVAR+4)/ZL*B(NVAR+1)
         IF(DUMAX.GT.100) THEN
           WRITE(*,1048)
           WRITE(KP,1048)
         END IF
      END IF
C
C  --- STOP PROGRAM IF NERR>=1  ----
  102 IF(NERR.GE.1) THEN
         WRITE(*,1010)
         WRITE(KP,1010)
         STOP
      END IF
C  WARNING FOR PARAMETER ESTIMATION
      IF(MODE.EQ.4.OR.MODE.EQ.6.OR.MODE.EQ.8) THEN
        IF(MIT.GE.1) THEN
         WRITE(*,1032)
         WRITE(KP,1032)
        END IF
      END IF
C ----- CHECK INDEX(I) FOR PARAMETER ESTIMATION --------------------
C  V CANNOT BE ESTIMATED WHEN DIMENSIONLESS TIME (NREDU=2,3)
C   V HAS TO BE SPECIFIED FOR DIMENSIONLESS TIME
  110 IF(MIT.GE.1.AND.NREDU.GE.2) THEN
        IF(INDEX(1).EQ.1) THEN
           WRITE(*,1012)
           WRITE(KP,1012)
           INDEX(1)=0
        END IF
      END IF
C
C    ----- DISPERSION COEFFICIENT FOR STOCHASTIC MODEL FOR MODE 3,4,8--------
C    ASSUME POSITIVE CORRELATION V&D FOR MODEL 3,4,8
C MODD; Index for a stochastic dispersion coefficient
C    MODD=2  SDLND=SDLNV (FIELD SCALE DISPERSIVITY=CONSTANT)
C    MODD=1 POSITIVE CORRELATION V&D.
C    MODD=0  D=CONSTANT
C    MODD=-1 NEGATIVE CORRELATION V&D.
C
      MODD=0
      MODK=0
      MD56=0
      MK34=0
      MAL8=0
      IF(MODE.LE.2) GOTO 70
      IF (MODE.EQ.3) THEN
        IF(B(NVAR+7).LT.1.D-07) THEN
          MODD=0
        ELSE IF(ABS(B(NVAR+5)-B(NVAR+7)).LT.1.D-05) THEN
          MODD=2
        ELSE
          MODD=1
        END IF
      END IF
C
      IF (MODE.EQ.4.OR.MODE.EQ.8) THEN
        IF(B(NVAR+9).LT.1.D-07) THEN
          MODD=0
        ELSE IF(DABS(B(NVAR+7)-B(NVAR+9)).LT.1.D-05) THEN
          MODD=2
        ELSE
          MODD=1
        END IF
      END IF
C
C  --- DISPERSION COEFFICIENT FOR STOCHASTIC MODEL FOR MODE 5,6--------
C MD56; Index for a stochastic dispersion coefficient for MODE=5,6
C    MD56=1  SDLND=SDLNV
C    MD56=0  SDLND AND SALNV ARE INDEPENSENT.
      IF (MODE.EQ.5) THEN
       IF(ABS(B(NVAR+5)-B(NVAR+7)).LT.1.D-05) THEN
          MD56=1
       END IF
      END IF
      IF (MODE.EQ.6) THEN
       IF(ABS(B(NVAR+7)-B(NVAR+9)).LT.1.D-05) THEN
          MD56=1
       END IF
      END IF
C
C  --- WHEN MODD=0, NO ESTIMATE FOR SDLNV (SIG.V)
      IF(MIT.GE.1) THEN
       IF (MODE.EQ.3) THEN
        IF(MODD.EQ.0) THEN
          INDEX(7)=0
        ELSE IF((INDEX(5)-INDEX(7)).EQ.0.AND.MODD.EQ.2) THEN
          INDEX(7)=0
        END IF
       END IF
       IF (MODE.EQ.5) THEN
        IF((INDEX(5)-INDEX(7)).EQ.0.AND.MD56.EQ.1) THEN
          INDEX(7)=0
        ELSE IF(INDEX(5)*INDEX(7).EQ.0) THEN
          MD56=0
        END IF
       END IF
C
       IF (MODE.EQ.4.OR.MODE.EQ.8) THEN
        IF(MODD.EQ.0) THEN
          INDEX(9)=0
        ELSE IF((INDEX(7)-INDEX(9)).EQ.0.AND.MODD.EQ.2) THEN
          INDEX(9)=0
        END IF
       END IF
C
       IF (MODE.EQ.6) THEN
         IF((INDEX(7)-INDEX(9)).EQ.0.AND.MD56.EQ.1) THEN
           INDEX(7)=0
        ELSE IF(INDEX(7)*INDEX(9).EQ.0) THEN
          MD56=0
         END IF
       END IF
      END IF
C    ----- DISTRIBUTION COEF FOR STOCHASTIC MODEL FOR MODE 5,6,8--------
C    ASSUME NEGATIVE CORRELATION V&Kd FOR MODE 5,6,8
C MODK; Index for a stochastic distribution coefficient
C    MODK=1 POSITIVE CORRELATION V&Kd.
C    MODK=0  R=CONSTANT
C    MODK=-1 NEGATIVE CORRELATION V&Kd.
C    MODK=-2  SDLNK=SDLNV
C
      IF (MODE.EQ.5) THEN
        IF(B(NVAR+6).LT.1.D-07) THEN
          MODK=0
        ELSE IF(ABS(B(NVAR+5)-B(NVAR+6)).LT.1.D-05) THEN
          MODK=-2
        ELSE
          MODK=-1
        END IF
      END IF
C
      IF (MODE.EQ.6.OR.MODE.EQ.8) THEN
        IF(B(NVAR+8).LT.1.D-07) THEN
          MODK=0
        ELSE IF(DABS(B(NVAR+7)-B(NVAR+8)).LT.1.D-05) THEN
          MODK=-2
        ELSE
          MODK=-1
        END IF
      END IF
C
C
C  --- DISTRIBUTION COEFFICIENT FOR MODE 3,4--------
C MK34; Index for a stochastic distribution coefficient for MODE=3,4
C    MK34=1  SDLNK=SDLNV
C    MK34=0  SDLNK AND SALNV ARE INDEPENSENT.
      IF (MODE.EQ.3) THEN
       IF(ABS(B(NVAR+5)-B(NVAR+6)).LT.1.D-05) THEN
          MK34=1
       END IF
      END IF
      IF (MODE.EQ.4) THEN
       IF(ABS(B(NVAR+7)-B(NVAR+8)).LT.1.D-05) THEN
          MK34=1
       END IF
      END IF
CC
      IF(MIT.GE.1) THEN
      IF (MODE.EQ.5) THEN
        IF(MODK.EQ.0) THEN
          INDEX(6)=0
        ELSE IF((INDEX(5)-INDEX(6)).EQ.0.AND.MODK.EQ.-2) THEN
          INDEX(6)=0
        END IF
       END IF
C
      IF (MODE.EQ.6.OR.MODE.EQ.8) THEN
        IF(MODK.EQ.0) THEN
          INDEX(8)=0
        ELSE IF((INDEX(7)-INDEX(8)).EQ.0.AND.MODK.EQ.-2) THEN
          INDEX(8)=0
        END IF
      END IF
       IF (MODE.EQ.3) THEN
        IF((INDEX(5)-INDEX(6)).EQ.0.AND.MK34.EQ.1) THEN
          INDEX(6)=0
        ELSE IF(INDEX(5)*INDEX(6).EQ.0) THEN
        MK34=0
        END IF
       END IF
       IF (MODE.EQ.6) THEN
        IF((INDEX(7)-INDEX(8)).EQ.0.AND.MK34.EQ.1) THEN
          INDEX(8)=0
        ELSE IF(INDEX(7)*INDEX(8).EQ.0) THEN
        MK34=0
        END IF
       END IF
      END IF
C
CC  --- ALPHA FOR MODE 8--------
C MK34; Index for a stochastic ALPHA for MODE=8
C    MAL8=1  SDLNAL=SDLNV
C    MAL8=0  SDLNAL AND SALNV ARE INDEPENSENT.
      IF(MODE.EQ.8) THEN
       IF(ABS(B(NVAR+7)-B(NVAR+10)).LT.1.D-05) THEN
          MAL8=1
       END IF
       IF(MIT.GE.1) THEN
        IF((INDEX(7)-INDEX(10)).EQ.0.AND.MAL8.EQ.1) THEN
          INDEX(10)=0
        ELSE IF(INDEX(7)*INDEX(10).EQ.0) THEN
          MAL8=0
        END IF
       END IF
      END IF
C----  NO ESTIMATE FOR CORREATION WHEN EITHER OF PARAMETERS
C                   ARE DETERMINISTIC
      IF(MIT.GE.1) THEN
       IF(MOD(MODE,2).EQ.1) THEN
          IF(B(NVAR+5).LT.1.D-07.AND.INDEX(5).EQ.0) THEN
           IF(INDEX(8).EQ.1) THEN
            WRITE(*,1050)
            WRITE(KP,1050)
            INDEX(8)=0
           END IF
          END IF
       ELSE
          IF(B(NVAR+7).LT.1.D-07.AND.INDEX(7).EQ.0) THEN
            IF(MODE.EQ.8) THEN
               IF(INDEX(11).EQ.1) THEN
                 WRITE(*,1050)
                 WRITE(KP,1050)
                 INDEX(11)=0
               END IF
            ELSE
               IF(INDEX(10).EQ.1) THEN
                 WRITE(*,1050)
                 WRITE(KP,1050)
                 INDEX(10)=0
               END IF
            END IF
          END IF
       END IF
C
       IF(MODE.EQ.3) THEN
          IF(B(NVAR+6).LT.1.D-07.AND.INDEX(6).EQ.0) THEN
           IF(INDEX(8).EQ.1) THEN
            WRITE(*,1052)
            WRITE(KP,1052)
            INDEX(8)=0
           END IF
          END IF
       ELSE IF(MODE.EQ.4) THEN
          IF(B(NVAR+8).LT.1.D-07.AND.INDEX(8).EQ.0) THEN
           IF(INDEX(10).EQ.1) THEN
            WRITE(*,1052)
            WRITE(KP,1052)
            INDEX(10)=0
           END IF
          END IF
       ELSE IF(MODE.EQ.5) THEN
          IF(B(NVAR+7).LT.1.D-07.AND.INDEX(7).EQ.0) THEN
           IF(INDEX(8).EQ.1) THEN
            WRITE(*,1054)
            WRITE(KP,1054)
            INDEX(8)=0
           END IF
          END IF
       ELSE IF(MODE.EQ.6) THEN
          IF(B(NVAR+9).LT.1.D-07.AND.INDEX(9).EQ.0) THEN
           IF(INDEX(10).EQ.1) THEN
            WRITE(*,1054)
            WRITE(KP,1054)
            INDEX(10)=0
           END IF
          END IF
       ELSE IF(MODE.EQ.8) THEN
          IF(B(NVAR+10).LT.1.D-07.AND.INDEX(10).EQ.0) THEN
           IF(INDEX(11).EQ.1) THEN
            WRITE(*,1056)
            WRITE(KP,1056)
            INDEX(11)=0
           END IF
          END IF
       END IF
      END IF
C
C ------------- ONE-SITE MODEL (Beta=1/R) -------------
   70 IF(MIT.GE.1) THEN
      IF(MODE.EQ.2.AND.MNEQ.EQ.1) THEN
         IF(B(NVAR+3).LT.1.0) THEN
           WRITE(*,1016) BI(3),B(NVAR+3)
           WRITE(KP,1016) BI(3),B(NVAR+3)
           STOP
         END IF
         DUMMY=ABS(B(NVAR+4)-1.0/B(NVAR+3))
         IF(DUMMY.GE.0.001) THEN
            WRITE(*,1015)
            WRITE(KP,1015)
            B(NVAR+4)=1.0/B(NVAR+3)
         END IF
         IF(INDEX(4).EQ.1) THEN
            WRITE(*,1014)
            WRITE(KP,1014)
            INDEX(4)=0
         END IF
       END IF
C       ---- DEGRADATION COEFFICIENTS (mu1 & mu2)  -----------------
        IF(MODE.EQ.2.AND.(INDEX(6)+INDEX(7)).EQ.0) GO TO 99
          M468=0
        IF(MODE.EQ.4.OR.MODE.EQ.6.OR.MODE.EQ.8) M468=1
        IF(M468.EQ.1.AND.(INDEX(5)+INDEX(6)).EQ.0) GO TO 99
        IF(MDEG.EQ.1) THEN
           IF(MODE.EQ.2) THEN
               INDEX(6)=1
               INDEX(7)=0
               B(NVAR+7)=(1.0-B(NVAR+4))/B(NVAR+4)*B(NVAR+6)
           END IF
C          < MODE=4,6,8, ONE-SITE MODEL >
           IF(M468.EQ.1) THEN
               INDEX(5)=1
               INDEX(6)=0
               B(NVAR+6)=(1.0-B(NVAR+3))*B(NVAR+5)
           END IF
        ELSE IF(MDEG.EQ.2) THEN
           IF(MODE.EQ.2) THEN
               INDEX(6)=1
               INDEX(7)=0
               IF(MNEQ.NE.0.AND.MNEQ.NE.3) B(NVAR+7)=0.0
               IF(MNEQ.EQ.0.OR.MNEQ.EQ.3) B(NVAR+7)=PHIIM/PHIM*B(NVAR+6)
           ELSE IF(M468.EQ.1) THEN
               INDEX(5)=1
               INDEX(6)=0
               B(NVAR+6)=0.0
           END IF
        ELSE IF(MDEG.EQ.3) THEN
           IF(MNEQ.EQ.1) THEN
             IF(MODE.EQ.2) THEN
                 INDEX(6)=0
                 INDEX(7)=1
                 B(NVAR+6)=0.0
             ELSE IF(M468.EQ.1) THEN
                 INDEX(5)=0
                 INDEX(6)=1
                 B(NVAR+5)=0.0
             END IF
           ELSE IF(MNEQ.EQ.2) THEN
             IF(MODE.EQ.2) THEN
               INDEX(6)=1
               INDEX(7)=0
               B(NVAR+7)=(1.0-B(NVAR+4))*B(NVAR+3)/(B(NVAR+4)*B(NVAR+3)
     &                    -1.0)*B(NVAR+6)
              END IF
           ELSE IF(MNEQ.EQ.0.OR.MNEQ.EQ.3) THEN
             IF(MODE.EQ.2) THEN
               INDEX(6)=1
               INDEX(7)=0
               B(NVAR+7)=((1.0-B(NVAR+4))*B(NVAR+3)-PHIIM)/PHIM
     &                *B(NVAR+6)
             END IF
           END IF
        END IF
C
C
      END IF
C  IF # OF EDIMATED PARAMETERS IS ZERO, MIT=0
      IF(INVERSE.EQ.1) THEN
        III=0
        DO 95 I=1,NVAR
          III=III+INDEX(I)
   95   CONTINUE
        IF(III.EQ.0) THEN
          WRITE(*,1058)
          WRITE(KP,1058)
          MIT=0
        END IF
      END IF
C
C  WARNING FOR PARAMETER ESTIMATION
   99    NP=0
      IF(MIT.GE.1) THEN
        DO 58 I=1,NVAR
          IF(INDEX(I).EQ.0) GO TO 58
               NP=NP+1
   58   CONTINUE
      IF(NP.GE.5) THEN
          WRITE(*,1030) NP
         WRITE(KP,1030) NP
      END IF
      END IF
 1000 FORMAT(5X,A6,4(1H.),F10.4,' <= 0, USE POSITIVE VALUE !')
 1001 FORMAT(5X,'RHO/THETA',2(1H.),F10.4,' < 0, DO NOT USE NEGATIVE VALU
     &E !')
 1002 FORMAT(5X,A6,4(1H.),F10.4,' < 0, DO NOT USE NEGATIVE VALUE !')
 1004 FORMAT(5X,A6,4(1H.),F10.4,'  Beta SHOULD BE  0 < Beta < 1')
 1006 FORMAT(5X,A6,4(1H.),F10.4,'  Rho SHOULD BE  -1 =< Rho =< 1')
 1007 FORMAT(5X,'  USE DETERMINISTIC MODEL (SD.v=SD.Kd=0).')
 1008 FORMAT(5X,'  USE DETERMINISTIC MODEL (SD.v=SD.D=0).')
 1009 FORMAT(5X,'  USE DETERMINISTIC MODEL (SD.v=SD.alp=0).')
 1010 FORMAT(//5X,' WARNING ! PLEASE CHECK TRANSPORT PARAMETERS. ')
 1012 FORMAT(//5X,' WARNING ! V MUST BE SPECIFIED FOR DIMENSION LESS TI
     &ME. ')
 1014 FORMAT(//5X,' WARNING ! ONE-SITE MODEL (Beta=1/R)',/,6X,'Beta CANN
     &OT BE ESTIMATED FOR ONE-SITE MODEL (MNEQ=1).')
 1015 FORMAT(//5X,' WARNING ! ONE-SITE MODEL (Beta=1/R) !')
 1016 FORMAT(5X,A6,4(1H.),F10.4,'  R SHOULD BE LARGER THAN 1')
 1020 FORMAT(5X,' WARNING ! IF SD.v=0, SD.D SHOULD ALSO BE ZERO.')
 1022 FORMAT(5X,' WARNING ! IF SD.v=0, SD.Kd SHOULD ALSO BE ZERO.')
 1024 FORMAT(5X,' WARNING ! ARBITRARY INPUT FUNCTION CANNOT BE USED FOR
     &THE STOCHASTIC ONE-SITE MODEL (MODE=4,6,8).')
 1025 FORMAT(5X,A6,4(1H.),' MAX CONSTRAINT >= MIN CONSTRAINT.')
 1026 FORMAT(5X,A6,4(1H.),' WARNING ! MIN CONSTRAINT SHOULD BE LARGER TH
     &AN ZERO.')
 1027 FORMAT(5X,' WARNING ! MAX CONSTRAINT FOR BETA SHOULD BE LESS THAN
     &ONE ! ')
 1028 FORMAT(5X,' WARNING ! MIN CONSTRAINTS ARE OUT OF RANGE.'/
     &10X, '(PHIM/R < BETA < (PHIM+R-1)/R FOR TWO REGION MODEL)')
 1029 FORMAT(5X,' WARNING ! MIN CONSTRAINT FOR BETA SHOULD BE LARGER THA
     &N 1/R')
 1030 FORMAT(5X,' WARNING ! YOU ESTIMATE',I5,'  PARAMETERS (GENERALLY TO
     &O MANY).')
 1031 FORMAT(5X,' WARNING ! MAX CONSTRAINTS ARE OUT OF RANGE.'/
     &10X, '(PHIM/R < BETA < (PHIM+R-1)/R FOR TWO REGION MODEL)')
 1032 FORMAT(5X,' WARNING ! THE STOCHASTIC ONE-SITE MODELS (MODE=4,6,8)'
     &,/,10X,'      CAN NOT BE USED FOR PARAMETER ESTIMATIONS.')
 1036 FORMAT(//5X,' WARNING ! IF <Kd>= 0, SD.Kd should be 0.  ')
 1037 FORMAT(5X,A6,4(1H.),' WARNING ! MIN CONSTRAINT SHOULD BE LARGER TH
     &AN -1.')
 1038 FORMAT(5X,A6,4(1H.),' WARNING ! MAN CONSTRAINT SHOULD BE SMALLER T
     &HAN 1.')
 1040 FORMAT(5X,' WARNING ! R SHOULD BE FIXED FOR MNEQ = 3')
 1042 FORMAT(5X,' WARNING ! OMEGA IS TOO LARGE ( > 100) FOR THE NONEQUIL
     &IBRIUM')
 1044 FORMAT(5X,' WARNING ! ALPHA IS TOO LARGE (OMEGA > 100) FOR THE NON
     &EQUILIBRIUM')
 1046 FORMAT(5X,' WARNING ! THE CONSTRAINT FOR OMEGA IS TOO LARGE (> 100
     &)')
 1048 FORMAT(5X,' WARNING ! THE CONSTRAINT FOR ALPHA IS TOO LARGE (> 100
     &)')
 1050 FORMAT(5X,'  NO ESTIMATION FOR RHO BECAUSE V IS DETERMINISTIC !')
 1052 FORMAT(5X,'  NO ESTIMATION FOR RHO BECAUSE Kd IS DETERMINISTIC !')
 1054 FORMAT(5X,'  NO ESTIMATION FOR RHO BECAUSE D IS DETERMINISTIC !')
 1056 FORMAT(5X,'  NO ESTIMATION FOR RHO BECAUSE ALPHA IS DETERMINISTIC
     &!')
 1058 FORMAT(5X,'  SET MIT = 0 BECAUSE OF NO FITTED PARAMETER !')
      RETURN
      END
C
C    -----------------------------------------------------------
      SUBROUTINE DATAOUT(TITLE1,TITLE2,ILMT,B,BMAX,BMIN,BI,KP,FILEIN)
C
C       PURPOSE: WRITE INPUT DATA TO FILE.
C
      IMPLICIT REAL*8(A-H,O-Z)
C
C
      CHARACTER TITLE1*60,TITLE2*60,BI(15)*6,FILEIN*15
      DIMENSION B(30),BMAX(15),BMIN(15)
      COMMON MODC,MM,A,BB,DA,CX,P,BETA,BETR,OMEGA,R,DMU1,DMU2,MODE,
     &   MCON,TT,ZZ,STOPER,LEVEL,ICHEB
      COMMON/MODAT/INDEX(15),INVERSE,NREDU,NVAR,ZL,MIT,MDEG,MASS,
     &     DUMTP(50),DUMGA1(50),MNEQ,ISKIP,PHIM,PHIIM
      COMMON/STOCH/MODD,MODK,MODS,MSTOCH,CORR,MCORR,MSD,SDLNK,
     & V,SDLNV,VMAX,VMIN,D,AVEY,SDLNY,YMAX,YMIN,ALPHA,SDLND
     & ,DK,RHOTH,MD56,MK34,MAL8
      COMMON/BOUN/MODB,NPULSE,TPULSE(50),PULSE(50),MASSST
      COMMON/INITI/MODI,NINI,CINI(50),ZINI(50)
      COMMON/PROD/MODP,NPRO1,NPRO2,GAMMA1(50),ZPRO1(50),GAMMA2(50),
     &        ZPRO2(50)
C
C ------ WRITE TITLE -----------------
      WRITE(*,800)
      IF(INVERSE.EQ.1.AND.MIT.GT.0) THEN
         WRITE(*,802)
      ELSE
         WRITE(*,804)
      END IF
      WRITE(*,806)
      WRITE(*,808) TITLE1
      WRITE(*,808) TITLE2
      WRITE(*,806)
      WRITE(*,809) FILEIN
      WRITE(*,810)
C
      WRITE(KP,800)
      IF(INVERSE.EQ.1.AND.MIT.GT.0) THEN
         WRITE(KP,802)
      ELSE
         WRITE(KP,804)
      END IF
      WRITE(KP,806)
      WRITE(KP,808) TITLE1
      WRITE(KP,808) TITLE2
      WRITE(KP,806)
      WRITE(KP,809) FILEIN
      WRITE(KP,810)
C    -------- WRITE MODEL DESCRIPTION -----------
      WRITE(KP,850)
      IF(MODE.GE.9) GOTO 17
      GO TO (11,12,13,14,15,16,17,18) MODE
   11   WRITE(KP,852)
        GO TO 20
   12   WRITE(KP,854)
        GO TO 20
   13   WRITE(KP,856)
        GO TO 20
   14   WRITE(KP,858)
        GO TO 20
   15   WRITE(KP,860)
        GO TO 20
   16   WRITE(KP,862)
        GO TO 20
   17   WRITE(KP,864)
        GO TO 20
   18   WRITE(KP,866)
C
C------  DISPERSION COEFF. FOR MODE=3,4,8
   20 IF(MODE.EQ.3.OR.MODE.EQ.4.OR.MODE.EQ.8) THEN
        IF(MODD.EQ.2) THEN
          WRITE(KP,1150)
        ELSE IF (MODD.EQ.0) THEN
          WRITE(KP,1152)
        ELSE IF (MODD.EQ.1) THEN
          WRITE(KP,1154)
        END IF
      END IF
C------  DISPERSION COEFF. FOR MODE=5,6
      IF(MODE.EQ.5.OR.MODE.EQ.6) THEN
        IF(MD56.EQ.1) WRITE(KP,1151)
      END IF
C
C------  RETARDATION FACTOR  ----------
      IF(MODE.GE.3) THEN
        WRITE(KP,1170) R,RHOTH
      END IF
C------  DISTRIBUTION COEFF. FOR MODE=5,6,8
      IF(MODE.EQ.5.OR.MODE.EQ.6.OR.MODE.EQ.8) THEN
        IF(MODK.EQ.-2) THEN
          WRITE(KP,1160)
        ELSE IF (MODK.EQ.0) THEN
          WRITE(KP,1162)
        ELSE IF (MODK.EQ.-1) THEN
          WRITE(KP,1164)
        END IF
      END IF

C------  DISTRIBUTION COEFF. FOR MODE=3,4
      IF(MODE.EQ.3.OR.MODE.EQ.4) THEN
        IF(MK34.EQ.1) WRITE(KP,1161)
      END IF
C------  ALPHA FOR MODE=8
      IF(MODE.EQ.8.AND.MAL8.EQ.1) WRITE(KP,1163)
C -----------------------------------------
      IF(MODC.LE.0.OR.MODC.GE.7) THEN
          WRITE(KP,880)
           GO TO 30
        END IF
        GO TO (21,22,23,24,25,26) MODC
   21   IF(MODE.LE.2) THEN
        WRITE(KP,882)
        ELSE
        WRITE(KP,884)
        END IF
        GO TO 30
   22   IF(MODE.LE.2) THEN
        WRITE(KP,882)
        ELSE
        WRITE(KP,886)
        END IF
        GO TO 30
   23   IF(MODE.LE.2) THEN
        WRITE(KP,888)
        ELSE
        WRITE(KP,890)
        END IF
        GO TO 30
   24   IF(MODE.LE.2) THEN
        WRITE(KP,892)
        ELSE
        WRITE(KP,894)
        END IF
        GO TO 30
   25   IF(MODE.LE.2) THEN
        WRITE(KP,896)
        ELSE
        WRITE(KP,898)
        END IF
        GO TO 30
   26   IF(MODE.LE.2) THEN
        WRITE(KP,900)
        ELSE
        WRITE(KP,902)
        END IF
C
   30   IF(NREDU.EQ.2) GOTO 32
        IF(NREDU.EQ.3) GOTO 33
C
        IF(MOD(MODE,2).EQ.1) THEN
        WRITE(KP,910)
        ELSE
        WRITE(KP,912)
        WRITE(KP,1005) ZL
         IF(MOD(MODE,2).EQ.0.AND.MODE.GE.4) THEN
         IF(MODC.NE.4.AND.MODC.NE.6) THEN
          IF(NREDU.EQ.1) WRITE (KP,918)
          IF(NREDU.EQ.0) WRITE (KP,920)
         END IF
         END IF

        END IF
        GO TO 35
   32   WRITE(KP,914)
        WRITE(KP,1005) ZL
        GO TO 35
   33   WRITE(KP,916)
        WRITE(KP,1005) ZL
C
   35 IF(INVERSE.LE.0) GO TO 38
      IF(INVERSE.EQ.1.AND.MIT.EQ.0) GO TO 38
       IF(ILMT.EQ.1) THEN
          WRITE(KP,1006)
          DO 36 I=1,NVAR
          IF(INDEX(I).EQ.0) THEN
            WRITE(KP,1013) BI(I),B(I+NVAR)
          ELSE
            WRITE(KP,1014) BI(I),B(I+NVAR),BMIN(I),BMAX(I)
          END IF
   36     CONTINUE
       ELSE
          WRITE(KP,1010)
          DO 37 I=1,NVAR
          IF(INDEX(I).EQ.0) THEN
            WRITE(KP,1013) BI(I),B(I+NVAR)
          ELSE
            WRITE(KP,1014) BI(I),B(I+NVAR)
          END IF
   37     CONTINUE
       END IF
       GO TO 40
   38 WRITE(KP,1009)
      DO 39 I=1,NVAR
   39 WRITE(KP,1008) BI(I),B(I+NVAR)
C
   40 WRITE(KP,1020)
C     ------- PARAMETERS FOR BOUNDARY VALUE PROBLEM ----------
      IF(MASS.EQ.1) THEN
        WRITE(KP,1064)
      END IF
      GO TO (41,42,43,44,45,46,47) MODB+1
   41 WRITE(KP,1065)
      GO TO 50
   42 IF(MODE.LE.2) THEN
       WRITE(KP,1067) PULSE(1)
      ELSE IF (MASSST.EQ.1) THEN
       WRITE(KP,1068) PULSE(1)
      ELSE
       WRITE(KP,1069) PULSE(1)
      END IF
      GO TO 50
   43 WRITE(KP,1070) PULSE(1)
      GO TO 50
   44 IF(MODE.LE.2) THEN
         WRITE(KP,1071) PULSE(1),TPULSE(2)
      ELSE IF (MASSST.EQ.0) THEN
C  MASSST=1 AMOUNT OF SOLUTE IS CONSTANT REGARDELESS IF V (pulse*<v>/v)
C        =0                  DEPENDS V
         WRITE(KP,1061) PULSE(1),TPULSE(2)
      ELSE
         WRITE(KP,1062) PULSE(1),TPULSE(2)
      END IF
      GO TO 50
   45 WRITE(KP,1072) NPULSE
      WRITE(KP,1073)
          DO 49 I=1, NPULSE
            WRITE(KP,1012) PULSE(I),TPULSE(I)
   49 CONTINUE
      GO TO 50
   46 WRITE(KP,1074)
      WRITE(KP,1075) PULSE(1),PULSE(2),TPULSE(1)
      GO TO 50
   47 WRITE(KP,1076)
C     ------- PARAMETERS FOR INITIAL VALUE PROBLEM ----------
   50   GO TO (51,52,53,54,56) MODI+1
   51   WRITE(KP,1078)
           GO TO 60
   52   WRITE(KP,1080) CINI(1)
           GO TO 60
   53   WRITE(KP,1081) NINI
           WRITE(KP,1082)
         DO 55 I=1,NINI
          WRITE(KP,1012) CINI(I),ZINI(I)
   55     CONTINUE
          GO TO 60
C
   54   WRITE(KP,1084)
        WRITE(KP,1085) CINI(1),CINI(2),ZINI(1)
          GOTO 60
   56   WRITE(KP,1086) CINI(2),ZINI(2),CINI(1)
C
C     ------- PARAMETERS FOR PRODUCTION VALUE PROBLEM  ----------
   60   GO TO (61,62,63,64) MODP+1
   61   WRITE(KP,1087)
        GO TO 95
   62   IF(MOD(MODE,2).EQ.1) THEN
          WRITE(KP,1088) GAMMA1(1)
        ELSE
          WRITE(KP,1089) GAMMA1(1),GAMMA2(1)
        END IF
        GO TO 95
   63   WRITE(KP,1091)
        WRITE(KP,1092) NPRO1
      DO 65 I=1,NPRO1
          WRITE(KP,1012) GAMMA1(I),ZPRO1(I)
   65     CONTINUE
      IF(MOD(MODE,2).EQ.1) GOTO 95
        WRITE(KP,1093) NPRO2
      DO 66 I=1,NPRO2
          WRITE(KP,1012) GAMMA2(I),ZPRO2(I)
   66     CONTINUE
          GO TO 95
C
   64   WRITE(KP,1094)
        WRITE(KP,1095) GAMMA1(1),GAMMA1(2),ZPRO1(1)
      IF(MOD(MODE,2).EQ.1) GOTO 95
        WRITE(KP,1096) GAMMA2(1),GAMMA2(2),ZPRO2(1)
C
C  ------ PARAMETER ESTIMATION MODE  -------
   95 IF(INVERSE.LE.0) GOTO 99
        WRITE(KP,1025)
        WRITE(KP,1098) MIT
C----  MASS ESTIMATION
      IF(MIT.GE.1.AND.MASS.EQ.1) THEN
         IF(MODB.EQ.1.AND.INDEX(NVAR).EQ.1) THEN
           WRITE(KP,1120)
           IF(ILMT.EQ.1) WRITE(KP,1121) BMIN(NVAR),BMAX(NVAR)
         ELSE IF(MODB.EQ.2.AND.INDEX(NVAR).EQ.1) THEN
           WRITE(KP,1122)
           IF(ILMT.EQ.1) WRITE(KP,1123) BMIN(NVAR),BMAX(NVAR)
         ELSE IF(MODB.EQ.3) THEN
           IF(INDEX(NVAR-1).EQ.1) THEN
               WRITE(KP,1122)
               IF(ILMT.EQ.1) WRITE(KP,1123) BMIN(NVAR-1),BMAX(NVAR-1)
           END IF
           IF(INDEX(NVAR).EQ.1) THEN
               WRITE(KP,1124)
               IF(ILMT.EQ.1) WRITE(KP,1125) BMIN(NVAR),BMAX(NVAR)
           END IF
         END IF
      END IF
C----   TYPE OF THE NONEQUILIBRIUM CDE & DEGRADATION CODE  -------
C
      IF(MOD(MODE,2).EQ.1) GOTO 99
      IF(MIT.GE.1) THEN
        IF(MNEQ.EQ.0) THEN
          WRITE(KP,1100)
        ELSE IF(MNEQ.EQ.1) THEN
          WRITE(KP,1102)
        ELSE IF(MNEQ.EQ.2) THEN
          WRITE(KP,1104)
        ELSE IF(MNEQ.EQ.3) THEN
          WRITE(KP,1105)
          WRITE(KP,1114) PHIM
       END IF
C
        IF(MODE.EQ.2.AND.(INDEX(6)+INDEX(7)).EQ.0) GO TO 99
        IF(MODE.EQ.4.AND.(INDEX(6)+INDEX(7)).EQ.0) GO TO 99
        IF(MDEG.EQ.0) THEN
          WRITE(KP,1106)
        ELSE IF(MDEG.EQ.1) THEN
          WRITE(KP,1108)
        ELSE IF(MDEG.EQ.2) THEN
          WRITE(KP,1110)
        ELSE IF(MDEG.EQ.3) THEN
          WRITE(KP,1112)
        END IF
        IF(MNEQ.EQ.0.AND.MDEG.GE.2) THEN
           WRITE(KP,1114) PHIM
        END IF
C
      END IF
C
  800 FORMAT('$',//5X,67(1H*),/5X,1H*,65X,1H*,/5X,1H*,5X,
     &'CXTFIT VERSION 2.1 (4/17/99)',32X,1H*,/5X,1H*,5X,
     &'ANALYTICAL SOLUTIONS FOR ONE-DIMENSIONAL CDE',16X,1H*)
  802 FORMAT(5X,1H*,5X,'NON-LINEAR LEAST-SQUARES ANALYSIS',27X,1H*)
  804 FORMAT(5X,1H*,5X,'DIRECT PROBLRM                        ',22X,1H*)
  806 FORMAT(5X,1H*,65X,1H*)
  808 FORMAT(5X,1H*,5X,A60,1H*)
  809 FORMAT(5X,1H*,5X,'DATA INPUT FILE:  ',A15,27X,1H*)
  810 FORMAT(5X,1H*,65X,1H*/5X,67(1H*)//)
  850 FORMAT(5X,'MODEL DESCRIPTION',/5X,17(1H=))
  852 FORMAT(8X,'DETERMINISTIC EQUILIBRIUM CDE (MODE=1)')
  854 FORMAT(8X,'DETERMINISTIC NONEQUILIBRIUM CDE (MODE=2)')
  856 FORMAT(8X,'STOCHASTIC EQUILIBRIUM CDE WITH VARIABLE V & Kd ',
     &'(MODE=3)')
  858 FORMAT(8X,'STOCHASTIC ONE-SITE NONEQUILIBRIUM CDE WITH VARIABLE V
     && Kd (MODE=4)')
  860 FORMAT(8X,'STOCHASTIC EQUILIBRIUM CDE WITH VARIABLE V & D (MODE=5)
     &')
  862 FORMAT(8X,'STOCHASTIC NONEQUILIBRIUM CDE WITH VARIABLE V & D (MODE
     &=6)')
  864 FORMAT(8X,'WARNING! PLEASE CHECK MODE.')
  866 FORMAT(8X,'STOCHASTIC NONEQUILIBRIUM CDE WITH VARIABLE V & ALPHA (
     &MODE=8)')
  880 FORMAT(8X,'WARNING! PLEASE CHECK CONCENTRATION MODE (MODC).')
  882 FORMAT(8X,'FLUX-AVERAGED CONCENTRATION')
  884 FORMAT(8X,'ENSEMBLE AVERAGE OF FLUX-AVERAGED CONCENTRATION <Cf>')
  886 FORMAT(8X,'FIELD-SCALE FLUX-AVERAGED CONCENTRATION <VCf>/<V>')
  888 FORMAT(8X,'RESIDENT CONCENTRATION (THIRD-TYPE INPUT)')
  890 FORMAT(8X,'ENSEMBLE AVERAGE OF RESIDENT CONCENTRATION (THIRD-TYPE
     &INPUT)')
  892 FORMAT(8X,'TOTAL RESIDENT CONCENTRATION (THIRD-TYPE INPUT)')
  894 FORMAT(8X,'ENSEMBLE AVERAGE OF TOTAL RESIDENT CONCENTRATION (THIRD
     &-TYPE INPUT)')
  896 FORMAT(8X,'RESIDENT CONCENTRATION (FIRST-TYPE INPUT)')
  898 FORMAT(8X,'ENSEMBLE AVERAGE OF RESIDENT CONCENTRATION (FIRST-TYPE
     &INPUT)')
  900 FORMAT(8X,'TOTAL RESIDENT CONCENTRATION (FIRST-TYPE INPUT)')
  902 FORMAT(8X,'ENSEMBLE AVERAGE OF TOTAL RESIDENT CONCENTRATION (FIRST
     &-TYPE INPUT)')
  910 FORMAT(8X,'REAL TIME (t), POSITION(x)',/,10X,'(D,V,mu, AND gamma A
     &RE ALSO DIMENSIONAL)')
  912 FORMAT(8X,'REAL TIME (t), POSITION(x)',/,10X,'(D AND V ARE DIMENSI
     &ONAL;',/,13X,' R,beta,omega,mu, AND gamma ARE DIMENSIONLESS)')
  914 FORMAT(8X,'REDUCED TIME (T), POSITION(Z)',/,10X,'(ALL PARAMETERS E
     &XCEPT D AND V ARE DIMENSIONLESS)')
  916 FORMAT(8X,'REDUCED TIME (T), DIMENSIONAL POSITION(Z)',/,10X,'(ALL
     &PARAMETERS EXCEPT D AND V ARE DIMENSIONLESS)')
  918 FORMAT(9X,'REDUCED MEAN PHASE 2 CONC. <C2>=<c2/Kd>')
  920 FORMAT(9X,'DIMENSIONAL MEAN PHASE 2 CONC. <c2>')
 1005 FORMAT(8X,'CHARACTERISTIC LENGTH = ',F9.4,/,10X,
     & 'FOR DIMENSIONLESS PARAMETERS')
 1006 FORMAT(//5X,'INITIAL VALUES OF COEFFICIENTS'/5X,30(1H=)/5X,
     &'NAME',8X,'INITIAL VALUE',2X,'FITTING',3X,'MIN VALUE ','    MAX VA
     &LUE')
 1008 FORMAT(5X,A6,4(1H.),3E14.4)
 1009 FORMAT(//5X,'INITIAL VALUES OF COEFFICIENTS'/5X,30(1H=)/5X,
     1'NAME',8X,'INITIAL VALUE ')
 1010 FORMAT(//5X,'INITIAL VALUES OF COEFFICIENTS'/5X,30(1H=)/5X,
     1'NAME',8X,'INITIAL VALUE ',2X,'FITTING')
 1012 FORMAT(7X,8F10.5)
 1013 FORMAT(5X,A6,4(1H.),E14.4,6X,'N')
 1014 FORMAT(5X,A6,4(1H.),E14.4,6X,'Y',1X,2E14.4)
 1020 FORMAT(//5X,'BOUNDARY, INITIAL, AND PRODUCTION CONDITIONS'/5X,
     & 43(1H=))
 1025 FORMAT(//5X,'PARAMETER ESTIMATION MODE'/5X,25(1H=))
 1064 FORMAT(5X,'  <INITIAL ESTIMATE OF B.C.>')
 1065 FORMAT(5X,'  SOLUTE FREE INPUT')
 1067 FORMAT(5X,'  DIRAC DELTA INPUT , MASS =',F9.4)
 1068 FORMAT(5X,'  DIRAC DELTA INPUT , MASS =',F9.4,/,
     &9X,'AMOUNT OF SOLUTE IS EVERYWHERE SAME REGARDLESS OF V')
 1069 FORMAT(5X,'  DIRAC DELTA INPUT , MASS =',F9.4,/
     &9X,'AMOUNT OF SOLUTE DEPENDS ON LOCAL V  (input=V delta/<V>)')
 1070 FORMAT(5X,'  STEP INPUT OF CONC. = ',F9.4)
 1071 FORMAT(5X,'  SINGLE PULSE OF CONC. = ',F9.4,' & DURATION = ',
     &F9.4)
 1061 FORMAT(5X,'  SINGLE PULSE OF CONC. = ',F9.4,' & DURATION = ',F9.4
     &,/,9X,'DURATION TIME IS CONSTANT FOE ALL TUBES')
 1062 FORMAT(5X,'  SINGLE PULSE OF CONC. = ',F9.4,' & Mo/THETA = ',
     &F9.4,/9X,'DURATION TIME VARIES INVERSELY WITH LOCAL V')
 1072 FORMAT(5X,'  MULTIPLE PULSE INPUT   No. OF PULSES =  ',I6)
 1073 FORMAT(5X,'  CONC.      Starting time')
 1074 FORMAT(5X,'  EXPONENTIAL INPUT')
 1075 FORMAT(5X,'    Cin(T) = ',F9.4,'  +',F9.4,' EXP ( - ',F9.4,' * ',
     &'T )')
 1076 FORMAT(5X,'  INPUT FUNCTION DEFINED IN SUBROUTINE CINPUT')
 1078 FORMAT(5X,'  SOLUTE FREE INITIAL CONDITION')
 1080 FORMAT(5X,'  CONSTANT INITIAL CONC. = ',F9.4)
 1081 FORMAT(5X,I5,' STEPS INITIAL CONDITION')
 1082 FORMAT(5X,'  CONC.    STARTING POSITION')
 1084 FORMAT(5X,'  EXPONENTIAL INITIAL CONDITION')
 1085 FORMAT(5X,'    C(Z,0) = ',F9.4,'  +',F9.4,' * EXP (  - ',F9.4,
     & ' * Z )')
 1086 FORMAT(5X,'  DELTA WITH CONSTANT INITIAL CONDITION',/,
     &5X,'    C(Z,0) = ',F9.4,' * DELTA( Z - ',F9.4,' ) +',F9.4)
 1087 FORMAT(5X,'  NO PRODUCTION TERM ')
 1088 FORMAT(5X,'  CONSTANT PRODUCTION   GAMMA = ',F9.4)
 1089 FORMAT(5X,'  CONSTANT PRODUCTION   GAMMA1 = ',F9.4,' GAMMA2 = '
     & ,F9.4)
 1091 FORMAT(5X,'  STEPWISE PRODUCTION')
 1092 FORMAT(8X,'GAMMA1 ',I5,' STEPS ')
 1093 FORMAT(/8X,'GAMMA2 ',I5,' STEPS ')
 1094 FORMAT(5X,'  EXPONENTIAL PRODUCTION')
 1095 FORMAT(5X,'  GAMMA1 = ',F9.4,' + ',F9.4,' * EXP ( - '
     & ,F9.4,' * Z )')
 1096 FORMAT(5X,'  GAMMA2 = ',F9.4,' + ',F9.4,' * EXP ( - '
     & ,F9.4,' * Z )')
 1098 FORMAT(7X,'MAXIMUM NUMBER OF ITERATIONS =',I5)
 1100 FORMAT(7X,'TWO-REGION PHYSICAL NONEQUILIBRIUM MODEL')
 1102 FORMAT(7X,'ONE-SITE CHEMICAL NONEQUILIBRIUM MODEL')
 1104 FORMAT(7X,'TWO-SITE CHEMICAL NONEQUILIBRIUM MODEL')
 1105 FORMAT(7X,'TWO-REGION PHYSICAL NONEQUILIBRIUM MODEL',
     & /10X,'INTERNAL CONSTRAINTS FOR BETA BASED ON R & PHIM')
 1106 FORMAT(7X,'SOLUTION AND SORBED PHASE DEGRADATION RATES ARE INDEPE
     &NDENT')
 1108 FORMAT(7X,'DEGRADATION EVERYWHERE THE SAME')
 1110 FORMAT(7X,'DEGRADATION ONLY IN THE LIQUID PHASE')
 1112 FORMAT(7X,'DEGRADATION ONLY IN THE SORBED PHASE')
 1114 FORMAT(9X,'PHIM =',F9.4)
 1120 FORMAT(7X,'MASS IS FITTED TO THE DATA')
 1121 FORMAT(7X,F9.4,' < MASS < ',F9.4)
 1122 FORMAT(7X,'INPUT CONCENTRATION, Cin, IS FITTED TO THE DATA')
 1123 FORMAT(7X,F9.4,' < Cin < ',F9.4)
 1124 FORMAT(7X,'DURATION TIME, T2, IS FITTED TO THE DATA')
 1125 FORMAT(7X,F9.4,' < T2 < ',F9.4)
 1150 FORMAT(9X,'SD.v=SD.D, FIELD SCALE DISPERSIVITY IS CONSTANT ( =<D>/
     &<V> )')
 1151 FORMAT(9X,'SD.v=SD.D, FOR THE JOINT PDF, f(V,D) ')
 1152 FORMAT(9X,'SD.D=0, LOCAL SCALE DISPERSION IS CONSTANT EVERYWHERE')
 1154 FORMAT(9X,'POSITIVE CORRELATION BETWEEN V & D')
 1160 FORMAT(9X,'SD.v=SD.Kd, NEGATIVE CORRELATION BETWEEN V & Kd')
 1161 FORMAT(9X,'SD.v=SD.Kd, FOR THE JOINT PDF, f(V,Kd)')
 1162 FORMAT(9X,'SD.Kd=0, LOCAL SCALE ADSORPTION IS CONSTANT EVERYWHERE'
     &)
 1163 FORMAT(9X,'SD.v=SD.alp, FOR THE JOINT PDF, f(V,alpha) ')
 1164 FORMAT(9X,'NEGATIVE CORRELATION BETWEEN V & Kd')
 1170 FORMAT(9X,' <R> = 1+RHO*<Kd>/THETA = ',F9.4,',  RHO/THETA = ',
     & F9.4)
C
   99 RETURN
      END
