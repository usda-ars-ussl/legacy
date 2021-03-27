C    -----------------------------------------------------------
      SUBROUTINE MODEL(BN,T,Z,CXT1,CXT2,VAR1,VAR2)
C
C     PURPOSE: ASSIGN COEFFICIENTS AND ROUTE EXECUTION FOR MODELS
C
      IMPLICIT REAL*8(A-H,O-Z)
C
C
      DIMENSION BN(30)
      COMMON MODC,MM,A,B,DA,CX,P,BETA,BETR,OMEGA,R,DMU1,DMU2,MODE,
     &       MCON,TT,ZZ,STOPER,LEVEL,ICHEB
      COMMON/MODAT/INDEX(15),INVERSE,NREDU,NVAR,ZL,MIT,MDEG,MASS,
     &      DUMTP(50),DUMGA1(50),MNEQ,ISKIP,PHIM,PHIIM
      COMMON/STOCH/MODD,MODK,MODS,MSTOCH,CORR,MCORR,MSD,SDLNK,
     & V,SDLNV,VMAX,VMIN,D,AVEY,SDLNY,YMAX,YMIN,ALPHA,SDLND
     & ,DK,RHOTH,MD56,MK34,MAL8
      COMMON/BOUN/MODB,NPULSE,TPULSE(50),PULSE(50),MASSST
      COMMON/INITI/MODI,NINI,CINI(50),ZINI(50)
      COMMON/PROD/MODP,NPRO1,NPRO2,GAMMA1(50),ZPRO1(50),GAMMA2(50),
     &        ZPRO2(50)
C
      VAR1=0.0
      VAR2=0.0
      CXT2=0.0
      TT=T
      ZZ=Z
C
C     ----UPDATE COEFFICIENT ARRAY----
      IF(INVERSE.EQ.1) THEN
        K=0
        NU1=NVAR+1
        NU2=NVAR*2
       DO 2 I=NU1,NU2
        IF(INDEX(I-NVAR).EQ.0) GO TO 2
        K=K+1
        BN(I)=BN(K)
    2  CONTINUE
      END IF
C
      IF(MOD(MODE,2).EQ.1) THEN
        BETA=1.D00
        OMEGA=0.D00
        DMU2=0.D00
      END IF
C--------CHANGE FROM ACTUAL TO REDUCED VARIABLE
       V=BN(NVAR+1)
C
C   NREDU=0; REAL t & x  C2=<C2>
C   NREDU=1; REAL t & x  C2=<C2/k> FOR THE STOCAHSTIC NONEQUILIBRIUM MODEL
C   NREDU=2; DIMENSIONLESS T & Z
C   NREDU=3; DIMENSIONLESS T & REAL x
C     (MODE=1,3,5) & NREDU=1; mu & gamma ARE DIMENSIONAL (SEE (2.5),(2.6)).
C         OTHERWISE, mu & gamma ARE DIMENSIONLESS AS DESCRIBED IN TABLE 2.1.
C
       IF(MODB.EQ.3.AND.MASS.EQ.1) DUMTP(2)=BN(NVAR*2)
c
      IF(NREDU.LE.1.OR.NREDU.EQ.3) ZZ=ZZ/ZL
       IF(NREDU.LE.1) THEN
          TT=TT*V/ZL
          IF(MODB.EQ.3.OR.MODB.EQ.4) THEN
            DO 21 I=1,NPULSE
              TPULSE(I)=DUMTP(I)*V/ZL
   21        CONTINUE
          END IF
          IF(MODB.EQ.5)  TPULSE(1)=DUMTP(1)/V*ZL

C     ---- CHANGE gamma TO DIMENSIONLESS VARIABLES FOR MODE=1,3,5
          IF(MOD(MODE,2).EQ.1) THEN
            IF(MODP.EQ.1.OR.MODP.EQ.2) THEN
             DO 23 I=1,NPRO1
              GAMMA1(I)=DUMGA1(I)*ZL/V
   23        CONTINUE
            END IF
          IF(MODP.EQ.4) THEN
              GAMMA1(1)=DUMGA1(1)*ZL/V
              GAMMA1(2)=DUMGA1(2)*ZL/V
          END IF
        END IF
       END IF
C --------------------------------------------------------
C  PARAMETER FOR TOTAL MASS
      IF(MASS.EQ.1) THEN
          IF(MODB.EQ.1.OR.MODB.EQ.2) PULSE(1)=BN(NVAR*2)
          IF(MODB.EQ.3) PULSE(1)=BN(NVAR*2-1)
      END IF
C
       GO TO (31,32,33,34,35,36,37,38) MODE
C
C -------  DETERMINISTIC CDE (MODE 1 AMD 2) ------------------
C
C----MODE 1 ----------
   31    D=BN(NVAR+2)
         P=V*ZL/D
         R=BN(NVAR+3)
         BETR=R
         DMU1=BN(NVAR+4)
          IF(NREDU.LE.1) THEN
            DMU1=DMU1*ZL/V
          END IF
C
      CALL   DETCDE(CXT1,CXT2)
      IF(INVERSE.EQ.1.AND.(MODC.EQ.4.OR.MODC.EQ.6)) THEN
         CXT1=R*CXT1
         CXT2=0.0
      END IF
      IF(DABS(CXT1).LT.1.D-30) CXT1=0.0D+0
      IF(DABS(CXT2).LT.1.D-30) CXT2=0.0D+0
C
      RETURN
C----MODE 2 --------------
   32    D=BN(NVAR+2)
         P=V*ZL/D
         R=BN(NVAR+3)
         BETA=BN(NVAR+4)
         OMEGA=BN(NVAR+5)
         DMU1=BN(NVAR+6)
         DMU2=BN(NVAR+7)
C      < ONE-SITE MODEL>
         IF(MNEQ.EQ.1.AND.MIT.GE.1) THEN
         BETA=1/R
         BN(NVAR+4)=BETA
         END IF
C      < DEGRADATION COEFFICIENTS >
        IF(MIT.GE.1.AND.MDEG.GE.1) THEN
          IF(MDEG.EQ.1) THEN
            DMU2=(1.0-BETA)/BETA*DMU1
            BN(NVAR+7)=DMU2
          ELSE IF(MDEG.EQ.2.AND.(MNEQ.EQ.0.OR.MNEQ.EQ.3)) THEN
            DMU2=PHIIM/PHIM*DMU1
            BN(NVAR+7)=DMU2
          ELSE IF(MDEG.EQ.3.AND.(MNEQ.EQ.0.OR.MNEQ.EQ.3)) THEN
            DMU2=DMU1/(BETA*R-PHIM)*((1.0-BETA)*R-PHIIM)
            BN(NVAR+7)=DMU2
          ELSE IF(MDEG.EQ.3.AND.MNEQ.EQ.2) THEN
            DMU2=DMU1/(BETA*R-1.0)*(1.0-BETA)*R
            BN(NVAR+7)=DMU2
          ELSE IF(MDEG.EQ.3.AND.MNEQ.EQ.1) THEN
            DMU1=0.0
            BN(NVAR+6)=DMU1
          ELSE IF(MDEG.EQ.2.AND.(MNEQ.EQ.2.OR.MNEQ.EQ.1)) THEN
            DMU2=0.0
            BN(NVAR+7)=DMU2
        END IF
      END IF
C
        IF(INVERSE.EQ.1.AND.(MODC.NE.4.AND.MODC.NE.6)) MCON=1
C
      CALL   DETCDE(CXT1,CXT2)
        IF(INVERSE.EQ.1.AND.(MODC.EQ.4.OR.MODC.EQ.6)) THEN
         CXT1=BETA*R*CXT1
         CXT2=(1-BETA)*R*CXT2
        END IF
C
      IF(DABS(CXT1).LT.1.D-30) CXT1=0.0D+0
      IF(DABS(CXT2).LT.1.D-30) CXT2=0.0D+0
      RETURN
C
C ------------- STOCHASTIC CDE  (MODE >= 3)  -----------------
C
C--- MODE 3 ---------
   33    D=BN(NVAR+2)
         P=V*ZL/D
         DK=BN(NVAR+3)
         R=1.0+RHOTH*DK
         BETR=R
         Y=DK
         DMU1=BN(NVAR+4)
           IF(NREDU.LE.1) THEN
             DMU1=DMU1*ZL/V
           END IF
         SDLNV=BN(NVAR+5)
         SDLNY=BN(NVAR+6)
         IF(MK34.EQ.1.AND.MIT.GE.1) THEN
           SDLNY=SDLNV
           BN(NVAR+6)=SDLNY
         END IF
C     <DISPERSION COEFFICEIENT>
         IF(MODD.EQ.2.AND.MIT.GE.1) THEN
           SDLND=SDLNV
           BN(NVAR+7)=SDLND
         ELSE
           SDLND=BN(NVAR+7)
         END IF
         CORR=BN(NVAR+8)
       GO TO 40
C------- MODE 4 --------------------------
   34    D=BN(NVAR+2)
         P=V*ZL/D
         DK=BN(NVAR+3)
         R=1.0+RHOTH*DK
         Y=DK
         BETA=1.0/R
         BETR=1.0
         OMEGA=BN(NVAR+4)
         DMU1=BN(NVAR+5)
         DMU2=BN(NVAR+6)
C      < DEGRADATION COEFFICIENTS >
        IF(MIT.GE.1) THEN
          IF(MDEG.EQ.1) THEN
            DMU2=(R-1.0)*DMU1
            BN(NVAR+6)=DMU2
          ELSE IF(MDEG.EQ.2) THEN
            DMU2=0
            BN(NVAR+6)=DMU2
          ELSE IF(MDEG.EQ.3) THEN
            DMU1=0.0
            BN(NVAR+5)=DMU1
          END IF
      END IF
C
         SDLNV=BN(NVAR+7)
         SDLNY=BN(NVAR+8)
         IF(MK34.EQ.1.AND.MIT.GE.1) THEN
           SDLNY=SDLNV
           BN(NVAR+8)=SDLNY
         END IF
C     <DISPERSION COEFFICEIENT>
         IF(MODD.EQ.2.AND.MIT.GE.1) THEN
           SDLND=SDLNV
           BN(NVAR+9)=SDLND
         ELSE
           SDLND=BN(NVAR+9)
         END IF
         CORR=BN(NVAR+10)
       GO TO 40
C----MODE 5 ----------------
   35    D=BN(NVAR+2)
         P=V*ZL/D
         DK=BN(NVAR+3)
         R=1.0+RHOTH*DK
         BETR=R
         Y=D
         DMU1=BN(NVAR+4)
           IF(NREDU.LE.1) THEN
             DMU1=DMU1*ZL/V
           END IF
         SDLNV=BN(NVAR+5)
         SDLNY=BN(NVAR+7)
         IF(MD56.EQ.1.AND.MIT.GE.1) THEN
           SDLNY=SDLNV
           BN(NVAR+7)=SDLNY
         END IF
C     <RETARDATION COEFFICIENT>
         IF(MODK.EQ.-2.AND.MIT.GE.1) THEN
           SDLNK=SDLNV
           BN(NVAR+6)=SDLNK
         ELSE
           SDLNK=BN(NVAR+6)
         END IF
         CORR=BN(NVAR+8)
       GO TO 40
C----MODE 6 ------------------
   36    D=BN(NVAR+2)
         P=V*ZL/D
         DK=BN(NVAR+3)
         R=1.0+RHOTH*DK
         Y=D
         BETA=1.0/R
         BETR=1.0
         OMEGA=BN(NVAR+4)
         DMU1=BN(NVAR+5)
         DMU2=BN(NVAR+6)
C      < DEGRADATION COEFFICIENTS >
        IF(MIT.GE.1) THEN
          IF(MDEG.EQ.1) THEN
            DMU2=(R-1.0)*DMU1
            BN(NVAR+6)=DMU2
          ELSE IF(MDEG.EQ.2) THEN
            DMU2=0
            BN(NVAR+6)=DMU2
          ELSE IF(MDEG.EQ.3) THEN
            DMU1=0.0
            BN(NVAR+5)=DMU1
          END IF
      END IF
C
         SDLNV=BN(NVAR+7)
         SDLNY=BN(NVAR+9)
         IF(MD56.EQ.1.AND.MIT.GE.1) THEN
           SDLNY=SDLNV
           BN(NVAR+9)=SDLNY
         END IF
C     <RETARDATION COEFFICIENT>
         IF(MODD.EQ.-2.AND.MIT.GE.1) THEN
           SDLNK=SDLNV
           BN(NVAR+8)=SDLNK
         ELSE
           SDLNK=BN(NVAR+8)
         END IF
         CORR=BN(NVAR+10)
       GO TO 40
   37 GO TO 40
C----MODE 8 -------------------
   38    D=BN(NVAR+2)
         P=V*ZL/D
         DK=BN(NVAR+3)
         R=1.0+RHOTH*DK
         ALPHA=BN(NVAR+4)
         Y=ALPHA
         OMEGA=ALPHA*(R-1.0)*ZL/V
         BETA=1.0/R
         BETR=1.0
         DMU1=BN(NVAR+5)
         DMU2=BN(NVAR+6)
C      < DEGRADATION COEFFICIENTS >
        IF(MIT.GE.1) THEN
          IF(MDEG.EQ.1) THEN
            DMU2=(R-1.0)*DMU1
            BN(NVAR+6)=DMU2
          ELSE IF(MDEG.EQ.2) THEN
            DMU2=0
            BN(NVAR+6)=DMU2
          ELSE IF(MDEG.EQ.3) THEN
            DMU1=0.0
            BN(NVAR+5)=DMU1
          END IF
      END IF
C
         SDLNV=BN(NVAR+7)
         SDLNY=BN(NVAR+10)
         IF(MAL8.EQ.1.AND.MIT.GE.1) THEN
           SDLNY=SDLNV
           BN(NVAR+10)=SDLNY
         END IF
C     <DISPERSION COEFFICEIENT>
         IF(MODD.EQ.2.AND.MIT.GE.1) THEN
           SDLND=SDLNV
           BN(NVAR+9)=SDLND
         ELSE
           SDLND=BN(NVAR+9)
         END IF
C     <RETARDATION COEFFICIENT>
         IF(MODD.EQ.-2.AND.MIT.GE.1) THEN
           SDLNK=SDLNV
           BN(NVAR+8)=SDLNK
         ELSE
           SDLNK=BN(NVAR+8)
         END IF
         CORR=BN(NVAR+11)
C
C -------------------------------------------
C MSTOCH;    Index for stochastic v and Y
C   MSTOCH=1: V; VARIABLE Y;CONSTANT
C   MSTOCH=2: V; CONSTANT Y;VARIABLE
C   MSTOCH=3: POSITIVE OR NEGATIVE CORRELATION BETWEEN  V & Y
C   MSTOCH=4: V; VARIABLE Y;VARIABLE
C
C ISKIP   Calculation control code for the
C         evaluation of the integral limits
C     ISKIP=0 ; CALCULATE INTEGRAL LIMITS
C     ISKIP=1 ;  SKIP THIS PART.
   40 IF (ISKIP.EQ.1) THEN
             GOTO 50
      ELSE IF (SDLNY.LT.1.D-07) THEN
             MSTOCH=1
             CORR=0.0
      ELSE IF (SDLNV.LT.1.D-07) THEN
             MSTOCH=2
             CORR=0.0
      ELSE
             MSTOCH=4
      END IF
C
C    MCORR=1 POSITIVE CORRELATION V&Y.
C    MCORR=-1 NEGATIVE CORRELATION V&Y.
C    MCORR=0 NO CORRELATION V&Y.
C    MCORR=2 BIVARIATE DISTRIBUTION V&Y.
      MCORR=0
      IF(MSTOCH.EQ.4) THEN
      IF(CORR.LT.-0.999999) THEN
          MCORR=-1
          MSTOCH=3
      ELSE IF(CORR.GT.0.999999) THEN
          MCORR=1
          MSTOCH=3
      ELSE IF(DABS(CORR).LT.1.D-20) THEN
          MCORR=0
      ELSE
          MCORR=2
      END IF
      END IF
C
C------- EVALUATE THE MAX AMD MIM  LIM IT FOR THE INTERGRATION  --------
C
      IF(MSTOCH.NE.2) THEN
        X=V
         SDLNX=SDLNV
         CALL LIMIT(X,SDLNX,XMIN,XMAX)
         VMIN=XMIN
         VMAX=XMAX
C        WRITE(*,*) 'V',V,SDLNV,VMIN,VMAX
      END IF
C     WRITE(*,*) 'MSTOCH',MSTOCH,MOD(MSTOCH,2)
      IF(MOD(MSTOCH,2).EQ.0) THEN
         X=Y
         SDLNX=SDLNY
         CALL LIMIT(X,SDLNX,XMIN,XMAX)
         YMIN=XMIN
         YMAX=XMAX
C         WRITE(*,*) 'Y',Y,SDLNY,YMIN,YMAX
      END IF
C
C   ADUJUST STOP CRITERIA FOR ROMBERG
      CALL CONST2(MODE,STOPER,LEVEL,MSTOCH,P,ICHEB)
C
C ------ CALULATE FILD SCALE CONCENTRATION ---------------
C
   50   CALL   STOCDE(CXT1,CXT2,VAR1,VAR2)
C
      ISKIP=1
      IF(DABS(CXT1).LT.1.D-30) CXT1=0.0D+0
      IF(DABS(CXT2).LT.1.D-30) CXT2=0.0D+0
      IF(DABS(VAR1).LT.1.D-30) VAR1=0.0D+0
      IF(DABS(VAR2).LT.1.D-30) VAR2=0.0D+0
C
      RETURN
      END
C-----------------------------------------------------
      SUBROUTINE DIRECT(T,Z,B,C1,C2,VAR1,VAR2,MPRINT,DT,DZ,NT,NZ,
     &  MAXOB,KP)
C
C        PURPOSE: CALCULATE DIRECT PROBLEMS
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION T(MAXOB),Z(MAXOB),C1(MAXOB),C2(MAXOB),
     & VAR1(MAXOB),VAR2(MAXOB),B(30)
      COMMON MODC,MM,A,BB,DA,CX,PEC,BETA,BETR,OMEGA,R,DMU1,DMU2,
     &   MODE,MCON,TT,ZZ,STOPER,LEVEL,ICHEB
      COMMON/MODAT/INDEX(15),INVERSE,NREDU,NVAR,ZL,MIT,MDEG,MASS,
     &         DUMTP(50),DUMGA1(50),MNEQ,ISKIP,PHIM,PHIIM
      COMMON/STOCH/MODD,MODK,MODS,MSTOCH,CORR,MCORR,MSD,SDLNK,
     & V,SDLNV,VMAX,VMIN,D,AVEY,SDLNY,YMAX,YMIN,ALPHA,SDLND
     & ,DK,RHOTH,MD56,MK34,MAL8
C
      WRITE(KP,*)
      WRITE(KP,*)
C
C    MPRINT=1 Conc. vs. Time
C          =2 Conc. vs. Depth
C
      IF(MPRINT.NE.2) GO TO 55
C----------------------------------------
C   CONCENTRATION vs. DEPTH
C
C
       DO 50 J=1,NT
       DO 51 I=1,NZ
C
         CALL MODEL(B,T(J),Z(I),C1(I),C2(I),VAR1(I),VAR2(I))
C
C
       IF(MOD(MODE,2).EQ.1.AND.MODE.GE.3) THEN
          IF(I.EQ.1.AND.J.EQ.1) THEN
            WRITE(*,1205)
          END IF
            WRITE(*,1202) Z(I),T(J),C1(I),VAR1(I)
       ELSE IF(MOD(MODE,2).EQ.0.AND.MODE.GE.4) THEN
          IF(I.EQ.1) THEN
            WRITE(*,1200)
            WRITE(*,1202) Z(I),T(J),C1(I),C2(I),VAR1(I),VAR2(I)
C             IF(MSTOCH.EQ.4) THEN
C                WRITE(KP,1200)
C                WRITE(KP,1202) Z(I),T(J),C1(I),C2(I),VAR1(I),VAR2(I)
C             END IF
         ELSE
C          IF(MSTOCH.EQ.4) THEN
C             WRITE(KP,1202) Z(I),T(J),C1(I),C2(I),VAR1(I),VAR2(I)
C           END IF
             WRITE(*,1202) Z(I),T(J),C1(I),C2(I),VAR1(I),VAR2(I)
         END IF
        END IF
   51  CONTINUE
C
      AREA1=0.0
      AREA2=0.0
      AMASS1=0.0
      AMASS2=0.0
      DO 60 I=1,NZ-1
       AREA1=AREA1+(C1(I)+C1(I+1))/2.0*DZ
       AREA2=AREA2+(C2(I)+C2(I+1))/2.0*DZ
   60 CONTINUE
C
      IF(MODE.GE.3) GO TO 29
C ****** DETERMINSTIC (MODE=1,2) ****************************
C ------ NONEQUILIBRIUM CDE (MODE=2)-------------------------
      IF(MODE.EQ.2) THEN
      DO 46 I=1,NZ-1
       AMASS1=AMASS1+BETA*R*(C1(I)+C1(I+1))/2.0*DZ
       AMASS2=AMASS2+(R-BETA*R)*(C2(I)+C2(I+1))/2.0*DZ
   46 CONTINUE
C    -- FLUX CONC.--
        IF(MODC.LE.2) THEN
         WRITE(KP,1000) T(J),AREA1,AREA2
         DO 41 I=1,NZ
          WRITE(KP,2020) Z(I), C1(I),C2(I)
   41    CONTINUE
        ELSE
C    -- RESIDENT AND TOTAL RESIDENT CONC.(SAME RESULT) (MODC>=3)--
         WRITE(KP,1100) T(J),AREA1,AREA2,AMASS1,AMASS2,AMASS1+AMASS2
         DO 61 I=1,NZ
          C1T=BETA*R*C1(I)
          C2T=(1-BETA)*R*C2(I)
          TCON=C1T+C2T
          WRITE(KP,2020) Z(I), C1(I),C2(I),C1T,C2T,TCON
   61   CONTINUE
        END IF
C
      ELSE IF (MODE.EQ.1) THEN
C ------ EQUILIBRIUM CDE (MODE=1)-------------------------
       IF(MODC.LE.2) THEN
C    -- FLUX CONC.--
        WRITE(KP,1001) T(J),AREA1
        DO 42 I=1,NZ
          WRITE(KP,2020) Z(I), C1(I)
   42   CONTINUE
       ELSE
C    -- RESIDENT AND TOTAL RESIDENT CONC.(SAME RESULT) (MODC>=3)--
        WRITE(KP,1101) T(J),AREA1,AREA1*R
        DO 62 I=1,NZ
          WRITE(KP,2020) Z(I),C1(I),R*C1(I)
   62   CONTINUE
       END IF
      END IF
      GO TO 50
C
C ****** STOCHASTIC (MODE>=3) ****************************
   29 IF(MODC.EQ.4.OR.MODC.EQ.6) THEN
C ------ TOTAL RESIDENT CONC  ---------------------------
          IF(INVERSE.EQ.-1) THEN
              WRITE(KP,1102) T(J),AREA1
          ELSE
              WRITE(KP,1103) T(J),AREA1
          END IF
             DO 34 I=1,NZ
           IF(INVERSE.EQ.-1) THEN
                WRITE(KP,2020) Z(I),C1(I)
           ELSE
                WRITE(KP,2020) Z(I),C1(I),VAR1(I)
           END IF
   34           CONTINUE
           GO TO 50
        ELSE
C  ----- NONEQUILIBRIUM FLUX AND RESIDENT CONC.  ---------
          IF(MOD(MODE,2).EQ.0) THEN
            IF(INVERSE.EQ.-1) THEN
              IF(NREDU.EQ.0) THEN
                 WRITE(KP,1109) T(J),AREA1,AREA2
               ELSE
                 WRITE(KP,1106) T(J),AREA1,AREA2
               END IF
            ELSE
              IF(NREDU.EQ.0) THEN
                 WRITE(KP,1108) T(J),AREA1,AREA2
              ELSE
                 WRITE(KP,1104) T(J),AREA1,AREA2
              END IF
            END IF
            DO 36 I=1,NZ
             IF(INVERSE.EQ.-1) THEN
               WRITE(KP,2020) Z(I), C1(I),C2(I)
             ELSE
               WRITE(KP,2020) Z(I), C1(I),C2(I),VAR1(I),VAR2(I)
             END IF
   36       CONTINUE
C
          ELSE
C  ----- EQUILIBRIUM FLUX AND RESIDENT CONC.  ---------
            IF(INVERSE.EQ.-1) THEN
              WRITE(KP,1107) T(J),AREA1
            ELSE
              WRITE(KP,1105) T(J),AREA1
            END IF
            DO 38 I=1,NZ
            IF(INVERSE.EQ.-1) THEN
             WRITE(KP,2020) Z(I),C1(I)
            ELSE
             WRITE(KP,2020) Z(I),C1(I),VAR1(I)
            END IF
   38       CONTINUE
          END IF
         END IF
   50 CONTINUE
      GO TO 99
C------------------------------------------------
C  CONCENTRATION vs. TIME
   55  DO 52 I=1,NZ
C
       DO 53 J=1,NT
C
        CALL MODEL(B,T(J),Z(I),C1(J),C2(J),VAR1(J),VAR2(J))
C
C
C       IF(MOD(MODE,2).EQ.1.AND.MODE.GE.3) THEN
       IF(MOD(MODE,2).EQ.1.AND.MODE.GE.2) THEN
          IF(I.EQ.1.AND.J.EQ.1) THEN
            WRITE(*,1205)
          END IF
            WRITE(*,1202) Z(I),T(J),C1(J),VAR1(J)
       ELSE IF(MOD(MODE,2).EQ.0.AND.MODE.GE.4) THEN
C       ELSE IF(MOD(MODE,2).EQ.0.AND.MODE.GE.2) THEN
         IF(J.EQ.1) THEN
          WRITE(*,1200)
          WRITE(*,1202) Z(I),T(J),C1(J),C2(J),VAR1(J),VAR2(J)
C          IF(MSTOCH.EQ.4) THEN
C            WRITE(KP,1200)
C            WRITE(KP,1202) Z(I),T(J),C1(J),C2(J),VAR1(J),VAR2(J)
C          END IF
         ELSE
          WRITE(*,1202) Z(I),T(J),C1(J),C2(J),VAR1(J),VAR2(J)
C          IF(MSTOCH.EQ.4) THEN
C            WRITE(KP,1202) Z(I),T(J),C1(J),C2(J),VAR1(J),VAR2(J)
C          END IF
         END IF
        END IF
   53  CONTINUE
C
C
C AREA1; O TH, AREA2; 1 st, AREA3; 2ND TIME MOMENT
      AREA11=0.0
      AREA12=0.0
      AREA13=0.0
      AREA21=0.0
      AREA22=0.0
      AREA23=0.0
      AREAT1=0.0
      AREAT2=0.0
      AREAT3=0.0
C   MOFORM = 1 ; WRITE 1 ST AND 2ND MOMENTS
      MOFORM=0
C TIME MOMENTS OF RESIDENT CONC. FOR THE DETERMINISTIC CDE
      IF(MODE.LE.2.AND.MODC.GE.3) THEN
        DO 63 J=1,NT-1
         AREA11=AREA11+(C1(J)+C1(J+1))/2.0*DT
         AREA12=AREA12+(T(J)*C1(J)+T(J)*C1(J+1))/2.0*DT
         AREA13=AREA13+(T(J)*T(J)*C1(J)+T(J)*T(J)*C1(J+1))/2.0*DT
         AREA21=AREA21+(C2(J)+C2(J+1))/2.0*DT
         AREA22=AREA22+(T(J)*C2(J)+T(J)*C2(J+1))/2.0*DT
         AREA23=AREA23+(T(J)*T(J)*C2(J)+T(J)*T(J)*C2(J+1))/2.0*DT
         C1J=BETA*R*C1(J)+(1-BETA)*R*C2(J)
         C1J1=BETA*R*C1(J+1)+(1-BETA)*R*C2(J+1)
         AREAT1=AREAT1+(C1J+C1J1)/2.0*DT
         AREAT2=AREAT2+(T(J)*C1J+T(J)*C1J1)/2.0*DT
         AREAT3=AREAT3+(T(J)*T(J)*C1J+T(J)*T(J)*C1J1)/2.0*DT
   63   CONTINUE
      ELSE
C TIME MOMENTS OF FLUX CONC. FOR THE DETERMINISTIC CDE
C  & THE STOCHASTIC CDE
      DO 65 J=1,NT-1
         AREA11=AREA11+(C1(J)+C1(J+1))/2.0*DT
         AREA12=AREA12+(T(J)*C1(J)+T(J)*C1(J+1))/2.0*DT
         AREA13=AREA13+(T(J)*T(J)*C1(J)+T(J)*T(J)*C1(J+1))/2.0*DT
         AREA21=AREA21+(C2(J)+C2(J+1))/2.0*DT
         AREA22=AREA22+(T(J)*C2(J)+T(J)*C2(J+1))/2.0*DT
         AREA23=AREA23+(T(J)*T(J)*C2(J)+T(J)*T(J)*C2(J+1))/2.0*DT
   65 CONTINUE
      END IF
C
      IF(MODE.GE.3) GO TO 75
      IF(MODE.EQ.2.AND.MODC.LE.2) THEN
C ****** DETERMINSTIC (MODE=1,2) ****************************
C ------ FLUX CONC. FOR NONEQUILIBRIUM CDE (MODE=2) ---------
        IF(MOFORM.EQ.1) THEN
        WRITE(KP,2110) Z(I),AREA11,AREA12,AREA13,AREA21,AREA22,AREA23
        ELSE
        WRITE(KP,1110) Z(I),AREA11,AREA21
        END IF
        DO 64 J=1,NT
         WRITE(KP,2020) T(J), C1(J), C2(J)
   64   CONTINUE
C
C ------ (TOTAL) RESIDENT CONC. FOR NONEQUILIBRIUM CDE (MODE=2) ---------
       ELSE IF(MODE.EQ.2.AND.MODC.GE.3) THEN
        IF(MOFORM.EQ.1) THEN
         WRITE(KP,2111) Z(I),AREA11,AREA12,AREA13,AREA21,AREA22,AREA23,
     & AREAT1,AREAT2,AREAT3
        ELSE
         WRITE(KP,1111) Z(I),AREA11,AREA21,AREAT1
        END IF
        DO 66 J=1,NT
         C1T=BETA*R*C1(J)
         C2T=(1-BETA)*R*C2(J)
         CT=C1T+C2T
         WRITE(KP,2020) T(J), C1(J), C2(J),C1T,C2T,CT
   66   CONTINUE
C
       ELSE IF(MODE.EQ.1.AND.MODC.LE.2) THEN
C ------ FLUX CONC. FOR EQUILIBRIUM CDE (MODE=1) ---------
        IF(MOFORM.EQ.1) THEN
        WRITE(KP,2112) Z(I),AREA11,AREA12,AREA13
        ELSE
        WRITE(KP,1112) Z(I),AREA11
        END IF
       DO 67 J=1,NT
          WRITE(KP,2020) T(J),C1(J)
   67    CONTINUE
       ELSE IF(MODE.EQ.1.AND.MODC.GE.3) THEN
C ------ (TOTAL) RESIDENT CONC. FOR EQUILIBRIUM CDE (MODE=1) ---------
        IF(MOFORM.EQ.1) THEN
        WRITE(KP,2113) Z(I),AREA11,AREA12,AREA13,AREAT1,AREAT2,AREAT3
        ELSE
        WRITE(KP,1113) Z(I),AREA11,AREAT1
        END IF
        DO 87 J=1,NT
          WRITE(KP,2020) T(J),C1(J),C1(J)*R
   87    CONTINUE
      END IF
C
      GO TO 52
C ****** STOCHASTIC (MODE>=3) ****************************
C ------ TOTAL RESIDENT CONC. ---------
C
   75  IF(MODC.EQ.4.OR.MODC.EQ.6) THEN
           IF(INVERSE.EQ.-1) THEN
            IF(MOFORM.EQ.1) THEN
              WRITE(KP,2117) Z(I),AREA11,AREA12,AREA13
             ELSE
             WRITE(KP,1117) Z(I),AREA11
            END IF
           ELSE
            IF(MOFORM.EQ.1) THEN
             WRITE(KP,2118) Z(I),AREA11,AREA12,AREA13
            ELSE
             WRITE(KP,1118) Z(I),AREA11
            END IF
           END IF
            DO 78 J=1,NT
           IF(INVERSE.EQ.-1) THEN
            WRITE(KP,2020) T(J),C1(J)
           ELSE
            WRITE(KP,2020) T(J),C1(J),VAR1(J)
           END IF
   78  CONTINUE
       GO TO 52
       END IF
C   ----- NONEQUILBRIUM (MODE=4,6,8) -------------------
      IF(MOD(MODE,2).EQ.0) THEN
           IF(INVERSE.EQ.-1) THEN
            IF(NREDU.EQ.0) THEN
             IF(MOFORM.EQ.1) THEN
               WRITE(KP,2115) Z(I),AREA11,AREA12,AREA13,AREA21,AREA22,
     &                        AREA23
              ELSE
               WRITE(KP,1115) Z(I),AREA11,AREA21
              END IF
            ELSE
             IF(MOFORM.EQ.1) THEN
              WRITE(KP,2210) Z(I),AREA11,AREA12,AREA13,AREA21,AREA22,
     &                       AREA23
             ELSE
               WRITE(KP,1210) Z(I),AREA11,AREA21
             END IF
            END IF
           ELSE
            IF(NREDU.EQ.0) THEN
             IF(MOFORM.EQ.1) THEN
              WRITE(KP,2114) Z(I),AREA11,AREA12,AREA13,AREA21,AREA22,
     &                       AREA23
             ELSE
              WRITE(KP,1114) Z(I),AREA11,AREA21
             END IF
            ELSE
             IF(MOFORM.EQ.1) THEN
              WRITE(KP,2211) Z(I),AREA11,AREA12,AREA13,AREA21,AREA22,
     &                       AREA23
             ELSE
              WRITE(KP,1211) Z(I),AREA11,AREA21
             END IF
            END IF
           END IF
        DO 76 J=1,NT
           IF(INVERSE.EQ.-1) THEN
                WRITE(KP,2020) T(J),C1(J),C2(J)
           ELSE
                WRITE(KP,2020) T(J),C1(J),C2(J),VAR1(J),VAR2(J)
           END IF
   76   CONTINUE
C   ----- EQUILBRIUM (MODE=6,8) -------------------
      ELSE
C
           IF(INVERSE.EQ.-1) THEN
             IF(MOFORM.EQ.1) THEN
              WRITE(KP,2119) Z(I),AREA11,AREA12,AREA13
             ELSE
              WRITE(KP,1119) Z(I),AREA11
             END IF
           ELSE
             IF(MOFORM.EQ.1) THEN
              WRITE(KP,2116) Z(I),AREA11,AREA12,AREA13
             ELSE
              WRITE(KP,1116) Z(I),AREA11
             END IF
           END IF
        DO 77 J=1,NT
           IF(INVERSE.EQ.-1) THEN
              WRITE(KP,2020) T(J), C1(J)
           ELSE
              WRITE(KP,2020) T(J), C1(J),VAR1(J)
           END IF
   77   CONTINUE
      END IF
C
   52  CONTINUE
      GO TO 99
C
C
 1000 FORMAT('$ T=',F9.5,'   (FLUX-AVERAGED CONC. VS. DEPTH)',
     &/,'$ Sum(C1*dZ) = ',F11.5,'   Sum(C2*dZ) = ',F11.5,
     &/'$    Z   ',8X,'C1',12X,'C2')
 1001 FORMAT('$ T=',F9.5,'   (FLUX-AVERAGED CONC. VS. DEPTH)',
     &/,'$ Sum(C*dZ) = ',F12.5,
     &/,'$    Z   ',8X,' C')
 1100 FORMAT('$ T=',F8.5,'   (RESIDENT AND TOTAL RESIDENT CONC. VS. DEPT
     &H)',/,'$ Sum(C1*dZ) = ',F11.5,'   Sum(C2*dZ) = ',F11.5,
     &/,'$ MASS  PHASE 1 = ',F11.5,'  PHASE 2 =',F11.5,'  TOTAL = ',F11.
     &5,/,'$    Z   ',10X,'C1',11X,'C2',9X,'BR*C1',6X,'(1-B)R*C2',8X,'Ct
     &')
 1101 FORMAT('$ T=',F9.5,'    (RESIDENT AND TOTAL RESIDENT CONC. VS. DEP
     &TH)',/,'$ Sum(C*dZ) = ',F11.5,'   Sum(Ct*dZ) = ',F11.5,' (TOTAL MA
     &SS)',/,'$    Z   ',5X,'  C  ',9X,'Ct (=R*C)')
 1102 FORMAT('  No results for variance (INVERSE = -1)',//,
     &'$ T=',F9.5,'   (TOTAL RESIDENT CONC. VS. DEPTH)',
     &/,'$ Sum(Ct*dZ) = ',F12.5,' (TOTAL MASS)',
     &/,'$    Z   ',6X,'<Ct>')
 1103 FORMAT('$ T=',F9.5,'    (TOTAL RESIDENT CONC. VS. DEPTH)',
     &/,'$ Sum(Ct*dZ) = ',F12.5,' (TOTAL MASS)',
     &/,'$    Z   ',6X,'<Ct>',10X,'VAR(Ct)')
 1104 FORMAT('$ T=',F9.5,'     (CONC. VS. DEPTH)',
     &/,'$ Sum(C1*dZ) = ',F11.5,'   Sum(C2*dZ) = ',F11.5,
     &/,'$    Z   ',6X,'<C1>',10X,'<C2>',10X, 'VAR(C1)',7X,'VAR(C2)')
 1105 FORMAT('$ T=',F9.5,'     (CONC. VS. DEPTH)',
     &/,'$ Sum(C*dZ) = ',F12.5,
     &/,'$    Z   ',6X,' <C>',10X,'VAR(C)')
 1106 FORMAT('  No results for variance (INVERSE = -1)',//,
     &'$ T=',F9.5,'     (CONC. VS. DEPTH)',
     &/,'$ Sum(C1*dZ) = ',F11.5,'    Sum(C2*dZ) = ',F11.5,
     &/,'$    Z   ',6X,'<C1>',10X,'<C2>')
 1107 FORMAT('  No results for variance (INVERSE = -1)',//,
     &'$ T=',F9.5,'     (CONC. VS. DEPTH)',
     &/,'$ Sum(C*dZ) = ',F12.5,
     &/,'$    Z   ',6X,' <C>')
 1108 FORMAT('$ T=',F9.5,'     (CONC. VS. DEPTH)',
     &/,'$ Sum(C1*dZ) = ',F11.5,'   Sum(s*dZ) = ',F11.5,
     &/,'$    Z   ',6X,'<C1>',11X,'<s>',10X, 'VAR(C1)',7X,'VAR(s)')
 1109 FORMAT('  No results for variance (INVERSE = -1)',//,
     &'$ T=',F9.5,'     (CONC. VS. DEPTH)',
     &/,'$ Sum(C1*dZ) = ',F11.5,'    Sum(s*dZ) = ',F11.5,
     &/,'$    Z   ',6X,'<C1>',11X,'<s>')
 1110 FORMAT('$ Z=',F8.4,'     (FLUX CONC. VS. TIME)',
     &/,'$ Sum(C1*dT)=',F12.4,', Sum(C2*dT)=',F12.4)
 1111 FORMAT('$ Z=',F8.4,'     (RESIDENT CONC. VS. TIME)',
     &/,'$ Sum(C1*dT)=',F11.4,', Sum(C2*dT)=',F11.4,', Sum(Ct*dT)=',F11.
     &4,/,'$',3X,'TIME ',9X,'C1',11X,'C2',9X,'BR*C1',7X,'(1-B)R*C2',7X,
     &'Ct')
 1112 FORMAT('$ Z=',F8.4,'     (FLUX CONC. VS. TIME)',
     &/,'$ Sum(C*dT)=',F12.4,/,'$',3X,'TIME ',9X,'C')
 1113 FORMAT('$ Z=',F8.4,'     (RESIDENT CONC. VS. TIME)',
     &/,'$ Sum(C*dT) =',F12.4,',  Sum(Ct*dT)=',F12.4,
     &/,'$',3X,'TIME ',9X,'C',9X,'Ct (=R*C)')
 1114 FORMAT('$ Z=',F8.4,'    (CONC. VS. TIME)'
     &/,'$ Sum(C1*dT)=',F12.4,', Sum(s*dT) =',F12.4,
     &/,'$',3X,'TIME ',7X,'<C1>',10X,' <s>',10X,'VAR(C1)',7X,' VAR(s)')
 1115 FORMAT('  No results for variance (INVERSE = -1)',//,
     &'$ Z=',F8.4,'     (CONC. VS. TIME)'
     &/,'$ Sum(C1*dT)=',F12.4,', Sum(s*dT) =',F12.4,
     &/,'$',3X,'TIME ',7X,'<C1>',10X,' <s>')
 1116 FORMAT('$ Z=',F8.4,'     (CONC. VS. TIME)'
     &/,'$ Sum(C*dT)=',F12.4,/,'$',3X,'TIME ',7X,'<C>',10X,'VAR(C)')
 1117 FORMAT('  No results for variance (INVERSE = -1)',//,
     &'$ Z=',F8.4,'     (TOTAL RESIDENT CONC. VS. TIME)'
     &/,'$ Sum(Ct*dT)=',F12.4,
     &/,'$',3X,'TIME ',7X,'<Ct>')
 1118 FORMAT('$ Z=',F8.4,'     (TOTAL RESIDENT CONC. VS. TIME)'
     &/,'$ Sum(Ct*dT)=',F12.4,
     &/,'$',3X,'TIME ',7X,'<Ct>',9X,'VAR(Ct)')
 1119 FORMAT('  No results for variance (INVERSE = -1)',//,
     &'$ Z=',F8.4,'     (CONC. VS. TIME)'
     &/,'$ Sum(C*dT)=',F12.4,/,'$',3X,'TIME ',7X,'<C>')
 1200 FORMAT(//,5X,'Z',9X,'T',9X,'C1',9X,'C2',10X,'VAR1',8X,'VAR2')
 1202 FORMAT(2F9.4,4E12.4)
 1205 FORMAT(//,5X,'Z',9X,'T',9X,'C1',9X,'VAR1')
 1210 FORMAT('  No results for variance (INVERSE = -1)',//,
     &'$ Z=',F8.4,'      (CONC. VS. TIME)'
     &/,'$ Sum(C1*dT)=',F12.4,', Sum(C2*dT)=',F12.4,
     &/,'$',3X,'TIME ',7X,'<C1>',10X,'<C2>')
 1211 FORMAT('$ Z=',F8.4,'      (CONC. VS. TIME)'
     &/,'$ Sum(C1*dT)=',F12.4,',  Sum(C2*dT)=',F12.4,
     &/,'$',3X,'TIME ',7X,'<C1>',10X,'<C2>',10X,'VAR(C1)',7X,'VAR(C2)')
 2020 FORMAT(F11.4,5E13.5)
 2110 FORMAT('$ Z=',F8.4,'     (FLUX CONC. VS. TIME)',
     &/,'$ Sum(C1*dT)=',F9.4,',  Sum(T*C1*dT)=',F9.4,
     &',  Sum(T*T*C1*dT)=',F11.4,
     &/,'$ Sum(C2*dT)=',F9.4,',  Sum(T*C2*dT)=',F9.4,
     &',  Sum(T*T*C2*dT)=',F11.4,
     &/,'$',3X,'TIME ',9X,'C1',11X,'C2')
 2111 FORMAT('$ Z=',F8.4,'     (RESIDENT CONC. VS. TIME)',
     &/,'$ Sum(C1*dT)=',F9.4,',  Sum(T*C1*dT)=',F9.4,
     &',  Sum(T*T*C1*dT)=',F11.4,
     &/,'$ Sum(C2*dT)=',F9.4,',  Sum(T*C2*dT)=',F9.4,
     &',  Sum(T*T*C2*dT)=',F11.4,
     &/,'$ Sum(Ct*dT)=',F9.4,',  Sum(T*Ct*dT)=',F9.4,
     &',  Sum(T*T*Ct*dT)=',F11.4,
     &/,'$',3X,'TIME ',9X,'C1',11X,'C2',9X,'BR*C1',7X,'(1-B)R*C2',7X,
     &'Ct')
 2112 FORMAT('$ Z=',F8.4,'     (RESIDENT CONC. VS. TIME)',
     &/,'$ Sum(C1*dT)=',F11.4,', Sum(C2*dT)=',F11.4,', Sum(Ct*dT)=',F11.
     &4,/,'$',3X,'TIME ',9X,'C1',11X,'C2',9X,'BR*C1',7X,'(1-B)R*C2',7X,
     &'Ct')
 2113 FORMAT('$ Z=',F8.4,'     (RESIDENT CONC. VS. TIME)',
     &/,'$ Sum(C*dT) =',F9.4,',  Sum(T*C*dT) =',F9.4,
     &',  Sum(T*T*C*dT) =',F11.4,
     &/,'$ Sum(Ct*dT)=',F9.4,',  Sum(T*Ct*dT)=',F9.4,
     &',  Sum(T*T*Ct*dT)=',F11.4,
     &/,'$',3X,'TIME ',9X,'C',9X,'Ct (=R*C)')
 2114 FORMAT('$ Z=',F8.4,'    (CONC. VS. TIME)'
     &/,'$ Sum(C1*dT)=',F9.4,',  Sum(T*C1*dT)=',F9.4,
     &',  Sum(T*T*C1*dT)=',F11.4,
     &/,'$ Sum(s*dT) =',F9.4,',  Sum(T*s*dT) =',F9.4,
     &',  Sum(T*T*s*dT) =',F11.4,
     &/,'$',3X,'TIME ',7X,'<C1>',10X,' <s>',10X,'VAR(C1)',7X,' VAR(s)')
 2115 FORMAT('  No results for variance (INVERSE = -1)',//,
     &'$ Z=',F8.4,'     (CONC. VS. TIME)'
     &/,'$ Sum(C1*dT)=',F9.4,',  Sum(T*C1*dT)=',F9.4,
     &',  Sum(T*T*C1*dT)=',F11.4,
     &/,'$ Sum(s*dT) =',F9.4,',  Sum(T*s*dT) =',F9.4,
     &',  Sum(T*T*s*dT) =',F11.4,
     &/,'$',3X,'TIME ',7X,'<C1>',10X,' <s>')
 2116 FORMAT('$ Z=',F8.4,'     (CONC. VS. TIME)'
     &/,'$ Sum(C*dT)=',F9.4,',  Sum(T*C*dT)=',F9.4,
     &',  Sum(T*T*C*dT)=',F11.4,
     &/,'$',3X,'TIME ',7X,'<C>',10X,'VAR(C)')
 2117 FORMAT('  No results for variance (INVERSE = -1)',//,
     &'$ Z=',F8.4,'     (TOTAL RESIDENT CONC. VS. TIME)'
     &/,'$ Sum(Ct*dT)=',F9.4,',  Sum(T*Ct*dT)=',F9.4,
     &',  Sum(T*T*Ct*dT)=',F11.4,
     &/,'$',3X,'TIME ',7X,'<Ct>')
 2118 FORMAT('$ Z=',F8.4,'     (TOTAL RESIDENT CONC. VS. TIME)'
     &/,'$ Sum(Ct*dT)=',F9.4,',  Sum(T*Ct*dT)=',F9.4,
     &',  Sum(T*T*Ct*dT)=',F11.4,
     &/,'$',3X,'TIME ',7X,'<Ct>',9X,'VAR(Ct)')
 2119 FORMAT('  No results for variance (INVERSE = -1)',//,
     &'$ Z=',F8.4,'     (CONC. VS. TIME)'
     &/,'$ Sum(C*dT)=',F9.4,',  Sum(T*C*dT)=',F9.4,
     &',  Sum(T*T*C*dT)=',F11.4,
     &/,'$',3X,'TIME ',7X,'<C>')
 2210 FORMAT('  No results for variance (INVERSE = -1)',//,
     &'$ Z=',F8.4,'      (CONC. VS. TIME)'
     &/,'$ Sum(C1*dT)=',F9.4,',  Sum(T*C1*dT)=',F9.4,
     &',  Sum(T*T*C1*dT)=',F11.4,
     &/,'$ Sum(C2*dT)=',F9.4,',  Sum(T*C2*dT)=',F9.4,
     &',  Sum(T*T*C2*dT)=',F11.4,
     &/,'$',3X,'TIME ',7X,'<C1>',10X,'<C2>')
 2211 FORMAT('$ Z=',F8.4,'      (CONC. VS. TIME)',
     &/,'$ Sum(C1*dT)=',F9.4,',  Sum(T*C1*dT)=',F9.4,
     &',  Sum(T*T*C1*dT)=',F11.4,
     &/,'$ Sum(C2*dT)=',F9.4,',  Sum(T*C2*dT)=',F9.4,
     &',  Sum(T*T*C2*dT)=',F11.4,
     &/,'$',3X,'TIME ',7X,'<C1>',10X,'<C2>',10X,'VAR(C1)',7X,'VAR(C2)')
   99 RETURN
      END
