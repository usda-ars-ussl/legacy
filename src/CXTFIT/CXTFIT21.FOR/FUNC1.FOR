C     ------------------------------------------------------------------
      REAL*8 FUNCTION CTRAN(TAU,MC)
C
C     PURPOSE: CALCULATE ARGUMENT IN INTEGRAL
C                     FOR DELTA INPUT (TRANSFER FUNCTION MODEL)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON MODC,MM,A,B,DA,CX,P,BETA,BETR,OMEGA,R,DMU1,DMU2,MODE,
     &       MCON,T,Z,STOPER,LEVEL,ICHEB
C
      CTRAN=0.0
      G=CC0(TAU)
      IF(G.LT.1.D-07) RETURN
C
      AA=OMEGA*TAU/BETR
      AT=A*TAU
      BT=B*(T-TAU)
      XII=2.0*DSQRT(AT*BT)
      BETA=BETR/R
C
      IF(MC.EQ.1) THEN
        CBI1=DSQRT(TAU/((1-BETA)*BETA*(T-TAU)))*OMEGA/R
        CTRAN=G*EXPBI1(XII,-AA-BT)*CBI1
      ELSE                                          
      
        CTRAN=OMEGA/(R-BETR)*G*EXPBI0(XII,-AA-BT)
      END IF
      RETURN
      END
C
C     ------------------------------------------------------------------
C
      REAL*8 FUNCTION CBJ(TAU,MC)
C
C     PURPOSE: CALCULATE ARGUMENT IN INTEGRAL FOR STEP INPUT
C               (SOLUTION USING GOLDSTEIN'S J-FUNCTION)
C     updated  April 25, 1997
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON MODC,MM,A,B,DA,CX,P,BETA,BETR,OMEGA,R,DMU1,DMU2,MODE,
     &       MCON,T,Z,STOPER,LEVEL,ICHEB
C
      CBJ=0.0
      G=CC0(TAU)
      IF(G.LT.1.D-07) RETURN
        AT=A*TAU
        BT=B*(T-TAU)
      IF(MC.EQ.1) THEN
        CBJ=G*GOLD(AT,BT)*DBEXP(-AT*DMU2/OMEGA)
        IF(CBJ.LT.1.D-07) CBJ=0.0
      ELSE
          IF(MODC.NE.3.AND.MODC.NE.4) THEN
             IF(Z.LT.1.D-07.AND.T.GT.0.D00) THEN
                CBJ=OMEGA/(R-BETR)*DBEXP(-BT)
                RETURN 
             END IF
           END IF
            CBJ=CX*G*(1.0-GOLD(BT,AT))*DBEXP(-AT*DMU2/OMEGA)
            IF(CBJ.LT.1.D-07) CBJ=0.0
      END IF
      RETURN
      END
C
C     ------------------------------------------------------------------
      REAL*8 FUNCTION CBAL(TAU,MC)
C
C     PURPOSE: CALCULATE ARGUMENT IN INTEGRAL FOR FOR STEP INPUT
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON MODC,MM,A,B,DA,CX,P,BETA,BETR,OMEGA,R,DMU1,DMU2,MODE,
     &       MCON,T,Z,STOPER,LEVEL,ICHEB
      CBAL=0.0
      G=CC4(TAU,DA)
      IF(G.LT.1.D-07) RETURN
C
      BETA=BETR/R
      AT=A*TAU
      BT=B*(T-TAU)
      XII=2.0*DSQRT(AT*BT)
C
      IF(MC.EQ.1) THEN
        CBI1=DSQRT(TAU/((1-BETA)*BETA*(T-TAU)))
        CBI0=CX/BETA
        CBAL=EXPBI0(XII,-AT-BT)*CBI0+EXPBI1(XII,-AT-BT)*CBI1
        CBAL=OMEGA/R*G*CBAL
      ELSE
        CBI1=DSQRT((1-BETA)*(T-TAU)/TAU/BETA)*CX
        CBAL=EXPBI0(XII,-AT-BT)+EXPBI1(XII,-AT-BT)*CBI1
        CBAL=OMEGA/(R-BETR)*G*CBAL
      END IF
C
      RETURN
      END
C
C     ------------------------------------------------------------------
      REAL*8 FUNCTION CBEXP(TAU,MC)
C
C     PURPOSE: CALCULATE ARGUMENT IN INTEGRAL FOR EXPONENTIAL INPUT
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON MODC,MM,A,B,DA,CX,P,BETA,BETR,OMEGA,R,DMU1,DMU2,MODE,
     &       MCON,T,Z,STOPER,LEVEL,ICHEB
      COMMON/BOUN/MODB,NPULSE,TPULSE(50),PULSE(50),MASSST
C
      CBEXP=0.0
      G=CC0(TAU)
      IF(G.LT.1.D-07) RETURN
C
      IF(MODE.EQ.1) THEN
      CBEXP=G*DBEXP(-TPULSE(1)*(T-TAU))
      RETURN
      END IF
C
      C1=-OMEGA*TAU/BETR+B*TAU
      C2=-OMEGA*TAU/BETR+TPULSE(1)*TAU
      C3=TPULSE(1)-B
      IF(MC.EQ.1) THEN
        CBEXP=G*(DBEXP(-TPULSE(1)*T+C2-A*B/C3*TAU)
     &   -DBEXP(-B*T+C1)*PHI1(TAU))
      ELSE
       CBEXP=G*(DBEXP(-B*T+C1-A*B*TAU/C3)/C3-DBEXP(-TPULSE(1)*T
     &  +C2-A*B*TAU/C3)/C3-DBEXP(-B*T+C1)*PHI2(TAU))
      END IF
C
      RETURN
      END
C
C     ------------------------------------------------------------------
      REAL*8 FUNCTION CBIN1(TAU)
C
C     PURPOSE: CALCULATE ARGUMENT IN INTEGRAL FOR ARBITARY INPUT 
C                GIVEN IN FUNCTION CINPUT 
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON MODC,MM,A,B,DA,CX,P,BETA,BETR,OMEGA,R,DMU1,DMU2,MODE,
     &       MCON,T,Z,STOPER,LEVEL,ICHEB
      EXTERNAL CTRAN
C
      TT=T
      A1=0.0
      CBIN1=0.0
      G=CC0(TAU)
      IF(G.LT.1.D-07.AND.MODE.EQ.1) RETURN
C
         C1=G*DBEXP(-OMEGA*TAU/BETR)
      IF(MODE.EQ.1) GOTO 10
           MC=1
           T=TAU
           CALL CHEBYCON(CTRAN,A1,0.D00,TAU,MC)
   10   T=TT
           CBIN1=(C1+A1)*CINPUT(T-TAU)
C
      RETURN
      END
C
C     ------------------------------------------------------------------
      REAL*8 FUNCTION CBIN2(TAU)
C
C     PURPOSE: CALCULATE ARGUMENT IN INTEGRAL FOR ARBITARY INPUT 
C                GIVEN IN FUNCTION CINPUT FOR NONEQUILIBRIUM PHASE
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON MODC,MM,A,B,DA,CX,P,BETA,BETR,OMEGA,R,DMU1,DMU2,MODE,
     &       MCON,T,Z,STOPER,LEVEL,ICHEB
      EXTERNAL CBIN1
C
      TT=T
      CBIN2=0.0
C
           T=TAU
           CALL CHEBY2(CBIN1,A1,0.D00,TAU)
           T=TT
           CBIN2=A1*DBEXP(-B*(T-TAU))
C
      RETURN
      END
C
C     ------------------------------------------------------------------
      REAL*8 FUNCTION CIVP(TAU,MC)
C
C     PURPOSE: CALCULATE ARGUMENT IN INITIAL VALUE PROBLEM
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON MODC,MM,A,B,DA,CX,P,BETA,BETR,OMEGA,R,DMU1,DMU2,MODE,
     &       MCON,T,Z,STOPER,LEVEL,ICHEB
      COMMON/INITI/MODI,NINI,CINI(50),ZINI(50)
C
C    
      BETA=BETR/R
      CIVP=0.0
      DG=DBEXP(-(DMU1+DA)/BETR*TAU)
C      
C      MODI=1   CONSTANT INITIAL CONCENTRATION
C      MODI=2   STEPWISE INITIAL CONCENTRATION
C      MODI=3   EXPONENTIAL INITIAL CONCENTRATION
C      MODI=4   DELTA(Z-Z0) & CONSTANT BACKGROUND
C 
      GO TO (10,10,25,40) MODI
   10 DO 15 I=1,NINI
      IF (I.EQ.1) THEN
        G=-CINI(I)*DG*(CC1(TAU)-1.0)
       ELSE
        G=G+(CINI(I-1)-CINI(I))*DG*(CC2(TAU,ZINI(I))-1.0)
       END IF
   15 CONTINUE
      GO TO 30
C
   25 G=CINI(1)*DG*(1.0-CC1(TAU))+CINI(2)*DG*CC3(TAU,ZINI(1))
      GO TO 30
C
   40 IF(MODC.LE.4.AND.ABS(ZINI(2)).LT.1.E-10) THEN
       G=DG*CC0(TAU)*BETA*R
      ELSE
       G=DG*CC5(TAU,ZINI(2))
      END IF
C
   30   AT=A*TAU
        BT=B*(T-TAU)
        XII=2.0*DSQRT(AT*BT)
C
      IF(MC.EQ.1) THEN
         CBI1=DSQRT(TAU/((1-BETA)*BETA*(T-TAU)))*OMEGA/R
         CIVP=EXPBI0(XII,-AT-BT)*OMEGA/BETR+EXPBI1(XII,-AT-BT)*CBI1
         CIVP=G*CIVP
      ELSE
         CBI1=DSQRT((1-BETA)*(T-TAU)/BETA/TAU)
         CIVP=EXPBI0(XII,-AT-BT)+EXPBI1(XII,-AT-BT)*CBI1
         CIVP=OMEGA/(R-BETR)*G*CIVP
      END IF
      RETURN
      END
C
C     ------------------------------------------------------------------
      REAL*8 FUNCTION C1PRO(TAU)
C
C     PURPOSE: CALCULATE ARGUMENT IN PRODUCTION TERM
C                    FOR EQUILIBRIUM CONCENTRATION
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON MODC,MM,A,B,DA,CX,P,BETA,BETR,OMEGA,R,DMU1,DMU2,MODE,
     &        MCON,T,Z,STOPER,LEVEL,ICHEB
      COMMON/PROD/MODP,NPRO1,NPRO2,GAMMA1(50),ZPRO1(50),GAMMA2(50),
     &   ZPRO2(50)
C
C          CX=OMEGA/(OMEGA+DMU2)
C
      C1PRO=0.0
      G=0.0
      H=0.0
      DG=DBEXP(-(DMU1+DA)/BETR*TAU)
C
C      MODP=0   ZERO PRODUCTION
C      MODP=1   CONSTANT PRODUCTION
C      MODP=2   STEPWISE PRODUCTION
C      MODP=3   EXPONENTIAL PRODUCTION
C
      GO TO (10,10,20) MODP
   10  IF(NPRO1.EQ.0) GO TO 17
       DO 15 I=1,NPRO1
       IF(I.EQ.1) THEN
         G=-GAMMA1(1)/BETR*DG*(CC1(TAU)-1.0)
        ELSE
         G=G+(GAMMA1(I-1)-GAMMA1(I))/BETR*DG*(CC2(TAU,ZPRO1(I))-1.0)
        END IF
   15  CONTINUE
C
   17  IF(MODE.EQ.1) GOTO 40
       IF(NPRO2.EQ.0) GO TO 18
       DO 16 I=1,NPRO2
       IF(I.EQ.1) THEN
         H=-CX*GAMMA2(1)/BETR*DG*(CC1(TAU)-1.0)
        ELSE
         H=H+CX*(GAMMA2(I-1)-GAMMA2(I))/BETR*DG*(CC2(TAU,ZPRO2(I))-1.0)
        END IF
   16  CONTINUE
C
   18   G=G+H
C
        GO TO 40
   20 G=(GAMMA1(1)*(1.0-CC1(TAU))+GAMMA1(2)*CC3(TAU,ZPRO1(1)))*DG/BETR
      IF(MODE.EQ.1) GOTO 40
      H=(GAMMA2(1)*(1.0-CC1(TAU))+GAMMA2(2)*CC3(TAU,ZPRO2(1)))
     &               *DG/BETR*CX
        G=G+H
C   40 IF(G.LT.1.D-07) RETURN
   40 CONTINUE
      IF(MODE.EQ.1) THEN
      C1PRO=G
      RETURN
      END IF
C
      AT=A*TAU
      BT=B*(T-TAU)
      XII=2.0*DSQRT(AT*BT)
      C1PRO=G*GOLD(AT,BT)-H*EXPBI0(XII,-AT-BT)
      RETURN
      END
C
C     ------------------------------------------------------------------
      REAL*8 FUNCTION C2PRO(TAU)
C
C     PURPOSE: CALCULATE ARGUMENT IN PRODUCTION TERM
C                    FOR NONEQUILIBRIUM CONCENTRATION
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON MODC,MM,A,B,DA,CX,P,BETA,BETR,OMEGA,R,DMU1,DMU2,MODE,
     &       MCON,T,Z,STOPER,LEVEL,ICHEB
      COMMON/PROD/MODP,NPRO1,NPRO2,GAMMA1(50),ZPRO1(50),GAMMA2(50),
     &  ZPRO2(50)
C
C          CX=OMEGA/(OMEGA+DMU2)
C
      C2PRO=0.0
      G=0.0
      H=0.0
      DG=DBEXP(-(DMU1+DA)/BETR*TAU)
C
C      MODP=0   ZERO PRODUCTION
C      MODP=1   CONSTANT PRODUCTION
C      MODP=2   STEPWISE PRODUCTION
C      MODP=3   EXPONENTIAL PRODUCTION
C
      GO TO (10,10,20) MODP
   10  IF(NPRO1.EQ.0) GO TO 17
       DO 15 I=1,NPRO1
       IF(I.EQ.1) THEN
         G=-GAMMA1(1)/BETR*DG*(CC1(TAU)-1.0)
        ELSE
         G=G+(GAMMA1(I-1)-GAMMA1(I))/BETR*DG*(CC2(TAU,ZPRO1(I))-1.0)
        END IF
   15  CONTINUE
C
   17  IF(NPRO2.EQ.0) GO TO 18
       DO 16 I=1,NPRO2
       IF(I.EQ.1) THEN
         H=-CX*GAMMA2(1)/BETR*DG*(CC1(TAU)-1.0)
        ELSE
         H=H+CX*(GAMMA2(I-1)-GAMMA2(I))/BETR*DG*(CC2(TAU,ZPRO2(I))-1.0)
        END IF
   16  CONTINUE
C
   18   G=G+H
C
        GO TO 40
   20 G=(GAMMA1(1)*(1.0-CC1(TAU))+GAMMA1(2)*CC3(TAU,ZPRO1(1)))*DG/BETR
      H=(GAMMA2(1)*(1.0-CC1(TAU))+GAMMA2(2)*CC3(TAU,ZPRO2(1)))
     &        *DG/BETR*CX
        G=G+H
C
C   40 IF(G.LT.1.D-07) RETURN
   40  CONTINUE
C
      AT=A*TAU
      BT=B*(T-TAU)
      BETA=BETR/R
      XII=2.0*DSQRT(AT*BT)
      CBI1=DSQRT(BETA/(1-BETA)*(T-TAU)/TAU)
      C2PRO=CX*G*(1.0-GOLD(BT,AT))-CBI1*H*EXPBI1(XII,-AT-BT)
C
      RETURN
      END
C
C     ------------------------------------------------------------------
      REAL*8 FUNCTION CC0(TAU)
C
C     PURPOSE: CALCULATE SOLUTIONS FOR DELTA INPUT
C              TRAVEL TIME DISTRIBUTION FOR EQUILIBRIUM CDE
C      (ARGUMENT IN BVP FOR STEP INPUT WITH GOLDSTEIN'S J-FUNCTION)
C      WHEN DMU1=0   CC0(TAU)=d(CC1(TAU))/dTAU 
C      IF NOT           CC0(TAU)=d(CC4(TAU))/dTAU (DUM=0 IN CC4)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON MODC,MM,A,B,DA,CX,P,BETA,BETR,OMEGA,R,DMU1,DMU2,MODE,
     &       MCON,T,Z,STOPER,LEVEL,ICHEB
      CC0=0.0
      DG=DBEXP(-DMU1/BETR*TAU)
      G1=DBEXP(P*(BETR*Z-TAU)*(TAU-BETR*Z)/(4.*BETR*TAU))
      IF(MODC.NE.3.AND.MODC.NE.4) GO TO 5
      G2=DSQRT(P/(BETR*TAU))
      CC0=DG*(0.56419*G2*G1-P/(2.*BETR)*EXF(P*Z,G2/2.*(BETR*Z+TAU)))
      GO TO 10
    5 CC0=DG*(Z/TAU)*DSQRT(P*BETR/(12.5664*TAU))*G1
   10 RETURN
      END
C
C     ------------------------------------------------------------------
      REAL*8 FUNCTION CC1(TAU)
C
C     PURPOSE: CALCULATE EQUILIBRIUM SOLUTIONS FOR STEP INPUT
C          SEE VAN GENUCHTEN & ALVES (1982) A(X,T) OF A1 & A2
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON MODC,MM,A,B,DA,CX,P,BETA,BETR,OMEGA,R,DMU1,DMU2,MODE,
     &       MCON,T,Z,STOPER,LEVEL,ICHEB
       PI=3.14159265
      BA=P/(4.0*BETR*TAU)
      RBA=DSQRT(BA)
      RBB=DSQRT(P*TAU/PI/BETR)
      G1=EXF(0.D00,RBA*(BETR*Z-TAU))
      G2=EXF(P*Z,RBA*(BETR*Z+TAU))
      IF(MODC.NE.3.AND.MODC.NE.4) GO TO 5
      G3=DBEXP(-BA*(BETR*Z-TAU)*(BETR*Z-TAU))
      CC1=G1/2.0+RBB*G3-(1.0+P*Z+P*TAU/BETR)*G2/2.0
      GO TO 10
    5 CC1=(G1+G2)/2.0
   10 RETURN
      END
C     ------------------------------------------------------------------
      REAL*8 FUNCTION CC2(TAU,Z1)
C
C     PURPOSE: CALCULATE ARGUMENT FOR STEPWISE INITIAL 
C                            & PRODUCTION PROFILES 
C          SEE VAN GENUCHTEN & ALVES (1982) A(X,T) OF A5 & A6
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON MODC,MM,A,B,DA,CX,P,BETA,BETR,OMEGA,R,DMU1,DMU2,MODE,
     &       MCON,T,Z,STOPER,LEVEL,ICHEB
C
C   MODC=1; flux conc. or area-averaged flux conc.
C   MODC=2; field-scale flux conc.
C   MODC=3; third-type resident conc.
C   MODC=4; third-type total resident concentration.
C   MODC=5; first-type resident conc.
C   MODC=6; first-type total resident conc.
c
       PI=3.14159265
      BA=P/(4.0*BETR*TAU)
      RBA=DSQRT(BA)
      RBB=DSQRT(P*TAU/PI/BETR)
      G1=EXF(0.D00,RBA*(BETR*(Z-Z1)-TAU))
      G2=EXF(P*Z,RBA*(BETR*(Z+Z1)+TAU))
      IF(MODC.EQ.3.OR.MODC.EQ.4) GO TO 5
        IF(MODC.EQ.5.OR.MODC.EQ.6) THEN
          CC2=(G1+G2)/2.0
         ELSE
           RBC=DSQRT(BETR/PI/P/TAU)/2.0
           A1=RBC*DBEXP(-BA*(BETR*(Z-Z1)-TAU)*(BETR*(Z-Z1)-TAU))
           G3=DBEXP(P*Z-BA*(BETR*(Z+Z1)+TAU)*(BETR*(Z+Z1)+TAU))
           A2=RBC*G3
           CC2=(G1+G2)/2.0+A1-A2
        END IF
        GO TO 10
    5 G3=DBEXP(P*Z-BA*(BETR*(Z+Z1)+TAU)*(BETR*(Z+Z1)+TAU))
        CC2=G1/2.0+RBB*G3-(1.0+P*(Z+Z1)+P*TAU/BETR)*G2/2.0
C
   10 RETURN
      END
C     ------------------------------------------------------------------
      REAL*8 FUNCTION CC3(TAU,Z1)
C
C     PURPOSE: CALCULATE ARGUMENT FOR EXPOTENTIAL 
C                   INITIAL AND PRODUCTION PROFILES
C          SEE VAN GENUCHTEN & ALVES (1982) B(X,T) OF A7 & A8
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON MODC,MM,A,B,DA,CX,P,BETA,BETR,OMEGA,R,DMU1,DMU2,MODE,
     &       MCON,T,Z,STOPER,LEVEL,ICHEB
C
C   MODC=1; flux conc. or area-averaged flux conc.
C   MODC=2; field-scale flux conc.
C   MODC=3; third-type resident conc.
C   MODC=4; third-type total resident concentration.
C   MODC=5; first-type resident conc.
C   MODC=6; first-type total resident conc.
C
      PI=3.14159265
      BA=P/(4.0*BETR*TAU)
      RBA=DSQRT(BA)
      RBB=DSQRT(P*TAU/PI/BETR)
      A1=DBEXP(Z1*Z1*TAU/BETR/P+Z1*TAU/BETR-Z1*Z)
      G1=EXF(0.D00,RBA*(BETR*Z-(1.0+2.0*Z1/P)*TAU))
      G2=EXF(P*Z+2.0*Z1*Z,RBA*(BETR*Z+(1.0+2.0*Z1/P)*TAU))
      IF(MODC.EQ.3.OR.MODC.EQ.4) GO TO 5
      IF(MODC.EQ.1.OR.MODC.EQ.2) THEN
        CC3=(1.0+Z1/P)*A1*(2.0-G1-G2)/2.0
      ELSE IF(MODC.EQ.5.OR.MODC.EQ.6) THEN
        CC3=A1*(2.0-G1-G2)/2.0
      ELSE
      END IF
      GO TO 10
    5  G3=EXF(P*Z,RBA*(BETR*Z+TAU))
       CC3=A1*(1.0-G1/2.0+(1+P/Z1)*G2/2.0)-P/Z1/2.0*G3
   10  RETURN
      END
C -------------------------------------------------
      REAL*8 FUNCTION CC4(TAU,DUM)
C
C     PURPOSE: ARGUMENT FOR STEP INPUT OR
C                CONSTANT PRODUCTION  
C    CC4=INT FROM O TO T CC0(TAU) d TAU (DUM not 0)
C          SEE VAN GENUCHTEN & ALVES (1982) B(X,T) OF C5 & C6
C                                         & A(X,T) OF B5 & B6
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON MODC,MM,A,B,DA,CX,P,BETA,BETR,OMEGA,R,DMU1,DMU2,MODE,
     &       MCON,T,Z,STOPER,LEVEL,ICHEB
       DCC4=0.0
       PI=3.14159265
      U=DSQRT(1.0+4.0*(DMU1+DUM)/P)
      BA=P/(4.0*BETR*TAU)
      RBA=DSQRT(BA)
      RBB=DSQRT(P*TAU/PI/BETR)
      G1=EXF(P*(1.0-U)*Z/2.0,RBA*(BETR*Z-U*TAU))
      G2=EXF(P*(1.0+U)*Z/2.0,RBA*(BETR*Z+U*TAU))
      IF (MODC.NE.3.AND.MODC.NE.4) GO TO 5
       IF(ABS(DMU1+DUM).LT.1.D-07) THEN
           G3=DBEXP(-BA*(BETR*Z-TAU)*(BETR*Z-TAU))
           CC4=G1/2.0+RBB*G3-(1.0+P*Z+P*TAU/BETR)*G2/2.0
       ELSE
           G3=EXF(P*Z-(DMU1+DUM)*TAU/BETR,RBA*(BETR*Z+TAU))
           CC4=1/(1+U)*G1+1/(1-U)*G2+P/2.0/(DMU1+DUM)*G3
       END IF
      GO TO 10
    5 CC4=(G1+G2)/2.0
   10 RETURN
      END
C     ------------------------------------------------------------------
      REAL*8 FUNCTION CC5(TAU,Z1)
C
C     PURPOSE: CALCULATE ARGUMENT FOR DELTA INITIAL CONDITION
C                  OR GENERAL IVP & PVP
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON MODC,MM,A,B,DA,CX,P,BETA,BETR,OMEGA,R,DMU1,DMU2,MODE,
     &       MCON,T,Z,STOPER,LEVEL,ICHEB
C
C   MODC=1; flux conc. or area-averaged flux conc.
C   MODC=2; field-scale flux conc.
C   MODC=3; third-type resident conc.
C   MODC=4; third-type total resident concentration.
C   MODC=5; first-type resident conc.
C   MODC=6; first-type total resident conc.
c
       PI=3.14159265
      BA=BETR*P/(4.0*PI*TAU)
      RBA=DSQRT(BA)
      BB=P/4.0/BETR/TAU
      G1=DBEXP(-BB*(BETR*(Z1-Z)+TAU)*(BETR*(Z1-Z)+TAU))
      G2=DBEXP(P*Z-BB*(BETR*(Z1+Z)+TAU)*(BETR*(Z1+Z)+TAU))
      IF(MODC.EQ.3.OR.MODC.EQ.4) GO TO 5
        IF(MODC.EQ.5.OR.MODC.EQ.6) THEN
          CC5=(G1-G2)*RBA
         ELSE
           BC=1.0-(BETR*(Z1-Z)*TAU)/2.0/TAU
           BD=1.0-(BETR*(Z1+Z)*TAU)/2.0/TAU
           CC5=(BC*G1-BD*G2)*RBA
        END IF
        GO TO 10
    5 RBB=DSQRT(BB)
      G3=EXF(P*Z,RBB*(BETR*(Z1+Z)+TAU))
        CC5=(G1+G2)*RBA-P/2.0*G3
C
   10 RETURN
      END
C     ------------------------------------------------------------------
      REAL*8 FUNCTION PROD0(TAU)
C
C     PURPOSE: CALCULATE ANALYTICAL SOLUTIONS OF CONSTANT PRODUCTION TERM
C         FOR EQUILIBRIUIM CDE IN CASE OF DMU1=0
C          SEE VAN GENUCHTEN & ALVES (1982) B(X,T) OF B5 & B6
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON MODC,MM,A,B,DA,CX,P,BETA,BETR,OMEGA,R,DMU1,DMU2,MODE,
     &       MCON,T,Z,STOPER,LEVEL,ICHEB
      COMMON/PROD/MODP,NPRO1,NPRO2,GAMMA1(50),ZPRO1(50),GAMMA2(50),
     &   ZPRO2(50)
       PI=3.14159265
      BA=P/(4.0*BETR*TAU)
      RBA=DSQRT(BA)
      BR1=(BETR*Z-TAU)/2.0
      BR2=(BETR*Z+TAU)/2.0
      G1=EXF(0.D00,RBA*(BETR*Z-TAU))
      G2=EXF(P*Z,RBA*(BETR*Z+TAU))
C
      IF(MODC.NE.3.AND.MODC.NE.4) GO TO 5
        G3=DBEXP(-BA*(BETR*Z-TAU)*(BETR*Z-TAU))
        RBB=DSQRT(P*TAU/PI/BETR/4.0)
        BRP=BETR/P
        PROD0=GAMMA1(1)/BETR*(TAU+(BR1+BRP/2.0)*G1-RBB*2.0*(BR2+BRP)*G3
     &     +(TAU/2.0-BRP/2.0+BR2*BR2/BRP)*G2)
         GO TO 10
    5  PROD0=GAMMA1(1)/BETR*(TAU+BR1*G1-BR2*G2)
C
   10  RETURN
      END
C
