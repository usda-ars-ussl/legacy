      SUBROUTINE CONST1(MAXTRY,STOPCR,GA,GD,DERL,STSQ,MM,ICHEB,NU1,NU2,
     &                  OMMAX,MIT)
C
C  PARAMETERS FOR CONTROL NUMERICAL EVALUATIONS
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION DERL(30)
C      
C MAXTRY    Maximum number of trial allowed within an iteration to find
C  new parameter values that decreases SSQ ( It is suggested that MAXTRY
C  be in the range 10 to 50; smaller values may reduce the run time but 
C  miss a convergent solution.).
        MAXTRY=50
C   ------------------------------------------------------------
C  STOPCR; Iteration criterion.  The curve-fitting process stops when the 
C   relative change in the ratio of all coefficients becomes less than STOPCR.
        STOPCR=0.0005
C Parameters for the Marquardt inversion method
C     GA:   GA*GD = INITIAL VALUE FOR FUDGE FACTOR 
C              (lambda OF EQ.(5.24), SEE ALSO Marquardt, 1963)
C     GD:   TRY AND ERROR FACTOR FOR GA
C         If a step succeeds, GA decreases to GA/GD.  If it fails,
C             GA increases to GA*GD.
C     DERL: INCREMENT TO EVALUATE VECTOR DERIVATIVES IN TERMS OF 
C             MODEL PARAMETERS
C           DERL FOR B(I)  = DERL(NU1+I)   SEE TABLE 7.4
C   CXTFIT version 1  
C      GA=0.02
C      GD=10.0
C      DERL=1.D-2
C       (GA=0.05)
C       (GD=20.0)
       DO 5 I=1,30
    5  DERL(I)=0.0
C
        GA=0.01
        GD=10.0
        DO 10 I=NU1,NU2
        DERL(I)=1.D-2
   10   CONTINUE
C STSQ: Stop criteria for the iteration based on the improvement 
C  of SSQ.  If the relative improvement in SSQ is less than STSQ 
C for three times, the iteration will stop.
        STSQ=1.E-6
C   -----------------------------------------------------------
C    MM; Initial number of integration points for Gauss Chebychev
        MM=100
C  IVHEB; Integration mode  for Gauss Chebychev.  If ICHEB=1, number of 
C   integration will increase until the result satisfies the error criteria. 
C  
        ICHEB=0
C  FOR INVERSE PROBLEM, USE A FIXED NUMBER OF INTEGRATION POINTS FOR 
C    THE NONEQUILIBRIUM CDE (SEE SUBROUTINE CHEBYCON)
      IF(MIT.GE.1) ICHEB=0
C
C  MAXIMUM CONSTRAINT FOR OMEGA
C   OMMAX = 100 IS RECOMMENDED WHEN L IS EQUAL TO THE OBSERVATION SCALE
        OMMAX=100.
C
        RETURN
        END
C
C ------------------------------------------------------------------
      SUBROUTINE CONST2(MODE,STOPER,LEVEL,MSTOCH,P,ICHEB)
C  
C  PARAMETER FOR LOG TRANSFORMED ROMBERG
      IMPLICIT REAL*8(A-H,O-Z)
C  STOPER;  Stop criteria for log-transformed Romberg integration  
C   (inial setting is 5.E-9 and the value will be adjusted internally 
C    depending on the Peclet number)
C  LEVEL; The maximum order for log-transformed Romberg integration
C  (the value will be adjusted internally depending on the Peclet number). 
C   
        STOPER=1.E-9
        LEVEL=12
      IF(MOD(MODE,2).EQ.1) THEN
       IF(P.LE.1.E+03) THEN
        STOPER=5.E-7
        LEVEL=11
       ELSE IF(P.GE.1.E+05) THEN
        STOPER=1.E-10
        LEVEL=15
       ELSE IF(P.GE.1.E+03.AND.P.LE.5.E+05) THEN
        STOPER=5.E-8
        LEVEL=12
      END IF
      END IF
      IF(MOD(MODE,2).EQ.0.AND.MSTOCH.EQ.4) THEN
      STOPER=0.00005
      LEVEL=9
      END IF
C STOCHASTIC NONEQUILIBRIUM MODEL 
      IF(MOD(MODE,2).EQ.0.AND.MSTOCH.EQ.4) ICHEB=0
C
      RETURN
      END
C
C--------------------------------------------------------------
      REAL*8 FUNCTION CINPUT(TAU)
C
C     PURPOSE: ARBITRARY FUNCTION DEFINED BY USER 
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/MODAT/INDEX(15),INVERSE,NREDU,NVAR,ZL,MIT,MDEG,MASS,
     &         DUMTP(50),DUMGA1(50),MNEQ,ISKIP,PHIM,PHIIM
      COMMON/STOCH/MODD,MODK,MODS,MSTOCH,CORR,MCORR,MSD,SDLNK,
     & V,SDLNV,VMAX,VMIN,D,AVEY,SDLNY,YMAX,YMIN,ALPHA,SDLND
     & ,DK,RHOTH,MD56,MK34,MAL8
C
C    
       IF(NREDU.LE.1) THEN 
         T=TAU/V*ZL
       ELSE
         T=TAU
       END IF
C
C  USE TIME DIMENSION FOR CINPUT ACCORDING TO NREDU
C   NREDU=0; REAL t & x 
C   NREDU=1; REAL t & x 
C   NREDU=2; DIMENSIONLESS T & Z
C   NREDU=3; DIMENSIONLESS T & REAL x
C
C       IF (T.LE.2) THEN
C        CINPUT=1.0
C       ELSE
C        CINPUT=0 
C       ENDIF
       IF (T.LE.9.1) THEN
        CINPUT=(-560.39*T) + 16528
       ELSE
        CINPUT=0 
       ENDIF
C  
       RETURN
       END
