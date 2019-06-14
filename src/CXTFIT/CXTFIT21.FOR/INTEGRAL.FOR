C       ------------------------------------------------------------------
      SUBROUTINE ROMB(FUNC,AREA,XMIN,XMAX)
C
C     PURPOSE: PERFORM ROMBERG INTEGRATION ON LOG-TRANSFORMED INTERVAL
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON MODC,MM,AA,BB,DA,CX,PEC,BETA,BETR,OMEGA,RE,DMU1,DMU2,
     &   MODE,MCON,TT,ZZ,STOPER,LEVEL,ICHEB
      DIMENSION R(15,15)
      AREA=0.0
      SSTOP=STOPER
      IF(XMIN.GE.XMAX)GOTO 30
      XLN1=DLOG(XMIN)
      XLN2=DLOG(XMAX)
      DX1=XLN2-XLN1
      SUM=0.5*XMIN*FUNC(XMIN)
      AREA=SUM*DX1
      R(1,1)=AREA
      DX2=DX1/2.
      INC=1
      DO 20 I=2,LEVEL
      XLN=XLN1+DX2
      DO 4 M=1,INC
      XE=DBEXP(XLN)
      SUM=SUM+XE*FUNC(XE)
    4 XLN=XLN+DX1
      AREA=SUM*DX2
      R(1,I)=AREA
      DO 17 J=2,I
      K=I+1-J
   17 R(J,K)=(4.**(J-1)*R(J-1,K+1)-R(J-1,K))/(4.**(J-1)-1)
      IF(R(I,1).GT.0.0) GO TO 18
C      IF(I.GT.4) GOTO 30
      GO TO 19
   18 ERROR=DABS((R(I,1)-R(I-1,2))/R(I,1))
C      WRITE(*,*) 'I,ERROR AREA',I,ERROR,AREA
      IF((ERROR.LT.STOPER).AND.(I.GT.4)) THEN
       GOTO 30
      END IF
      IF((ERROR.LT.STOPER*100.0).AND.(DABS(AREA).LT.1.D-10)) GOTO 30
   19 DX1=DX1/2.
      DX2=DX2/2.
   20 INC=INC*2
   30 STOPER=SSTOP
      RETURN
      END
C
C       ------------------------------------------------------------------
      SUBROUTINE ROMB2(FUNC,AREA,XMIN,XMAX)
C
C     PURPOSE: PERFORM ROMBERG INTEGRATION ON LOG-TRANSFORMED INTERVAL
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON MODC,MM,AA,BB,DA,CX,PEC,BETA,BETR,OMEGA,RE,DMU1,DMU2,
     &   MODE,MCON,TT,ZZ,STOPER,LEVEL,ICHEB
      DIMENSION R(15,15)
      AREA=0.0
      SSTOP=STOPER
      IF(XMIN.GE.XMAX)GOTO 30
      XLN1=DLOG(XMIN)
      XLN2=DLOG(XMAX)
      DX1=XLN2-XLN1
      SUM=0.5*XMIN*FUNC(XMIN)
      AREA=SUM*DX1
      R(1,1)=AREA
      DX2=DX1/2.
      INC=1
      DO 20 I=2,LEVEL
      XLN=XLN1+DX2
      DO 4 M=1,INC
      XE=DBEXP(XLN)
      SUM=SUM+XE*FUNC(XE)
    4 XLN=XLN+DX1
      AREA=SUM*DX2
      R(1,I)=AREA
      DO 17 J=2,I
      K=I+1-J
   17 R(J,K)=(4.**(J-1)*R(J-1,K+1)-R(J-1,K))/(4.**(J-1)-1)
      IF(R(I,1).GT.0.0) GO TO 18
      IF(I.GT.4) GOTO 30
      GO TO 19
   18 ERROR=DABS((R(I,1)-R(I-1,2))/R(I,1))
      IF((ERROR.LT.STOPER).AND.(I.GT.4)) THEN
       GOTO 30
      END IF
      IF((ERROR.LT.STOPER*100.0).AND.(DABS(AREA).LT.1.D-10)) GOTO 30
   19 DX1=DX1/2.
      DX2=DX2/2.
   20 INC=INC*2
   30 STOPER=SSTOP
      RETURN
      END
C
C   ------------------------------------------------
      SUBROUTINE CHEBY(FUNC,AREA,A,B)
C
C       PURPOSE: PERFORM INTEGRATION OF F(X) BETWEEN A AMD B
C             USING M-POINT GAUSS-CHEBYSHEV QUANDRATURE FORMULA
C             THIS SUBROUTINE IS TO EVALUATE INTEGRALS FOR NONEQUILIBRIUM CDE
C              ICEBY=0 USE FIXED NO. OF INTEHRATION POINTS
C              ICHEB=1 INCREASE NO. OF POINTS USING SOTOP CRITERIA
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON MODC,M,AA,BB,DA,CX,PEC,BETA,BETR,OMEGA,RE,DMU1,DMU2,
     &   MODE,MCON,TT,ZZ,STOPER,LEVEL,ICHEB
C
       AREA=0.0
       SUM=0.0
       PI=3.14159265
       IF(ICHEB.NE.1) THEN
         DO 5 I=1,M
           Z1=DCOS(DFLOAT(2*(I-1)+1)*PI/DFLOAT(2*M))
           X1=(Z1*(B-A)+B+A)/2.0
    5    SUM=SUM+FUNC(X1)*DSQRT(1.0-Z1*Z1)
           AREA=(B-A)*PI*SUM/DFLOAT(2*M)
           RETURN
       END IF
       MM=M
       STOPCH=1.E-3
       AREA1=0.0
       AREA2=0.0
       DO 20 J=1,10
       SUM=0.0
       DO 10 I=1,M
        Z1=DCOS(DFLOAT(2*(I-1)+1)*PI/DFLOAT(2*M))
        X1=(Z1*(B-A)+B+A)/2.0
   10 SUM=SUM+FUNC(X1)*DSQRT(1.0-Z1*Z1)
      AREA2=(B-A)*PI*SUM/DFLOAT(2*M)
      IF(ABS(AREA2).LT.1E-10) RETURN 
      IF(J.EQ.1) THEN
         AREA1=AREA2
         GOTO 15
      END IF
      ERROR=ABS(AREA2-AREA1)/AREA2
C      WRITE(*,*) ' J M ERROR',J,M,ERROR
C      WRITE(7,*) ' J M ERROR',J,M,ERROR
      IF(ERROR.LT.STOPCH.OR.M.GE.6400) THEN
        AREA=AREA2
        GOTO 30
      ELSE
        AREA1=AREA2
      END IF
   15 M=M*2
   20 CONTINUE
   30 M=MM
      RETURN
      END
C   ------------------------------------------------
      SUBROUTINE CHEBY2(FUNC,AREA,A,B)
C
C       PURPOSE: PERFORM INTEGRATION OF F(X) BETWEEN A AMD B
C             USING M-POINT GAUSS-CHEBYSHEV QUANDRATURE FORMULA
C             THIS SUBROUTINE IS TO EVALUATE INTEGRALS FOR NONEQUILIBRIUM CDE
C              ICEBY=0 USE FIXED NO. OF INTEHRATION POINTS
C              ICHEB=1 INCREASE NO. OF POINTS USING SOTOP CRITERIA
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON MODC,M,AA,BB,DA,CX,PEC,BETA,BETR,OMEGA,RE,DMU1,DMU2,
     &   MODE,MCON,TT,ZZ,STOPER,LEVEL,ICHEB
C
       AREA=0.0
       SUM=0.0
       PI=3.14159265
       IF(ICHEB.NE.1) THEN
         DO 5 I=1,M
           Z1=DCOS(DFLOAT(2*(I-1)+1)*PI/DFLOAT(2*M))
           X1=(Z1*(B-A)+B+A)/2.0
    5    SUM=SUM+FUNC(X1)*DSQRT(1.0-Z1*Z1)
           AREA=(B-A)*PI*SUM/DFLOAT(2*M)
           RETURN
       END IF
       MM=M
       STOPCH=1.E-3
       AREA1=0.0
       AREA2=0.0
       DO 20 J=1,10
       SUM=0.0
       DO 10 I=1,M
        Z1=DCOS(DFLOAT(2*(I-1)+1)*PI/DFLOAT(2*M))
        X1=(Z1*(B-A)+B+A)/2.0
   10 SUM=SUM+FUNC(X1)*DSQRT(1.0-Z1*Z1)
      AREA2=(B-A)*PI*SUM/DFLOAT(2*M)
      IF(ABS(AREA2).LT.1E-10) RETURN 
      IF(J.EQ.1) THEN
         AREA1=AREA2
         GOTO 15
      END IF
      ERROR=ABS(AREA2-AREA1)/AREA2
C      WRITE(*,*) ' J M ERROR',J,M,ERROR
C      WRITE(7,*) ' J M ERROR',J,M,ERROR
      IF(ERROR.LT.STOPCH.OR.M.GE.6400) THEN
        AREA=AREA2
        GOTO 30
      ELSE
        AREA1=AREA2
      END IF
   15 M=M*2
   20 CONTINUE
   30 M=MM
      RETURN
      END
C   ------------------------------------------------
      SUBROUTINE CHEBYCON(FUNC,AREA,A,B,MC)
C
C       PURPOSE: PERFORM INTEGRATION OF F(X,MC) BETWEEN A AMD B
C             USING M-POINT GAUSS-CHEBYSHEV QUANDRATURE FORMULA
C             THIS SUBROUTINE IS TO EVALUATE INTEGRALS FOR NONEQUILIBRIUM CDE
C              ICEBY=0 USE FIXED NO. OF INTEHRATION POINTS
C              ICHEB=1 INCREASE NO. OF POINTS USING SOTOP CRITERIA
C              MC: CONCENTRATION MODE; =1 PHASE 1; =2 PHASE 2 
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON MODC,M,AA,BB,DA,CX,PEC,BETA,BETR,OMEGA,RE,DMU1,DMU2,
     &   MODE,MCON,TT,ZZ,STOPER,LEVEL,ICHEB
C
       AREA=0.0
       SUM=0.0
       PI=3.14159265
       IF(ICHEB.NE.1) THEN
         DO 5 I=1,M
           Z1=DCOS(DFLOAT(2*(I-1)+1)*PI/DFLOAT(2*M))
           X1=(Z1*(B-A)+B+A)/2.0
    5    SUM=SUM+FUNC(X1,MC)*DSQRT(1.0-Z1*Z1)
           AREA=(B-A)*PI*SUM/DFLOAT(2*M)
           RETURN
       END IF
       MM=M
       STOPCH=1.E-3
       AREA1=0.0
       AREA2=0.0
       DO 20 J=1,10
       SUM=0.0
       DO 10 I=1,M
        Z1=DCOS(DFLOAT(2*(I-1)+1)*PI/DFLOAT(2*M))
        X1=(Z1*(B-A)+B+A)/2.0
   10 SUM=SUM+FUNC(X1,MC)*DSQRT(1.0-Z1*Z1)
      AREA2=(B-A)*PI*SUM/DFLOAT(2*M)
      IF(ABS(AREA2).LT.1E-10) RETURN 
      IF(J.EQ.1) THEN
         AREA1=AREA2
         GOTO 15
      END IF
      ERROR=ABS(AREA2-AREA1)/AREA2
c        WRITE(*,*) ' J M ERROR',J,M,ERROR
c        WRITE(7,*) ' J M ERROR',J,M,ERROR
      IF(ERROR.LT.STOPCH.OR.M.GE.6400) THEN
        AREA=AREA2
        GOTO 30
      ELSE
        AREA1=AREA2
      END IF
   15 M=M*2
   20 CONTINUE
   30 M=MM
      RETURN
      END
C   ------------------------------------------------
      SUBROUTINE CHEBYLOG(FUNC,AREA,AA,BB)
C
C       PURPOSE: PERFORM INTEGRATION OF F(X) BETWEEN A AMD B
C             USING M-POINT GAUSS-CHEBYSHEV QUANDRATURE FORMULA
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON MODC,M,AAA,BBB,DA,CX,PEC,BETA,BETR,OMEGA,RE,DMU1,DMU2,
     &   MODE,MCON,TT,ZZ,STOPER,LEVEL,ICHEB
C
       PI=3.14159265
       AREA=0.0
       SUM=0.0
       A=DLOG(AA)
       B=DLOG(BB)
      DO 10 I=1,M
        Z1=DCOS(DFLOAT(2*(I-1)+1)*PI/DFLOAT(2*M))
        X1=(Z1*(B-A)+B+A)/2.0
        DX1=DBEXP(X1)
   10 SUM=SUM+DX1*FUNC(DX1)*DSQRT(1.0-Z1*Z1)
      AREA=(B-A)*PI*SUM/DFLOAT(2*M)
      RETURN
      END
C   ------------------------------------------------
      SUBROUTINE CHEBYLOG2(FUNC,AREA,AA,BB)
C
C       PURPOSE: PERFORM INTEGRATION OF F(X) BETWEEN A AMD B
C             USING M-POINT GAUSS-CHEBYSHEV QUANDRATURE FORMULA
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON MODC,M,AAA,BBB,DA,CX,PEC,BETA,BETR,OMEGA,RE,DMU1,DMU2,
     &   MODE,MCON,TT,ZZ,STOPER,LEVEL,ICHEB
C
       PI=3.14159265
       AREA=0.0
       SUM=0.0
       A=DLOG(AA)
       B=DLOG(BB)
      DO 10 I=1,M
        Z1=DCOS(DFLOAT(2*(I-1)+1)*PI/DFLOAT(2*M))
        X1=(Z1*(B-A)+B+A)/2.0
        DX1=DBEXP(X1)
   10 SUM=SUM+DX1*FUNC(DX1)*DSQRT(1.0-Z1*Z1)
      AREA=(B-A)*PI*SUM/DFLOAT(2*M)
      RETURN
      END
