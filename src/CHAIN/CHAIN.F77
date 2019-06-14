C
C     ******************************************************************
C     *                                                                *
C     *       MULTI-ION SOLUTE TRANSPORT                   CHAIN       *
C     *                                                                *
C     *       FOUR SPECIES                                             *
C     *       PREFERENTIAL RELEASE                                     *
C     *                                                                *
C     *                                                                *
C     ******************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON R(4),DONE(4),RLAM(4),D,V,X,T,KSURF
      DIMENSION ZERO(4),GAMMA(4),A(10),B(10)
      CHARACTER*60 TITLE
C
C     ----- OPEN I/O FILES -----
      OPEN(6,FILE='CHAIN.OUT')
      CLOSE(6,STATUS='DELETE')
      OPEN(5,FILE='CHAIN.IN')
      OPEN(6,FILE='CHAIN.OUT',STATUS='NEW')
C
C     ----- READ NUMBER OF CASES CONSIDERED -----
      READ(5,1000) NC,KP
      IF(KP.NE.6) KP=1
      DO 14 NCASE=1,NC
      READ(5,1001) TITLE
C
C     ----- READ INPUT PARAMETERS -----
      READ(5,1000) NS,NR,NB,KSURF,KDB,KPR,V,D,WC,TP
      READ(5,1002) (R(I),I=1,4),(DONE(I),I=1,4)
      READ(5,1002) (ZERO(I),I=1,4),(GAMMA(I),I=1,4)
      READ(5,1002) XI,DX,XM,TI,DT,TM
C
C     ----------
      DO 1 I=1,4
      ZERO(I)=NB*ZERO(I)
    1 RLAM(I)=NB*DONE(I)+GAMMA(I)
      IF(NB.EQ.0) GO TO 3
C
C     ----- CALCULATE BATEMAN EQUATIONS IF NB=1 -----
      Q=V*WC
      A(1)=ZERO(1)
      A(2)=DONE(1)*A(1)/(RLAM(2)-RLAM(1))
      A(3)=ZERO(2)-A(2)
      A(4)=DONE(2)*A(2)/(RLAM(3)-RLAM(1))
      A(5)=DONE(2)*A(3)/(RLAM(3)-RLAM(2))
      A(6)=ZERO(3)-A(4)-A(5)
      A(7)=DONE(3)*A(4)/(RLAM(4)-RLAM(1))
      A(8)=DONE(3)*A(5)/(RLAM(4)-RLAM(2))
      A(9)=DONE(3)*A(6)/(RLAM(4)-RLAM(3))
      A(10)=ZERO(4)-A(7)-A(8)-A(9)
      K=0
      DO 2 I=1,4
      DO 2 J=1,I
      K=K+1
    2 B(K)=A(K)*GAMMA(I)/Q
      GO TO 4
    3 READ(5,1002) (B(I),I=1,10)
C
C     ----- WRITE INPUT PARAMETERS -----
    4 IF(NS.LT.0) GO TO 14
      WRITE(KP,1003)
      IF(KSURF.EQ.0) WRITE(KP,1004)
      IF(KSURF.NE.0) WRITE(KP,1005)
      WRITE(KP,1006) TITLE
      WRITE(KP,1007) NS,V,D,NR,WC,TP,NB
      WRITE(KP,1008) (I,ZERO(I),R(I),DONE(I),GAMMA(I),RLAM(I),I=1,NS)
      WRITE(KP,1009) (I,B(I),I=1,10)
C
C     ----- MULTIPLY DONE WITH R IF NR=1 -----
      IF(NR.EQ.0) GO TO 6
      DO 5 I=1,4
    5 DONE(I)=R(I)*DONE(I)
C
C     ----- CALCULATE VARIOUS CONSTANTS -----
    6 D12=DONE(1)-DONE(2)
      D13=DONE(1)-DONE(3)
      D23=DONE(2)-DONE(3)
      D14=DONE(1)-DONE(4)
      D24=DONE(2)-DONE(4)
      D34=DONE(3)-DONE(4)
      R12=R(1)-R(2)
      R13=R(1)-R(3)
      R23=R(2)-R(3)
      R14=R(1)-R(4)
      R24=R(2)-R(4)
      R34=R(3)-R(4)
      B123=-DONE(1)*R23+DONE(2)*R13-DONE(3)*R12
      B124=-DONE(1)*R24+DONE(2)*R14-DONE(4)*R12
      B134=-DONE(1)*R34+DONE(3)*R14-DONE(4)*R13
      B234=-DONE(2)*R34+DONE(3)*R24-DONE(4)*R23
      DD=DONE(1)*DONE(2)
      DP=DONE(2)*DONE(3)
      DDD=DD*DONE(3)
      C21=DONE(1)*B(1)/(D12-R12*RLAM(1))
      C31=DONE(2)*B(2)/(D23-R23*RLAM(1))
      C32=DONE(2)*B(3)/(D23-R23*RLAM(2))
      C331=DD*B(1)*R12/(B123*(R12*RLAM(1)-D12))
      C332=DD*B(1)*R13/(B123*(R13*RLAM(1)-D13))
      C333=DD*B(1)*R23/(B123*(R23*RLAM(1)-D23))
      C41=DONE(3)*B(4)/(D34-R34*RLAM(1))
      C42=DONE(3)*B(5)/(D34-R34*RLAM(2))
      C43=DONE(3)*B(6)/(D34-R34*RLAM(3))
      C441=DP*B(2)*R34/(B234*(R34*RLAM(1)-D34))
      C442=DP*B(2)*R23/(B234*(R23*RLAM(1)-D23))
      C443=DP*B(2)*R24/(B234*(R24*RLAM(1)-D24))
      C451=DP*B(3)*R34/(B234*(R34*RLAM(2)-D34))
      C452=DP*B(3)*R23/(B234*(R23*RLAM(2)-D23))
      C453=DP*B(3)*R24/(B234*(R24*RLAM(2)-D24))
      C461=DDD*B(1)*R14**2/(B134*B124*(R14*RLAM(1)-D14))
      C462=DDD*B(1)*R24**2/(B124*B234*(R24*RLAM(1)-D24))
      C463=DDD*B(1)*R34**2/(B134*B234*(R34*RLAM(1)-D34))
      C464=DDD*B(1)*R13**2/(B123*B134*(R13*RLAM(1)-D13))
      C465=DDD*B(1)*R23**2/(B123*B234*(R23*RLAM(1)-D23))
      C466=DDD*B(1)*R12**2/(B123*B124*(R12*RLAM(1)-D12))
C
C     ----- PRINT COEFFICIENTS IF KDB.GT.0 -----
      IF(KDB.GT.0) WRITE(KP,1010) C21,C31,C32,C331,C332,C333,C41,C42,C43
     1,C441,C442,C443,C451,C452,C453,C461,C462,C463,C464,C465,C466
      F12=0.
      F13=0.
      F23=0.
      F32=0.
      F14=0.
      F24=0.
      F42=0.
      F34=0.
      F43=0.
      F44=0.
      E=0.
      P12=DABS(1000.*R12/R(1))
      P13=DABS(1000.*R13/R(1))
      P23=DABS(1000.*R23/R(2))
      P14=DABS(1000.*R14/R(1))
      P24=DABS(1000.*R24/R(2))
      P34=DABS(1000.*R34/R(3))
      IF(DX.LE.0.) DX=1.0
      IF(DT.LE.0.) DT=1.0
      IMAX=(XM+DX-XI)/DX
      JMAX=(TM+DT-TI)/DT
C
C     ----- DYNAMIC PART OF PROGRAM -----
      DO 12 JJ=1,JMAX
      IF(IMAX.GE.JJ) WRITE(KP,1011)
      TIME=TI+(JJ-1)*DT
      DO 12 II=1,IMAX
      X=XI+(II-1)*DX
      VVO=0.0
      IF(X.LE.0.) GO TO 7
      VVO=V*TIME/X
    7 DO 10 KK=1,2
      C1=0.0
      C2=0.0
      C3=0.0
      C4=0.0
      E1=1.0
      E2=1.0
      E3=1.0
      E4=1.0
      T=TIME+(1-KK)*TP
      IF(T.LE.0.) GO TO 10
      IF(KK.EQ.1) GO TO 8
      E1=EXF(-RLAM(1)*TP,E)
      E2=EXF(-RLAM(2)*TP,E)
      E3=EXF(-RLAM(3)*TP,E)
      E4=EXF(-RLAM(4)*TP,E)
    8 F110=F(1,0)*E1
C
C     ----- FIRST TRACER -----
      C1=B(1)*F110
      IF(NS.EQ.1) GO TO 9
C
C     ----- SECOND TRACER -----
      F220=F(2,0)*E2
      F210=F(2,-1)*E1
      IF(P12.GT.1.) F12=F(1,2)*E1
      C2=B(2)*F210+B(3)*F220+C21*(F210-F110-F12)
      IF(NS.EQ.2) GO TO 9
C
C     ----- THIRD TRACER -----
      F310=F(3,-1)*E1
      F320=F(3,-2)*E2
      F330=F(3,0)*E3
      IF(P13.GT.1.) F13=F(1,3)*E1
      IF(P23.GT.1.) F23=F(2,3)
      F32=F23*E2
      F23=F23*E1
      IF(KDB.EQ.2) WRITE(KP,1012) F110,F210,F220,F12,F310,F320,F330,F13,F
     123,F32
      C3=B(4)*F310+B(5)*F320+B(6)*F330+C31*(F310-F210-F23)+C32*(F320-F22
     10-F32)+C331*(F210-F110-F12)+C332*(F13-F310+F110)+C333*(F310-F210-F
     223)
      IF(NS.EQ.3) GO TO 9
C
C     ----- FOURTH TRACER -----
      F410=F(4,-1)*E1
      F420=F(4,-2)*E2
      F430=F(4,-3)*E3
      F440=F(4,0)*E4
      IF(P14.GT.1.) F14=F(1,4)*E1
      IF(P24.GT.1.) F24=F(2,4)
      IF(P34.GT.1.) F34=F(3,4)
      F42=F24*E2
      F24=F24*E1
      F43=F34*E3
      F44=F34*E2
      F34=F34*E1
      IF(KDB.EQ.2) WRITE(KP,1013) F410,F420,F430,F440,F14,F24,F42,F34,F4
     13,F44
      C4=B(7)*F410+B(8)*F420+B(9)*F430+B(10)*F440+C41*(F410-F310-F34)+C4
     12*(F420-F320-F44)+C43*(F430-F330-F43)+C441*(F410-F310-F34)+C442*(F
     2310-F210-F23)+C443*(F24-F410+F210)+C451*(F420-F320-F44)+C452*(F320
     3-F220-F32)+C453*(F42-F420+F220)+C461*(F110-F410+F14)+C462*(F410-F2
     410-F24)+C463*(F34-F410+F310)+C464*(F310-F110-F13)+C465*(F23-F310+F
     5210)+C466*(F12-F210+F110)
    9 IF(KK.EQ.2) GO TO 10
      CONC1=C1
      CONC2=C2
      CONC3=C3
      CONC4=C4
   10 CONTINUE
      CONC1=CONC1-C1
      CONC2=CONC2-C2
      CONC3=CONC3-C3
      CONC4=CONC4-C4
      IF(KPR.EQ.0) WRITE(KP,1014) X,TIME,VVO,CONC1,CONC2,CONC3,CONC4
   12 IF(KPR.EQ.1) WRITE(KP,1015) X,TIME,VVO,CONC1,CONC2,CONC3,CONC4
   14 CONTINUE
      CLOSE(5)
      CLOSE(6)
C
C     ----------
 1000 FORMAT(6I5,5F10.0)
 1001 FORMAT(A60)
 1002 FORMAT(8F10.0)
 1003 FORMAT(//5X,67(1H*)/5X,1H*,65X,1H*/5X,1H*,5X,'MULTI-ION SOLUTE TRA
     1NSPORT',24X,'CHAIN',5X,1H*/5X,1H*,65X,1H*/5X,1H*,5X,'FOUR SPECIES'
     2,48X,1H*/5X,1H*,5X,'PREFERENTIAL RELEASE',40X,1H*/5X,1H*,65X,1H*)
 1004 FORMAT(5X,1H*,5X,'FIRST-TYPE BOUNDARY CONDITION',31X,1H*)
 1005 FORMAT(5X,1H*,5X,'THIRD-TYPE BOUNDARY CONDITION',31X,1H*)
 1006 FORMAT(5X,1H*,65X,1H*/5X,1H*,A60,5X,1H*/5X,1H*,65X,1H*/5X,67(1H*))
 1007 FORMAT(//5X,'INPUT PARAMETERS'/5X,16(1H=)//5X,'NS =',I5,10X,'V ='
     1,F12.4,10X,'D =',F12.4/5X,'NR =',I5,10X,'WC =',F11.4,10X,'TP =',F
     211.4/5X,'NB =',I5)
 1008 FORMAT(///5X,'TRACER',6X,'ZERO',8X,'R',10X,'DONE',10X,'GAMMA',10X
     1,'RLAM'//(7X,I2,D15.5,F9.3,3F14.8))
 1009 FORMAT(///5X,'B(I)-COEFFICIENTS'/5X,17(1H=)//6X,4('I',6X,'B(I)'
     1,7X)/(5X,4(I2,D14.6,2X)))
 1010 FORMAT(//5X,'C21= ',D13.5,6X,'C31= ',D13.5,6X,'C32= ',D13.5/5X,
     1'C331= ',D12.5,6X,'C332= ',D12.5,6X,'C333= ',D12.5/5X,'C41= ',
     2D13.5,6X,'C42= ',D13.5,6X,'C43= ',D13.5/5X,'C441= ',D12.5,6X,
     3'C442= ',D12.5,6X,'C443= ',D12.5/5X,'C451= ',D12.5,6X,'C452= ',
     4D12.5,6X,'C453= ',D12.5/5X,'C461= ',D12.5,6X,'C462= ',D12.5,6X,
     5'C463= ',D12.5/5X,'C464= ',D12.5,6X,'C465= ',D12.5,6X,'C466= ',
     6D12.5)
 1011 FORMAT(///5X,'DISTANCE',2X,'TIME',3X,'PORE VOLUME',2X,14(1H-),'CON
     1CENTRATION',13(1H-)/7X,'(X)',6X,'(T)',7X,'(VVO)',6X,'(1)',8X,'(2)'
     2,8X,'(3)',8X,'(4)')
 1012 FORMAT(//5X,'F110= ',D12.5,4X,'F210= ',D12.5,4X,'F220= ',D12.5/5X,
     1'F12= ',D13.5,4X,'F310= ',D12.5,4X,'F320= ',D12.5/9X,'F330= ',D12.
     25,4X,'F13= ',D13.5,4X,'F23= ',D13.5/5X,'F32= ',D13.5)
 1013 FORMAT(//5X,'F410= ',D12.5,4X,'F420= ',D12.5,4X,'F430= ',D12.5/5X,
     1'F440= ',D12.5,4X,'F14= ',D13.5,4X,'F24= ',D13.5/5X,'F42= ',D13.5,
     24X,'F34= ',D13.5,4X,'F43= ',D13.5/5X,'F44= ',D13.5)
 1014 FORMAT(2X,2F9.2,1X,F10.3,4F11.5)
 1015 FORMAT(2X,2F9.2,1X,F10.3,4D11.4)
C
C     ----- CLOSE FILES -----
      CLOSE(5)
      CLOSE(6)
      STOP
      END
C
C     ------------------------------------------------------------------
      FUNCTION F(I,JJ)
C
C     PURPOSE: TO CALCULATE FIJK-TERMS
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON R(4),DONE(4),RLAM(4),D,V,X,T,KSURF
C
C     ----------
      J=IABS(JJ)
      IF(JJ.LT.0) A=RLAM(J)
      IF(JJ.EQ.0) A=RLAM(I)
      IF(JJ.GT.0) A=(DONE(I)-DONE(J))/(R(I)-R(J))
      U2=V*V+4.*D*(DONE(I)-A*R(I))
      S1=2.*DSQRT(D*R(I)*T)
      U=DSQRT(DABS(U2))
      VT=V*T
      UT=U*T
      E=0.0
C
C     ----- CASE FOR POSITIVE VALUES OF U2 -----
      IF(U2.LT.0.) GO TO 28
      AM=0.5*(V-U)*X/D-A*T
      BM1=(R(I)*X-UT)/S1
      BP1=(R(I)*X+UT)/S1
      AP=0.5*(U+V)*X/D-A*T
      IF(KSURF.EQ.0) GO TO 20
C
C     ----- CASE FOR THIRD-TYPE BOUNDARY CONDITION -----
      CP1=(R(I)*X+VT)/S1
      ALP=.3275911
      CM1=(R(I)*X-VT)/S1
      DUV=DABS(1.-U/V)
      IF((BP1.GT.4.).AND.(CP1.GT.4.)) GO TO 4
      IF((BP1.LT.2.5).AND.(CP1.LT.2.5)) GO TO 2
      IF(DUV.GT.0.1) GO TO 1
      IF(CP1.GT.3.) GO TO 4
      GO TO 2
    1 E1=2./((U/V)**2-1.)
      G1=V*X/D-DONE(I)*T/R(I)
      F1=V/(U-V)*EXF(AP,BP1)-E1*EXF(G1,CP1)
      GO TO 6
    2 D1=-CM1**2-DONE(I)*T/R(I)
      BETA=ALP*V*DSQRT(T/(D*R(I)))
      T1=1./(1.+ALP*BP1)
      T2=1./(1.+ALP*CP1)
      C1=T1*T2*(1.+ALP*CM1)
      C2=T1*(C1-BETA*T2**2)
      C3=T1*(C2-BETA*T2**3)
      C4=T1*(C3-BETA*T2**4)
      C5=T1*(C4-BETA*T2**5)
      F1=V*(.2548296*C1-.2844967*C2+1.421414*C3-1.453152*C4+1.061405*C5)
     1*EXF(D1,E)/(V+U)
      GO TO 6
    4 B=DONE(I)*T/R(I)+CM1**2
      H=1./(BP1*(R(I)*X+VT))
      WT=H*VT
      HT=H*UT
      B1=WT+2.*HT
      B2=WT*B1+2.*HT**2
      B3=WT*B2+2.*HT**3
      B4=WT*B3+2.*HT**4
      B5=WT*B4+2.*HT**5
      B6=WT*B5+2.*HT**6
      B7=WT*B6+2.*HT**7
      B8=WT*B7+2.*HT**8
      Z1=H*R(I)*X
      Z2=Z1**2
      Z3=Z1*Z2
      Z4=Z1*Z3
      Z5=Z1*Z4
      Z6=Z1*Z5
      Z7=Z1*Z6
      Z8=Z1*Z7
      Z9=Z1*Z8
      A1=Z1-WT
      A2=Z3-WT*(B2+3.*(B1*Z1+Z2))
      A3=Z5-WT*(B4+5.*(B3*Z1+Z4)+10.*(B2*Z2+B1*Z3))
      A4=Z7-WT*(B6+7.*(B5*Z1+Z6)+21.*(B4*Z2+B1*Z5)+35.*(B3*Z3+B2*Z4))
      A5=Z9-WT*(B8+9.*(B7*Z1+Z8)+36.*(B6*Z2+B1*Z7)+84.*(B5*Z3+B2*Z6)+126
     1.*(B4*Z4+B3*Z5))
      F1=.5641896*V*EXF(-B,E)*(A1-.5*A2+.75*A3-1.875*A4+6.5625*A5)/(U+V)
    6 IF(JJ.GT.0) GO TO 8
      F=V/(V+U)*EXF(AM,BM1)-F1
      RETURN
    8 S2=2.*DSQRT(D*R(J)*T)
      BM2=(R(J)*X-UT)/S2
      BP2=(R(J)*X+UT)/S2
      CP2=(R(J)*X+VT)/S2
      CM2=(R(J)*X-VT)/S2
      IF((BP2.GT.4.).AND.(CP2.GT.4.)) GO TO 14
      IF((BP2.LT.2.5).AND.(CP2.LT.2.5)) GO TO 10
      IF(DUV.GT.0.1) GO TO 9
      IF(CP2.GT.3.) GO TO 14
      GO TO 10
    9 E2=2./((U/V)**2-1.)
      G2=V*X/D-DONE(J)*T/R(J)
      AP=0.5*(V+U)*X/D-A*T
      F2=V/(U-V)*EXF(AP,BP2)-E2*EXF(G2,CP2)
      GO TO 18
   10 D2=-CM2**2-DONE(J)*T/R(J)
      BETA=ALP*V*DSQRT(T/(D*R(J)))
      T1=1./(1.+ALP*BP2)
      T2=1./(1.+ALP*CP2)
      C1=T1*T2*(1.+ALP*CM2)
      C2=T1*(C1-BETA*T2**2)
      C3=T1*(C2-BETA*T2**3)
      C4=T1*(C3-BETA*T2**4)
      C5=T1*(C4-BETA*T2**5)
      F2=V*(.2548296*C1-.2844967*C2+1.421414*C3-1.453152*C4+1.061405*C5)
     1*EXF(D2,E)/(V+U)
      GO TO 16
   14 B=DONE(J)*T/R(J)+CM2**2
      H=1./(BP2*(R(J)*X+VT))
      WT=H*VT
      HT=H*UT
      B1=WT+2.*HT
      B2=WT*B1+2.*HT**2
      B3=WT*B2+2.*HT**3
      B4=WT*B3+2.*HT**4
      B5=WT*B4+2.*HT**5
      B6=WT*B5+2.*HT**6
      B7=WT*B6+2.*HT**7
      B8=WT*B7+2.*HT**8
      Z1=H*R(J)*X
      Z2=Z1**2
      Z3=Z1*Z2
      Z4=Z1*Z3
      Z5=Z1*Z4
      Z6=Z1*Z5
      Z7=Z1*Z6
      Z8=Z1*Z7
      Z9=Z1*Z8
      A1=Z1-WT
      A2=Z3-WT*(B2+3.*(B1*Z1+Z2))
      A3=Z5-WT*(B4+5.*(B3*Z1+Z4)+10.*(B2*Z2+B1*Z3))
      A4=Z7-WT*(B6+7.*(B5*Z1+Z6)+21.*(B4*Z2+B1*Z5)+35.*(B3*Z3+B2*Z4))
      A5=Z9-WT*(B8+9.*(B7*Z1+Z8)+36.*(B6*Z2+B1*Z7)+84.*(B5*Z3+B2*Z6)+126
     1.*(B4*Z4+B3*Z5))
      F2=.5641896*V*EXF(-B,E)*(A1-.5*A2+.75*A3-1.875*A4+6.5625*A5)/(U+V)
   16 IF((BM1.GT.0.).OR.(BM2.GT.0.)) GO TO 18
      Q=-BM2
      BM2=-BM1
      BM1=Q
   18 F=V/(V+U)*(EXF(AM,BM2)-EXF(AM,BM1))+F1-F2
      RETURN
C
C     ----- CASE FOR FIRST-TYPE BOUNDARY CONDITION -----
   20 IF(JJ.GT.0) GO TO 22
      F=0.5*(EXF(AM,BM1)+EXF(AP,BP1))
      RETURN
   22 S2=2.*DSQRT(D*R(J)*T)
      BM2=(R(J)*X-UT)/S2
      BP2=(R(J)*X+UT)/S2
      IF((BM1.GT.0.).OR.(BM2.GT.0.)) GO TO 24
      Q=-BM2
      BM2=-BM1
      BM1=Q
   24 F=0.5*(EXF(AM,BM2)-EXF(AM,BM1)+EXF(AP,BP2)-EXF(AP,BP1))
      RETURN
C
C     ----- CASE FOR NEGATIVE VALUES OF U2 -----
   28 AM=0.5*V*X/D-A*T
      ARG=0.5*U*X/D
      CS=DCOS(ARG)
      SN=DSIN(ARG)
      A1=R(I)*X/S1
      A2=UT/S1
      CALL CEXF(A1,A2,AM,X1,Y1)
      IF(KSURF.EQ.0) GO TO 32
      CP1=(R(I)*X+VT)/S1
      CM1=V*X/D-DONE(I)*T/R(I)
      C1=2.*V/(V*V-U2)
      C2=0.5*V*V/(D*(DONE(I)-A*R(I)))
      IF(JJ.GT.0) GO TO 30
      F=C1*(CS*(X1*V-Y1*U)-SN*(X1*U+Y1*V))+C2*EXF(CM1,CP1)
      RETURN
   30 CONTINUE
      S2=2.*DSQRT(D*R(J)*T)
      CM2=V*X/D-DONE(J)*T/R(J)
      C3=0.5*V*V/(D*(DONE(J)-A*R(J)))
      CP2=(R(J)*X+VT)/S2
      A1=R(J)*X/S2
      A2=UT/S2
      CALL CEXF(A1,A2,AM,X2,Y2)
      F=C1*(CS*(V*(X2-X1)-U*(Y2-Y1))-SN*(U*(X2-X1)+V*(Y2-Y1)))+C3*EXF(CM
     12,CP2)-C2*EXF(CM1,CP1)
      RETURN
   32 IF(JJ.GT.0) GO TO 34
      F=CS*X1-SN*Y1
      RETURN
   34 S2=2.*DSQRT(D*R(J)*T)
      A1=R(J)*X/S2
      A2=UT/S2
      CALL CEXF(A1,A2,AM,X2,Y2)
      F=CS*(X2-X1)-SN*(Y2-Y1)
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
      IF(DABS(B).GT.1.D-10) GO TO 1
      EXF=DEXP(A)
      RETURN
    1 C=A-B*B
      IF((DABS(C).GT.170.).AND.(B.GT.0.)) RETURN
      IF(C.LT.-170.) GO TO 4
      X=DABS(B)
      IF(X.GT.2.5) GO TO 2
      T=1.D0/(1.D0+.3275911D0*X)
      Y=T*(.254829592D0-T*(.284496736D0-T*(1.421413741D0-T*(1.453152027D
     10-1.061405429D0*T))))
      GO TO 3
    2 Y=.564189584D0/(X+.5D0/(X+1.D0/(X+1.5D0/(X+2.D0/(X+2.5D0/(X+3.D0/(
     1X+1.D0)))))))
    3 EXF=Y*DEXP(C)
    4 IF(B.LT.0.0) EXF=2.D0*DEXP(A)-EXF
      RETURN
      END
C
C     ------------------------------------------------------------------
      SUBROUTINE CEXF(A,B,Z,U,V)
C
C     COMPLEX ERFC-FUNCTION: U+IV=EXP(Z)ERFC(A+IB)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION W(10),H(10)
      DATA W/.4622437,.2866755,.1090172,.02481052,.3243773D-2,
     1.2283386D-3,.7802556D-5,.1086069D-6,.4399341D-9,.2229394D-12/,
     2H/.2453407,.7374737,1.234076,1.738538,2.254974,
     32.788806,3.347855,3.944764,4.603682,5.387481/
C
C     ----------
      X=DABS(A)
      Y=DABS(B)
      U=0.0
      V=0.0
      E=0.0
      XYZ=Y*Y+Z-X*X
      IF(DABS(XYZ).GT.170.)RETURN
      COS=DCOS(2.*X*Y)
      SIN=DSIN(2.*X*Y)
      IF((X+Y).LT.6.) GO TO 2
      T=0.0
      DO 1 K=1,10
      T=T+W(K)*((X/((Y-H(K))**2+X*X))+(X/((Y+H(K))**2+X*X)))
    1 V=V+W(K)*((Y-H(K))/((Y-H(K))**2+X*X)+(Y+H(K))/((Y+H(K))**2+X*X))
      U=.3183099*DEXP(XYZ)*(T*COS-V*SIN)
      V=.3183099*DEXP(XYZ)*(-T*SIN-V*COS)
      IF(X.LE.0.)U=DEXP(DMIN1(Z,1.7D2))
      GO TO 8
    2 IF(X.GT.2.5) GO TO 3
      T=1./(1.+.3275911*X)
      U=T*(.2548296-T*(.2844967-T*(1.421414-T*(1.453152-1.061405*T))))
      GO TO 4
    3 U=.5641896/(X+.5/(X+1./(X+1.5/(X+2./(X+2.5/(X+3./(X+1.)))))))
    4 IF(Y.LE.0.) GO TO 7
      IF(X.LE.0.) V=-.3183099*Y
      IF(X.LE.0.) GO TO 5
      U=U-.1591549*(1.D0-COS)/X
      V=V-.1591549*SIN/X
    5 NT=12.+2.*Y
      DO 6 I=1,NT
      P=I
      ARG=P*Y
      F1=X*(DEXP(ARG)+DEXP(-ARG))
      F2=0.5*P*(DEXP(ARG)-DEXP(-ARG))
      EX=.6366198*DEXP(-0.25*P*P)/(4.*X*X+P*P)
      U=U-EX*(2.*X-F1*COS+F2*SIN)
    6 V=V-EX*(F1*SIN+F2*COS)
      V=V*DEXP(Z-X*X)
    7 U=U*DEXP(Z-X*X)
    8 IF(B.LT.0.) V=-V
      IF(A.LT.0.) U=2.*EXF(Z,E)-U
      RETURN
      END
