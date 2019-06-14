C     ********************************************************************
C     *                                                                  *
C     *                           N3DADE                                 *
C     *                                                                  *
C     *         ANALYTICAL SOLUTIONS FOR NONEQUILIBRIUM SOLUTE           *
C     *          TRANSPORT IN 3-D SEMI-INFINITE POROUS MEDIA             *
C     *                 JOURNAL OF HYDROLOGY 151:193-228                 *
C     *                                                                  *
C     *                 FEIKE J. LEIJ and NOBUO TORIDE                   *
C     *                                                                  *   
C     *                    U.S. SALINITY LABORATORY                      *
C     *                    450 W. BIG SPRINGS ROAD                       *
C     *                    RIVERSIDE, CA 92507-4617                      *
C     *                    U.S.A.                                        *
C     *                                                                  *           
C     *                      PHONE (909) 369-4851                        *
C     *                      FAX   (909) 342-4964                        *
C     *                                                                  *   
C     *                         MODIFIED: 3/98                           *
C     ********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (ZIP=1.0D-10,ZERO=0.0D+00,ONE=1.0D+00)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/BVP/ NBP,T0B,DMB,GPB(10),TPB(10),GEB,ALB,JB  
      COMMON/IVP/ NIP,X0I,DMI,FPI(10),XPI(10),FEI,ALI,JI
      COMMON/PVP/ NPP,OLAM(10),TLAM(10),XPP(10),X0P,OLAME,TLAME,ALP,JP
      COMMON/TRANS/ Y0,Z0,NYP,YP(5),NZP,ZP(5),ALY,ALZ
      DIMENSION GP(5),TP(5),FP(5),XP(5),XD(4),XLOW(4),XUP(4),NXST(4),DX(
     &4)
      CHARACTER TITLE*60,VN(4)*5,FILENAME*15
      DATA VN/'t','x','y','z'/
C
C     --- OPEN I/O FILES ---
      WRITE(*, 150)
      READ(*, 160) FILENAME
      IF (FILENAME.EQ.'  ') FILENAME='N3DADE.IN'
      OPEN(5, FILE=FILENAME, STATUS='OLD')
C
      WRITE(*,170)
      READ(*,160) FILENAME
      IF (FILENAME.EQ.'  ') FILENAME='N3DADE.OUT'
      OPEN(7, FILE=FILENAME, STATUS='UNKNOWN')
C
C     --- CLEAR SCREEN AND WRITE TITLE ---
      WRITE(*,190) FILENAME
      READ(5,200) NC,NGC
C
      DO 120 NCASE=1,NC
C
      READ(5,201) TITLE
      READ(5,202) V,R,DCX,DCY,DCZ
      READ(5,202) BETA,OMEGA,U1,U2
C                     
      READ(5,204) UL,CO,MP,MC
      IF(MC.EQ.1) THEN
       WRITE(7,300) TITLE
      ELSE
       WRITE(7,299) TITLE
      ENDIF 
      WRITE(7,210) V,R,DCX,DCY,DCZ,BETA,OMEGA,U1,U2
      PX=V*UL/DCX
      PY=V*UL/DCY
      PZ=V*UL/DCZ
      WRITE(7,301) UL,CO
      OMU1=(OMEGA+U1)/BETA/R
      OMU2=(OMEGA+U2)/(ONE-BETA)/R
      OM1=OMEGA/(ONE-BETA)/R
      OM2=OMEGA/BETA/R
      OMOM=OM1*OM2
      U1R=U1/BETA/R
      U2R=OMEGA*U2/(OMEGA+U2)/BETA/R
      IF(MP.EQ.2) GOTO 3
      IF(MP.EQ.3) GOTO 11
C
C     --- BOUNDARY VALUE PROBLEM ---
      READ(5,201) TITLE
      READ(5,200) NP,NBP,NYP,NZP
      WRITE(7,211) NP          
C     --- DIRAC ---
      IF(NP.GE.1 .AND. NP.LE.9 .OR. NP.EQ.28 .OR. NP.EQ.29) THEN
       READ(5,202) T0,DMB
       WRITE(7,212) T0,DMB
       DMB=DMB/CO/UL
       T0B=V*T0/UL
C     --- HEAVISIDE ---
      ELSE IF(NP.GE.10.AND.NP.LE.18.OR.NP.EQ.30.OR.NP.EQ.31)THEN
       WRITE(7,213)
       TP(1)=0.
       TPB(1)=0. 
       READ(5,202) (GP(I),TP(I+1),I=1,NBP)
        DO 2 I=1,NBP
         WRITE(7,214) I,GP(I),TP(I),TP(I+1)
         GPB(I)=GP(I)/CO
    2    TPB(I+1)=V*TP(I+1)/UL
C     --- EXPONENTIAL ---
      ELSE IF(NP.GE.19.AND.NP.LE.27.OR.NP.EQ.32.OR.NP.EQ.33)THEN
       READ(5,202) GE,ALB
       WRITE(7,215) GE,ALB
       GEB=GE/CO
      ENDIF
      GOTO 13 
C
C     --- INITIAL VALUE PROBLEM ---
    3 READ(5,201) TITLE
      READ(5,200) NP,NIP,NYP,NZP
      WRITE(7,221) NP
C     --- DIRAC ---
      IF(NP.GE.1 .AND. NP.LE.9 .OR. NP.EQ.28 .OR. NP.EQ.29) THEN
       READ(5,202) X0,DMI
       WRITE(7,222) X0,DMI
       DMI=DMI/CO/UL
       X0I=X0/UL
C     --- HEAVISIDE ---
      ELSE IF(NP.GE.10.AND.NP.LE.18.OR.NP.EQ.30.OR.NP.EQ.31)THEN
       WRITE(7,223)
       READ(5,202) XP(1),(FP(I),XP(I+1),I=1,NIP)
       XPI(1)=XP(1)/UL
        DO 4 I=1,NIP
        WRITE(7,214) I,FP(I),XP(I),XP(I+1)
        FPI(I)=FP(I)/CO
    4   XPI(I+1)=XP(I+1)/UL
C     --- EXPONENTIAL ---
      ELSE IF(NP.GE.19.AND.NP.LE.27.OR.NP.EQ.32.OR.NP.EQ.33)THEN
       READ(5,202) FE,ALI
       WRITE(7,224) FE,ALI
       FEI=FE/CO
C     --- spherical HEAVISIDE ---
      ELSE IF(NP.EQ.34)THEN
       READ(5,202) X0,(FP(I),XP(I),I=1,NIP)
       X0I=X0/UL
       I=1
       FPI(I)=FP(I)/CO
       XPI(I)=XP(I)/UL
       WRITE(7,225) X0
       WRITE(7,214) I,FP(I),DUM,XP(I)
       DO 5 I=2,NIP
        WRITE(7,214) I,FP(I),XP(I-1),XP(I)
        FPI(I)=FP(I)/CO
    5   XPI(I)=XP(I)/UL
C     --- spherical EXPONENTIAL ---
      ELSE IF(NP.EQ.35)THEN
       READ(5,202) X0,FE,ALI,BINT
       WRITE(7,226) X0,FE,ALI
       X0I=X0/UL
       BINT=BINT/UL
       FEI=FE/CO
      ENDIF        
      GOTO 13
C
C     --- PRODUCTION VALUE PROBLEM ---
   11 READ(5,201) TITLE
      READ(5,200) NP,NPP,NYP,NZP
      WRITE(7,231) NP
C     --- HEAVISIDE ---
      IF(NP.EQ.14.OR.NP.EQ.15.OR.NP.EQ.17.OR.NP.EQ.18.OR.NP.EQ
     &.30.OR.NP.EQ.31) THEN
       WRITE(7,232)
       READ(5,202) XP(1),(OLAM(I),TLAM(I),XP(I+1),I=1,NPP)
       XPP(1)=XP(1)/UL
        DO 6 I=1,NPP
         WRITE(7,214) I,OLAM(I),TLAM(I),XP(I),XP(I+1)
    6    XPP(I+1)=XP(I+1)/UL
C     --- EXPONENTIAL ---
      ELSE IF(NP.EQ.23.OR.NP.EQ.24.OR.NP.EQ.26.OR.NP.EQ.27.OR.NP.EQ.32
     & .OR.NP.EQ.33) THEN
       READ(5,202) OLAME,TLAME,ALP
       WRITE(7,233) OLAME,TLAME,ALP
C     --- spherical HEAVISIDE ---
      ELSE IF(NP.EQ.34)THEN
       READ(5,202) X0,(OLAM(I),TLAM(I),XP(I),I=1,NPP)
       DUM=0.
       I=1
       X0P=X0/UL
       XPP(1)=XP(1)/UL
       WRITE(7,234) X0
       WRITE(7,214) I,OLAM(I),TLAM(I),DUM,XP(I)
       DO 7 I=2,NPP
    7   WRITE(7,214) I,OLAM(I),TLAM(I),XP(I-1),XP(I)
C     --- spherical EXPONENTIAL ---
      ELSE IF(NP.EQ.35)THEN
       READ(5,202) X0,OLAME,TLAME,ALP,BINT
       WRITE(7,235) X0,OLAME,TLAME,ALP
       X0P=X0/UL                               
       BINT=BINT/UL
      ENDIF
C
C     --- TRANSVERSAL CONDITIONS ---
C     --- DIRAC/DIRAC ---
   13 IF(NP.EQ.1 .OR. NP.EQ.10 .OR. NP.EQ.19) THEN
       READ(5,202) Y0
       WRITE(7,241) Y0
       READ(5,202) Z0
       WRITE(7,242) Z0
       Y0=Y0/UL
       Z0=Z0/UL              
C     --- DIRAC/HEAVISIDE ---
      ELSE IF(NP.EQ.2 .OR. NP.EQ.11 .OR. NP.EQ.20) THEN
       READ(5,202) Y0
       READ(5,202) (ZP(I),ZP(I+1),I=1,2*NZP,2)
       WRITE(7,241) Y0
       WRITE(7,244) 
       DO 14 I=1,2*NZP,2
   14   WRITE(7,214) (I+1)/2,ZP(I),ZP(I+1)
       DO 814 I=1,2*NZP
  814   ZP(I)=ZP(I)/UL       
       Y0=Y0/UL
C     --- DIRAC/EXPONENTIAL ---
      ELSE IF(NP.EQ.3 .OR. NP.EQ.12 .OR. NP.EQ.21) THEN
       READ(5,202) Y0
       READ(5,202) ALZ,Z0
       WRITE(7,241) Y0
       WRITE(7,246) ALZ,Z0
       Y0=Y0/UL             
       Z0=Z0/UL
C     --- HEAVISIDE/DIRAC ---
      ELSE IF(NP.EQ.4 .OR. NP.EQ.13 .OR. NP.EQ.22) THEN
       READ(5,202) (YP(I),YP(I+1),I=1,2*NYP,2)
       READ(5,202) Z0
       WRITE(7,243) 
       DO 15 I=1,2*NYP,2
   15   WRITE(7,214) (I+1)/2,YP(I),YP(I+1)
       DO 815 I=1,2*NYP
  815   YP(I)=YP(I)/UL        
       WRITE(7,242) Z0
       Z0=Z0/UL
C     --- HEAVISIDE/HEAVISIDE ---
      ELSE IF(NP.EQ.5 .OR. NP.EQ.14 .OR. NP.EQ.23) THEN
       READ(5,202) (YP(I),YP(I+1),I=1,2*NYP,2)
       READ(5,202) (ZP(I),ZP(I+1),I=1,2*NZP,2)
       WRITE(7,243)
       DO 16 I=1,2*NYP,2
   16   WRITE(7,214) (I+1)/2,YP(I),YP(I+1)
       DO 816 I=1,2*NYP
  816   YP(I)=YP(I)/UL        
       WRITE(7,244)
       DO 17 I=1,2*NZP,2   
   17   WRITE(7,214) (I+1)/2,ZP(I),ZP(I+1)
       DO 817 I=1,2*NZP 
  817   ZP(I)=ZP(I)/UL        
C     --- HEAVISIDE/EXPONENTIAL ---
      ELSE IF(NP.EQ.6 .OR. NP.EQ.15 .OR. NP.EQ.24) THEN
       READ(5,202) (YP(I),YP(I+1),I=1,2*NYP,2)
       READ(5,202) ALZ,Z0
       WRITE(7,243)  
       DO 18 I=1,2*NYP,2
   18   WRITE(7,214) (I+1)/2,YP(I),YP(I+1)
       DO 818 I=1,2*NYP
  818   YP(I)=YP(I)/UL        
       WRITE(7,246) ALZ,Z0
       Z0=Z0/UL
C     --- EXPONENTIAL/DIRAC ---
      ELSE IF(NP.EQ.7 .OR. NP.EQ.16 .OR. NP.EQ.25) THEN
       READ(5,202) ALY,Y0
       READ(5,202) Z0
       WRITE(7,245) ALY,Y0       
       WRITE(7,242) Z0
       Y0=Y0/UL             
       Z0=Z0/UL
C     --- EXPONENTIAL/HEAVISIDE ---
      ELSE IF(NP.EQ.8 .OR. NP.EQ.17 .OR. NP.EQ.26) THEN
       READ(5,202) ALY,Y0
       WRITE(7,245) ALY,Y0       
       READ(5,202) (ZP(I),ZP(I+1),I=1,2*NZP,2)
       WRITE(7,244)
       DO 19 I=1,2*NZP,2
   19   WRITE(7,214) (I+1)/2,ZP(I),ZP(I+1)
       DO 819 I=1,2*NZP
  819   ZP(I)=ZP(I)/UL        
       Y0=Y0/UL 
C     --- EXPONENTIAL/EXPONENTIAL ---
      ELSE IF(NP.EQ.9 .OR. NP.EQ.18 .OR. NP.EQ.27) THEN
       READ(5,202) ALY,Y0
       WRITE(7,245) ALY,Y0       
       Y0=Y0/UL  
       READ(5,202) ALZ,Z0
       WRITE(7,246) ALZ,Z0
       Z0=Z0/UL
C     --- circular HEAVISIDE ---
      ELSE IF(NP.EQ.28.OR. NP.EQ.30 .OR. NP.EQ.32) THEN
       READ(5,202) (YP(I),YP(I+1),I=1,2*NYP,2)
       WRITE(7,247)
       DO 20 I=1,2*NYP,2
   20   WRITE(7,214) (I+1)/2,YP(I),YP(I+1)
       DO 820 I=1,2*NYP
  820   YP(I)=YP(I)/UL        
C     --- circular EXPONENTIAL ---
      ELSE IF(NP.EQ.29.OR. NP.EQ.31 .OR. NP.EQ.33) THEN
       READ(5,202) ALY,BINT
       BINT=BINT/UL
       WRITE(7,248) ALY
      ENDIF
C
C     --- CALCULATE NODE SPACING ---  
      WRITE(7,338)
      DO 21 I=1,4
       READ(5,202) DX(I),XLOW(I),XUP(I)
   21  WRITE(7,340) VN(I),DX(I),XLOW(I),XUP(I)
C       
       WRITE(7,342) CHAR(12)
C
      DO 23 I=1,4
      IF (ABS(DX(I)).LE.ZIP) THEN
        NXST(I)=1
      ELSE
        NXST(I)=ZIP+(XUP(I)+DX(I)-XLOW(I))/DX(I)
      ENDIF
   23 CONTINUE
C
      NUM=0
      DO 22 I=1,NXST(1)
       XD(1)=XLOW(1)+(I-1)*DX(1)
       X(1)=XD(1)*V/UL
       DO 22 J=1,NXST(2)
        XD(2)=XLOW(2)+(J-1)*DX(2)
        X(2)=XD(2)/UL
        DO 22 K=1,NXST(3)
         XD(3)=XLOW(3)+(K-1)*DX(3)   
         X(3)=XD(3)/UL
         DO 22 L=1,NXST(4)
          XD(4)=XLOW(4)+(L-1)*DX(4)   
          X(4)=XD(4)/UL
          NUM=NUM+1
          Y1=ZERO
          Y2=ZERO     
          CALL MODELS(Y1,Y2)
       YT=BETA*R*Y1+(ONE-BETA)*R*Y2       
       YD1=Y1*CO
       YD2=Y2*CO
       YDT=YT*CO
       WRITE(7,346) NUM,YD1,Y1,YD2,Y2,YDT,YT,(XD(JJ),JJ=1,4)
   22 CONTINUE
  120 CONTINUE
C
C     --- END OF PROBLEM ---
  150 FORMAT(' Enter input file name (default = N3DADE.IN)')
  160 FORMAT(A15)
  170 FORMAT(' Enter output file name (default = N3DADE.OUT)')
  190 FORMAT(////////////////5X,67(1H*)/5X,1H*,65X,1H*/5X,1H*,29X,'N3DAD
     $E',30X,1H*/5X,1H*,65X,1H*/5X,1H*,4X,'PLEASE WAIT - The concentrati
     $ons are being calculated and',4X,1H*/5x,1H*,18X,'written to ',A15,
     $21X,1H*/,5X,1H*,65X,1H*/5X,67(1H*)///)
  200 FORMAT(5I5)
  201 FORMAT(A60)
  202 FORMAT(8F10.2)
  204 FORMAT(2F10.2,2I10)
  210 FORMAT(//5X,'MODEL PARAMETERS'/5X,16(1H=)//5X,'NAME',12X,'VALUE'//
     $5X,'dimensional'/5X,'V',10X,F10.2/5X,'R',10X,F10.2/5X,'Dx',9X,F10.
     &2/5X,'Dy',9X,F10.2/5X,'Dz',9X,F10.2//5X,'dimensionless'/5X,'BETA',
     &7X,F10.2/5X,'OMEGA',6X,F10.2/5X,'MU1',$8X,F
     $10.2/5X,'MU2',8X,F10.2)
  211 FORMAT(//5X,'BOUNDARY VALUE PROBLEM ',i5/5X,22(1H=))
  212 FORMAT(5X, 'DIRAC INPUT FUNCTION'/7X,'t0=',F10.2,8X,'m=',F10.2)
  213 FORMAT(5X, 'HEAVISIDE INPUT FUNCTION'/9X,'i',9X,'g',7X,'tin',6X,
     &'tend')
  214 FORMAT(5X,I5,6(F10.2))
  215 FORMAT(5X, 'EXPONENTIAL INPUT FUNCTION'/7X,'ge=',F10.2,4X,'alpha='
     &,F10.2)
  221 FORMAT(//5X,'INITIAL VALUE PROBLEM ',i5/5X,21(1H=))
  222 FORMAT(/5X,'DIRAC INITIAL DISTRIBUTION'/7X,'x0=',F10.2,8X,'m=',F10
     $.2)
  223 FORMAT(/5X,'HEAVISIDE INITIAL DISTRIBUTION'/9X,'i',9X,'f',7X,'xin'
     &,6X,'xend')
  224 FORMAT(/5X,'EXPONENTIAL INITIAL DISTRIBUTION'/7X,'fe=',F10.2,4X,'a
     &lpha=',F10.2)
  225 FORMAT(/5X,'HEAVISIDE INITIAL DISTRIBUTION FOR SPHERICAL GEOMETRY  
     &   x0=',F10.2/9X,'i',9X,'f',7X,'rin',6X,'rend')
  226 FORMAT(/5X,'EXPONENTIAL INITIAL DISTRIBUTION FOR SPHERICAL GEOMETR
     &Y    x0=',F10.2/7X,'fe=',F10.2,4X,'alpha=',F10.2)
  231 FORMAT(//5X,'PRODUCTION VALUE PROBLEM ',i5/5X,24(1H=))
  232 FORMAT(5X,  'HEAVISIDE PRODUCTION PROFILE'/9X,'i',6X,'lam1',6X,'la
     &m2',7X,'xin',6X,'xend')
  233 FORMAT(5X,  'EXPONENTIAL PRODUCTION PROFILE'/4X,'lam1e=',F10.2,4X,
     &'lam2e=',F10.2,4x,'alpha=',F10.2)
  234 FORMAT(5X, 'HEAVISIDE PRODUCTION PROFILE FOR SPHERICAL GEOMETRY  
     &   x0=',F10.2/9X,'i',6X,'lam1',6X,'lam2',7X,'rin',6X,'rend')
  235 FORMAT(5X, 'EXPONENTIAL INITIAL DISTRIBUTION FOR SPHERICAL GEOMETR
     &Y    x0=',F10.2/4X,'lam1e=',F10.2,4X,'lam2e=',F10.2,4X,'alpha=',
     &F10.2)
  241 FORMAT(/5X, 'DIRAC DISTRIBUTION FOR Y-DIRECTION'/7X,'y0=',F10.2)
  242 FORMAT(/5X, 'DIRAC DISTRIBUTION FOR Z-DIRECTION'/7X,'z0=',F10.2)
  243 FORMAT(/5X,'HEAVISIDE DISTRIBUTION FOR Y-DIRECTION'/9X,'i',7X,'yin
     &',6X,'yend')
  244 FORMAT(/5X,'HEAVISIDE DISTRIBUTION FOR Z-DIRECTION'/9X,'i',7X,'zin
     &',6X,'zend')
  245 FORMAT(/5X,'EXPONENTIAL DISTRIBUTION FOR Y-DIRECTION'/4X,'alpha=',
     &F10.2,4X,'y0=',F10.2)
  246 FORMAT(/5X,'EXPONENTIAL DISTRIBUTION FOR Z-DIRECTION'/4X,'alpha=',
     &F10.2,4X,'z0=',F10.2)
  247 FORMAT(/5X,'HEAVISIDE DISTRIBUTION FOR CIRCULAR GEOMETRY'/9X,'i',7
     &X,'rin',6X,'rend')
  248 FORMAT(/5X,'EXPONENTIAL DISTRIBUTION FOR CIRCULAR GEOMETRY'/4X,'al
     &pha=',F10.2)
  299 FORMAT(//5X,67(1H*)/5X,1H*,65X,1H*/5X,1H*,29X,'N3DADE',30X,1H*/5X,
     $1H*,65X,1H*/5X,1H*,4X,'ANALYTICAL SOLUTIONS FOR NONEQUILIBRIUM SOL
     $UTE TRANSPORT',5X,1H*/5X,1H*,16X,'IN 3-D SEMI-INFINITE POROUS MEDI
     $A',16X,1H*/5X,1H*,65X,1H*/5X,1H*,5X,'VOLUME-AVERAGED CONCENTRATION
     $S',30X,1H*/5X,1H*,5X,A60,1H*/5X,1H*,65X,1H*/5X,67(1H*))
  300 FORMAT(//5X,67(1H*)/5X,1H*,65X,1H*/5X,1H*,29X,'N3DADE',30X,1H*/5X,
     $1H*,65X,1H*/5X,1H*,4X,'ANALYTICAL SOLUTIONS FOR NONEQUILIBRIUM SOL
     $UTE TRANSPORT',5X,1H*/5X,1H*,16X,'IN 3-D SEMI-INFINITE POROUS MEDI
     $A',16X,1H*/5X,1H*,65X,1H*/5X,1H*,5X,'FLUX-AVERAGED CONCENTRATIONS'
     $,32X,1H*/5X,1H*,5X,A60,1H*/5X,1H*,65X,1H*/5X,67(1H*))
  301 FORMAT(//5X,'CONSTANTS'/5X,9(1H=)/5X,'L  =',F8.2/5X,'Co =',F8.2)
  338 FORMAT(//5X,'INDEPENDENT VARIABLES'/5X,21(1H=)/5X,'VARIABLE',3X,'I
     $NCREMENT',8X,'MINIMUM',8X,'MAXIMUM')
  340 FORMAT(5X,A5,3(5X,F10.3))
  342 FORMAT(/A1/'PREDICTED CONCENTRATIONS'/24(1H=)//6X,'EQUILIBRIUM',7X
     $,'NONEQUILIBRIUM',7X,'TOTAL'/1X,'NO',4X,'C',6X,'C/Co',6X,'C',6X,'C
     $/Co',6X,'Ct',5x,'Ct/Co',8X,'t',7X,'x',7X,'y',7X,'z')
  346 FORMAT(I4,X,3(E8.2,E8.2,X),4F8.2)
C
      CLOSE(5)
      CLOSE(7)
C
      STOP
      END
C
C     ------------
      SUBROUTINE MODELS(Y1,Y2)
C
C     PURPOSE: TO CALCULATE EQUILIBRIUM AND NONEQUILIBRIUM CONCENTRATION
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/BVP/ NBP,T0B,DMB,GPB(10),TPB(10),GEB,ALB,JB  
      COMMON/IVP/ NIP,X0I,DMI,FPI(10),XPI(10),FEI,ALI,JI
      COMMON/PVP/ NPP,OLAM(10),TLAM(10),XPP(10),X0P,OLAME,TLAME,ALP,JP
      COMMON/TRANS/ Y0,Z0,NYP,YP(5),NZP,ZP(5),ALY,ALZ
      PARAMETER(ZIP=1.0D-8,ZERO=0.0D+00,ONE=1.0D+00,TWO=2.0D+00,
     $ FOUR=4.0D+00,EIGHT=8.0D+00)
      EXTERNAL DI1B,DI2B,HE1B,HE2B,EX1B,EX2B,
     $DI1I,DI2I,HE1I,HE2I,EX1I,EX2I,H1S1I,H1S2I,H2S1I,E1S1I,E1S2I,E2S1I,
     $HE1P,HE2P,EX1P,EX2P,H1S1P,H2SP,E1S1P,E2SP 
      PI=FOUR*DATAN(ONE)
C
      Y1=0.
      Y2=0. 
      IF(MP.EQ.2) GOTO 3
      IF(MP.EQ.3) GOTO 11
C
C     --- BOUNDARY VALUE PROBLEM ---
C     --- DIRAC ---
      IF(NP.GE.1 .AND. NP.LE.9 .OR. NP.EQ.28 .OR. NP.EQ.29) THEN
       TT=X(1)-T0B            
       IF(TT.LE.ZERO) GOTO 103
       Y1=DEXP(-OM2*TT)*GAMY(TT)*GAMZ(TT)*GB(X(2),TT)
       CALL CHEB10(DI1B,RES,ZERO,TT,NGC)
       Y1=DMB*(Y1+RES)
       CALL CHEB10(DI2B,RES,ZERO,TT,NGC)
       Y2=DMB*OM1*RES
C  
C     --- HEAVISIDE ---                                               
      ELSE IF(NP.GE.10.AND.NP.LE.18.OR.NP.EQ.30.OR.NP.EQ.31)THEN
      IF(X(1).LE.TPB(1)) GOTO 103
      DO 31 JB=1,NBP
       IF(X(1).GT.TPB(JB)) THEN
        TLOW=TPB(JB)
        TUP=X(1)
        CALL CHEBHB(HE1B,RES1,TLOW,TUP,NGC)
        CALL CHEBHB(HE2B,RES2,TLOW,TUP,NGC)
        Y1=Y1+GPB(JB)*RES1
        Y2=Y2+GPB(JB)*RES2
        IF(X(1).GT.TPB(JB+1)) THEN
	   TLOW=TPB(JB)
	   TUP=X(1)-TPB(JB+1)
         CALL CHEBHB(HE1B,RES1,TLOW,TUP,NGC)
         CALL CHEBHB(HE2B,RES2,TLOW,TUP,NGC)
         Y1=Y1-GPB(JB)*RES1
         Y2=Y2-GPB(JB)*RES2
        ENDIF
       ELSE
        GOTO 32
       ENDIF   
   31 CONTINUE
   32 Y2=OMEGA*Y2/(OMEGA+U2)
C
C     --- EXPONENTIAL ---
      ELSE IF (NP.GE.19.AND.NP.LE.27.OR.NP.EQ.32.OR.NP.EQ.33)THEN
       CALL CHEB10(EX1B,RES,ZERO,X(1),NGC)
       Y1=GEB*RES
       CALL CHEB10(EX2B,RES,ZERO,X(1),NGC) 
       Y2=OM1*GEB*RES 
      ENDIF
      GOTO 103
C
C     --- INITIAL VALUE PROBLEM ---
C     --- DIRAC ---
    3 IF(NP.GE.1 .AND. NP.LE.9 .OR. NP.EQ.28 .OR. NP.EQ.29) THEN
       YI=DEXP(-OM2*X(1))*GAMY(X(1))*GAMZ(X(1))*G(X0I,X(1),X(2))
       CALL CHEB10(DI1I,RES,ZERO,X(1),NGC)
       Y1=DMI*(YI+RES)
C
       CALL CHEB10(DI2I,RES,ZERO,X(1),NGC)
       Y2=DMI*OM1*RES
C  
C     --- HEAVISIDE ---                                         
      ELSE IF(NP.GE.10.AND.NP.LE.18.OR.NP.EQ.30.OR.NP.EQ.31)THEN
       DO 34 JI=1,NIP
        if(x(1).le.zip) then
         if(x(2).ge.xpi(ji) .and. x(2).lt.xpi(ji+1)) y1=fpi(ji)
     &    *YTRANS(X(3))*ZTRANS(X(4))    
        else   
         Y1I=DEXP(-OM2*X(1))*(G1(X(2),X(1),XPI(JI))
     &     -G1(X(2),X(1),XPI(JI+1)))*GAMY(X(1))*GAMZ(X(1))
         CALL CHEB10(HE1I,RES1,ZERO,X(1),NGC)
         Y1=Y1+FPI(JI)*(Y1I+RES1)
        endif 
        Y2I=0.
        IF(X(2).GE.XPI(JI).AND.X(2).LT.XPI(JI+1)) Y2I=YTRANS(X(3))
     &   *ZTRANS(X(4))*DEXP(-OMU2*X(1))
        CALL CHEB10(HE2I,RES2,ZERO,X(1),NGC)
        Y2=Y2+FPI(JI)*(Y2I+OM1*RES2)
   34  CONTINUE
C 
C     --- EXPONENTIAL ---
      ELSE IF(NP.GE.19.AND.NP.LE.27.OR.NP.EQ.32.OR.NP.EQ.33)THEN
        if(x(1).le.zip) then
         y1=fei*dexp(-ali*x(2))*YTRANS(X(3))*ZTRANS(X(4))    
        else   
         Y1=YTRANS(X(3))*ZTRANS(X(4))*DEXP(-OM2*X(1))*G2(X(2),X(1),ALI)
         CALL CHEB10(EX1I,RES1,ZERO,X(1),NGC)
         Y1=FEI*(Y1+RES1)
        endif 
       Y2=YTRANS(X(3))*ZTRANS(X(4))*DEXP(-ALI*X(2)-OMU2*X(1))
       CALL CHEB10(EX2I,RES2,ZERO,X(1),NGC) 
       Y2=FEI*(Y2+OM1*RES2)
C
C     --- spherical HEAVISIDE ---
      ELSE IF(NP.EQ.34)THEN
       RS=SQRT((X(2)-X0I)**2+X(3)**2+X(4)**2)
       JI=1                 
       RI=XPI(JI)
       if(x(1).le.zip) then
        if(rs.lt.xpi(ji)) then
         y1=fpi(ji)
        else
         do 77 i=2,nip
          if(rs.lt.xpi(i) .and. rs.ge.xpi(i-1)) y1=fpi(i)
   77    continue
        endif
       else         
        BOT=X0I-RI
        TOP=X0I+RI
        CALL CHEB11(H1S1I,RI,RES1,BOT,TOP,NGC)
        CALL CHEB11(H1S2I,RI,RES2,ZERO,X(1),NGC)
        Y1=FPI(JI)*(DEXP(-OM2*X(1))*RES1+RES2) 
       endif
       Y2I=0       
       IF(RS.LT.XPI(JI)) Y2I=DEXP(-OMU2*X(1))               
       CALL CHEB11(H2S1I,RI,RES,ZERO,X(1),NGC)
       Y2=FPI(JI)*(Y2I+OM1*RES)
       DO 35 JI=2,NIP
        RI=XPI(JI)
        BOT=X0I-RI
        TOP=X0I+RI
        CALL CHEB11(H1S1I,RI,RES1,BOT,TOP,NGC)
        CALL CHEB11(H1S2I,RI,RES2,ZERO,X(1),NGC)
        Y1=Y1+FPI(JI)*(DEXP(-OM2*X(1))*RES1+RES2)
        Y2I=0.
        IF(RS.GE.XPI(JI-1) .AND. RS.LT.XPI(JI)) Y2I=DEXP(-OMU2*X(1))
        CALL CHEB11(H2S1I,RI,RES,ZERO,X(1),NGC)
        Y2=FPI(JI)*(Y2I+OM1*RES)
        RI=XPI(JI-1)
        BOT=X0I-RI
        TOP=X0I+RI
        CALL CHEB11(H1S1I,RI,RES1,BOT,TOP,NGC)
        CALL CHEB11(H1S2I,RI,RES2,ZERO,X(1),NGC)
        Y1=Y1-FPI(JI)*(DEXP(-OM2*X(1))*RES1+RES2)
        CALL CHEB11(H2S1I,RI,RES,ZERO,X(1),NGC)
        Y2=Y2-FPI(JI)*OM1*RES
   35  CONTINUE
C
C     --- spherical EXPONENTIAL ---
      ELSE IF(NP.EQ.35) THEN
       if(x(1).lt.zip) then
        Y1=FEI*DEXP(-ALI*SQRT((X(2)-X0I)**2+X(3)**2+X(4)**2))
       else
        bot=zero
        top=x0i+bint
        CALL CHEB10(E1S1I,res1,bot,top,NGC)
       endif
        CALL CHEB10(E1S2I,RES2,zero,x(1),NGC) 
        Y1=FEI*(res1+RES2)
        CALL CHEB10(E2S1I,RES,zero,top,NGC) 
        Y2=FEI*(DEXP(-OMU2*X(1)-ALI*SQRT((X(2)-X0I)**2+X(3)**2+X(4)**2))
     &   +OM1*RES)  
c 
      ENDIF   
      GOTO 103 
C          
C     --- PRODUCTION VALUE PROBLEM ---
C     --- HEAVISIDE ---
   11 IF(X(1).LT.ZIP) GOTO 103
      IF(NP.EQ.14.OR.NP.EQ.15.OR.NP.EQ.17.OR.NP.EQ.18.OR.NP.EQ.30
     & .OR.NP.EQ.31) THEN
       DO 36 JP=1,NPP
        CALL CHEB10(HE1P,RES,ZERO,X(1),NGC)
        Y1=Y1+RES
        Y2P=0.
        IF(X(2).GE.XPP(JP).AND.X(2).LT.XPP(JP+1)) Y2P=YTRANS(X(3))
     &   *ZTRANS(X(4))*TLAM(JP)*(1-DEXP(-OMU2*X(1)))
        CALL CHEB10(HE2P,RES,ZERO,X(1),NGC)
        Y2=Y2P+OM2*RES
   36  CONTINUE
C
       Y1=Y1/BETA/R
       Y2=Y2/(OMEGA+U2)
C
C     --- EXPONENTIAL ---
      ELSE IF(NP.EQ.23.OR.NP.EQ.24.OR.NP.EQ.26.OR.NP.EQ.27.OR.NP.EQ.32
     & .OR.NP.EQ.33) THEN
       CALL CHEB10(EX1P,RES,ZERO,X(1),NGC)
       Y1=RES/BETA/R
       Y2=YTRANS(X(3))*ZTRANS(X(4))*TLAME*DEXP(-ALP*X(2))
     &   *(1-DEXP(-OMU2*X(1))) 
       CALL CHEB10(EX2P,RES,ZERO,X(1),NGC)
       Y2=(Y2+OM2*RES)/(OMEGA+U2)         
C
C     --- spherical HEAVISIDE ---
      ELSE IF(NP.EQ.34)THEN           
       RS=SQRT((X(2)-X0P)**2+X(3)**2+X(4)**2)
       JP=1
       RI=XPP(JP)
       CALL CHEB11(H1S1P,RI,RES,ZERO,X(1),NGC)
       Y1=RES/BETA/R
       IF(RS.LT.XPP(JP)) Y2P=TLAM(JP)*(1-DEXP(-OMU2*X(1)))
       CALL CHEB11(H2SP,RI,RES,ZERO,X(1),NGC)
       Y2=(Y2P+OM2*RES)/(OMEGA+U2)
       DO 37 JP=2,NPP
        RI=XPP(JP)
        CALL CHEB11(H1S1P,RI,RES,ZERO,X(1),NGC)
        Y1=Y1+RES/BETA/R
        Y2P=0.
        IF(RS.LT.XPP(JP).AND.RS.GE.XPP(JP-1))
     &   Y2P=TLAM(JP)*(1-DEXP(-OMU2*X(1)))
        CALL CHEB11(H2SP,RI,RES,ZERO,X(1),NGC)
        Y2=Y2+(Y2P+OM2*RES)/(OMEGA+U2)
        RI=XPP(JP-1)
        CALL CHEB11(H1S1P,RI,RES,ZERO,X(1),NGC)
        Y1=Y1-RES/BETA/R
        CALL CHEB11(H2SP,RI,RES,ZERO,X(1),NGC)
        Y2=Y2-OM2*RES/(OMEGA+U2)
   37  CONTINUE
C     
C     --- spherical EXPONENTIAL ---
      ELSE IF(NP.EQ.35) THEN
       CALL CHEB10(E1S1P,RES,ZERO,X(1),NGC)
       Y1=RES/BETA/R
       CALL CHEB10(E2SP,RES,ZERO,X(1),NGC)
       Y2=TLAME*DEXP(-ALP*SQRT((X(2)-X0P)**2+X(3)**2+X(4)**2))
     &   *(1-DEXP(-OMU2*X(1)))+OM2*RES/(OMEGA+U2)
C     
      ENDIF       
  103 CONTINUE
C
      RETURN
      END
C                           
C     AUXILIARY FUNCTIONS FOR LONG. EQ. AND NONEQ. PROBLEMS
C     -----------------------------------------------------
      FUNCTION DI1B(TAU)
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/BVP/ NBP,T0B,DMB,GPB(10),TPB(10),GEB,ALB,JB  
      PARAMETER (ZIP=1.0D-10,TWO=2.0D+00)
C
      TT=X(1)-T0B 
      DUM=TT-TAU
      IF(DUM.LT.ZIP) DUM=ZIP
      P=TWO*DSQRT(DABS(OMOM*DUM*TAU))
      Q=-OM2*TAU-OMU2*DUM
      DI1B=DSQRT(OMOM*TAU/DUM)*EXPBI1(P,Q)*GAMY(TAU)*GAMZ(TAU)
     &*GB(X(2),TAU)
C
      RETURN
      END
C
C     ----
      FUNCTION DI2B(TAU)
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/BVP/ NBP,T0B,DMB,GPB(10),TPB(10),GEB,ALB,JB  
      PARAMETER (TWO=2.0D+00)
C
      TT=X(1)-T0B              
      P=TWO*DSQRT(DABS(OMOM*(TT-TAU)*TAU))
      Q=-OMU2*(TT-TAU)-OM2*TAU
      DI2B=EXPBI0(P,Q)*GB(X(2),TAU)*GAMY(TAU)*GAMZ(TAU)
      RETURN
      END
C
C     ----
      FUNCTION HE1B(TAU,TUP)      
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/BVP/ NBP,T0B,DMB,GPB(10),TPB(10),GEB,ALB,JB  
C
      P=OMEGA*OMEGA*TAU/(OMEGA+U2)/BETA/R
      Q=OMU2*(TUP-TAU)
      HE1B=GOLD(P,Q)*DEXP(-U2R*TAU)*GAMY(TAU)*GAMZ(TAU)*GB(X(2),TAU)
C     
      RETURN
      END
C
C     ----
      FUNCTION HE2B(TAU,TUP)
C                  
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/BVP/ NBP,T0B,DMB,GPB(10),TPB(10),GEB,ALB,JB  
      PARAMETER (ONE=1.0D+00)
C
      P=OMEGA*OMEGA*TAU/(OMEGA+U2)/BETA/R
      Q=OMU2*(TUP-TAU)                                              
      HE2B=(ONE-GOLD(Q,P))*DEXP(-U2R*TAU)*GAMY(TAU)*GAMZ(TAU)
     &*GB(X(2),TAU)
C
      RETURN
      END
C
C     ----
      FUNCTION EX1B(TAU)
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/BVP/ NBP,T0B,DMB,GPB(10),TPB(10),GEB,ALB,JB  
      PARAMETER (ZERO=0.0D+00,ZIP=1.0D-10,ONE=1.0D+00,TWO=2.0D+00)
C
      TT=X(1)-TAU                                     
      U=ALB-OMU2
      IF(ABS(U).LE.ZIP) U=ZIP
      EX1B=DEXP(-ALB*TT-OM2*TAU)*(DEXP(-OMOM*TAU/U)
     &    -DEXP(U*TT)*PHI1(TAU))*GAMY(TAU)*GAMZ(TAU)*GB(X(2),TAU)
C
      RETURN
      END
C
C     ----
      FUNCTION EX2B(TAU)
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/BVP/ NBP,T0B,DMB,GPB(10),TPB(10),GEB,ALB,JB  
      PARAMETER (ZERO=0.0D+00,ZIP=1.0D-10,ONE=1.0D+00,TWO=2.0D+00)
C
      TT=X(1)-TAU
      U=ALB-OMU2               
      IF(ABS(U).LE.ZIP) U=ZIP
      EX2B=DEXP(-OM2*TAU-OMU2*TT)*(DEXP(-OMOM*TAU/U)/U
     &    *(1-DEXP(-U*TT))-PHI2(TAU))*GAMY(TAU)*GAMZ(TAU)*GB(X(2),TAU)
C
      RETURN
      END
C
C     ----------------
      FUNCTION DI1I(TAU)
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/IVP/ NIP,X0I,DMI,FPI(10),XPI(10),FEI,ALI,JI
      PARAMETER (ZIP=1.0D-10,TWO=2.0D+00)
C     
      TT=X(1)-TAU
      IF(TT.LE.ZIP) TT=ZIP
      P=TWO*DSQRT(OMOM*TT*TAU)
      Q=-OM2*TAU-OMU2*TT
      DI1I=(OM2*EXPBI0(P,Q)+DSQRT(OMOM*TAU/TT)*EXPBI1(P,Q))
     $     *GAMY(TAU)*GAMZ(TAU)*G(X0I,TAU,X(2))
C
      RETURN
      END
C
C     ----
      FUNCTION DI2I(TAU)     
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/IVP/ NIP,X0I,DMI,FPI(10),XPI(10),FEI,ALI,JI
      PARAMETER (ONE=1.0D+00,TWO=2.0D+00,ZIP=1.0D-10)
C
      TT=X(1)-TAU
      IF(TT.LE.ZIP) TT=ZIP
      DUM=TAU
      IF(TAU.LE.ZIP) DUM=ZIP
      P=TWO*DSQRT(OMOM*TT*TAU)
      Q=-OM2*TAU-OMU2*TT
      DI2I=(EXPBI0(P,Q)+DSQRT(OM2*TT/OM1/DUM)*EXPBI1(P,Q))
     $     *GAMY(TAU)*GAMZ(TAU)*G(X0I,TAU,X(2))
C
      RETURN
      END
C
C     ----
      FUNCTION HE1I(TAU)
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/IVP/ NIP,X0I,DMI,FPI(10),XPI(10),FEI,ALI,JI
      PARAMETER (ZIP=1.0D-10,TWO=2.0D+00)
C
      TT=X(1)-TAU
      IF(TT.LE.ZIP) TT=ZIP
      P=TWO*DSQRT(OMOM*TT*TAU)
      Q=-OM2*TAU-OMU2*TT
      HE1I=(OM2*EXPBI0(P,Q)+DSQRT(OMOM*TAU/TT)*EXPBI1(P,Q))
     $*(G1(X(2),TAU,XPI(JI))-G1(X(2),TAU,XPI(JI+1)))*GAMY(TAU)*GAMZ(TAU)
C
      RETURN
      END
C
C     ----
      FUNCTION HE2I(TAU)
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/IVP/ NIP,X0I,DMI,FPI(10),XPI(10),FEI,ALI,JI
      PARAMETER (ZIP=1.0D-10,TWO=2.0D+00)
C
      TT=X(1)-TAU
      IF(TT.LE.ZIP) TT=ZIP                  
      DUM=TAU
      IF(TAU.LE.ZIP) DUM=ZIP
      P=TWO*DSQRT(OMOM*TT*TAU)
      Q=-OM2*TAU-OMU2*TT                    
      HE2I=(EXPBI0(P,Q)+DSQRT(OM2*TT/OM1/DUM)*EXPBI1(P,Q))
     $*(G1(X(2),TAU,XPI(JI))-G1(X(2),TAU,XPI(JI+1)))*GAMY(TAU)*GAMZ(TAU)
C
      RETURN
      END
C
C     ----
      FUNCTION EX1I(TAU)
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/IVP/ NIP,X0I,DMI,FPI(10),XPI(10),FEI,ALI,JI
      PARAMETER (ZERO=0.0D+00,ZIP=1.0D-10,ONE=1.0D+00,TWO=2.0D+00)
C
      TT=X(1)-TAU
      IF(TT.LE.ZIP) TT=ZIP
      DUM=TAU
      IF(TAU.LE.ZIP) DUM=ZIP
      P=TWO*DSQRT(OMOM*TT*TAU)
      Q=-OM2*TAU-OMU2*TT
      EX1I=(OM2*EXPBI0(P,Q)+DSQRT(OMOM*dum/TT)*EXPBI1(P,Q))
     $     *GAMY(TAU)*GAMZ(TAU)*G2(X(2),TAU,ALI)
C
      RETURN
      END
C
C     ----
      FUNCTION EX2I(TAU)
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/IVP/ NIP,X0I,DMI,FPI(10),XPI(10),FEI,ALI,JI
      PARAMETER (ZERO=0.0D+00,ZIP=1.0D-10,ONE=1.0D+00,TWO=2.0D+00)
C
      TT=X(1)-TAU
      IF(TT.LE.ZIP) TT=ZIP                  
      DUM=TAU
      IF(TAU.LE.ZIP) DUM=ZIP
      P=TWO*DSQRT(OMOM*TT*TAU)
      Q=-OM2*TAU-OMU2*TT
      EX2I=(EXPBI0(P,Q)+DSQRT(OM2*TT/OM1/DUM)*EXPBI1(P,Q))
     $     *GAMY(TAU)*GAMZ(TAU)*G2(X(2),TAU,ALI)
      RETURN
      END
C
C     ----
      FUNCTION H1S1I(XI,RI)
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/IVP/ NIP,X0I,DMI,FPI(10),XPI(10),FEI,ALI,JI
      PARAMETER (ZERO=0.0D+00,ZIP=1.0D-10,ONE=1.0D+00,TWO=2.0D+00)
C     
      H1S1I=GAM6(XI,X(3),X(4),X(1),RI,X0I)*G(XI,X(1),X(2))
C
      RETURN
      END
C
C     ----
      FUNCTION H1S2I(TAU,RI)
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/IVP/ NIP,X0I,DMI,FPI(10),XPI(10),FEI,ALI,JI
      PARAMETER (ZERO=0.0D+00,ZIP=1.0D-10,ONE=1.0D+00,TWO=2.0D+00)
      EXTERNAL H1S3I
C
      TT=X(1)-TAU
      IF(TT.LE.ZIP) TT=ZIP
      P=TWO*DSQRT(OMOM*TAU*TT)
      Q=-OM2*TAU-OMU2*TT
      BOT=X0I-RI
      TOP=X0I+RI
      CALL CHEB21(H1S3I,TAU,RI,RES,bot,top,NGC)
      H1S2I=(OM2*EXPBI0(P,Q)+DSQRT(OMOM*TAU/TT)*EXPBI1(P,Q))*RES
C
      RETURN
      END                       
C
C     ----
      FUNCTION H1S3I(TAU,XI,RI)
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/IVP/ NIP,X0I,DMI,FPI(10),XPI(10),FEI,ALI,JI
      PARAMETER (ZERO=0.0D+00,ZIP=1.0D-10,ONE=1.0D+00,TWO=2.0D+00)
C
      H1S3I=GAM6(XI,X(3),X(4),TAU,RI,X0I)*G(XI,TAU,X(2))      
C
      RETURN
      END
C
C     ----
      FUNCTION H2S1I(TAU,RI)
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/IVP/ NIP,X0I,DMI,FPI(10),XPI(10),FEI,ALI,JI
      PARAMETER (ZERO=0.0D+00,ZIP=1.0D-10,ONE=1.0D+00,TWO=2.0D+00)
      EXTERNAL H2S2I
C
      BOT=X0I-RI
      TOP=X0I+RI
      CALL CHEB21(H2S2I,TAU,RI,RES,bot,Top,NGC)
      H2S1I=RES
C
      RETURN
      END                       
C
C     ----
      FUNCTION H2S2I(TAU,XI,RI)
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/IVP/ NIP,X0I,DMI,FPI(10),XPI(10),FEI,ALI,JI
      PARAMETER (ZERO=0.0D+00,ZIP=1.0D-10,ONE=1.0D+00,TWO=2.0D+00)
C
      H2S2I=GAM6(XI,X(3),X(4),TAU,RI,X0I)*G(XI,TAU,X(2))      
C
      RETURN
      END
C                              
C     ----
      FUNCTION E1S1I(XI)
C      
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/IVP/ NIP,X0I,DMI,FPI(10),XPI(10),FEI,ALI,JI
      PARAMETER (ZERO=0.0D+00,ZIP=1.0D-10,ONE=1.0D+00,TWO=2.0D+00)
C
      E1S1I=dexp(-om2*x(1))
     $      *GAM7(XI,X(3),X(4),X(1),ALI,X0I)*G(XI,X(1),X(2))      
C
      RETURN
      END                                                           
C                                                                 
C     ----
      FUNCTION E1S2I(XI)
C      
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/IVP/ NIP,X0I,DMI,FPI(10),XPI(10),FEI,ALI,JI
      PARAMETER (ZERO=0.0D+00,ZIP=1.0D-10,ONE=1.0D+00,TWO=2.0D+00)
      EXTERNAL E1S3I
C
      CALL CHEB20(E1S3I,XI,RES,ZERO,X(1),NGC)
      E1S2I=RES
C
      RETURN
      END
C                              
C     ----
      FUNCTION E1S3I(XI,TAU)
C      
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/IVP/ NIP,X0I,DMI,FPI(10),XPI(10),FEI,ALI,JI
      PARAMETER (ZERO=0.0D+00,ZIP=1.0D-10,ONE=1.0D+00,TWO=2.0D+00)
C     
      TT=X(1)-TAU 
      P=TWO*DSQRT(OMOM*TAU*TT)
      Q=-OM2*TAU-OMU2*TT
      E1S3I=(OM2*EXPBI0(P,Q)+DSQRT(OMOM*TAU/TT)*EXPBI1(P,Q))
     &     *GAM7(XI,X(3),X(4),TAU,ALI,X0I)*G(XI,TAU,X(2))
C
      RETURN
      END
C
C     ----
      FUNCTION E2S1I(XI)
C      
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/IVP/ NIP,X0I,DMI,FPI(10),XPI(10),FEI,ALI,JI
      PARAMETER (ZERO=0.0D+00,ZIP=1.0D-10,ONE=1.0D+00,TWO=2.0D+00)
      EXTERNAL E2S2I
C
      CALL CHEB20(E2S2I,XI,RES,ZERO,X(1),NGC)
      E2S1I=RES
C
      RETURN
      END
C                              
C     ----
      FUNCTION E2S2I(XI,TAU)
C      
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/IVP/ NIP,X0I,DMI,FPI(10),XPI(10),FEI,ALI,JI
      PARAMETER (ZERO=0.0D+00,ZIP=1.0D-10,ONE=1.0D+00,TWO=2.0D+00)
C                                         
      IF(TT.LE.ZIP) TT=ZIP                  
      DUM=TAU
      TT=X(1)-TAU 
      P=TWO*DSQRT(OMOM*TAU*TT)
      Q=-OM2*TAU-OMU2*TT
      E2S2I=(EXPBI0(P,Q)+DSQRT(OM2*TT/OM1/DUM)*EXPBI1(P,Q))
     &     *GAM7(XI,X(3),X(4),TAU,ALI,X0I)*G(XI,TAU,X(2))
C
      RETURN
      END
C
C     ----
      FUNCTION HE1P(TAU)
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/PVP/ NPP,OLAM(10),TLAM(10),XPP(10),X0P,OLAME,TLAME,ALP,JP
      PARAMETER (ZERO=0.0D+00,ZIP=1.0D-10,ONE=1.0D+00,TWO=2.0D+00)
C
      DUM=TAU
      IF(TAU.LE.ZIP) DUM=ZIP
      TT=X(1)-TAU                          
      P=OMEGA*OMEGA*TAU/(OMEGA+U2)/BETA/R
      Q=OMU2*TT
C
      HE1P=DEXP(-U2R*TAU)*(OLAM(JP)*GOLD(P,Q)
     &    +OMEGA*TLAM(JP)*(ONE-GOLD(Q,P))/(OMEGA+U2))
     &    *(G1(X(2),TAU,XPP(JP))-G1(X(2),TAU,XPP(JP+1)))
     &    *GAMY(TAU)*GAMZ(TAU)
C
      RETURN
      END
C
C     ----
      FUNCTION HE2P(TAU)
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/PVP/ NPP,OLAM(10),TLAM(10),XPP(10),X0P,OLAME,TLAME,ALP,JP
      PARAMETER (ZERO=0.0D+00,ZIP=1.0D-10,ONE=1.0D+00,TWO=2.0D+00)
C
      TT=X(1)-TAU
      IF(TT.LE.ZIP) TT=ZIP                                
      P=OMEGA*OMEGA*TAU/(OMEGA+U2)/BETA/R
      Q=OMU2*TT
      ARG1=TWO*DSQRT(P*Q)
      ARG2=-P-Q
      HE2P=DEXP(-U2R*TAU)*(OLAM(JP)*(GOLD(Q,P)-EXPBI0(arg1,arg2))
     &    +OMEGA*TLAM(JP)*(1.-GOLD(Q,P)+SQRT(Q/P)*EXPBI1(ARG1,arg2))
     &    /(OMEGA+U2))
     &    *(G1(X(2),TAU,XPP(JP))-G1(X(2),TAU,XPP(JP+1)))
     &    *GAMY(TAU)*GAMZ(TAU)
C
      RETURN
      END
C
C     ----
      FUNCTION EX1P(TAU)
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/PVP/ NPP,OLAM(10),TLAM(10),XPP(10),X0P,OLAME,TLAME,ALP,JP
      PARAMETER (ZIP=1.0D-10,ONE=1.0D+00,TWO=2.0D+00)
C
      TT=X(1)-TAU 
      IF(TT.LE.ZIP) TT=ZIP                 
      P=OMEGA*OMEGA*TAU/(OMEGA+U2)/BETA/R
      Q=OMU2*TT
      EX1P=DEXP(-U2R*TAU)*(OLAME*GOLD(P,Q)
     &    +OMEGA*TLAME*(ONE-GOLD(Q,P))/(OMEGA+U2))
     &    *G2(X(2),TAU,ALP)*GAMY(TAU)*GAMZ(TAU)
C
      RETURN
      END
C
C     ----
      FUNCTION EX2P(TAU)
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/PVP/ NPP,OLAM(10),TLAM(10),XPP(10),X0P,OLAME,TLAME,ALP,JP
      PARAMETER (ZIP=1.0D-10,ZERO=0.0D+00,ONE=1.0D+00,TWO=2.0D+00)
C
      TT=X(1)-TAU
      IF(TT.LE.ZIP) TT=ZIP                                
      P=OMEGA*OMEGA*TAU/(OMEGA+U2)/BETA/R
      Q=OMU2*TT
      ARG1=TWO*DSQRT(P*Q)
      ARG2=-P-Q
      EX2P=DEXP(-U2R*TAU)*(OLAME*(GOLD(Q,P)-EXPBI0(arg1,arg2))
     &    +OMEGA*TLAME*(1.-GOLD(Q,P)+SQRT(Q/P)*EXPBI1(ARG1,arg2))
     &    /(OMEGA+U2))*G2(X(2),TAU,ALP)*GAMY(TAU)*GAMZ(TAU)
C
      RETURN
      END
C
C     ----
      FUNCTION H1S1P(TAU,RI)
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/PVP/ NPP,OLAM(10),TLAM(10),XPP(10),X0P,OLAME,TLAME,ALP,JP
      PARAMETER (ZIP=1.0D-10,ONE=1.0D+00,TWO=2.0D+00)
      EXTERNAL H1S2P
C
      TT=X(1)-TAU                          
      IF(TT.LE.ZIP) TT=ZIP                 
      P=OMEGA*OMEGA*TAU/(OMEGA+U2)/BETA/R
      Q=OMU2*TT
C                       
      BOT=X0P-RI
      TOP=X0P+RI
      CALL CHEB21(H1S2P,TAU,RI,RES,BOT,TOP,NGC)
C
      H1S1P=DEXP(-U2R*TAU)*(OLAM(JP)*GOLD(P,Q)
     &    +OMEGA*TLAM(JP)*(ONE-GOLD(Q,P))/(OMEGA+U2))*RES
C
      RETURN
      END
C
C     ----
      FUNCTION H1S2P(TAU,XI,RI)
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/PVP/ NPP,OLAM(10),TLAM(10),XPP(10),X0P,OLAME,TLAME,ALP,JP
      PARAMETER (ZERO=0.0D+00,ONE=1.0D+00,TWO=2.0D+00)
C
      H1S2P=GAM6(XI,X(3),X(4),TAU,RI,X0P)*G(XI,TAU,X(2))
C
      RETURN
      END
C
C     ----
      FUNCTION H2SP(TAU,RI)
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/PVP/ NPP,OLAM(10),TLAM(10),XPP(10),X0P,OLAME,TLAME,ALP,JP
      PARAMETER (ZERO=0.0D+00,ZIP=1.0D-10,ONE=1.0D+00,TWO=2.0D+00)
      EXTERNAL H1S2P
C                                              
      BOT=X0P-RI
      TOP=X0P+RI
      CALL CHEB21(H1S2P,TAU,RI,RES,BOT,TOP,NGC)
      TT=X(1)-TAU
      IF(TT.LE.ZIP) TT=ZIP                                
      P=OMEGA*OMEGA*TAU/(OMEGA+U2)/BETA/R
      Q=OMU2*TT
      ARG1=TWO*DSQRT(P*Q)
      ARG2=-P-Q
      H2SP=RES*DEXP(-U2R*TAU)*(OLAM(JP)*(GOLD(Q,P)-EXPBI0(arg1,arg2))
     &     +OMEGA*TLAM(JP)*(1.-GOLD(Q,P)+SQRT(Q/P)*EXPBI1(ARG1,arg2))
     &     /(OMEGA+U2))
C
      RETURN
      END
C
C     ----
      FUNCTION E1S1P(TAU)
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/PVP/ NPP,OLAM(10),TLAM(10),XPP(10),X0P,OLAME,TLAME,ALP,JP
      PARAMETER (ZERO=0.0D+00,ONE=1.0D+00)
      EXTERNAL E1S2P
C                   
      BOT=ZERO
      TOP=BINT
      CALL CHEB20(E1S2P,TAU,RES,BOT,TOP,NGC)
      TT=X(1)-TAU 
      P=OMEGA*OMEGA*TAU/(OMEGA+U2)/BETA/R
      Q=OMU2*TT
      E1S1P=DEXP(-U2R*TAU)*RES*(OLAME*GOLD(P,Q)
     &    +OMEGA*TLAME*(ONE-GOLD(Q,P))/(OMEGA+U2))
C
      RETURN
      END
C
C     ----
      FUNCTION E1S2P(TAU,XI)
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/PVP/ NPP,OLAM(10),TLAM(10),XPP(10),X0P,OLAME,TLAME,ALP,JP
      PARAMETER (ZERO=0.0D+00,ONE=1.0D+00)
C                   
      E1S2P=GAM7(XI,X(3),X(4),TAU,ALP,X0P)*G(XI,TAU,X(2))
C
      RETURN
      END
C
C     ----
      FUNCTION E2SP(TAU)
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/PVP/ NPP,OLAM(10),TLAM(10),XPP(10),X0P,OLAME,TLAME,ALP,JP
      PARAMETER (ZERO=0.0D+00,ZIP=1.0D-10,ONE=1.0D+00,TWO=2.0D+00)
      EXTERNAL E1S1P
C                                              
      BOT=ZERO
      TOP=BINT
      CALL CHEB20(E1S1P,TAU,RES,BOT,TOP,NGC)
      TT=X(1)-TAU
      IF(TT.LE.ZIP) TT=ZIP                                
      P=OMEGA*OMEGA*TAU/(OMEGA+U2)/BETA/R
      Q=OMU2*TT
      ARG1=TWO*DSQRT(P*Q)
      ARG2=-P-Q       
      E2SP=DEXP(-U2R*TAU)*RES*(OLAME*(GOLD(Q,P)-EXPBI0(arg1,arg2))
     &    +OMEGA*TLAME*(1.-GOLD(Q,P)+SQRT(Q/P)*EXPBI1(ARG1,arg2))
     &    /(OMEGA+U2))
C
      RETURN
      END
C
C     LONGITUDINAL AUXILIARY EXPRESSIONS (Table 2)
C     ----------------------------------
      REAL*8 FUNCTION GB(XX,T)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      PARAMETER(ZIP=1.0D-10,ONE=1.0D+00,TWO=2.0D+00,FOUR=4.0D+00)
      PI=FOUR*DATAN(ONE)
C
      DUM=T
      IF(DUM.LT.ZIP) DUM=ZIP
      IF (MC.EQ.1) THEN
       GB=XX*DSQRT(BETA*R*PX/FOUR/PI/DUM/DUM/DUM)
     $   *DEXP(-U1R*T-PX*(BETA*R*XX-T)**2/FOUR/BETA/R/DUM)
      ELSE
       GB=DEXP(-U1R*T)*(DSQRT(PX/PI/BETA/R/DUM)
     $   *DEXP(-PX*(BETA*R*XX-T)**2/FOUR/BETA/R/DUM)
     $   -PX/TWO/BETA/R*DEXP(PX*XX)*ERFC((BETA*R*XX+T)/
     $    DSQRT(FOUR*BETA*R*DUM/PX)))
      ENDIF
C
      RETURN
      END
C
C     ------
      REAL*8 FUNCTION G(XI,T,XX)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      PARAMETER(ZIP=1.0D-10,HALF=0.5D+00,ONE=1.0D+00,TWO=2.0D+00,
     $FOUR=4.0D+00)
      PI=FOUR*DATAN(ONE)
C
      DUM=T
      IF(DUM.LT.ZIP) DUM=ZIP
C   
      IF (MC.EQ.1) THEN
       G=DEXP(-U1R*T)
     &    *DSQRT(BETA*R*PX/FOUR/PI/DUM)*((ONE-(BETA*R*(XI-XX)+T)/
     $    (TWO*DUM))*DEXP(-PX*(BETA*R*(XI-XX)+T)**2/FOUR/BETA/R/DUM)
     $    -(ONE-(BETA*R*(XI+XX)+T)/TWO/DUM)*DEXP(PX*XX-PX*(BETA*R
     $    *(XI+XX)+T)**2/FOUR/BETA/R/DUM))
      ELSE
       G=(DSQRT(BETA*R*PX/FOUR/PI/DUM)*DEXP(PX*XX-PX*(BETA*R*(XI+
     $    XX)+T)**2/FOUR/BETA/R/DUM)-HALF*PX*DEXP(PX*XX)*ERFC((BETA*
     $    R*(XI+XX)+T)/DSQRT(FOUR*BETA*R*DUM/PX))+DSQRT(BETA*R*PX/
     $    FOUR/PI/DUM)*DEXP(-PX*(BETA*R*(XI-XX)+T)**2/FOUR/BETA/R/
     $    DUM))*DEXP(-U1R*T)
      ENDIF
C
      RETURN
      END
C
C     ------
      REAL*8 FUNCTION G1(XX,T,XXI)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      PARAMETER(ZIP=1.0D-10,HALF=0.5D+00,ONE=1.0D+00,FOUR=4.0D+00)
      PI=FOUR*DATAN(ONE)
C
      DUM=T
      IF(DUM.LT.ZIP) DUM=ZIP
C   
      IF (MC.EQ.1) THEN
       G1=DEXP(-U1R*T)*(ONE-HALF*ERFC((BETA*R*(XX-XXI)-T)/DSQRT(FOUR
     $  *BETA*R*DUM/PX))-HALF*DEXP(PX*XX)*ERFC((BETA*R*(XX+XXI)+T)
     $  /DSQRT(FOUR*BETA*R*DUM/PX))-DSQRT(BETA*R/FOUR/PI/DUM/PX)
     $  *DEXP(-PX*(BETA*R*(XX-XXI)-T)**2/FOUR/BETA/R/DUM)
     $  +DSQRT(BETA*R/FOUR/PI/DUM/PX)*DEXP(PX*XX)
     $  *DEXP(-PX*(BETA*R*(XX+XXI)+T)**2/FOUR/BETA/R/DUM))
      ELSE
       G1=(ONE-HALF*ERFC((BETA*R*(XX-XXI)-T)/DSQRT(FOUR*BETA*R*DUM/PX))
     $    +HALF*(ONE+PX*(XX+XXI)+PX*T/BETA/R)*DEXP(PX*XX)*ERFC((
     $    BETA*R*(XX+XXI)+T)/DSQRT(FOUR*BETA*R*DUM/PX))-DSQRT(PX*T/PI
     $    /BETA/R)*DEXP(PX*XX-PX*(BETA*R*(XX+XXI)+T)**2/FOUR/BETA
     $    /R/DUM))*DEXP(-U1R*T)
      ENDIF
C
      RETURN
      END
C
C     ------
      REAL*8 FUNCTION G2(XX,T,AL)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      PARAMETER(ZIP=1.0D-10,HALF=0.5D+00,ONE=1.0D+00,TWO=2.0D+00,
     $FOUR=4.0D+00)
      PI=FOUR*DATAN(ONE)
C
      DUM=T
      IF(DUM.LT.ZIP) DUM=ZIP
C   
      IF (MC.EQ.1) THEN
       G2=(ONE+AL/PX)*DEXP((AL+AL**2/PX-U1)*T/BETA/R)*
     $   (DEXP(-AL*XX)*(ONE-HALF*ERFC((BETA*R*XX-T*(ONE+TWO*AL/
     $   PX))/DSQRT(FOUR*BETA*R*DUM/PX)))+HALF*DEXP((AL+PX)*XX)*ERFC
     $   ((BETA*R*XX+T*(ONE+TWO*AL/PX))/DSQRT(FOUR*BETA*R*DUM/PX)))
      ELSE
       G2=DEXP((AL+AL**2/PX-U1)*T/(BETA*R))*(DEXP(-AL*XX)*(ONE-
     $   HALF*ERFC((BETA*R*XX-T*(ONE+TWO*AL/PX))/DSQRT(FOUR*BETA*R*
     $   DUM/PX)))+HALF*(ONE+PX/AL)*DEXP((AL+PX)*XX)*ERFC((BETA*R*
     $   XX+T*(ONE+TWO*AL/PX))/DSQRT(FOUR*BETA*R*DUM/PX)))-(HALF*PX/
     $   AL)*EXP(PX*XX-U1*T/BETA/R)*ERFC((BETA*R*XX+T)/DSQRT(FOUR*BETA*R
     $   *DUM/PX))
      ENDIF
C
      RETURN
      END
C
C     TRANSVERSE AUXILIARY EXPRESSIONS (Table 3)
C     --------------------------------
      REAL*8 FUNCTION GAMY(T)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/TRANS/ Y0,Z0,NYP,YP(5),NZP,ZP(5),ALY,ALZ
      PARAMETER (ZIP=1.0D-10,ONE=1.0D+00,FOUR=4.0D+00)
      EXTERNAL GAM4,GAM5
C
      PI=FOUR*DATAN(ONE)
C
      IF(T.LT.ZIP) T=ZIP                     
C
C     DIRAC
      IF(NP.EQ.1.OR.NP.EQ.2.OR.NP.EQ.3.OR.NP.EQ.10.OR.NP.EQ.11.OR.
     &   NP.EQ.12.OR.NP.EQ.19.OR.NP.EQ.20.OR.NP.EQ.21) THEN       
       GAMY=DSQRT(BETA*R*PY/FOUR/PI/T)*DEXP(-BETA*R*PY*(X(3)-Y0)**2
     &   /FOUR/T)
C
C     HEAVISIDE
      ELSE IF (NP.EQ.4.OR.NP.EQ.5.OR.NP.EQ.6.OR.NP.EQ.13.OR.NP.EQ.14.OR.
     &   NP.EQ.15.OR.NP.EQ.22.OR.NP.EQ.23.OR.NP.EQ.24) THEN
       SUMY=0.
       DO 11 I=1,2*NYP,2
        SUMY=SUMY+0.5*ERFC(DSQRT(BETA*R*PY/FOUR/T)*(X(3)-YP(I+1)))
     &      -0.5*ERFC(DSQRT(BETA*R*PY/FOUR/T)*(X(3)-YP(I))) 
   11  CONTINUE
       GAMY=SUMY
C
C     EXPONENTIAL
      ELSE IF (NP.EQ.7.OR.NP.EQ.8.OR.NP.EQ.9.OR.NP.EQ.16.OR.NP.EQ.17.OR.
     &   NP.EQ.18.OR.NP.EQ.25.OR.NP.EQ.26.OR.NP.EQ.27) THEN
       GAMY=0.5*DEXP(ALY*(Y0-X(3))+ALY*ALY*T/BETA/R/PY)
     &     *ERFC((BETA*R*PY*(Y0-X(3))+2*ALY*T)/sqrt(4.*BETA*R*PY*T))
     &     +0.5*DEXP(ALY*(X(3)-Y0)+ALY*ALY*T/BETA/R/PY)
     &     *ERFC((BETA*R*PY*(X(3)-Y0)+2*ALY*T)/sqrt(4.*BETA*R*PY*T))
C        
C     circular HEAVISIDE
      ELSE IF (NP.EQ.28.OR.NP.EQ.30.OR.NP.EQ.32) THEN
       SUMY=0.
       DO 12 I=1,2*NYP,2
        IF(YP(I).LE.ZIP) THEN
         RES=0.
      ELSE
         BOT=X(3)-YP(I)
         TOP=X(3)+YP(I)
         CALL CHEB12(GAM4,T,YP(I),RES,BOT,TOP,NGC)
        ENDIF  
        SUMY=RES+SUMY
        BOT=X(3)-YP(I+1)
        TOP=X(3)+YP(I+1)
        CALL CHEB12(GAM4,T,YP(I+1),RES,BOT,TOP,NGC)
        SUMY=RES-SUMY
   12  CONTINUE
       GAMY=SUMY
C
C     circular EXPONENTIAL 
      ELSE IF (NP.EQ.29.OR.NP.EQ.31.OR.NP.EQ.33) THEN
       BOT=-BINT
       TOP=BINT
       CALL CHEB11(GAM5,T,RES,BOT,TOP,NGC)
       GAMY=RES
      ENDIF
C
      RETURN
      END
C                                        
C     ------------
      REAL*8 FUNCTION GAMZ(T)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/TRANS/ Y0,Z0,NYP,YP(5),NZP,ZP(5),ALY,ALZ
      PARAMETER (ZIP=1.0D-10,ONE=1.0D+00,FOUR=4.0D+00)
      PI=FOUR*DATAN(ONE)
C
      IF(T.LT.ZIP) T=ZIP                     
C
C     DIRAC
      IF(NP.EQ.1.OR.NP.EQ.10.OR.NP.EQ.19.OR.NP.EQ.4.OR.NP.EQ.13.OR.
     &   NP.EQ.22.OR.NP.EQ.7.OR.NP.EQ.16.OR.NP.EQ.25) THEN       
       GAMZ=DSQRT(BETA*R*PZ/FOUR/PI/T)*DEXP(-BETA*R*PZ*(X(4)-Z0)**2
     &   /FOUR/T)
C
C     HEAVISIDE
      ELSE IF(NP.EQ.2.OR.NP.EQ.11.OR.NP.EQ.20.OR.NP.EQ.5.OR.NP.EQ.14.OR.
     &   NP.EQ.23.OR.NP.EQ.8.OR.NP.EQ.17.OR.NP.EQ.26) THEN       
       SUMZ=0.
       DO 11 I=1,2*NZP,2
        SUMZ=SUMZ+0.5*ERFC(DSQRT(BETA*R*PZ/FOUR/T)*(X(4)-ZP(I+1)))
     &      -0.5*ERFC(DSQRT(BETA*R*PZ/FOUR/T)*(X(4)-ZP(I))) 
   11  CONTINUE
       GAMZ=SUMZ
C
C     EXPONENTIAL
      ELSE IF(NP.EQ.3.OR.NP.EQ.12.OR.NP.EQ.21.OR.NP.EQ.6.OR.NP.EQ.15.OR.
     &   NP.EQ.24.OR.NP.EQ.9.OR.NP.EQ.18.OR.NP.EQ.27) THEN       
       GAMZ=0.5*DEXP(ALZ*(Z0-X(4))+ALZ*ALZ*T/BETA/R/PZ)
     &     *ERFC((BETA*R*PZ*(Z0-X(4))+2*ALZ*T)/sqrt(4.*BETA*R*PZ*T))
     &     +0.5*DEXP(ALZ*(X(4)-Z0)+ALZ*ALZ*T/BETA/R/PZ)
     &     *ERFC((BETA*R*PZ*(X(4)-Z0)+2*ALZ*T)/sqrt(4.*BETA*R*PZ*T))
C        
C     circular
      ELSE IF (NP.GE.28) THEN
       GAMZ=1.0
      ENDIF
C
      RETURN
      END
C                                        
C     functions
C     ---------
      REAL*8 FUNCTION YTRANS(Y)
C
C     PURPOSE: DETERMINE DISTRIBUTION IN Y- OR RADIAL DIRECTION
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/TRANS/ Y0,Z0,NYP,YP(5),NZP,ZP(5),ALY,ALZ
      PARAMETER (ZERO=0.0D+00)
C 
      YTRANS=0.
C
C     HEAVISIDE
      IF (NP.EQ.4.OR.NP.EQ.5.OR.NP.EQ.6.OR.NP.EQ.13.OR.NP.EQ.14.OR.
     &NP.EQ.15.OR.NP.EQ.22.OR.NP.EQ.23.OR.NP.EQ.24
     &.OR. NP.EQ.28 .OR. NP.EQ.30 .OR. NP.EQ.32) THEN
      DO 11 I=1,2*NYP,2
       IF(Y.GE.YP(I) .AND. Y.LE.YP(I+1)) YTRANS=1.
   11 CONTINUE
C
C     EXPONENTIAL
      ELSE IF (NP.EQ.7.OR.NP.EQ.8.OR.NP.EQ.9.OR.NP.EQ.16.OR.NP.EQ.17.OR.
     &NP.EQ.18.OR.NP.EQ.25.OR.NP.EQ.26.OR.NP.EQ.27) THEN
       YTRANS=DEXP(-ALY*DABS(Y-Y0))
C
C     circular EXPONENTIAL
      ELSE IF (NP.EQ.29 .OR. NP.EQ.31 .OR. NP.EQ.33) THEN
       YTRANS=DEXP(-ALY*Y)
      ENDIF
C       
      RETURN
      END
C
C     ---------
      FUNCTION ZTRANS(Z)
C
C     PURPOSE: DETERMINE DISTRIBUTION IN Z-DIRECTION
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/TRANS/ Y0,Z0,NYP,YP(5),NZP,ZP(5),ALY,ALZ
      PARAMETER (ZERO=0.0D+00)
C 
      ZTRANS=0.
C
C     HEAVISIDE
      IF (NP.EQ.2.OR.NP.EQ.11.OR.NP.EQ.20.OR.NP.EQ.5.OR.NP.EQ.14.OR.
     & NP.EQ.23.OR.NP.EQ.8.OR.NP.EQ.17.OR.NP.EQ.26) THEN
       DO 11 I=1,2*NZP,2
        IF(Z.GE.ZP(I) .AND. Z.LE.ZP(I+1)) ZTRANS=1.
   11  CONTINUE
C
C     circular HEAVISIDE
      ELSE IF (NP.EQ.28 .OR. NP.EQ.30 .OR. NP.EQ.32) THEN
       ZTRANS=1.
C
C     EXPONENTIAL
      ELSE IF (NP.EQ.3.OR.NP.EQ.12.OR.NP.EQ.21.OR.NP.EQ.6.OR.NP.EQ.15
     & .OR.NP.EQ.24.OR.NP.EQ.9.OR.NP.EQ.18.OR.NP.EQ.27) THEN
       ZTRANS=DEXP(-ALZ*DABS(Z-Z0))
C
C     circular EXPONENTIAL
      ELSE IF (NP.EQ.29 .OR. NP.EQ.31 .OR. NP.EQ.33) THEN
       ZTRANS=1.
      ENDIF
C
      RETURN
      END
C
C     ---------
      REAL*8 FUNCTION GAM4(PHI,T,RI)
C 
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/TRANS/ Y0,Z0,NYP,YP(5),NZP,ZP(5),ALY,ALZ
      PARAMETER(ZIP=1.0D-10,ONE=1.0D+00,TWO=2.0D+00,FOUR=4.0D+00)
      PI=FOUR*DATAN(ONE)
C
      IF(T.LT.ZIP) T=ZIP
C
      GAM4=DSQRT(BETA*R*PY/PI/T)*DEXP(-BETA*R*PY*PHI*PHI/FOUR/T)
     $      *(ERFC(DSQRT(BETA*R*PZ/FOUR/T)*(X(4)-DSQRT(DABS(RI*RI-(X(3)
     $      -PHI)**2))))-ERFC(DSQRT(BETA*R*PZ/FOUR/T)*(X(4)+DSQRT(DABS
     $      (RI*RI-(X(3)-PHI)**2)))))/FOUR
C
      RETURN
      END
C
C     ----------
      REAL*8 FUNCTION GAM5(PHI,T)
C 
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/TRANS/ Y0,Z0,NYP,YP(5),NZP,ZP(5),ALY,ALZ
      PARAMETER(ONE=1.0D+00,FOUR=4.0D+00)
      EXTERNAL GAM51
C      
      PI=FOUR*DATAN(ONE)
C
      BOT=-BINT
      TOP=BINT
C 
      CALL CHEB21(GAM51,PHI,T,RES,BOT,TOP,NGC)
      GAM5=RES
C
      RETURN
      END
C
C     ----------
      REAL*8 FUNCTION GAM51(PHI,OME,T)
C 
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/TRANS/ Y0,Z0,NYP,YP(5),NZP,ZP(5),ALY,ALZ
      PARAMETER(ZIP=1.0D-10,ONE=1.0D+00,TWO=2.0D+00,FOUR=4.0D+00)
      PI=FOUR*DATAN(ONE)
C
      IF(T.LT.ZIP) T=ZIP
C
      GAM51=(BETA*R*DSQRT(PY*PZ)/PI/T/FOUR)
     & *DEXP(-ALY*DSQRT((X(3)-PHI)**2+(X(4)-OME)**2)
     & -BETA*R*(PY*PHI*PHI+PZ*OME*OME)/FOUR/T)
C 
      RETURN
      END
C
C     ----------
      REAL*8 FUNCTION GAM6(XI,Y,Z,T,RI,X0)
C 
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/TRANS/ Y0,Z0,NYP,YP(5),NZP,ZP(5),ALY,ALZ
      PARAMETER(ZIP=1.0D-10,ONE=1.0D+00,TWO=2.0D+00,FOUR=4.0D+00)
      EXTERNAL GAM61
      PI=FOUR*DATAN(ONE)
C     
      IF(T.LE.ZIP) T=ZIP
      BOT=X(3)-SQRT(RI*RI-(XI-X0)**2)
      TOP=X(3)+SQRT(RI*RI-(XI-X0)**2)
C 
      CALL CHEB25(GAM61,XI,Y,Z,T,RI,X0,RES,BOT,TOP,NGC)
      GAM6=RES
C
      RETURN
      END
C             
C     ----------
      REAL*8 FUNCTION GAM61(XI,PHI,Y,Z,T,RI,X0)
C 
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/TRANS/ Y0,Z0,NYP,YP(5),NZP,ZP(5),ALY,ALZ
      PARAMETER(ZIP=1.0D-10,ONE=1.0D+00,TWO=2.0D+00,FOUR=4.0D+00)
      PI=FOUR*DATAN(ONE)
C
      IF(T.LE.ZIP) T=ZIP
      GAM61=DSQRT(BETA*R*PY/PI/T)*DEXP(-BETA*R*PY*PHI*PHI/FOUR/T)
     $      *(ERFC(DSQRT(BETA*R*PZ/FOUR/T)*(Z-DSQRT(DABS(RI*RI
     $      -(XI-X0)**2-(Y-PHI)**2))))-ERFC(DSQRT(BETA*R*PZ/FOUR/T)
     $      *(Z+DSQRT(DABS(RI*RI-(XI-X0)**2-(Y-PHI)**2)))))/FOUR
C 
      RETURN
      END
C
C     ----------
      REAL*8 FUNCTION GAM7(XI,Y,Z,T,AL,X0)
C 
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/TRANS/ Y0,Z0,NYP,YP(5),NZP,ZP(5),ALY,ALZ
      PARAMETER(ZERO=0.0D+00,ONE=1.0D+00,FOUR=4.0D+00)
      EXTERNAL GAM71
C
      PI=FOUR*DATAN(ONE)          
c
      GAM7=0.  
      BOT=-BINT
      TOP=bint
C 
      CALL CHEB25(GAM71,XI,Y,Z,T,AL,X0,RES,BOT,TOP,NGC)
      GAM7=RES+GAM7
C
      RETURN
      END
C
C     ----------
      REAL*8 FUNCTION GAM71(XI,PHI,Y,Z,T,AL,X0)
C 
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/TRANS/ Y0,Z0,NYP,YP(5),NZP,ZP(5),ALY,ALZ
      PARAMETER(ZERO=0.0D+00,ONE=1.0D+00,FOUR=4.0D+00)
      EXTERNAL GAM72
C
      GAM71=0.  
      BOT=-BINT
      TOP=bint
C 
      CALL CHEB35(GAM72,XI,PHI,Y,Z,T,AL,X0,RES,BOT,TOP,NGC)
      GAM71=RES+GAM71
C
      RETURN
      END
C
C     ----------
      REAL*8 FUNCTION GAM72(XI,PHI,OME,Y,Z,T,AL,X0)
C 
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/TRANS/ Y0,Z0,NYP,YP(5),NZP,ZP(5),ALY,ALZ
      PARAMETER(ZIP=1.0D-10,ONE=1.0D+00,TWO=2.0D+00,FOUR=4.0D+00)
      PI=FOUR*DATAN(ONE)
c
      if(t.lt.zip) t=zip
      GAM72=(BETA*R*DSQRT(PY*PZ)/PI/t/FOUR)
     & *DEXP(-AL*DSQRT((XI-X0)**2+(Y-PHI)**2+(Z-OME)**2)
     & -BETA*R*(PY*PHI*PHI+PZ*OME*OME)/FOUR/t)
C 
      RETURN
      END
      
C
C     INTEGRATION ROUTINES
C     --------------------
      SUBROUTINE CHEB10(FUNC,AREA,A,B,M)
C
C     PURPOSE: PERFORM INTEGRATION OF F(X) BETWEEN A AND B
C              USING M-POINT GAUSS-CHEBYSHEV QUADRATURE FORMULA
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(ZERO=0.0D+00,ONE=1.0D+00,TWO=2.0D+00,FOUR=4.0D+00)
      PI=FOUR*DATAN(ONE)
      AREA=ZERO
      SUM=ZERO
      DO 10 I=1,M
        Z1=DCOS(DFLOAT(2*(I-1)+1)*PI/DFLOAT(2*M))
        VAR1=(Z1*(B-A)+B+A)/TWO
   10 SUM=SUM+FUNC(VAR1)*DSQRT(ONE-Z1*Z1)
      AREA=(B-A)*PI*SUM/DFLOAT(2*M)
C
      RETURN
      END
C
C     --------------------
      SUBROUTINE CHEBHB(FUNC,AREA,A,B,M)
C
C     PURPOSE: PERFORM INTEGRATION OF F(X) BETWEEN A AND B
C              USING M-POINT GAUSS-CHEBYSHEV QUADRATURE FORMULA
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(ZERO=0.0D+00,ONE=1.0D+00,TWO=2.0D+00,FOUR=4.0D+00)
      PI=FOUR*DATAN(ONE)
      AREA=ZERO
      SUM=ZERO
      DO 10 I=1,M
        Z1=DCOS(DFLOAT(2*(I-1)+1)*PI/DFLOAT(2*M))
        VAR1=(Z1*(B-A)+B+A)/TWO
   10 SUM=SUM+FUNC(VAR1,B)*DSQRT(ONE-Z1*Z1)
      AREA=(B-A)*PI*SUM/DFLOAT(2*M)
C
      RETURN
      END
C
C     ----------
      SUBROUTINE CHEB11(FUNC,CON1,AREA,A,B,M)
C
C     PURPOSE: PERFORM INTEGRATION OF F(X) BETWEEN A AND B
C              USING M-POINT GAUSS-CHEBYSHEV QUADRATURE FORMULA
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(ZERO=0.0D+00,ONE=1.0D+00,TWO=2.0D+00,FOUR=4.0D+00)
      PI=FOUR*DATAN(ONE)
      AREA=ZERO
      SUM=ZERO
      DO 10 I=1,M
        Z1=DCOS(DFLOAT(2*(I-1)+1)*PI/DFLOAT(2*M))
        VAR1=(Z1*(B-A)+B+A)/TWO
   10 SUM=SUM+FUNC(VAR1,CON1)*DSQRT(ONE-Z1*Z1)
      AREA=(B-A)*PI*SUM/DFLOAT(2*M)
C
      RETURN
      END
C
C     ----------
      SUBROUTINE CHEB12(FUNC,CON1,CON2,AREA,A,B,M)
C
C     PURPOSE: PERFORM INTEGRATION OF F(X) BETWEEN A AND B
C              USING M-POINT GAUSS-CHEBYSHEV QUADRATURE FORMULA
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(ZERO=0.0D+00,ONE=1.0D+00,TWO=2.0D+00,FOUR=4.0D+00)
      PI=FOUR*DATAN(ONE)
      AREA=ZERO
      SUM=ZERO
      DO 10 I=1,M
        Z1=DCOS(DFLOAT(2*(I-1)+1)*PI/DFLOAT(2*M))
        VAR1=(Z1*(B-A)+B+A)/TWO
   10 SUM=SUM+FUNC(VAR1,CON1,CON2)*DSQRT(ONE-Z1*Z1)
      AREA=(B-A)*PI*SUM/DFLOAT(2*M)
C
      RETURN
      END
C
C     ----------
      SUBROUTINE CHEB20(FUNC,VAR1,AREA,A,B,M)
C
C     PURPOSE: PERFORM INTEGRATION OF F(X) BETWEEN A AND B
C              USING M-POINT GAUSS-CHEBYSHEV QUADRATURE FORMULA
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(ZERO=0.0D+00,ONE=1.0D+00,TWO=2.0D+00,FOUR=4.0D+00)
      PI=FOUR*DATAN(ONE)
      AREA=ZERO
      SUM=ZERO
      DO 10 I=1,M
        Z1=DCOS(DFLOAT(2*(I-1)+1)*PI/DFLOAT(2*M))
        VAR2=(Z1*(B-A)+B+A)/TWO
   10 SUM=SUM+FUNC(VAR1,VAR2)*DSQRT(ONE-Z1*Z1)
      AREA=(B-A)*PI*SUM/DFLOAT(2*M)
C
      RETURN
      END       
C
C    ----------
      SUBROUTINE CHEB21(FUNC,VAR1,CON1,AREA,A,B,M)
C
C     PURPOSE: PERFORM INTEGRATION OF F(X) BETWEEN A AND B
C              USING M-POINT GAUSS-CHEBYSHEV QUADRATURE FORMULA
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(ZERO=0.0D+00,ONE=1.0D+00,TWO=2.0D+00,FOUR=4.0D+00)
      PI=FOUR*DATAN(ONE)
      AREA=ZERO
      SUM=ZERO
      DO 10 I=1,M
        Z1=DCOS(DFLOAT(2*(I-1)+1)*PI/DFLOAT(2*M))
        VAR2=(Z1*(B-A)+B+A)/TWO
   10 SUM=SUM+FUNC(VAR1,VAR2,CON1)*DSQRT(ONE-Z1*Z1)
      AREA=(B-A)*PI*SUM/DFLOAT(2*M)
C
      RETURN
      END       
C                                        
C     ----------
      SUBROUTINE CHEB25(FUNC,VAR1,CON1,CON2,CON3,CON4,CON5,AREA,A,B,M)
C
C     PURPOSE: PERFORM INTEGRATION OF F(X) BETWEEN A AND B
C              USING M-POINT GAUSS-CHEBYSHEV QUADRATURE FORMULA
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(ZERO=0.0D+00,ONE=1.0D+00,TWO=2.0D+00,FOUR=4.0D+00)
      PI=FOUR*DATAN(ONE)
      AREA=ZERO
      SUM=ZERO
      DO 10 I=1,M
        Z1=DCOS(DFLOAT(2*(I-1)+1)*PI/DFLOAT(2*M))
        VAR2=(Z1*(B-A)+B+A)/TWO
   10 SUM=SUM+FUNC(VAR1,VAR2,CON1,CON2,CON3,CON4,CON5)*DSQRT(ONE-Z1*Z1)
      AREA=(B-A)*PI*SUM/DFLOAT(2*M)
C
      RETURN
      END       
C                                        
C     ----------
      SUBROUTINE CHEB35(FUNC,VAR1,VAR2,CON1,CON2,CON3,CON4,CON5,AREA,A,B
     &,M)
C
C     PURPOSE: PERFORM INTEGRATION OF F(X) BETWEEN A AND B
C              USING M-POINT GAUSS-CHEBYSHEV QUADRATURE FORMULA
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(ZERO=0.0D+00,ONE=1.0D+00,TWO=2.0D+00,FOUR=4.0D+00)
      PI=FOUR*DATAN(ONE)
      AREA=ZERO
      SUM=ZERO
      DO 10 I=1,M
        Z1=DCOS(DFLOAT(2*(I-1)+1)*PI/DFLOAT(2*M))
        VAR3=(Z1*(B-A)+B+A)/TWO
   10 SUM=SUM+FUNC(VAR1,VAR2,VAR3,CON1,CON2,CON3,CON4,CON5)
     & *DSQRT(ONE-Z1*Z1)
      AREA=(B-A)*PI*SUM/DFLOAT(2*M)
C
      RETURN
      END       
      FUNCTION ERFC(A)
C
C     ----------
C     PURPOSE: TO CALCULATE ERFC(A)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      ERFC=0.0
      C=-A*A
      IF((DABS(C).GT.170.).AND.(A.GE.0.)) RETURN
      IF(C.LT.-170.) GO TO 3
      X=DABS(A)
      IF(X.GT.3.0) GO TO 1
      T=1./(1.+.3275911*X)
      Y=T*(.2548296-T*(.2844967-T*(1.421414-T*(1.453152-1.061405*T))))
      GO TO 2
    1 Y=.5641896/(X+.5/(X+1./(X+1.5/(X+2./(X+2.5/(X+1.))))))
    2 ERFC=Y*DEXP(C)
    3 IF(A.LT.0.0) ERFC=2.-ERFC
      RETURN
      END
C
C     ----
      FUNCTION EXPBI0(X,Z)
      IMPLICIT REAL*8(A-Z)
C
C     PURPOSE: TO CALCULATE PRODUCT OF EXP(Z) AND I0(X) (ZERO-ORDER
C              MODIFIED BESSEL FUNCTION OF THE FIRST KIND)
C
      DATA P1,P2,P3,P4,P5,P6,P7/1.0D0,3.5156229D0,3.0899424D0,
     $     1.2067492D0,0.2659732D0,0.360768D-1,0.45813D-2/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1,
     $     0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1,
     $     0.2635537D-1,-0.1647633D-1,0.392377D-2/
C
      IF (DABS(X).LT.3.75) THEN
        Y=(X/3.75)**2
        EXPBI0=DEXP(Z)*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
        AX=DABS(X)
        Y=3.75/AX 
       EXPBI0=DEXP(AX+Z)/DSQRT(AX)*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*
     $         (Q6+Y*(Q7+Y*(Q8+Y*Q9))))))))
      ENDIF
      RETURN
      END
C
C     ----
      FUNCTION EXPBI1(X,Z)
      IMPLICIT REAL*8(A-Z)
C
C     PURPOSE: TO CALCULATE PRODUCT OF EXP(Z) AND I1(X) (FIRST-ORDER
C              MODIFIED BESSEL FUNCTION OF THE FIRST KIND)
C
      DATA P1,P2,P3,P4,P5,P6,P7/0.5D0,0.87890594D0,0.51498869D0,
     $     0.15084984D0,0.2658733D-1,0.301532D-2,0.32411D-3/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,-0.3988024D-1,
     $     -0.3662018D-2,0.163801D-2,-0.1031555D-1,0.2282967D-1,
     $     -0.2895312D-1,0.1787654D-1,-0.420059D-2/
C
      IF (DABS(X).LT.3.75) THEN
        Y=(X/3.75)**2
        EXPBI1=DEXP(Z)*X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
        AX=DABS(X)
        Y=3.75/AX 
        EXPBI1=DEXP(AX+Z)/DSQRT(AX)*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*
     $         (Q6+Y*(Q7+Y*(Q8+Y*Q9))))))))
        IF(X.LT.0.) EXPBI1=-EXPBI1
      ENDIF
      RETURN
      END
C
C     ---
      FUNCTION GOLD(X,Y)
C
C     PURPOSE: TO CALCULATE GOLDSTEIN'S J-FUNCTION J(X,Y)
C
      IMPLICIT REAL*8(A-H,O-Z)
      GOLD=0.0
      BF=0.0
      E=2.*DSQRT(DMAX1(1.D-35,X*Y))
      Z=X+Y-E
      IF(Z.GT.17.) GO TO 8
      IF(E.GT.1.D-15) GO TO 2
      GOLD=DEXP(-X)
      RETURN
    2 A=DMAX1(X,Y)
      B=DMIN1(X,Y)
      NT=11.+2.*B+0.3*A
      IF(NT.GT.25) GO TO 6
      I=0
      IF(X.LT.Y) I=1
      GXY=1.+I*(B-1.)
      GXYO=GXY
      GX=1.0
      GY=GXY
      GZ=1.0
      DO 4 K=1,NT
      GX=GX*A/K
      GY=GY*B/(K+I)
      GZ=GZ+GX
      GXY=GXY+GY*GZ
      IF((GXY-GXYO)/GXY.LT.1.D-08) GO TO 5
    4 GXYO=GXY
    5 GOLD=GXY*DEXP(-X-Y)
      GO TO 8
    6 DA=DSQRT(A)
      DB=DSQRT(B)
      P=3.75/E
      B0=(.3989423+P*(.01328592+P*(.00225319-P*(.00157565-P*(.00916281-P
     $*(.02057706-P*(.02635537-P*(.01647633-.00392377*P))))))))/DSQRT(E)
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
     $*(1.8-3.3*P*Z)+97.45313*P)))))
      GOLD=0.5*BF+(.3535534*(DA+DB)*ERF+SUM)*BF/(B0*DSQRT(E))
    8 IF(X.LT.Y) GOLD=1.+BF-GOLD
      RETURN
      END
C
C     ---
      FUNCTION PHI1(TAU)
C
C     PURPOSE: TO CALCULATE SERIES IN EXPONENTIAL BVP
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/BVP/ NBP,T0B,DMB,GPB(10),TPB(10),GEB,ALB,JB  
      PARAMETER(ZERO=0.0D+00,ZIP=1.0D-10,ONE=1.0D+00,TWO=2.0D+00,
     & FOUR=4.0D+00)
C                  
      U=ALB-OMU2
      IF(ABS(U).LE.ZIP) U=ZIP
      PHI1=ZERO
      FN=ONE
      DO 12 N=1,25
       SUMK=ZERO
       K=1
       CK=ONE
       FNK=FN
       FN=FN*N 
    9  CK=-CK
       SUMK=SUMK+CK*(X(1)-TAU)**(N-K)/(U**K)/FNK
       IF (N.EQ.K) GOTO 12
       FNK=FNK/(N-K)
       K=K+1
       GOTO 9
   12  PHI1=PHI1+SUMK*((OMOM*TAU)**N)/FN
C  
       RETURN
       END
C
C     ---
      FUNCTION PHI2(TAU)
C
C     PURPOSE: TO CALCULATE SERIES IN EXPONENTIAL BVP
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON X(4),R,PX,PY,PZ,BETA,OMEGA,U1,U2,MP,MC,NGC,NP,BINT
      COMMON OMU1,OMU2,OM1,OM2,OMOM,U1R,U2R
      COMMON/BVP/ NBP,T0B,DMB,GPB(10),TPB(10),GEB,ALB,JB  
      PARAMETER(ZERO=0.0D+00,ZIP=1.0D-10,ONE=1.0D+00,TWO=2.0D+00,
     & FOUR=4.0D+00)
C
      U=ALB-OMU2
      IF(ABS(U).LE.ZIP) U=ZIP
      PHI2=ZERO
      FN=ONE
      DO 12 N=1,25
       SUMK=ZERO
       K=1
       CK=ONE
       FNK=FN
       FN=FN*N 
    9  CK=-CK
       SUMK=SUMK+CK*(X(1)-TAU)**(N-K+1)/(U**K)/FNK/(N-K+1)
       IF (N.EQ.K) GOTO 12
       FNK=FNK/(N-K)
       K=K+1
       GOTO 9
   12  PHI2=PHI2+SUMK*((OMOM*TAU)**N)/FN
C  
       RETURN
       END