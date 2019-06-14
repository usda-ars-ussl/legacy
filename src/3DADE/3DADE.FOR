C     ********************************************************************	
C     *                                                                  *
C     *                            3DADE                                 *
C     *                                                                  *
C     *      ANALYTICAL SOLUTIONS FOR SOLUTE TRANSPORT IN 3-D SEMI-      *
C     *       INFINITE POROUS MEDIA (WRR 27(10):2719-2733, 1991)         *
C     *        WITH FITTING BY NONLINEAR LEAST SQUARES ANALYSIS          *
C     *                                                                  *
C     *                FEIKE J. LEIJ AND SCOTT BRADFORD                  *
C     *                                                                  * 
C     *                    U.S. SALINITY LABORATORY                      *
C     *                      4500 GLENWOOD DRIVE                         *
C     *                      RIVERSIDE, CA 92501                         *
C     *                                                                  *           
C     *                      PHONE (714) 369-4846                        *
C     *                      FAX   (714) 369-4818                        *
C     *                                                                  *   
C     *                         VERSION  12/93                           *
C     ********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL PROD
      PARAMETER (NCON=250,ZIP=1.0D-08)
      COMMON NPAR,C,MODEL,NIV,INDEX
      DIMENSION Y(NCON),X(4,NCON),F(NCON),R(NCON),DELZ(NCON,8),B(16),
     $TH(16),P(8),PHI(8),Q(8),LSORT(NCON),TB(16),A(8,8),D(8,8),INDEX(8),
     $DX(4),XLOW(4),XUP(4),NXST(4),C(5),E(8),DEV(NCON)
      CHARACTER TITLE*60,BI(8)*10,TBI(8)*10,CN(5)*5,VN(4)*5,FILENAME*15,
     $JUNK*60
      DATA TBI/'V','R','To','MU','LAM','Dx',' ','Dz'/,
     $STOPCR/0.001/,MAXTRY/50/
      DATA VN/'t','x',' ','z'/,CN/'Co','a','b','x1','x2'/
C
C     --- OPEN I/O FILES ---
      WRITE(*, 150)
      READ(*, 160) FILENAME
      IF (FILENAME.EQ.'  ') FILENAME='3DADE.IN'
      OPEN(5, FILE=FILENAME, STATUS='OLD')
C
      WRITE(*, 170)
      READ(*,160) FILENAME
      IF (FILENAME.EQ.'  ') FILENAME='3DADE.OUT'
      OPEN(7, FILE=FILENAME, STATUS='UNKNOWN')
C
C     --- READ MODEL ---
      READ(5,201) JUNK
      READ(5,203) MODEL,MIT,IND,NGC 
C
C     --- INITIALIZATIONS ---			
      DO 2 I=1,8
    2  BI(I)=TBI(I)
      DO 3 I=1,4,3
       DX(I)=0.0
       XLOW(I)=0.0
       XUP(I)=1.0
    3  NXST(I)=1.0
      IF ((MODEL.GE.1.AND.MODEL.LE.6).OR.MODEL.EQ.11.OR.MODEL.EQ.12)THEN
       NIV=4
       NPAR=16
       BI(7)='Dy'         
       VN(3)='y'
      ELSE IF ((MODEL.GE.7.AND.MODEL.LE.10).OR.MODEL.EQ.13) THEN
       NIV=3
       NPAR=15
       BI(7)='Dr'                          
       VN(3)='r'
      ENDIF
      PROD=.TRUE.
C
C     --- READ INPUT FILE ---
      READ(5,201) TITLE
      READ(5,201) JUNK
      READ(5,202) (B(I),I=9,13)
      READ(5,203) (INDEX(I),I=1,5)
      READ(5,201) JUNK
      READ(5,202) (B(I),I=14,NPAR)
      READ(5,203) (INDEX(I),I=6,NPAR-8)
      READ(5,201) JUNK
      READ(5,202) (C(I),I=1,5)
C
C     --- DIRECT PROBLEM ---
      IF (IND.EQ.0) THEN
       NOB=1
       DO 4 I=1,NPAR-8
   4    INDEX(I)=0
C
       DO 5 I=1,4
        IF ((MODEL.GT.10.AND.I.EQ.1).OR.
     $  (((MODEL.GE.7.AND.MODEL.LE.10).OR.MODEL.EQ.13).AND.I.EQ.4)) THEN
         READ(5,201) JUNK
         READ(5,205)
         GOTO 5
        ENDIF
        READ(5,201) JUNK
        READ(5,202) DX(I),XLOW(I),XUP(I)
    5  CONTINUE
       GOTO 8
      ENDIF
C
C     --- INDIRECT PROBLEM --- 
      IF (MODEL.GT.10) INDEX(2)=0
      IF (MODEL.EQ.5 .OR. MODEL.EQ.6 .OR. MODEL.GT.8) INDEX(3)=0
      IF (MODEL.GT.10) INDEX(4)=0
      IF (MODEL.GT.10) INDEX(5)=0
      IF (MODEL.GT.10) INDEX(6)=0
C			   
      NOB=0
      DO 6 I=1,NCON
       IF(MODEL.GE.11) THEN	    
        READ(5,*,END=8) Y(I),(X(J,I),J=2,NIV)
       ELSE    
        READ(5,*,END=8) Y(I),(X(J,I),J=1,NIV)
       ENDIF
    6 NOB=NOB+1 
C
C     --- ERROR MESSAGES ---
    8 IF (MODEL.GT.13.OR.MODEL.LT.1) 
     $ STOP 'LINE 2 INPUT FILE : MODEL NUMBER 1-13'
      IF (MIT.LT.0)
     $ STOP 'LINE 2 INPUT FILE : MIT >= 0'
      IF (MODEL.LE.10.AND.B(10).LE.0.0) STOP 'LINE 4 INPUT FILE : R > 0'
      IF (B(9).LT.0.0) STOP 'LINE 4 INPUT FILE : V >=0'
      IF (MODEL.LE.4 .OR. MODEL.EQ.7 .OR. MODEL.EQ.8) THEN
       IF (B(11).LT.0.0) STOP 'LINE 4 INPUT FILE : To >= 0'
      ENDIF
      IF ((MODEL.GE.1.AND.MODEL.LE.6).AND.
     $ (B(14).LE.0.0.OR.B(15).LE.0.0.OR.B(16).LE.0.0))
     $  STOP 'LINE 6 INPUT FILE : Dx, Dy, Dz > 0'
      IF ((MODEL.GE.7.AND.MODEL.LE.10).AND.
     $ (B(14).LE.0.0.OR.B(15).LE.0.0))                  
     $  STOP 'LINE 6 INPUT FILE : Dx, Dr > 0'        
      IF ((MODEL.EQ.11.OR.MODEL.EQ.12).AND.
     $ (B(15).LE.0.0.OR.B(16).LE.0.0))
     $  STOP 'LINE 6 INPUT FILE : Dy, Dz > 0'
      IF (MODEL.EQ.13.AND.B(15).LE.0.0)                  
     $  STOP 'LINE 6 INPUT FILE : Dr > 0'        
      IF (C(1).LT.0.0)
     $  STOP 'LINE 8 INPUT FILE : Co >=0'
      IF ((MODEL.EQ.3.OR.MODEL.EQ.4.OR.MODEL.EQ.12) .AND.
     $  (C(2).LT.0.0 .OR. C(3).LT.0.0))
     $  STOP 'LINE 8 INPUT FILE : a,b >= 0'
      IF ((MODEL.EQ.5.OR.MODEL.EQ.6) .AND.
     $  (C(2).LT.0.0.OR.C(3).LT.0.0.OR.C(4).LT.0.0.OR.C(5).LT.0.0))
     $  STOP 'LINE 8 INPUT FILE : a,b,x1,x2 >= 0'
      IF ((MODEL.EQ.7.OR.MODEL.EQ.8.OR.MODEL.EQ.13) .AND. C(2).LT.0.0)
     $  STOP 'LINE 8 INPUT FILE : a >= 0'
      IF ((MODEL.EQ.9.OR.MODEL.EQ.10) .AND.
     $  (C(2).LT.0.0.OR.C(4).LT.0.0.OR.C(5).LT.0.0))
     $  STOP 'LINE 8 INPUT FILE : a,x1,x2 >= 0'
      IF ((MODEL.EQ.5.OR.MODEL.EQ.6.OR.MODEL.EQ.9.OR.MODEL.EQ.10)
     $  .AND. C(4).GT.C(5))
     $  STOP 'LINE 8 INPUT FILE : x1 <= x2'
      IF(IND.EQ.1) GOTO 11
      DO 9 I=1,2
    9 IF (XLOW(I).LT.0.0.OR.XUP(I).LT.0.0) 
     $  STOP 'INPUT ERROR GRID SYSTEM : x,t >= 0'
      DO 10 I=1,4
      IF (DX(I).LT.0.0)
     $  STOP 'INPUT ERROR GRID SYSTEM : INCREMENTS MUST BE POSITIVE'
   10 IF (XLOW(I).GT.XUP(I)) 
     $  STOP 'INPUT ERROR GRID SYSTEM : UPPER AND LOWER LIMITS'
      IF (((MODEL.GE.7.AND.MODEL.LE.10).OR.MODEL.EQ.13).AND.
     $XLOW(3).LT.0.0)
     $  STOP 'INPUT ERROR GRID SYSTEM : r >= 0'
      IF(IND.EQ.0) GOTO 14
   11 IF (INDEX(1).EQ.1.AND.INDEX(2).EQ.1.AND.INDEX(6).EQ.1)
     $  STOP 'EITHER v, R, OR Dx SHOULD BE FIXED'
      NINDEX = 0
      DO 12 I=1,8
      IF (INDEX(I).EQ.1.AND.ABS(B(I+8)).LE.ZIP)
     $  STOP 'INITIAL ESTIMATE OF FITTED VARIABLE CANNOT BE ZERO'
      IF (INDEX(I).EQ.1) NINDEX=NINDEX+1
   12 CONTINUE
      IF (NINDEX.EQ.1) STOP 'AT LEAST TWO VARIABLES MUST BE FITTED'
       DO 13 I=1,NOB
   13   IF (Y(I).LT.0.0.OR.X(1,I).LT.0.0.OR.X(2,I).LT.0.0)
     $   STOP 'INPUT ERROR OBSERVATIONS - C/Co, t, AND x >= 0.0'
C
C     --- CLEAR SCREEN AND WRITE TITLE ---
   14 WRITE(*,180)
      IF (IND.EQ.0) WRITE(*,190) FILENAME
      IF (IND.EQ.1) WRITE(*,191) FILENAME 
      WRITE(*,195)
C
C     --- START MAIN PROGRAM ---
      WRITE(7,300)
      WRITE(7,302) TITLE
C
C     --- WRITE MODEL TYPE ---
      IF (MODEL.EQ.1)  WRITE(7,304)
      IF (MODEL.EQ.2)  WRITE(7,306)
      IF (MODEL.EQ.3)  WRITE(7,308)
      IF (MODEL.EQ.4)  WRITE(7,310)
      IF (MODEL.EQ.5)  WRITE(7,312)
      IF (MODEL.EQ.6)  WRITE(7,314)
      IF (MODEL.EQ.7)  WRITE(7,316)
      IF (MODEL.EQ.8)  WRITE(7,318)
      IF (MODEL.EQ.9)  WRITE(7,320)
      IF (MODEL.EQ.10) WRITE(7,322)
      IF (MODEL.EQ.11) WRITE(7,324)
      IF (MODEL.EQ.12) WRITE(7,326)
      IF (MODEL.EQ.13) WRITE(7,328)
C  
C     --- WRITE PARAMETER VALUES AND MODEL CONSTANTS ---
      IF (IND.EQ.0) WRITE(7,330)
      IF (IND.EQ.1) WRITE(7,331)
      DO 15 I=9,NPAR
       IF ((MODEL.EQ.5.OR.MODEL.EQ.6) .AND. I.EQ.11) GOTO 15
       IF ((MODEL.EQ.9.OR.MODEL.EQ.10) .AND. I.EQ.11) GOTO 15
       IF ((MODEL.EQ.11.OR.MODEL.EQ.12).AND.(I.GT.9.AND.I.LT.15))GOTO 15
       IF (MODEL.EQ.13 .AND. (I.GT.9.AND.I.LT.15)) GOTO 15
       IF (IND.EQ.0) WRITE(7,336) BI(I-8),B(I)
       IF (IND.EQ.1) THEN
	IF(INDEX(I-8).EQ.0) THEN
         WRITE(7,332) BI(I-8),B(I),'-'
	ELSE
         WRITE(7,332) BI(I-8),B(I),'+'
        ENDIF
       ENDIF   
   15 CONTINUE
C
      WRITE(7,334)
      DO 16 I=1,5
       IF ((MODEL.LE.2.OR.MODEL.EQ.11).AND.I.GT.1) GOTO 16  
       IF ((MODEL.EQ.3.OR.MODEL.EQ.4.OR.MODEL.EQ.12).AND.I.GT.3) GOTO 16  
       IF ((MODEL.EQ.7.OR.MODEL.EQ.8.OR.MODEL.GE.13).AND.I.GT.2) GOTO 16
       IF ((MODEL.EQ.9.OR.MODEL.EQ.10).AND.I.EQ.3) GOTO 16
       WRITE(7,336) CN(I),C(I)
   16 CONTINUE   
C                   
C     --- CALCULATE NODE SPACING ---  
      IF (IND.EQ.0) THEN
      WRITE(7,338)
      DO 17 J=1,NIV
       I=J
       IF (MODEL.GT.10) I=J+1
       IF (ABS(DX(I)).LE.ZIP) THEN
        NXST(I)=1
       ELSE
        NXST(I)=ZIP+(XUP(I)+DX(I)-XLOW(I))/DX(I)
       ENDIF
       WRITE(7,340) VN(I),DX(I),XLOW(I),XUP(I) 
       IF (MODEL.GT.10.AND.I.EQ.NIV) GOTO 18
   17 CONTINUE
C
   18 IF (MODEL.LE.6.) THEN
       WRITE(7,342) CHAR(12)
      ELSE IF (MODEL.GE.7 .AND. MODEL.LE.10) THEN
       WRITE(7,343) CHAR(12)
      ELSE IF (MODEL.EQ.11 .OR. MODEL.EQ.12) THEN
       WRITE(7,344) CHAR(12)
      ELSE
       WRITE(7,345) CHAR(12)
      ENDIF
C
      NUM=0
      DO 19 I=1,NXST(1)
       X(1,1)=XLOW(1)+(I-1)*DX(1)
       DO 19 J=1,NXST(2)
        X(2,1)=XLOW(2)+(J-1)*DX(2)
        PROD=.TRUE.
        DO 19 K=1,NXST(3)
         X(3,1)=XLOW(3)+(K-1)*DX(3)
         DO 19 L=1,NXST(4)
          X(4,1)=XLOW(4)+(L-1)*DX(4)
          NUM=NUM+1
          CALL MODELS(B,Y,NOB,NCON,X,NGC,PROD)
          PROD=.FALSE.
          IF (C(1).GT.0.0) YY=Y(1)/C(1)  
       IF (MODEL.LE.10) WRITE(7,346) NUM,Y(1),YY,(X(JJ,1), JJ=1,NIV)
   19  IF (MODEL.GT.10) WRITE(7,346) NUM,Y(1),YY,(X(JJ,1), JJ=2,NIV) 
       GO TO 119
      ENDIF
C
C     --- ELSE WRITE EXPERIMENTAL DATA ---
      II=1
      IF (MODEL.LE.6.) THEN
       WRITE(7,347) CHAR(12)
      ELSE IF (MODEL.GE.7 .AND. MODEL.LE.10) THEN
       WRITE(7,348) CHAR(12)
      ELSE IF (MODEL.EQ.11 .OR. MODEL.EQ.12) THEN
       WRITE(7,349) CHAR(12)
      ELSE
       WRITE(7,350) CHAR(12)
      ENDIF
      DO 24 I=1,NOB
       IF (MODEL.LE.10) WRITE(7,346) I,Y(I),(X(J,I), J=1,NIV)
   24  IF (MODEL.GT.10) WRITE(7,346) I,Y(I),(X(J,I), J=2,NIV) 
C
C     --- BEGIN FITTING ---
      NP=0
      DO 26 I=9,NPAR
       TB(I)=B(I)
       IF (INDEX(I-8).EQ.0) GO TO 26
       NP=NP+1
       BI(NP)=BI(I-8)
       B(NP)=B(I)
       TB(NP)=B(I)
       TH(NP)=B(NP)
   26 TH(I)=B(I)
C
C     ------
      GA=0.02
      NIT=0
      CALL MODELS(TH,F,NOB,NCON,X,NGC,PROD)
      SSQ=0.
      DO 32 I=1,NOB
       R(I)=Y(I)-F(I)
   32  SSQ=SSQ+R(I)*R(I)
      WRITE(7,357) (BI(J),J=1,NP)
      WRITE(7,358) NIT,SSQ,(B(I),I=1,NP)
C
C     --- BEGIN OF ITERATION ---
   34 NIT=NIT+1
      NTRIAL=0 
      GA=0.1*GA
      DO 38 J=1,NP
      TEMP=TH(J)
      TH(J)=1.01*TH(J)
      Q(J)=0
      CALL MODELS(TH,DELZ(1,J),NOB,NCON,X,NGC,PROD)
      DO 36 I=1,NOB
       DELZ(I,J)=DELZ(I,J)-F(I)
   36  Q(J)=Q(J)+DELZ(I,J)*R(I)
      Q(J)=100.*Q(J)/TH(J)
C
C     --- Q=XT*R (STEEPEST DESCENT) ---
   38 TH(J)=TEMP
      DO 44 I=1,NP
      DO 42 J=1,I
      SUM=0
      DO 40 K=1,NOB
   40  SUM=SUM+DELZ(K,I)*DELZ(K,J)
      D(I,J)=10000.*SUM/(TH(I)*TH(J))
   42 D(J,I)=D(I,J)
   44 E(I)=DSQRT(D(I,I))
   50 DO 52 I=1,NP
       DO 52 J=1,NP 
   52   A(I,J)=D(I,J)/(E(I) * E(J))
C
C     --- A IS THE SCALED MOMENT MATRIX ---
      DO 54 I=1,NP
       P(I)=Q(I)/E(I)
       PHI(I)=P(I)
   54  A(I,I)=A(I,I)+GA
      CALL MATINV(A,NP,P)
C
C     --- P/E IS THE CORRECTION VECTOR ---
      STEP=1.0
   56 DO 58 I=1,NP
   58  TB(I)=P(I)*STEP/E(I)+TH(I)
      DO 62 I=1,NP
       IF (TH(I)*TB(I)) 66,66,62
   62 CONTINUE
      SUMB=0
      CALL MODELS(TB,F,NOB,NCON,X,NGC,PROD)
      DO 64 I=1,NOB
       R(I)=Y(I)-F(I)
   64  SUMB=SUMB+R(I)*R(I)
   66 SUM1=0.0
      SUM2=0.0
      SUM3=0.0
      DO 68 I=1,NP
       SUM1=SUM1+P(I)*PHI(I)
       SUM2=SUM2+P(I)*P(I)
   68  SUM3=SUM3+PHI(I)*PHI(I)
      ARG=SUM1/DSQRT(SUM2*SUM3)
      ARG1=0.
      IF (NP.GT.1) ARG1=DSQRT(1.-ARG*ARG)
      ANGLE=57.29578*DATAN2(ARG1,ARG)
C
C     ------
      DO 72 I=1,NP
       IF (TH(I)*TB(I)) 74,74,72
   72 CONTINUE
      NTRIAL=NTRIAL+1
      IF(NTRIAL.GT.MAXTRY) GO TO 95
      IF (SUMB/SSQ-1.0) 80,80,74
   74 IF (ANGLE-30.0) 76,76,78
   76 STEP=0.5*STEP
      GOTO 56
   78 GA=10.*GA
      GOTO 50
C
C     --- PRINT COEFFICIENTS AFTER EACH ITERATION ---
   80 CONTINUE
      DO 82 I=1,NP
   82  TH(I)=TB(I)
      WRITE(7,358) NIT,SUMB,(TH(I),I=1,NP)
      DO 86 I=1,NP
       IF(DABS(P(I)*STEP/E(I))/(1.0D-20+DABS(TH(I)))-STOPCR) 86,86,94
   86 CONTINUE
      GOTO 96
   94 SSQ=SUMB        
      IF (NIT.LT.MIT) GOTO 34
      IF (NIT.EQ.MIT) WRITE(7,353) NIT
      GOTO 96
   95 WRITE(7,354) MAXTRY
C
C     --- END OF ITERATION LOOP ---
   96 CONTINUE
      CALL MATINV(D,NP,P)
C
C     --- WRITE CORRELATION MATRIX ---
      DO 98 I=1,NP
   98  E(I)=DSQRT(DMAX1(D(I,I),1.D-20))
      WRITE(7,359) (BI(I),I=1,NP)
      DO 102 I=1,NP
       DO 100 J=1,I
  100   A(J,I)=D(J,I)/(E(I)*E(J))
  102  WRITE(7,360) BI(I),(A(J,I),J=1,I)
C
C     ----- CALCULATE r2 OF FITTED VS OBSERVED CONCENTRATION VALUES -----
      AVC=0.
      AVCDEV=0.
      IF(NOB.GT.0) THEN
       DO 103 I=1,NOB
        DEV(I)=Y(I)-F(I)
        AVC=AVC+Y(I)
  103   AVCDEV=AVCDEV+DEV(I)
       AVC=AVC/NOB
       AVCDEV=AVCDEV/NOB
C
       SUM1=0.
       SUM2=0.
       VARI=0.
       DO 104 I=1,NOB
        SUM1=SUM1+DEV(I)**2
        SUM2=SUM2+(Y(I)-AVC)**2
  104  VARI=VARI+(DEV(I)-AVCDEV)**2
       RSQ=1-SUM1/SUM2
       WRITE(7,355) RSQ     
C
       IF(NOB.GT.NP) THEN
        UMSE=SUM1/(NOB-NP)
        WRITE(7,366) UMSE
       ENDIF
       IF(NOB.GT.0) THEN
        VARI=100*VARI/SUM1
        BIAS=100*NOB*AVCDEV**2/SUM1
        WRITE(7,367) BIAS,VARI
       ENDIF		       
      ENDIF
C       
C     --- CALCULATE 95% CONFIDENCE INTERVAL ---
      Z=1./FLOAT(NOB-NP)
      SDEV=DSQRT(Z*SUMB)
      TVAR=1.96+Z*(2.3779+Z*(2.7135+Z*(3.187936+2.466666*Z**2)))
      WRITE(7,361)
      DO 108 I=1,NP
       SECOEF=E(I)*SDEV
       TVALUE=TH(I)/SECOEF
       TSEC=TVAR*SECOEF
       TMCOE=TH(I)-TSEC
       TPCOE=TH(I)+TSEC
  108  WRITE(7,362) I,BI(I),TH(I),SECOEF,TVALUE,TMCOE,TPCOE
C
C     --- PREPARE FINAL OUTPUT ---
      LSORT(1)=1
      DO 116 J=2,NOB
       TEMP=R(J)
       K=J-1
       DO 111 L=1,K
        LL=LSORT(L)
        IF(TEMP-R(LL)) 112,112,111
  111  CONTINUE
       LSORT(J)=J
       GO TO 116
  112  KK=J
  113  KK=KK-1
       LSORT(KK+1)=LSORT(KK)
       IF(KK-L) 115,115,113
  115  LSORT(L)=J
  116 CONTINUE
      WRITE(7,363)
      DO 118 I=1,NOB
      J=LSORT(NOB+1-I)
  118 WRITE(7,364) I,Y(I),F(I),R(I),J,Y(J),F(J),R(J)
  119 WRITE(7,365)
C
C     --- END OF PROBLEM ---
  150 FORMAT(' Enter input file name (default = 3DADE.IN)')
  160 FORMAT(A15)
  170 FORMAT(' Enter output file name (default = 3DADE.OUT)')
  180 FORMAT(////////////////)
  190 FORMAT(//5X,67(1H*)/5X,1H*,65X,1H*/5X,1H*,30X,'3DADE',30X,1H*/5X,1
     $H*,65X,1H*/5X,1H*,4X,'PLEASE WAIT - The concentrations are being c
     $alculated and',4X,1H*/5x,1H*,18X,'written to ',A15,21X,1H*/,5X,1H*
     $,65X,1H*)
  191 FORMAT(//5X,67(1H*)/5X,1H*,65X,1H*/5X,1H*,30X,'3DADE',30X,1H*/5X,1
     $H*,65X,1H*/5X,1H*,2X,'PLEASE WAIT - The observed concentrations ar
     $e being fitted and',1X,1H*/5x,1H*,16X,'written to ',A15,23X,1H*/,5
     $X,1H*,65X,1H*)
  195 FORMAT(5X,67(1H*),///) 
  200 FORMAT(5I5)
  201 FORMAT(A60)
  202 FORMAT(7F10.2)
  203 FORMAT(7I10)
  205 FORMAT()
  300 FORMAT(//5X,67(1H*)/5X,1H*,65X,1H*/5X,1H*,30X,'3DADE',30X,1H*/5X,
     $1H*,65X,1H*/5X,1H*,6X,'ANALYTICAL SOLUTIONS FOR SOLUTE TRANSPORT I
     $N 3-D SEMI-',5X,1H*/5X,1H*,8X,'INFINITE POROUS MEDIA (WRR 27(10):2
     $719-2733, 1991)',7X,1H*/5X,1H*,9X,'WITH FITTING BY NONLINEAR LEAST
     $ SQUARES ANALYSIS',8X,1H*/5X,1H*,65X,1H*/5X,1H*,65X,1H*)
  302 FORMAT(5X,1H*,5X,A60,1H*/5X,1H*,65X,1H*)
  304 FORMAT(5X,1H*,5X,'MODEL 1                      FIRST-TYPE INLET CO
     $NDITION',5X,1H*/5X,1H*,5X,'SCENARIO 1: SOLUTE APPLICATION FROM A S
     $EMI-INFINITE',9X,1H*/5X,1H*,17X,'QUADRANT OF AN INFINITE SOIL SURF
     $ACE',12X,1H*/5X,1H*,65X,1H*/5X,67(1H*)) 
  306 FORMAT(5X,1H*,5X,'MODEL 2                      THIRD-TYPE INLET CO
     $NDITION',5X,1H*/5X,1H*,5X,'SCENARIO 1: SOLUTE APPLICATION FROM A S
     $EMI-INFINITE',9X,1H*/5X,1H*,17X,'QUADRANT OF AN INFINITE SOIL SURF
     $ACE',12X,1H*/5X,1H*,65X,1H*/5X,67(1H*))
  308 FORMAT(5X,1H*,5X,'MODEL 3                      FIRST-TYPE INLET CO
     $NDITION',5X,1H*/5X,1H*,5X,'SCENARIO 2: SOLUTE APPLICATION FROM A R
     $ECTANGULAR',11X,1H*/5X,1H*,17X,'SOURCE AT THE SOIL SURFACE',22X,1H
     $*/5X,1H*,65X,1H*/5X,67(1H*)) 
  310 FORMAT(5X,1H*,5X,'MODEL 4                      THIRD-TYPE INLET CO
     $NDITION',5X,1H*/5X,1H*,5X,'SCENARIO 2: SOLUTE APPLICATION FROM A R
     $ECTANGULAR',11X,1H*/5X,1H*,17X,'SOURCE AT THE SOIL SURFACE',22X,1H
     $*/5X,1H*,65X,1H*/5X,67(1H*)) 
  312 FORMAT(5X,1H*,5X,'MODEL 5                      FIRST-TYPE INLET CO
     $NDITION',5X,1H*/5X,1H*,5X,'SCENARIO 3: SOLUTE INITIALLY UNIFORMLY 
     $DISTRIBUTED IN A',5X,1H*/5X,1H*,17X,'PARALLELEPIPEDAL REGION OF TH
     $E SOIL',13X,1H*/5X,1H*,65X,1H*/5X,67(1H*)) 
  314 FORMAT(5X,1H*,5X,'MODEL 6                      THIRD-TYPE INLET CO
     $NDITION',5X,1H*/5X,1H*,5X,'SCENARIO 3: SOLUTE INITIALLY UNIFORMLY 
     $DISTRIBUTED IN A',5X,1H*/5X,1H*,17X,'PARALLELEPIPEDAL REGION OF TH
     $E SOIL',13X,1H*/5X,1H*,65X,1H*/5X,67(1H*)) 
  316 FORMAT(5X,1H*,5X,'MODEL 7                      FIRST-TYPE INLET CO
     $NDITION',5X,1H*/5X,1H*,5X,'SCENARIO 4: SOLUTE APPLICATION FROM A C
     $IRCULAR',14X,1H*/5X,1H*,17X,'SOURCE AT THE SOIL SURFACE',22X,1H*/5
     $X,1H*,65X,1H*/5X,67(1H*)) 
  318 FORMAT(5X,1H*,5X,'MODEL 8                      THIRD-TYPE INLET CO
     $NDITION',5X,1H*/5X,1H*,5X,'SCENARIO 4: SOLUTE APPLICATION FROM A C
     $IRCULAR',14X,1H*/5X,1H*,17X,'SOURCE AT THE SOIL SURFACE',22X,1H*/5
     $X,1H*,65X,1H*/5X,67(1H*)) 
  320 FORMAT(5X,1H*,5X,'MODEL 9                      FIRST-TYPE INLET CO
     $NDITION',5X,1H*/5X,1H*,5X,'SCENARIO 5: SOLUTE INITIALLY UNIFORMLY 
     $DISTRIBUTED IN A',5X,1H*/5X,1H*,17X,'CYLINDRICAL REGION OF THE SOI
     $L',18X,1H*/5X,1H*,65X,1H*/5X,67(1H*)) 
  322 FORMAT(5X,1H*,5X,'MODEL 10                     THIRD-TYPE INLET CO
     $NDITION',5X,1H*/5X,1H*,5X,'SCENARIO 5: SOLUTE INITIALLY UNIFORMLY 
     $DISTRIBUTED IN A',5X,1H*/5X,1H*,17X,'CYLINDRICAL REGION OF THE SOI
     $L',18X,1H*/5X,1H*,65X,1H*/5X,67(1H*)) 
  324 FORMAT(5X,1H*,5X,'MODEL 11                                   STEAD
     $Y STATE',5X,1H*/5X,1H*,5X,'SCENARIO 1: SOLUTE APPLICATION FROM A S
     $EMI-INFINITE',9X,1H*/5X,1H*,17X,'QUADRANT OF AN INFINITE SOIL SURF
     $ACE',12X,1H*/5X,1H*,65X,1H*/5X,67(1H*)) 
  326 FORMAT(5X,1H*,5X,'MODEL 12                                   STEAD
     $Y STATE',5X,1H*/5X,1H*,5X,'SCENARIO 2: SOLUTE APPLICATION FROM A R
     $ECTANGULAR',11X,1H*/5X,1H*,17X,'SOURCE AT THE SOIL SURFACE',22X,1H
     $*/5X,1H*,65X,1H*/5X,67(1H*)) 
  328 FORMAT(5X,1H*,5X,'MODEL 13                                   STEAD
     $Y STATE',5X,1H*/5X,1H*,5X,'SCENARIO 4: SOLUTE APPLICATION FROM A C
     $IRCULAR',14X,1H*/5X,1H*,17X,'SOURCE AT THE SOIL SURFACE',22X,1H*/5
     $X,1H*,65X,1H*/5X,67(1H*)) 
  330 FORMAT(//5X,'MODEL PARAMETERS'/5X,16(1H=)/5X,'NAME',11X,'VALUE')
  331 FORMAT(//5X,'INITIAL ESTIMATES OF MODEL PARAMETERS'/5X,37(1H=)/5X,
     $'NAME',11X,'VALUE',7X,'FITTING')
  332 FORMAT(5X,A5,5X,F10.2,10X,A1)
  334 FORMAT(//5X,'CONSTANTS'/5X,9(1H=)/5X,'NAME',11X,'VALUE')
  336 FORMAT(5X,A5,5X,F10.2)
  338 FORMAT(//5X,'INDEPENDENT VARIABLES'/5X,21(1H=)/5X,'VARIABLE',3X,'I
     $NCREMENT',8X,'MINIMUM',8X,'MAXIMUM')
  340 FORMAT(5X,A5,5X,F10.3,5X,F10.3,5X,F10.3)
  342 FORMAT(/A1/5X,'PREDICTED CONCENTRATIONS'/5X,24(1H=)/8X,'NO',8X,'C'
     $,9X,'C/Co',8X,'t',9X,'x',9X,'y',9X,'z')
  343 FORMAT(/A1/5X,'PREDICTED CONCENTRATIONS'/5X,24(1H=)/8X,'NO',8X,'C'
     $,9X,'C/Co',8X,'t',9X,'x',9X,'r')
  344 FORMAT(/A1/5X,'PREDICTED CONCENTRATIONS'/5X,24(1H=)/8X,'NO',8X,'C'
     $,9X,'C/Co',8X,'x',9X,'y',9X,'z')
  345 FORMAT(/A1/5X,'PREDICTED CONCENTRATIONS'/5X,24(1H=)/8X,'NO',8X,'C'
     $,9X,'C/Co',8X,'x',9X,'r')
  346 FORMAT(5X,I5,3X,2F10.4,5F10.2)
  347 FORMAT(/A1/5X,'OBSERVED CONCENTRATIONS'/5X,23(1H=)/8X,'NO',8X,'C',
     $9X,'t',11X,'x',9X,'y',9X,'z',9X)
  348 FORMAT(/A1/5X,'OBSERVED CONCENTRATIONS'/5X,23(1H=)/8X,'NO',8X,'C',
     $9X,'t',11X,'x',9X,'r',9X)
  349 FORMAT(/A1/5X,'OBSERVED CONCENTRATIONS'/5X,23(1H=)/8X,'NO',8X,'C',
     $9X,'x',11X,'y',9X,'z',9X)
  350 FORMAT(/A1/5X,'OBSERVED CONCENTRATIONS'/5X,23(1H=)/8X,'NO',8X,'C',
     $9X,'x',11X,'r',9X)
  353 FORMAT(/5X,'CONVERGENCE CRITERION NOT MET IN ',I2,' ITERATIONS') 
  354 FORMAT(/5X,'NO FURTHER DECREASE IN SSQ OBTAINED AFTER ',I2,' TRIAL
     $S')
  355 FORMAT(//5X,'r^2 FOR REGRESSION OF OBSERVED VS. PREDICTED VALUES =
     $',F9.6/5X,51(1H=))
  357 FORMAT(//5X,'RESULTS FOR EACH ITERATION STEP'/5X,31(1H=)/5X,'STEP'
     $,6X,'SSQ',5X,7(3X,A3,2X))
  358 FORMAT(5X,I3,3X,F10.5,2X,7(F7.3,1X))
  359 FORMAT(//5X,'CORRELATION MATRIX'/5X,18(1H=)/8X,10(2X,A3,6X))
  360 FORMAT(3X,A3,10(2X,F7.4,2X))
  361 FORMAT(///5X,'NON-LINEAR LEAST SQUARES ANALYSIS, FINAL RESULTS'/5X
     $,48(1H=)//54X,'95% CONFIDENCE LIMITS'/5X,'VAR',5X,'NAME',3X,'VALUE
     $',6X,'S.E.COEFF.',3X,'T-VALUE',5X,'LOWER',8X,'UPPER')
  362 FORMAT(5X,I2,7X,A3,F10.5,1X,F11.4,1X,F11.2,1X,F10.3,1X,F12.3)
  363 FORMAT(//5X,4(1H-),'ORDERED BY COMPUTER INPUT',5(1H-), 3X,8(1H-),'
     $ORDERED BY RESIDUALS',6(1H-)/11X,3X,'CONCENTRATION',4X,'RESI-',8X,
     $8X,'CONCENTRATION',4X,'RESI-'/6X,'NO',2X,5X,'OBS ',2X,'FITTED',4X,
     $'DUAL',9X,'NO',1X,6X,'OBS ',2X,'FITTED',4X,'DUAL')
  364 FORMAT(5X,I3,3X,3F8.3,8X,I3,3X,3F8.3)
  365 FORMAT(///5X,'END OF PROBLEM'/5X,14(1H=))
  366 FORMAT(/5X,'MEAN SQUARE FOR ERROR (MSE) =',F11.6/5X,27(1H=))
  367 FORMAT(/5X,'BIAS (%)       =',F8.2,5X,'VARIANCE (%)      =',F8.2)
C
      CLOSE(5)
      CLOSE(7)
C
      STOP
      END
C
C     ------------
      SUBROUTINE MATINV(A,NP,B)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(8,8),B(16),INDEX(8,2)
      DO 2 J=1,5
    2 INDEX(J,1)=0
      I=0
    4 AMAX=-1.0
      DO 11 J=1,NP
      IF(INDEX(J,1)) 11,6,11
    6 DO 10 K=1,NP
      IF(INDEX(K,1)) 10,8,10
    8 P=DABS(A(J,K))
      IF(P.LE.AMAX) GO TO 10
      IR=J
      IC=K
      AMAX=P
   10 CONTINUE
   11 CONTINUE
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
C     ------------
      SUBROUTINE MODELS(BB,YY,NOB,NCON,XX,M,PROD)
C
C     purpose: to calculate YY(X1,X2,X3,...)
C
      IMPLICIT REAL*8 (A-H,L,O-Z)
      LOGICAL PROD 
      COMMON NPAR,CC,MODEL,NIV,INDEX
      COMMON/EVAL/ T,X,Y,Z,V,R,TO,RMU,RLAM,DX,DY,DZ,CO,A,B,X1,X2,NGC
      PARAMETER(ZERO=0.0D+00,ZIP=1.0D-08,ONE=1.0D+00,TWO=2.0D+00,
     $FOUR=4.0D+00,EIGHT=8.0D+00)
      DIMENSION XX(4,NCON),YY(NCON),BB(16),CC(5),INDEX(8)
      EXTERNAL CF1,CT1,CF2,CT2,CF4,LAM2,LAM2C,CT4,LAM4,LAM4C,CF5,CT5,CS3 		 
      PI=FOUR*DATAN(ONE)
C            			
      NGC=M
      K=0
      DO 5 I=9,NPAR
       IF (INDEX(I-8).EQ.0) GOTO 5
       K=K+1
       BB(I)=BB(K)
    5 CONTINUE
      V=BB(9)
      R=BB(10)
      TO=BB(11)
      RMU=BB(12)
      RLAM=BB(13)
      DX=BB(14)
      DY=BB(15)
      DZ=BB(16)
      CO=CC(1)
      A=CC(2)
      B=CC(3)
      X1=CC(4)
      X2=CC(5)
      DO 85 I=1,NOB
       T=XX(1,I)
       X=XX(2,I)
       Y=XX(3,I)
       Z=XX(4,I)	     
C
C     --- INLET AND INITIAL CONCENTRATIONS ---
       ZIP1=ABS(X-X1)
       ZIP2=ABS(X-X2)  
       ZIPA=ABS(ABS(Y)-A)
       ZIPB=ABS(ABS(Z)-B)
       GOTO (25,25,35,35,45,45,55,55,65,65,25,35,55) MODEL
   25  IF (X.LT.ZIP .AND. MODEL.NE.2) THEN
        YY(I)=ZERO
        IF (Y.LT.0.0 .AND. Z.LT.0.0) YY(I)=CO
        IF ((ABS(Y).LT.ZIP .AND. Z.LT.0.0) .OR. (Y.LT.0.0
     $   .AND. ABS(Z).LT.ZIP)) YY(I)=CO/TWO
        IF (ABS(Y).LT.ZIP .AND. ABS(Z).LT.ZIP) YY(I)=CO/FOUR
        IF (MODEL.NE.11 .AND. T.GE.TO) YY(I)=ZERO
       ELSE IF (MODEL.NE.11 .AND. T.LT.ZIP) THEN
        YY(I)=ZERO
       ELSE
        GOTO 75
       ENDIF                                      
       GOTO 85
C
   35  IF (X.LT.ZIP .AND. MODEL.NE.4) THEN
	YY(I)=ZERO
        IF (ABS(Y).LT.A .AND. ABS(Z).LT.B) YY(I)=CO 
        IF ((ZIPA.LE.ZIP .AND. ABS(Z).LT.B) .OR. (ABS
     $   (Y).LT.A .AND. ZIPB.LE.ZIP)) YY(I)=CO/TWO
        IF (ZIPA.LE.ZIP .AND.ZIPB.LE.ZIP) YY(I)=CO/FOUR 
        IF (MODEL.NE.12. AND. T.GE.TO) YY(I)=ZERO
       ELSE IF (MODEL.NE.12 .AND. T.LT.ZIP) THEN
        YY(I)=ZERO
       ELSE
        GOTO 75 
       ENDIF
       GOTO 85							      
C
   45  IF (T.LT.ZIP) THEN
        IF (ZIP1.LE.ZIP .OR. ZIP2.LE.ZIP) THEN
         IF (ZIPA.LE.ZIP .AND. ZIPB.LE.ZIP) THEN
          YY(I)=CO/EIGHT				  
         ELSE IF ((ABS(Y).LT.A .AND. ZIPB.LE.ZIP) .OR. (ZIPA
     $    .LE.ZIP .AND. ABS(Z).LT.B)) THEN
          YY(I)=CO/FOUR
         ELSE IF (ABS(Y).LT.A .AND. ABS(Z).LT.B) THEN
          YY(I)=CO/TWO
         ELSE
          YY(I)=ZERO
         ENDIF
        ELSE IF (X.GT.X1 .AND. X.LT.X2) THEN
         IF (ZIPA.LE.ZIP .AND. ZIPB.LE.ZIP) THEN
          YY(I)=CO/FOUR
         ELSE IF (ZIPA.LE.ZIP .AND. ABS(Z).LT.B) THEN
          YY(I)=CO/TWO
         ELSE IF (ABS(Y).LT.A .AND. ZIPB.LE.ZIP) THEN
          YY(I)=CO/TWO
         ELSE IF (ABS(Y).LT.A .AND. ABS(Z).LT.B) THEN
          YY(I)=CO
         ELSE IF (ABS(Y).GT.A .OR. ABS(Z).GT.B) THEN
          YY(I)=ZERO
         ENDIF 
        ELSE IF (X.LT.X1 .OR. X.GT.X2) THEN
         YY(I)=ZERO
        ENDIF
       ELSE
        GOTO 75 
       ENDIF
       GOTO 85
C 
  55   IF (X.LT.ZIP .AND. MODEL.NE.8) THEN
        YY(I)=ZERO
        IF (Y.LT.A) YY(I)=CO 
        IF (ZIPA.LE.ZIP) YY(I)=CO/TWO
        IF (MODEL.NE.13 .AND. T.GE.TO) YY(I)=ZERO
       ELSE IF (T.LT.ZIP .AND. MODEL.NE.13) THEN
        YY(I)=ZERO
       ELSE
        GOTO 75
       ENDIF
       GOTO 85
C
   65  IF (T.LT.ZIP) THEN
        IF (ZIP1.LE.ZIP .OR. ZIP2.LE.ZIP) THEN
         IF (ZIPA.LE.ZIP) YY(I)=CO/FOUR
         IF (Y.LT.A) YY(I)=CO/TWO
         IF (Y.GT.A) YY(I)=ZERO 
        ELSE IF (X.GT.X1 .AND. X.LT.X2) THEN
         IF (ZIPA.LE.ZIP) YY(I)=CO/TWO
         IF (Y.LT.A) YY(I)=CO
         IF (Y.GT.A) YY(I)=ZERO
        ELSE IF (X.LT.X1 .OR. X.GT.X2) THEN
         YY(I)=ZERO
        ENDIF
       ELSE
        GOTO 75
       ENDIF
       GOTO 85
C
   75  IF (MODEL.EQ.1) THEN
         CALL CHEBY(CF1,RES1,ZERO,T)
         IF(T.GT.TO) THEN
           PT=T-TO
           CALL CHEBY(CF1,RESP,ZERO,PT)
           RES1=RES1-RESP
         ENDIF
         IF (RLAM.GT.0. .AND. PROD)  CALL CHEBY(LAM2,RES2,ZERO,T)
         YY(I)=CO*RES1/FOUR+RLAM*RES2/TWO/R
       ELSE IF (MODEL.EQ.2) THEN
         CALL CHEBY(CT1,RES1,ZERO,T)
         IF(T.GT.TO) THEN
           PT=T-TO
           CALL CHEBY(CT1,RESP,ZERO,PT)
           RES1=RES1-RESP
         ENDIF
         IF (RLAM.GT.0..AND. PROD)  CALL CHEBY(LAM4,RES2,ZERO,T)
         YY(I)=CO*V*RES1/R/FOUR+RLAM*RES2/TWO/R
       ELSE IF (MODEL.EQ.3) THEN
         CALL CHEBY(CF2,RES1,ZERO,T)
         IF(T.GT.TO) THEN
           PT=T-TO
           CALL CHEBY(CF2,RESP,ZERO,PT)
           RES1=RES1-RESP
         ENDIF
         IF (RLAM.GT.0. .AND. PROD)  CALL CHEBY(LAM2,RES2,ZERO,T)
         YY(I)=CO*RES1/FOUR+RLAM*RES2/TWO/R
       ELSE IF (MODEL.EQ.4) THEN
         CALL CHEBY(CT2,RES1,ZERO,T)
         IF(T.GT.TO) THEN
           PT=T-TO
           CALL CHEBY(CT2,RESP,ZERO,PT)
           RES1=RES1-RESP
         ENDIF
         IF (RLAM.GT.0. .AND. PROD)  CALL CHEBY(LAM4,RES2,ZERO,T)
         YY(I)=CO*V*RES1/R/FOUR+RLAM*RES2/TWO/R
       ELSE IF (MODEL.EQ.5) THEN
         IF (RLAM.GT.0. .AND. PROD)  CALL CHEBY(LAM2,RES2,ZERO,T)
         YY(I)=CO*LAM5(T)*GAMMA2(T)/EIGHT+RLAM*RES2/TWO/R
       ELSE IF (MODEL.EQ.6) THEN
         IF (RLAM.GT.0. .AND. PROD)  CALL CHEBY(LAM4,RES2,ZERO,T)
         Y1=CO*LAM6(T)*GAMMA2(T)/EIGHT
         Y2=RLAM*RES2/TWO/R
         YY(I)=y1+y2
       ELSE IF (MODEL.EQ.7) THEN
         CALL CHEBY(CF4,RES1,ZERO,T)
         IF(T.GT.TO) THEN
           PT=T-TO
           CALL CHEBY(CF4,RESP,ZERO,PT)
           RES1=RES1-RESP
         ENDIF
         IF (RLAM.GT.0. .AND. PROD)  CALL CHEBY(LAM2,RES2,ZERO,T)
         YY(I)=CO*RES1/TWO+RLAM*RES2/TWO/R
       ELSE IF (MODEL.EQ.8) THEN
         CALL CHEBY(CT4,RES1,ZERO,T)
         IF(T.GT.TO) THEN
           PT=T-TO
           CALL CHEBY(CT4,RESP,ZERO,PT)
           RES1=RES1-RESP
         ENDIF
         IF (RLAM.GT.0..AND. PROD) CALL CHEBY(LAM4,RES2,ZERO,T)
         YY(I)=CO*RES1*V/TWO/R+RLAM*RES2/TWO/R
       ELSE IF (MODEL.EQ.9) THEN
         CALL CHEBY(CF5,RES1,ZERO,A)
         IF (RLAM.GT.0. .AND. PROD) CALL CHEBY(LAM2,RES2,ZERO,T)
         YY(I)=CO*RES1/FOUR+RLAM*RES2/TWO/R
       ELSE IF (MODEL.EQ.10) THEN
         CALL CHEBY(CT5,RES1,ZERO,A)
         IF (RLAM.GT.0..AND. PROD) CALL CHEBY(LAM4,RES2,ZERO,T)
         YY(I)=CO*RES1/FOUR+RLAM*RES2/TWO/R
       ELSE IF (MODEL.EQ.11) THEN
         YY(I)=CO*EXF(ZERO,Y/DSQRT(FOUR*DY*X/V))
     $      *EXF(ZERO,Z/DSQRT(FOUR*DZ*X/V))/FOUR
       ELSE IF (MODEL.EQ.12) THEN
         YY(I)=CO/FOUR*(EXF(ZERO,(Y-A)/DSQRT(FOUR*DY*X/V))
     $       -EXF(ZERO,(Y+A)/DSQRT(FOUR*DY*X/V)))
     $       *(EXF(ZERO,(Z-B)/DSQRT(FOUR*DZ*X/V))
     $       -EXF(ZERO,(Z+B)/DSQRT(FOUR*DZ*X/V)))
       ELSE IF (MODEL.EQ.13) THEN
         CALL CHEBY(CS3,RES,ZERO,A)
         YY(I)=CO*RES
       ENDIF
C
   85 CONTINUE
C
      RETURN
      END
C
C     ----------
      REAL*8 FUNCTION CF1(T)
      IMPLICIT REAL*8(A-M,O-Z)
      COMMON/EVAL/TT,X,Y,Z,V,R,TO,RMU,RLAM,DX,DY,DZ,CO,A,B,X1,X2,NGC
      PARAMETER(ZERO=0.0D+00,ZIP=1.0D-08,FOUR=4.0D+00)
      PI=4.0D+00*DATAN(1.0D+00)
C
      IF(T.LE.ZIP) T=ZIP
      CF1=(X/T)*DSQRT(R/(FOUR*PI*DX*T))
     $ *DEXP(-RMU*T/R-(R*X-V*T)**2/FOUR/R/DX/T)
     $ *EXF(ZERO,Y/DSQRT(FOUR*DY*T/R))*EXF(ZERO,Z/DSQRT(FOUR*DZ*T/R))
C
      RETURN
      END
C
C     ----------
      REAL*8 FUNCTION CT1(T)
      IMPLICIT REAL*8(A-M,O-Z)
      COMMON/EVAL/TT,X,Y,Z,V,R,TO,RMU,RLAM,DX,DY,DZ,CO,A,B,X1,X2,NGC
      PARAMETER(ZERO=0.0D+00,ZIP=1.0D-08,TWO=2.0D+00,FOUR=4.0D+00)
      PI=4.0D+00*DATAN(1.0D+00)
C
      IF(T.LE.ZIP) T=ZIP
      CT1=(DSQRT(R/PI/DX/T)*DEXP(-RMU*T/R-(R*X-V*T)**2/FOUR/R/DX/T)
     $   -V/TWO/DX*EXF(V*X/DX-RMU*T/R,(R*X+V*T)/DSQRT(FOUR*R*DX*T)))
     $   *EXF(ZERO,Y/DSQRT(FOUR*DY*T/R))*EXF(ZERO,Z/DSQRT(FOUR*DZ*T/R))
C
      RETURN
      END
C
C     ----------
      REAL*8 FUNCTION CF2(T)
      IMPLICIT REAL*8(A-M,O-Z)
      COMMON/EVAL/TT,X,Y,Z,V,R,TO,RMU,RLAM,DX,DY,DZ,CO,A,B,X1,X2,NGC
      PARAMETER(ZERO=0.0D+00,ZIP=1.0D-08,FOUR=4.0D+00)
      PI=4.0D+00*DATAN(1.0D+00)
C
      IF(T.LE.ZIP) T=ZIP		 
      CF2=(X/T)*DSQRT(R/(FOUR*PI*DX*T))
     $    *DEXP(-RMU*T/R-(R*X-V*T)**2/FOUR/R/DX/T)
     $    *(EXF(ZERO,(Y-A)/DSQRT(FOUR*DY*T/R))
     $    -EXF(ZERO,(Y+A)/DSQRT(FOUR*DY*T/R)))
     $    *(EXF(ZERO,(Z-B)/DSQRT(FOUR*DZ*T/R))
     $    -EXF(ZERO,(Z+B)/DSQRT(FOUR*DZ*T/R)))
C			 
      RETURN
      END
C
C     ----------
      REAL*8 FUNCTION CT2(T)
      IMPLICIT REAL*8(A-M,O-Z)
      COMMON/EVAL/TT,X,Y,Z,V,R,TO,RMU,RLAM,DX,DY,DZ,CO,A,B,X1,X2,NGC
      PARAMETER(ZERO=0.0D+00,ZIP=1.0D-08,TWO=2.0D+00,FOUR=4.0D+00)
      PI=4.0D+00*DATAN(1.0D+00)
C
      IF(T.LE.ZIP) T=ZIP		 
      CT2=(DSQRT(R/PI/DX/T)*DEXP(-RMU*T/R-(R*X-V*T)**2/FOUR/R/DX/T)
     $    -V/TWO/DX*EXF(V*X/DX-RMU*T/R,(R*X+V*T)/DSQRT(FOUR*R*DX*T)))
     $    *(EXF(ZERO,(Y-A)/DSQRT(FOUR*DY*T/R))
     $    -EXF(ZERO,(Y+A)/DSQRT(FOUR*DY*T/R)))
     $    *(EXF(ZERO,(Z-B)/DSQRT(FOUR*DZ*T/R))
     $    -EXF(ZERO,(Z+B)/DSQRT(FOUR*DZ*T/R)))
C
      RETURN
      END
C
C     ----------
      REAL*8 FUNCTION CF4(T)
      IMPLICIT REAL*8(A-M,O-Z)
      COMMON/EVAL/TT,X,Y,Z,V,R,TO,RMU,RLAM,DX,DY,DZ,CO,A,B,X1,X2,NGC
      PARAMETER(ZERO=0.0D+00,FOUR=4.0D+00)
      EXTERNAL CF41
      PI=4.0D+00*DATAN(1.0D+00)
C
      CALL CHEONE(CF41,T,RES,ZERO,A)
      CF4=RES
      RETURN
      END
C
C     ----------
      REAL*8 FUNCTION CF41(RHO,T)
      IMPLICIT REAL*8(A-M,O-Z)
      COMMON/EVAL/TT,X,Y,Z,V,R,TO,RMU,RLAM,DX,DY,DZ,CO,A,B,X1,X2,NGC
      PARAMETER(ZERO=0.0D+00,ZIP=1.0D-08,TWO=2.0D+00,FOUR=4.0D+00)
      PI=4.0D+00*DATAN(1.0D+00)
C
      IF(T.LE.ZIP) T=ZIP
      CF41=(X/T)*DSQRT(R/(FOUR*PI*DX*T))
     $    *DEXP(-RMU*T/R-(R*X-V*T)**2/FOUR/R/DX/T)
     $    *RHO*R/DY/T*EXPBI0(R*Y*RHO/TWO/DY/T,
     $    -(R*(Y**2+RHO**2)/FOUR/DY/T))
C
      RETURN
      END
C
C     ----------
      REAL*8 FUNCTION CT4(T)
      IMPLICIT REAL*8(A-M,O-Z)
      COMMON/EVAL/TT,X,Y,Z,V,R,TO,RMU,RLAM,DX,DY,DZ,CO,A,B,X1,X2,NGC
      PARAMETER(ZERO=0.0D+00,FOUR=4.0D+00)
      EXTERNAL CT41
      PI=4.0D+00*DATAN(1.0D+00)
C
      CALL CHEONE(CT41,T,RES,ZERO,A)
      CT4=RES
      RETURN
      END
C
C     ----------
      REAL*8 FUNCTION CT41(RHO,T)
      IMPLICIT REAL*8(A-M,O-Z)
      COMMON/EVAL/TT,X,Y,Z,V,R,TO,RMU,RLAM,DX,DY,DZ,CO,A,B,X1,X2,NGC
      PARAMETER(ZERO=0.0D+00,ZIP=1.0D-08,TWO=2.0D+00,FOUR=4.0D+00)
      PI=4.0D+00*DATAN(1.0D+00)
C					 
      IF(T.LE.ZIP) T=ZIP
      CT41=(DSQRT(R/PI/DX/T)*DEXP(-RMU*T/R-(R*X-V*T)**2/FOUR/R/DX/T)
     $     -V/TWO/DX*EXF(V*X/DX-RMU*T/R,(R*X+V*T)/DSQRT(FOUR*R*DX*T)))
     $     *RHO*R/DY/T*EXPBI0(R*Y*RHO/TWO/DY/T,
     $     -(R*(Y**2+RHO**2)/FOUR/DY/T))
C
      RETURN
      END
C
C     ----------
      REAL*8 FUNCTION CF5(RHO)
      IMPLICIT REAL*8(A-M,O-Z)
      COMMON/EVAL/ T,X,Y,Z,V,R,TO,RMU,RLAM,DX,DY,DZ,CO,A,B,X1,X2,NGC
      PARAMETER(ZERO=0.0D+00,TWO=2.0D+00,FOUR=4.0D+00)
      PI=4.0D+00*DATAN(1.0D+00)
      CF5=(EXF(-RMU*T/R,(R*(X-X2)-V*T)/DSQRT(FOUR*R*DX*T))
     $    -EXF(-RMU*T/R,(R*(X-X1)-V*T)/DSQRT(FOUR*R*DX*T))
     $    +EXF(V*X/DX-RMU*T/R,(R*(X+X2)+V*T)/DSQRT(FOUR*R*DX*T))
     $    -EXF(V*X/DX-RMU*T/R,(R*(X+X1)+V*T)/DSQRT(FOUR*R*DX*T)))
     $    *RHO*R/DY/T*EXPBI0(R*Y*RHO/TWO/DY/T,
     $    -(R*(Y**2+RHO**2)/FOUR/DY/T))
C
      RETURN
      END
C
C     ----------
      REAL*8 FUNCTION CT5(RHO)
      IMPLICIT REAL*8(A-M,O-Z)
      COMMON/EVAL/ T,X,Y,Z,V,R,TO,RMU,RLAM,DX,DY,DZ,CO,A,B,X1,X2,NGC
      PARAMETER(ZERO=0.0D+00,ONE=1.0D+00,TWO=2.0D+00,FOUR=4.0D+00)
      PI=4.0D+00*DATAN(1.0D+00)
      CT5=((ONE+V/DX*(X+X1)+V**2*T/R/DX)
     $     *EXF(V*X/DX-RMU*T/R,(R*(X+X1)+V*T)/DSQRT(FOUR*R*DX*T))
     $     -(ONE+V/DX*(X+X2)+V**2*T/R/DX)
     $     *EXF(V*X/DX-RMU*T/R,(R*(X+X2)+V*T)/DSQRT(FOUR*R*DX*T))
     $     +EXF(-RMU*T/R,(R*(X-X2)-V*T)/DSQRT(FOUR*R*DX*T))
     $     -EXF(-RMU*T/R,(R*(X-X1)-V*T)/DSQRT(FOUR*R*DX*T))
     $     +DSQRT(FOUR*V**2*T/PI/R/DX)
     $     *(DEXP(V*X/DX-RMU*T/R-(R*(X+X2)+V*T)**2/FOUR/R/DX/T)
     $     -DEXP(V*X/DX-RMU*T/R-(R*(X+X1)+V*T)**2/FOUR/R/DX/T)))
     $     *RHO*R/DY/T*EXPBI0(R*Y*RHO/TWO/DY/T,
     $     -(R*(Y**2+RHO**2)/FOUR/DY/T))
C
      RETURN
      END
C
C     ----------
      REAL*8 FUNCTION CS3(RHO)
      IMPLICIT REAL*8(A-M,O-Z)
      COMMON/EVAL/ T,X,Y,Z,V,R,TO,RMU,RLAM,DX,DY,DZ,CO,A,B,X1,X2,NGC
      PARAMETER(ZERO=0.0D+00,TWO=2.0D+00,FOUR=4.0D+00)
C
      CS3=RHO*V/DY/X/TWO*EXPBI0(V*Y*RHO/TWO/DY/X,
     $   -(V*(Y**2+RHO**2)/FOUR/DY/X))
C
      RETURN
      END
C
C     ----------	 
      REAL*8 FUNCTION GAMMA2(T)
      IMPLICIT REAL*8(A-M,O-Z)
      COMMON/EVAL/TT,X,Y,Z,V,R,TO,RMU,RLAM,DX,DY,DZ,CO,A,B,X1,X2,NGC
      PARAMETER(ZERO=0.0D+00,FOUR=4.0D+00)
C
      GAMMA2=(EXF(ZERO,(Y-A)/DSQRT(FOUR*DY*T/R))
     $       -EXF(ZERO,(Y+A)/DSQRT(FOUR*DY*T/R)))
     $       *(EXF(ZERO,(Z-B)/DSQRT(FOUR*DZ*T/R))
     $       -EXF(ZERO,(Z+B)/DSQRT(FOUR*DZ*T/R)))
      RETURN
      END
C
C     ----------
      REAL*8 FUNCTION LAM2(T)
      IMPLICIT REAL*8(A-M,O-Z)
      COMMON/EVAL/TT,X,Y,Z,V,R,TO,RMU,RLAM,DX,DY,DZ,CO,A,B,X1,X2,NGC
      PARAMETER(ZERO=0.0D+00,ZIP=1.0D-08,FOUR=4.0D+00)
      PI=4.0D+00*DATAN(1.0D+00)
C
      LAM2=EXF(-RMU*T/R,(V*T-R*X)/DSQRT(FOUR*R*DX*T))
     $     -EXF(V*X/DX-RMU*T/R,(R*X+V*T)/DSQRT(FOUR*R*DX*T))
      RETURN
      END
C			  
C     ----------
      REAL*8 FUNCTION LAM4(T)
      IMPLICIT REAL*8(A-M,O-Z)
      COMMON/EVAL/TT,X,Y,Z,V,R,TO,RMU,RLAM,DX,DY,DZ,CO,A,B,X1,X2,NGC
      PARAMETER(ZERO=0.0D+00,ONE=1.0D+00,TWO=2.0D+00,FOUR=4.0D+00)
      PI=4.0D+00*DATAN(1.0D+00)
C
      LAM4=EXF(-RMU*T/R,(V*T-R*X)/DSQRT(FOUR*R*DX*T))+
     $     (ONE+V*X/DX+V**2*T/R/DX)*EXF(V*X/DX-RMU*T/R,(R*X+V*T)/
     $     DSQRT(FOUR*R*DX*T))-DSQRT(FOUR*V**2*T/PI/R/DX)*
     $     DEXP(-RMU*T/R-(R*X-V*T)**2/FOUR/R/DX/T)
      RETURN
      END
C
C     ----------
      REAL*8 FUNCTION LAM5(T)
      IMPLICIT REAL*8(A-M,O-Z)
      COMMON/EVAL/TT,X,Y,Z,V,R,TO,RMU,RLAM,DX,DY,DZ,CO,A,B,X1,X2,NGC
      PARAMETER(FOUR=4.0D+00)
C
      LAM5=EXF(-RMU*T/R,(R*(X-X2)-V*T)/DSQRT(FOUR*R*DX*T))
     $    -EXF(-RMU*T/R,(R*(X-X1)-V*T)/DSQRT(FOUR*R*DX*T))
     $    +EXF(V*X/DX-RMU*T/R,(R*(X+X2)+V*T)/DSQRT(FOUR*R*DX*T))
     $    -EXF(V*X/DX-RMU*T/R,(R*(X+X1)+V*T)/DSQRT(FOUR*R*DX*T))
      RETURN
      END
C
C     ----------
      REAL*8 FUNCTION LAM6(T)
      IMPLICIT REAL*8(A-M,O-Z)
      COMMON/EVAL/TT,X,Y,Z,V,R,TO,RMU,RLAM,DX,DY,DZ,CO,A,B,X1,X2,NGC
      PARAMETER(ZERO=0.0D+00,ONE=1.0D+00,TWO=2.0D+00,FOUR=4.0D+00)
      PI=4.0D+00*DATAN(1.0D+00)
C
      LAM6=(ONE+V/DX*(X+X1)+V**2*T/R/DX)
     $     *EXF(V*X/DX-RMU*T/R,(R*(X+X1)+V*T)/DSQRT(FOUR*R*DX*T))
     $     -(ONE+V/DX*(X+X2)+V**2*T/R/DX)
     $     *EXF(V*X/DX-RMU*T/R,(R*(X+X2)+V*T)/DSQRT(FOUR*R*DX*T))
     $     +EXF(-RMU*T/R,(R*(X-X2)-V*T)/DSQRT(FOUR*R*DX*T))
     $     -EXF(-RMU*T/R,(R*(X-X1)-V*T)/DSQRT(FOUR*R*DX*T))
     $     +DSQRT(FOUR*V**2*T/PI/R/DX)
     $     *(DEXP(V*X/DX-RMU*T/R-(R*(X+X2)+V*T)**2/FOUR/R/DX/T)
     $     -DEXP(V*X/DX-RMU*T/R-(R*(X+X1)+V*T)**2/FOUR/R/DX/T))
      RETURN
      END
C						 
C     --------------
      FUNCTION EXF(A,B)
C
C     PURPOSE:  TO CALCULATE EXP(A) ERFC(B)
C
      IMPLICIT REAL*8(A-M,O-Z)
      EXF=0.D00
      IF((DABS(A).GT.100.).AND.(B.LE.0.)) RETURN
      C=A-B*B
      IF((DABS(C).GT.100.).AND.(B.GE.0.)) RETURN
      IF(C.LT.-100.) GO TO 3
      X=DABS(B)
      IF(X.GT.3.0) GO TO 1
      T=1./(1.+.3275911*X)
      Y=T*(.2548296-T*(.2844967-T*(1.421414-T*(1.453152-1.061405*T))))
      GO TO 2
    1 Y=.5641896/(X+.5/(X+1./(X+1.5/(X+2./(X+2.5/X+1.)))))
    2 EXF=Y*DEXP(C)
    3 IF(B.LT.0.0) EXF=2.*DEXP(A)-EXF
C
      RETURN
      END
C
C     ----------
      REAL*8 FUNCTION EXPBI0(X,Z)
      IMPLICIT REAL*8(A-M,O-Z)
C
C     PURPOSE: TO CALCULATE EXP(Z)*Io(X) FOR ANY REAL X AND Z
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
     $        (Q6+Y*(Q7+Y*(Q8+Y*Q9))))))))
      ENDIF
      RETURN
      END
C
C   ------------------------------------------------
      SUBROUTINE CHEBY(FUNC,AREA,A,B)
C
C     PURPOSE: PERFORM INTEGRATION OF F(X) BETWEEN A AND B
C              USING M-POINT GAUSS-CHEBYSHEV QUADRATURE FORMULA
      IMPLICIT REAL*8 (A-M,O-Z)
      COMMON/EVAL/ T,X,Y,Z,V,R,TO,RMU,RLAM,DX,DY,DZ,CO,AA,BB,X1,X2,NGC
      PARAMETER(ZERO=0.0D+00,ONE=1.0D+00,TWO=2.0D+00,FOUR=4.0D+00)
      PI=FOUR*DATAN(ONE)
      AREA=ZERO
      SUM=ZERO
      DO 10 I=1,NGC
        Z1=DCOS(DFLOAT(2*(I-1)+1)*PI/DFLOAT(2*NGC))
        TAU=(Z1*(B-A)+B+A)/TWO
   10 SUM=SUM+FUNC(TAU)*DSQRT(ONE-Z1*Z1)
      AREA=(B-A)*PI*SUM/DFLOAT(2*NGC)
C
      RETURN
      END
C
C   ------------------------------------------------
      SUBROUTINE CHEONE(FUNC,T,AREA,A,B)
C
C     PURPOSE: PERFORM FIRST INTEGRATION OF F(X,T) BETWEEN A AND B
C              USING M-POINT GAUSS-CHEBYSHEV QUADRATURE FORMULA
      IMPLICIT REAL*8 (A-M,O-Z)
      COMMON/EVAL/TT,X,Y,Z,V,R,TO,RMU,RLAM,DX,DY,DZ,CO,AA,BB,X1,X2,NGC
      PARAMETER(ZERO=0.0D+00,ONE=1.0D+00,TWO=2.0D+00,FOUR=4.0D+00)
      PI=FOUR*DATAN(ONE)
      AREA=ZERO
      SUM=ZERO
      DO 10 I=1,NGC
        Z1=DCOS(DFLOAT(2*(I-1)+1)*PI/DFLOAT(2*NGC))
        RHO=(Z1*(B-A)+B+A)/TWO
   10 SUM=SUM+FUNC(RHO,T)*DSQRT(ONE-Z1*Z1)
      AREA=(B-A)*PI*SUM/DFLOAT(2*NGC)
C
      RETURN
      END
