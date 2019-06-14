*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
*                                                                      *
*     CHAIN_2D - Numerical model of two-dimensional variably saturated *
*                water flow, heat transport, and transport of solutes  *
*                involved in sequential first-order decay reactions.   *
*                Conjugate gradient solver for symmetric matrix        *
*                ORTHOMIN solver for asymmetric matrix                 *
*                version 1.1                                           *
*                                                                      *
*     Designed by J.Simunek and M. Th. van Genuchten (1994)            *
*     Based on model SWMS_2D (Simunek, Vogel, and van Genuchten, 1992) *
*                                                                      *
*                                         Last modified: October, 1994 *
*                                                                      *
*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*

      program CHAIN_2D

      parameter (NumNPD=12000,
     !           NumElD=12000,
     !           MBandD=20,
     !           NumBPD=300,
     !           NSeepD=2,
     !           NumSPD=50,
     !           NDrD=2,
     !           NElDrD=8,
     !           NMatD=15,
     !           NTabD=100,
     !           NumKD=6,
     !           NSD=6,
     !           NObsD=4,
     !           MNorth=4)

      double precision A,B,RTime1,RTime2
      double precision A1,B1,VRV,RES,RQI,RQ,QQ,QI,RQIDOT,ECNVRG,
     !                 RCNVRG,ACNVRG
      logical lWat,lChem,SinkF,qGWLF,FreeD,AtmInF,ShortF,SeepF,CheckF,
     !        FluxF,DrainF,Explic,lUpW,lTemp,lWDep,lLinear,lEquil,lArtD,
     !        lOrt
      integer PLevel,ALevel,TLevel

      dimension A(MBandD,NumNPD),B(NumNPD),Kode(NumNPD),Q(NumNPD),
     !  hNew(NumNPD),hTemp(NumNPD),hOld(NumNPD),ConSat(NMatD),F(NumNPD),
     !  hTab(NTabD),ConTab(NTabD,NMatD),CapTab(NTabD,NMatD),Con(NumNPD),
     !  Cap(NumNPD),x(NumNPD),y(NumNPD),MatNum(NumNPD),LayNum(NumElD),
     !  KX(NumElD,4),KXB(NumBPD),TPrint(50),Par(10,NMatD),Width(NumBPD),
     !  ConAxx(NumElD),ConAzz(NumElD),ConAxz(NumElD),SWidth(NumKD),
     !  NP(NSeepD,NumSPD),NSP(NSeepD),hSat(NMatD),WatIn(NumElD),
     !  Axz(NumNPD),Bxz(NumNPD),Dxz(NumNPD),thR(NMatD),thSat(NMatD),
     !  TheTab(NTabD,NMatD),ThNew(NumNPD),ThOld(NumNPD),ListNE(NumNPD),
     !  Sink(NumNPD),POptm(NMatD),Beta(NumNPD),DS(NumNPD),cNew(NumNPD),
     !  VxO(NumNPD),VzO(NumNPD),VxN(NumNPD),VzN(NumNPD),cTemp(NumNPD),
     !  Dispxx(NumNPD),Dispzz(NumNPD),Dispxz(NumNPD),CumQ(NumKD),
     !  vMean(NumKD),hMean(NumKD),Conc(NSD,NumNPD),WeTab(3,2*NumElD),
     !  cBound(NSD,9),Ac(NumNPD),Qc(NSD,NumNPD),Fc(NumNPD),Gc(NumNPD),
     !  KodCB(NumBPD),cPrevO(NumNPD),SMean(NSD,NumKD),ChemS(NSD,NumKD),
     !  ChPar(NSD*16+4,NMatD),SolIn(NumElD),CumCh0(NSD),CumCh1(NSD),
     !  CumChR(NSD),cPrec(NSD),cht(NSD),crt(NSD),cCumA(NSD),cCumT(NSD),
     !  cVolI(NSD),TDep(NSD*16+4),TempN(NumNPD),TempO(NumNPD),TBound(6),
     !  TPar(10,NMatD),KodTB(NumBPD),SorbN(NumNPD),Sorb(NSD,NumNPD),
     !  Node(NObsD),ND(NDrD),NED(NDrD),EfDim(2,NDrD),KElDr(NDrD,NElDrD),
     !  ConO(NumNPD),
     !  B1(NumNPD),IAD(MBandD,NumNPD),IADN(NumNPD),IADD(NumNPD),
     !  A1(MBandD,NumNPD),RES(NumNPD),VRV(NumNPD),RQI(NumNPD,MNorth),
     !  RQ(NumNPD),QQ(NumNPD),RQIDOT(MNorth),QI(NumNPD,MNorth)

      open(30,file='CHAIN_2D.IN\Selector.in',  status='old')
      open(32,file='CHAIN_2D.IN\Grid.in',      status='old')
      open(50,file='CHAIN_2D.OUT\Check.out',   status='unknown')
      open(70,file='CHAIN_2D.OUT\Run_Inf.out', status='unknown')
      open(71,file='CHAIN_2D.OUT\v_Mean.out',  status='unknown')
      open(72,file='CHAIN_2D.OUT\A_Level.out', status='unknown')
      open(73,file='CHAIN_2D.OUT\Q.out',       status='unknown')
      open(75,file='CHAIN_2D.OUT\h.out',       status='unknown')
      open(76,file='CHAIN_2D.OUT\th.out',      status='unknown')
      open(77,file='CHAIN_2D.OUT\h_Mean.out',  status='unknown')
      open(78,file='CHAIN_2D.OUT\Cum_Q.out',   status='unknown')
      open(79,file='CHAIN_2D.OUT\Boundary.out',status='unknown')
      open(80,file='CHAIN_2D.OUT\Balance.out', status='unknown')
      open(81,file='CHAIN_2D.OUT\vz.out',      status='unknown')
      open(82,file='CHAIN_2D.OUT\vx.out',      status='unknown')
      open(92,file='CHAIN_2D.OUT\ObsNod.out',  status='unknown')
            

      data SinkF  ,qGWLF  ,NTab,ItCum,IterW,IterC,TLevel,ALevel,PLevel
     !    /.false.,.false.,100 , 0   , 0   , 0   ,  1  ,  1    , 1   /
     !
     !     tInit,CumQ    ,Sink     ,CumQrT,CumQrR,CumQvR,rRoot,rTop
     !    / 0.  ,NumKD*0.,NumNPD*0.,0.    ,0.    ,0.    ,  0. , 0./
     !
     !     CumCh0,CumCh1,CumChR,dtMaxC,wCumA,cCumA ,Explic
     !    /NSD*0.,NSD*0.,NSD*0.,1.e+30,  0. ,NSD*0.,.false./

      data ECNVRG,ACNVRG,RCNVRG,MaxItO
     !    /1.0d-6,1.0d-6,1.0d-6, 200 /

* --- Reading of the input files and initial calculations --------------

      call BasInf (KAT,MaxIt,TolTh,TolH,lWat,lChem,AtmInF,ShortF,SeepF,
     !             CheckF,FluxF,FreeD,DrainF,lTemp,lWDep,lEquil)
      call NodInf (NumNP,NumEl,IJ,NumBP,NumNPD,NumElD,NumBPD,NumKD,NS,
     !             NSD,NObs,NObsD,Kode,Q,Conc,hNew,hOld,hTemp,x,y,
     !             MatNum,Beta,Axz,Bxz,Dxz,CheckF,TempN,TempO,Sorb,
     !             lEquil)
      call ElemIn (NumEl,NumElD,NumNP,KX,LayNum,ConAxx,ConAzz,ConAxz,
     !             CheckF,ListNE,IJ,MBand,MBandD,lChem,lOrt)
      call GeomIn (NumKD,NumNP,NumBP,NObs,NObsD,NS,SWidth,Width,Kode,
     !             KXB,rLen,Node)
      if(lOrt) call IADMake(KX,NumNP,NumEl,NumElD,MBandD,IAD,IADN,IADD)
      close(32)
      call MatIn  (NMatD,NMat,NLay,Par,hTab(1),hTab(NTab))
      call GenMat (NTab,NTabD,NMat,thR,hSat,Par,hTab,ConTab,CapTab,
     !             ConSat,TheTab,thSat)
      call SetMat (NumNP,NTab,NTabD,NMat,hTab,ConTab,CapTab,hNew,hOld,
     !             MatNum,Par,Con,Cap,ConSat,Axz,Bxz,Dxz,hSat,hTemp,
     !             Explic,TheTab,thSat,thR,ThOld,lWDep,TempN,IterW)
      if(lWat.and.(lChem.or.lTemp))
     !  call Veloc(KAT,NumNP,NumEl,NumElD,hOld,x,y,KX,ListNE,Con,ConAxx,
     !             ConAzz,ConAxz,VxO,VzO)
      if(AtmInF) then
        open(31,file='CHAIN_2D.IN\Atmosph.in', status='old')
        call AtmIn (GWL0L,SinkF,qGWLF,tInit,tMax,Aqh,Bqh,hCritS,MaxAL)
        call SetAtm(tAtm,rTop,rRoot,hCritA,Width,KXB,NumBP,Kode,hNew,Q,
     !              NumNP,GWL0L,qGWLF,FreeD,lTemp,lChem,cPrec,cht,crt,
     !              Th3,Th4,NS)
      end if
      call TmIn    (tInit,tMax,tAtm,tOld,dt,dtMax,dMul,dMul2,dtMin,
     !              TPrint,t,dtOpt,dtOld,AtmInF)
      if(SinkF) then
        call SinkIn(NMat,NumEl,NumNP,NumElD,KAT,KX,x,y,P0,POptm,P2H,P2L,
     !              P3,r2H,r2L,Beta)
        call SetSnk(NumNP,NMat,MatNum,hNew,rRoot,Sink,P0,POptm,P2H,P2L,
     !              P3,r2H,r2L,Beta,rLen)
      end if
      if(SeepF)
     !  call SeepIn(NSeepD,NumSPD,NSeep,NSP,NP)
      if(DrainF)
     !  call DrainIn(NDr,NDrD,NElDrD,NumEl,ND,NED,KElDr,EfDim,ConAxx,
     !               ConAxz,ConAzz)
      if(lChem) then
        call ChemIn(NMat,NumBP,NS,NSD,NumKD,NLevel,cBound,ChPar,epsi,
     !              tPulse,KodCB,ChemS,TDep,lUpW,cTolA,cTolR,MaxItC,
     !              lLinear,lEquil,lArtD,PeCr)
        open(74,file='CHAIN_2D.OUT\Conc.out'   ,status='unknown')
        open(84,file='CHAIN_2D.OUT\Solute1.out',status='unknown')
        if(NS.ge.2)
     !    open(85,file='CHAIN_2D.OUT\Solute2.out',status='unknown')
        if(NS.ge.3)
     !    open(86,file='CHAIN_2D.OUT\Solute3.out',status='unknown')
        if(NS.ge.4)
     !    open(87,file='CHAIN_2D.OUT\Solute4.out',status='unknown')
        if(NS.ge.5)
     !    open(88,file='CHAIN_2D.OUT\Solute5.out',status='unknown')
        if(NS.ge.6)
     !    open(89,file='CHAIN_2D.OUT\Solute6.out',status='unknown')
        do 14 j=1,NS
          call cOut(j,NumNP,Conc,x,y,tInit,IJ,NSD)
          if(.not.lEquil) call cOut(j,NumNP,Sorb,x,y,tInit,IJ,NSD)
14      continue
      end if
      if(lTemp) then
        call TempIn(NMat,NumBP,TPar,KodTB,TBound,Ampl,tPeriod)
        open(83,file='CHAIN_2D.OUT\Temp.out',    status='unknown')
        call TOut(NumNP,TempN,x,y,tInit,IJ)
      end if
      close(30)
      close(50)

      call hOut  (hNew,x,y,NumNP,tInit,IJ)
      call thOut (ThOld,x,y,NumNP,tInit,IJ)
      call SubReg(NumEl,NumElD,NumNP,NMat,hNew,ThOld,ThOld,x,y,MatNum,
     !            LayNum,KX,KAT,tInit,dt,NLay,0,lWat,lChem,Conc,ChPar,
     !            wCumA,wCumT,cCumA,cCumT,wVolI,cVolI,WatIn,SolIn,thSat,
     !            lTemp,TempN,TPar,TDep,NS,NSD,Sorb,lLinear,lEquil)
      if(NObs.gt.0)
     !  call ObsNod(tInit,NumNP,NS,NSD,NObs,NObsD,Node,hNew,ThOld,Conc,
     !              Sorb,TempN)

      write(*,*)'beginning of numerical solution'
      call getdat(i,i,iday)
      call gettim(ihours,mins,isecs,i)
      Rtime1=iday*24.*60.*60.+ihours*60.*60.+mins*60.+isecs
c      call seconds_since_1980@(RTime1)

* --- Beginning of time loop -------------------------------------------

11    continue

*     Calculate water flow
      if(lWat.or.TLevel.eq.1)
     !call WatFlow(NumNP,NumEl,NumElD,NTab,NTabD,MBand,MBandD,NMat,
     !             NSeep,NSeepD,NumSPD,NSP,NP,NumBP,ItCum,MaxIt,IterW,
     !             Kode,KAT,t,dt,dtMin,dtOpt,dtOld,tOld,hCritA,hCritS,
     !             TolTh,TolH,rLen,Width,rTop,vMeanR,hMeanR,AtmInf,
     !             SinkF,SeepF,qGWLF,FreeD,Par,hTab,ConTab,CapTab,
     !             TheTab,hNew,hOld,hTemp,thR,thSat,ThNew,ThOld,MatNum,
     !             Con,Cap,ConSat,Axz,Bxz,Dxz,hSat,A,B,Q,F,x,y,KX,Sink,
     !             DS,Beta,ConAxx,ConAzz,ConAxz,KXB,Explic,GWL0L,Aqh,
     !             Bqh,lWat,lWDep,TempN,lOrt,lChem,lTemp,TLevel,VxO,VzO,
     !             VxN,VzN,ListNE,rRoot,P0,POptm,P2H,P2L,P3,r2H,r2L,
     !             DrainF,ND,NDr,NDrD,ConO,
     !             A1,B1,NumNPD,IAD,IADN,IADD,VRV,RES,RQI,RQ,QQ,QI,
     !             RQIDOT,ECNVRG,RCNVRG,ACNVRG,MNorth,MaxItO)

*     Calculate heat transport
      if(lTemp) 
     !  call Temper(NumNP,NumEl,NumElD,MBand,MBandD,NMat,NumBP,t,A,B,Q,
     !              x,y,KX,KAT,dt,DS,MatNum,VxO,VzO,VxN,VzN,Dispxx,
     !              Dispzz,Dispxz,TPar,ThNew,ThOld,Ac,TempO,TempN,
     !              TBound,KodTB,KXB,Kode,Th3,Th4,Ampl,tPeriod,lOrt,
     !              A1,B1,NumNPD,IAD,IADN,IADD,VRV,RES,RQI,RQ,QQ,QI,
     !              RQIDOT,ECNVRG,RCNVRG,ACNVRG,MNorth,MaxItO)

*     Calculate solute transport
      if(lChem)
     !  call Solute(NumNP,NumEl,NumElD,MBand,MBandD,NMat,NS,NSD,t,Kode,
     !              A,B,Q,F,x,y,KX,KAT,dt,DS,Sink,MatNum,VxO,VzO,VxN,
     !              VzN,Dispxx,Dispzz,Dispxz,ChPar,ThNew,ThOld,thSat,Ac,
     !              Fc,Gc,Qc,Conc,cBound,tPulse,NumBP,KodCB,KXB,NLevel,
     !              cPrec,crt,cht,lUpW,WeTab,epsi,CumCh0,CumCh1,CumChR,
     !              dtMaxC,Peclet,Courant,cPrevO,TempO,TempN,TDep,Width,
     !              cTolA,cTolR,MaxItC,IterC,cTemp,cNew,Sorb,SorbN,
     !              hTemp,Cap,lLinear,lEquil,lArtD,PeCr,lOrt,
     !              A1,B1,NumNPD,IAD,IADN,IADD,VRV,RES,RQI,RQ,QQ,QI,
     !              RQIDOT,ECNVRG,RCNVRG,ACNVRG,MNorth,MaxItO)

* --- Printing of output information -----------------------------------
*     T-Level information
      call TLInf   (NumNP,NumBP,Kode,Q,hNew,CumQ,Width,SWidth,KXB,t,dt,
     !              TLevel,ShortF,TPrint(PLevel),IterW,IterC,ItCum,rTop,
     !              rRoot,vMeanR,hMeanT,hMeanR,hMeanG,AtmInF,SinkF,
     !              CumQrT,CumQrR,CumQvR,NumKD,hMean,vMean,lWat,lChem,
     !              rLen,Peclet,Courant,wCumT,wCumA)
      if(lChem)
     ! call SolInf(1,NumNP,NS,NSD,Kode,Qc,t,dt,IterC,TLevel,ShortF,
     !             TPrint(PLevel),NumKD,SMean,ChemS,CumCh0,CumCh1,
     !             CumChR,cCumA,cCumT,lWat)
      if(NObs.gt.0)
     !  call ObsNod(t,NumNP,NS,NSD,NObs,NObsD,Node,hNew,ThNew,Conc,Sorb,
     !              TempN)
  
*     P-Level information
      if(abs(TPrint(PLevel)-t).lt.0.001*dt) then
        if(lWat.or.(.not.lWat.and.PLevel.eq.1)) then
          call hOut (hNew,x,y,NumNP,t,IJ)
          call thOut(ThNew,x,y,NumNP,t,IJ)
          if(FluxF) then
            if(.not.lChem.and..not.lTemp)
     !        call Veloc(KAT,NumNP,NumEl,NumElD,hNew,x,y,KX,ListNE,Con,
     !                   ConAxx,ConAzz,ConAxz,VxN,VzN)
            call FlxOut(VxN,VzN,x,y,NumNP,t,IJ)
            call QOut  (Q,x,y,NumNP,t,IJ)
          endif
        end if
        call SubReg(NumEl,NumElD,NumNP,NMat,hNew,ThOld,ThNew,x,y,MatNum,
     !              LayNum,KX,KAT,t,dt,NLay,PLevel,lWat,lChem,Conc,
     !              ChPar,wCumA,wCumT,cCumA,cCumT,wVolI,cVolI,WatIn,
     !              SolIn,thSat,lTemp,TempN,TPar,TDep,NS,NSD,Sorb,
     !              lLinear,lEquil)
        call BouOut(NumNP,NumBP,NS,NSD,t,hNew,ThNew,Q,Width,KXB,Kode,x,
     !              y,TempN,Conc)
        if(lChem) then
          do 15 j=1,NS
            call cOut(j,NumNP,Conc,x,y,t,IJ,NSD)
            if(.not.lEquil) call cOut(j,NumNP,Sorb,x,y,t,IJ,NSD)
15        continue
        end if
        if(lTemp)
     !    call TOut(NumNP,TempN,x,y,t,IJ)
        PLevel=PLevel+1
      endif

*     A-level information
      if(abs(t-tAtm).le.0.001*dt.and.AtmInF) then
        if(lWat)
     !    call ALInf (t,CumQ,hMeanT,hMeanR,hMeanG,ALevel,CumQrT,CumQrR,
     !                CumQvR,NumKD)
        if(lChem)
     !    call SolInf(0,NumNP,NS,NSD,Kode,Qc,t,dt,IterC,ALevel,ShortF,
     !                0.,NumKD,SMean,ChemS,CumCh0,CumCh1,CumChR,cCumA,
     !                cCumT,lWat)
        if(ALevel.lt.MaxAL) then
          call SetAtm(tAtm,rTop,rRoot,hCritA,Width,KXB,NumBP,Kode,hNew,
     !                Q,NumNP,GWL0L,qGWLF,FreeD,lTemp,lChem,cPrec,cht,
     !                crt,Th3,Th4,NS)
          ALevel=ALevel+1
        endif
      endif

*     Root extraction
      if(SinkF)
     !  call SetSnk(NumNP,NMat,MatNum,hNew,rRoot,Sink,P0,POptm,P2H,P2L,
     !              P3,r2H,r2L,Beta,rLen)

*     Time governing
      if(abs(t-tMax).le.0.001*dt) then
        call getdat(i,i,iday)
        call gettim(ihours,mins,isecs,i)
        Rtime2=iday*24.*60.*60.+ihours*60.*60.+mins*60.+isecs
c        call seconds_since_1980@(RTime2)
        write(70,*)
        write(70,*) 'Real time [sec]',Rtime2-RTime1
        write( *,*) 'Real time [sec]',Rtime2-RTime1
        stop
      end if
      tOld=t
      dtOld=dt
      Iter=max0(IterW,IterC)
      call TmCont(dt,dtMax,dtOpt,dMul,dMul2,dtMin,Iter,TPrint(PLevel),
     !            tAtm,t,tMax,dtMaxC)
      TLevel=TLevel+1
      t=t+dt
      goto 11

* --- End of time loop -------------------------------------------------

      end

*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||