*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
*                                                                      *
*     SWMS_2D  - Numerical model of two-dimensional flow and solute    *
*                transport in a variably saturated porous medium       *
*                Conjugate gradient solver for symmetric matrix        *
*                ORTHOMIN solver for asymmetric matrix                 *
*                version 1.22                                          *
*                                                                      *
*     Updated by J.Simunek (1994)                                      *
*     Based on model SWMS_2D (Simunek et al., 1992)                    *
*                                                                      *
*                                         Last modified: January, 1996 *
*                                                                      *
*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*

      program SWMS_2D

      parameter (NumNPD=5000,
     !           NumElD=6000,
     !           NumBPD=250,
     !           MBandD=20,
     !           NSeepD=2,
     !           NumSPD=50,
     !           NDrD=2,
     !           NElDrD=8,
     !           NMatD=20,
     !           NTabD=100,
     !           NumKD=6,
     !           NObsD=4,
     !           MNorth=4)

      double precision A,B
c      double precision RTime1,RTime2
      double precision A1,B1,VRV,RES,RQI,RQ,QQ,QI,RQIDOT,ECNVRG,
     !                 RCNVRG,ACNVRG
      logical lWat,lChem,SinkF,qGWLF,AtmInF,ShortF,SeepF,CheckF,FluxF,
     !        Explic,lUpW,FreeD,DrainF,lArtD,lOrt,lMinStep
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
     !  Sink(NumNPD),POptm(NMatD),Beta(NumNPD),DS(NumNPD),CumQ(NumKD),
     !  vMean(NumKD),hMean(NumKD),KodCB(NumBPD),Qc(NumNPD),Vx(NumNPD),
     !  Vz(NUmNPD),ChPar(10,NMatD),Dispxx(NumNPD),Dispzz(NumNPD),
     !  Dispxz(NumNPD),cBound(6),Ac(NumNPD),Fc(NumNPD),SolIn(NumElD),
     !  Conc(NumNPD),SMean(NumKD),ChemS(NumKD),WeTab(3,2*NumElD),
     !  Gc(NumNPD),ND(NDrD),NED(NDrD),EfDim(2,NDrD),KElDr(NDrD,NElDrD),
     !  ConO(NumNPD),Node(NObsD),
     !  B1(NumNPD),IAD(MBandD,NumNPD),IADN(NumNPD),IADD(NumNPD),
     !  A1(MBandD,NumNPD),RES(NumNPD),VRV(NumNPD),RQI(NumNPD,MNorth),
     !  RQ(NumNPD),QQ(NumNPD),RQIDOT(MNorth),QI(NumNPD,MNorth)

      open(30,file='SWMS_2D.IN\Selector.in',  status='old')
      open(32,file='SWMS_2D.IN\Grid.in',      status='old')
      open(50,file='SWMS_2D.OUT\Check.out',   status='unknown')
      open(71,file='SWMS_2D.OUT\v_Mean.out',  status='unknown')
      open(72,file='SWMS_2D.OUT\A_Level.out', status='unknown')
      open(75,file='SWMS_2D.OUT\h.out',       status='unknown')
      open(76,file='SWMS_2D.OUT\th.out',      status='unknown')
      open(77,file='SWMS_2D.OUT\h_Mean.out',  status='unknown')
      open(78,file='SWMS_2D.OUT\Cum_Q.out',   status='unknown')
      open(79,file='SWMS_2D.OUT\Boundary.out',status='unknown')
      open(80,file='SWMS_2D.OUT\Balance.out', status='unknown')
      open(81,file='SWMS_2D.OUT\vz.out',      status='unknown')
      open(82,file='SWMS_2D.OUT\vx.out',      status='unknown')
      open(92,file='SWMS_2D.OUT\ObsNod.out',  status='unknown')

      data SinkF  ,qGWLF  ,tInit,NTab,ItCum,Iter,TLevel,ALevel,PLevel
     !    /.false.,.false., 0.  ,100 , 0   , 0  ,   1  ,  1   ,   1  /
     !
     !     CumQ    ,Sink     ,CumQrT,CumQrR,CumQvR,ChemS   ,rRoot,rTop
     !    /NumKD*0.,NumNPD*0.,0.    ,0.    ,0.    ,NumKD*0.,  0. , 0./
     !
     !     CumCh0,CumCh1,CumChR,dtMaxC,wCumA,cCumA,Explic ,lMinStep
     !    /   0. ,   0. ,   0. ,1.e+30,  0. ,  0. ,.false.,.false./

      data ECNVRG,ACNVRG,RCNVRG,MaxItO
     !    /1.0d-8,1.0d-8,1.0d-8, 200 /

* --- Reading of the input files and initial calculations --------------

      call BasInf (KAT,MaxIt,TolTh,TolH,lWat,lChem,AtmInF,ShortF,SeepF,
     !             CheckF,FluxF,FreeD,DrainF)
      call NodInf (NumNP,NumEl,IJ,NumBP,NumNPD,NumElD,NumBPD,NumKD,NObs,
     !             NObsD,Kode,Q,Conc,hNew,hOld,hTemp,x,y,MatNum,Beta,
     !             Axz,Bxz,Dxz,CheckF)
      call ElemIn (NumEl,NumElD,NumNP,KX,LayNum,ConAxx,ConAzz,ConAxz,
     !             CheckF,ListNE,IJ,MBand,MBandD,lChem,lOrt)
      call GeomIn (NumKD,NumNP,NumBP,NObs,NObsD,SWidth,Width,Kode,KXB,
     !             rLen,Node)
      call IADMake(KX,NumNP,NumEl,NumElD,MBandD,IAD,IADN,IADD)
      close(32)
      call MatIn  (NMatD,NMat,NLay,Par,hTab(1),hTab(NTab))
      call GenMat (NTab,NTabD,NMat,thR,hSat,Par,hTab,ConTab,CapTab,
     !             ConSat,TheTab,thSat)
      call SetMat (NumNP,NTab,NTabD,NMat,hTab,ConTab,CapTab,hNew,hOld,
     !             MatNum,Par,Con,Cap,ConSat,Axz,Bxz,Dxz,hSat,hTemp,
     !             Explic,TheTab,thSat,thR,ThOld)
      if(AtmInF) then
        open(31,file='SWMS_2D.IN\Atmosph.in', status='old')
        call AtmIn (GWL0L,SinkF,qGWLF,tInit,tMax,Aqh,Bqh,hCritS,MaxAL)
        call SetAtm(tAtm,rTop,rRoot,hCritA,Width,KXB,NumBP,Kode,hNew,Q,
     !              NumNP,GWL0L,qGWLF,FreeD,cPrec,cht,crt,lMinStep)
      end if
      call TmIn    (tInit,tMax,tAtm,tOld,dt,dtMax,dMul,dMul2,dtMin,
     !              TPrint,t,dtOpt,dtOld,AtmInF)
      dtInit=dt
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
        call ChemIn(NMat,NumBP,cBound,ChPar,epsi,tPulse,KodCB,NLevel,
     !              lUpW,lArtD,PeCr)
        if(lWat)
     !  call ChInit(NumNP,NumEl,NumElD,NMat,x,y,KX,MatNum,NLevel,Con,
     !              hNew,Sink,cBound(5),Vx,Vz,ConAxx,ConAzz,ConAxz,
     !              Dispxx,Dispzz,Dispxz,ChPar,ThOld,thSat,Conc,Fc,Gc,
     !              ListNE,lUpW,WeTab,dt,dtMaxC,Peclet,Courant,KAT,
     !              lArtD,PeCr,ConO)
        open(83,file='SWMS_2D.OUT\Conc.out'  ,status='unknown')
        open(74,file='SWMS_2D.OUT\Solute.out',status='unknown')
        call cOut  (NumNP,Conc,x,y,tInit,IJ)
      end if
      close(30)
      close(50)
      open(73,file='SWMS_2D.OUT\Q.out',       status='unknown')
      open(70,file='SWMS_2D.OUT\Run_Inf.out', status='unknown')

      call hOut  (hNew,x,y,NumNP,tInit,IJ)
      call thOut (ThOld,x,y,NumNP,tInit,IJ)
      call SubReg(NumEl,NumElD,NumNP,NMat,hNew,ThOld,ThOld,x,y,MatNum,
     !            LayNum,KX,KAT,tInit,dt,NLay,0,lWat,lChem,Conc,ChPar,
     !            wCumA,wCumT,cCumA,cCumT,wVolI,cVolI,WatIn,SolIn)
      if(NObs.gt.0) call ObsNod(tInit,NumNP,NObs,NObsD,Node,hNew,ThOld,
     !                          Conc)

      write(*,*)'beginning of numerical solution'
c      call getdat(i,i,iday)
c      call gettim(ihours,mins,isecs,i)
c      Rtime1=iday*24.*60.*60.+ihours*60.*60.+mins*60.+isecs

* --- Beginning of time loop -------------------------------------------

11    continue

*     Calculate water flow
      if(lWat.or.TLevel.eq.1)
     !  call WatFlow(NumNP,NumEl,NumElD,NTab,NTabD,MBand,MBandD,NMat,
     !              NSeep,NSeepD,NumSPD,NSP,NP,NumBP,ItCum,MaxIt,Iter,
     !              Kode,KAT,t,dt,dtMin,dtOpt,dtOld,tOld,hCritA,hCritS,
     !              TolTh,TolH,rLen,Width,rTop,vMeanR,hMeanR,AtmInf,
     !              SinkF,SeepF,qGWLF,FreeD,Par,hTab,ConTab,CapTab,
     !              TheTab,hNew,hOld,hTemp,thR,thSat,ThNew,ThOld,MatNum,
     !              Con,Cap,ConSat,Axz,Bxz,Dxz,hSat,A,B,Q,F,x,y,KX,Sink,
     !              DS,Beta,ConAxx,ConAzz,ConAxz,KXB,Explic,GWL0L,Aqh,
     !              Bqh,lWat,TLevel,lOrt,DrainF,ND,NDr,NDrD,rRoot,P0,
     !              POptm,P2H,P2L,P3,r2H,r2L,ConO,
     !              A1,B1,NumNPD,IAD,IADN,IADD,VRV,RES,RQI,RQ,QQ,QI,
     !              RQIDOT,ECNVRG,RCNVRG,ACNVRG,MNorth,MaxItO)
      if(.not.lWat.and.TLevel.eq.1) then
        if(lChem) then
          call ChInit(NumNP,NumEl,NumElD,NMat,x,y,KX,MatNum,NLevel,Con,
     !                hNew,Sink,cBound(5),Vx,Vz,ConAxx,ConAzz,ConAxz,
     !                Dispxx,Dispzz,Dispxz,ChPar,ThNew,thSat,Conc,Fc,Gc,
     !                ListNE,lUpW,WeTab,dt,dtMaxC,Peclet,Courant,KAT,
     !                lArtD,PeCr,ConO)
        else
          call Veloc(KAT,NumNP,NumEl,NumElD,hNew,x,y,KX,ListNE,Con,
     !               ConAxx,ConAzz,ConAxz,Vx,Vz)
        end if
        Iter=1
      end if

*     Calculate solute transport
      if(lChem)
     !  call Solute(NumNP,NumEl,NumElD,MBand,MBandD,NMat,t,Kode,A,B,Q,
     !              hNew,hOld,F,x,y,KX,KAT,dt,DS,Sink,MatNum,Con,ConO,
     !              ConAxx,ConAzz,ConAxz,Vx,Vz,Dispxx,Dispzz,Dispxz,
     !              ChPar,ThNew,ThOld,thSat,Ac,Fc,Gc,Qc,Conc,ListNE,
     !              cBound,tPulse,NumBP,KodCB,KXB,NLevel,cPrec,crt,cht,
     !              lWat,lUpW,WeTab,epsi,CumCh0,CumCh1,CumChR,dtMaxC,
     !              Peclet,Courant,lArtD,PeCr,lOrt,
     !              A1,B1,NumNPD,IAD,IADN,IADD,VRV,RES,RQI,RQ,QQ,QI,
     !              RQIDOT,ECNVRG,RCNVRG,ACNVRG,MNorth,MaxItO)

*     T-Level information
      call TLInf   (NumNP,NumBP,Kode,Q,hNew,CumQ,Width,SWidth,KXB,t,dt,
     !              TLevel,ShortF,TPrint(PLevel),Iter,ItCum,rTop,rRoot,
     !              vMeanR,hMeanT,hMeanR,hMeanG,AtmInF,SinkF,CumQrT,
     !              CumQrR,CumQvR,NumKD,hMean,vMean,lWat,lChem,rLen,
     !              Peclet,Courant,wCumT,wCumA)
      if(lChem)
     !  call SolInf(NumNP,Kode,Qc,t,dt,TLevel,ShortF,TPrint(PLevel),
     !              NumKD,SMean,ChemS,CumCh0,CumCh1,CumChR,cCumA,cCumT,
     !              lWat)
      if(NObs.gt.0) call ObsNod(t,NumNP,NObs,NObsD,Node,hNew,ThNew,Conc)

*     P-Level information
      if(abs(TPrint(PLevel)-t).lt.0.001*dt.or.
     !   (.not.lWat.and..not.lChem)) then
        if(lWat.or.(.not.lWat.and.PLevel.eq.1)) then
          call hOut (hNew,x,y,NumNP,t,IJ)
          call thOut(ThNew,x,y,NumNP,t,IJ)
          if(FluxF) then
            if(.not.lChem)
     !        call Veloc(KAT,NumNP,NumEl,NumElD,hNew,x,y,KX,ListNE,Con,
     !                   ConAxx,ConAzz,ConAxz,Vx,Vz)
            call FlxOut(Vx,Vz,x,y,NumNP,t,IJ)
            call QOut  (Q,x,y,NumNP,t,IJ)
          end if
        end if
        call SubReg(NumEl,NumElD,NumNP,NMat,hNew,ThOld,ThNew,x,y,MatNum,
     !              LayNum,KX,KAT,t,dt,NLay,PLevel,lWat,lChem,Conc,
     !              ChPar,wCumA,wCumT,cCumA,cCumT,wVolI,cVolI,WatIn,
     !              SolIn)
        call BouOut(NumNP,NumBP,t,hNew,ThNew,Q,Width,KXB,Kode,x,y,Conc)
        if(lChem)
     !    call cOut(NumNP,Conc,x,y,t,IJ)
        PLevel=PLevel+1
      end if

*     A-level information
      if(abs(t-tAtm).le.0.001*dt.and.AtmInF) then
        if(lWat)
     !    call ALInf (t,CumQ,hMeanT,hMeanR,hMeanG,ALevel,CumQrT,CumQrR,
     !                CumQvR,NumKD)
        if(ALevel.lt.MaxAL) then
          call SetAtm(tAtm,rTop,rRoot,hCritA,Width,KXB,NumBP,Kode,hNew,
     !                Q,NumNP,GWL0L,qGWLF,FreeD,cPrec,cht,crt,lMinStep)
          ALevel=ALevel+1
        end if
      end if

*     Root extraction
*      if(SinkF)
*     !  call SetSnk(NumNP,NMat,MatNum,hNew,rRoot,Sink,P0,POptm,P2H,P2L,
*     !              P3,r2H,r2L,Beta,rLen)

*     Time governing
      if(abs(t-tMax).le.0.001*dt.or.(.not.lWat.and..not.lChem)) then
c        call getdat(i,i,iday)
c        call gettim(ihours,mins,isecs,i)
c        Rtime2=iday*24.*60.*60.+ihours*60.*60.+mins*60.+isecs
c        write(70,*)
c        write(70,*) 'Real time [sec]',Rtime2-RTime1
c        write( *,*) 'Real time [sec]',Rtime2-RTime1
        stop
      end if
      tOld=t
      dtOld=dt
      call TmCont(dt,dtMax,dtOpt,dMul,dMul2,dtMin,Iter,TPrint(PLevel),
     !            tAtm,t,tMax,dtMaxC,lMinStep,dtInit)
      TLevel=TLevel+1
      t=t+dt
      goto 11

* --- end of time loop -------------------------------------------------

      end

*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||