*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
*                                                                      *
*     SWMS_3D  - Numerical model of three-dimensional flow and solute  *
*                transport in a variably saturated porous medium       *
*                Conjugate gradient solver for symmetric matrix        *
*                ORTHOMIN solver for asymmetric matrix                 *
*                version 1.0                                           *
*                                                                      *
*     Designed by J.Simunek (1994)                                     *
*     Based on model SWMS_2D (Simunek et al., 1992)                    *
*                                                                      *
*                                         Last modified: August, 1994  *
*                                                                      *
*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*

      program SWMS_3D

      parameter (NumNPD=12500,
     !           NumElD=11000,
     !           NumBPD=750,
     !           MBandD=20,
     !           NSeepD=1,
     !           NumSPD=250,
     !           NDrD=1,
     !           NElDrD=8,
     !           NMatD=5,
     !           NTabD=100,
     !           NumKD=6,
     !           NObsD=5,
     !           MNorth=4)

      double precision A,B,RTime1,RTime2
      double precision A1,B1,VRV,RES,RQI,RQ,QQ,QI,RQIDOT,ECNVRG,
     !                 RCNVRG,ACNVRG
      logical lWat,lChem,SinkF,qGWLF,AtmInF,ShortF,SeepF,CheckF,FluxF,
     !        Explic,lUpW,FreeD,DrainF,lArtD,lOrt
      integer PLevel,ALevel,TLevel

      dimension A(MBandD,NumNPD),B(NumNPD),Kode(NumNPD),Par(10,NMatD),
     !  hNew(NumNPD),hTemp(NumNPD),hOld(NumNPD),ConSat(NMatD),F(NumNPD),
     !  hTab(NTabD),ConTab(NTabD,NMatD),CapTab(NTabD,NMatD),Con(NumNPD),
     !  Cap(NumNPD),x(NumNPD),y(NumNPD),z(NumNPD),MatNum(NumNPD),
     !  LayNum(NumElD),KX(NumElD,9),KXB(NumBPD),TPrint(50),Q(NumNPD),
     !  ConAxx(NumElD),ConAyy(NumElD),ConAzz(NumElD),ConAxy(NumElD),
     !  ConAxz(NumElD),ConAyz(NumElD),Width(NumBPD),SWidth(NumKD),
     !  NP(NSeepD,NumSPD),NSP(NSeepD),hSat(NMatD),WatIn(NumElD),
     !  Axz(NumNPD),Bxz(NumNPD),Dxz(NumNPD),thR(NMatD),thSat(NMatD),
     !  TheTab(NTabD,NMatD),ThNew(NumNPD),ThOld(NumNPD),ListNE(NumNPD),
     !  Sink(NumNPD),POptm(NMatD),Beta(NumNPD),DS(NumNPD),CumQ(NumKD),
     !  vMean(NumKD),hMean(NumKD),KodCB(NumBPD),Vx(NumNPD),Vy(NumNPD),
     !  Vz(NUmNPD),ChPar(10,NMatD),Dispxx(NumNPD),Dispyy(NumNPD),
     !  Dispzz(NumNPD),Dispxy(NumNPD),Dispxz(NumNPD),Dispyz(NumNPD),
     !  Qc(NumNPD),cBound(12),Ac(NumNPD),Fc(NumNPD),SolIn(NumElD),
     !  Conc(NumNPD),SMean(NumKD),ChemS(NumKD),WeTab(6,5*NumElD),
     !  Gc(NumNPD),ND(NDrD),NED(NDrD),EfDim(2,NDrD),KElDr(NDrD,NElDrD),
     !  ConO(NumNPD),Node(NObsD),KNoDr(NDrD,NElDrD),
     !  B1(NumNPD),IAD(MBandD,NumNPD),IADN(NumNPD),IADD(NumNPD),
     !  A1(MBandD,NumNPD),RES(NumNPD),VRV(NumNPD),RQI(NumNPD,MNorth),
     !  RQ(NumNPD),QQ(NumNPD),RQIDOT(MNorth),QI(NumNPD,MNorth)

      open(30,file='SWMS_3D.IN\Selector.in',  status='old')
      open(32,file='SWMS_3D.IN\Grid.in',      status='old')
      open(50,file='SWMS_3D.OUT\Check.out',   status='unknown')
      open(71,file='SWMS_3D.OUT\v_Mean.out',  status='unknown')
      open(72,file='SWMS_3D.OUT\A_Level.out', status='unknown')
      open(75,file='SWMS_3D.OUT\h.out',       status='unknown')
      open(76,file='SWMS_3D.OUT\th.out',      status='unknown')
      open(77,file='SWMS_3D.OUT\h_Mean.out',  status='unknown')
      open(78,file='SWMS_3D.OUT\Cum_Q.out',   status='unknown')
      open(79,file='SWMS_3D.OUT\Boundary.out',status='unknown')
      open(80,file='SWMS_3D.OUT\Balance.out', status='unknown')
      open(81,file='SWMS_3D.OUT\vz.out',      status='unknown')
      open(82,file='SWMS_3D.OUT\vx.out',      status='unknown')
      open(84,file='SWMS_3D.OUT\vy.out',      status='unknown')
      open(92,file='SWMS_3D.OUT\ObsNod.out',  status='unknown')

      data SinkF  ,qGWLF  ,tInit,NTab,ItCum,Iter,TLevel,ALevel,PLevel
     !    /.false.,.false., 0.  ,100 , 0   , 0  ,   1  ,  1   ,   1  /
     !
     !     CumQ    ,Sink     ,CumQrT,CumQrR,CumQvR,ChemS   ,rRoot,rTop
     !    /NumKD*0.,NumNPD*0.,0.    ,0.    ,0.    ,NumKD*0.,  0. , 0./
     !
     !     CumCh0,CumCh1,CumChR,dtMaxC,wCumA,cCumA,Explic
     !    /   0. ,   0. ,   0. ,1.e+30,  0. ,  0. ,.false./

      data ECNVRG,ACNVRG,RCNVRG,MaxItO
     !    /1.0d-6,1.0d-6,1.0d-6, 200 /

* --- Reading of the input files and initial calculations --------------

      call BasInf (MaxIt,TolTh,TolH,lWat,lChem,AtmInF,ShortF,SeepF,
     !             CheckF,FluxF,FreeD,DrainF)
      call NodInf (NumNP,NumEl,IJ,NumBP,NumNPD,NumElD,NumBPD,NumKD,NObs,
     !             NObsD,Kode,Q,Conc,hNew,hOld,hTemp,x,y,z,MatNum,Beta,
     !             Axz,Bxz,Dxz,CheckF)
      call ElemIn (NumEl,NumElD,NumNP,KX,LayNum,ConAxx,ConAyy,ConAzz,
     !             ConAxy,ConAxz,ConAyz,CheckF,ListNE,MBand,MBandD,x,y,
     !             z,lChem,lOrt)
      call GeomIn (NumKD,NumNP,NumBP,NObs,SWidth,Width,Kode,KXB,rLen,
     !             Node)
      if(lOrt)
     !  call IADMake(KX,NumNP,NumEl,NumElD,MBandD,IAD,IADN,IADD)
      close(32)
      call MatIn  (NMatD,NMat,NLay,Par,hTab(1),hTab(NTab))
      call GenMat (NTab,NTabD,NMat,thR,hSat,Par,hTab,ConTab,CapTab,
     !             ConSat,TheTab,thSat)
      call SetMat (NumNP,NTab,NTabD,NMat,hTab,ConTab,CapTab,hNew,hOld,
     !             MatNum,Par,Con,Cap,ConSat,Axz,Bxz,Dxz,hSat,hTemp,
     !             Explic,TheTab,thSat,thR,ThOld)
      if(AtmInF) then
        open(31,file='SWMS_3D.IN\Atmosph.in', status='old')
        call AtmIn (GWL0L,SinkF,qGWLF,tInit,tMax,Aqh,Bqh,hCritS,MaxAL)
        call SetAtm(tAtm,rTop,rRoot,hCritA,Width,KXB,NumBP,Kode,hNew,Q,
     !              NumNP,GWL0L,qGWLF,FreeD,cPrec,cht,crt)
      end if
      call TmIn    (tInit,tMax,tAtm,tOld,dt,dtMax,dMul,dMul2,dtMin,
     !              TPrint,t,dtOpt,dtOld,AtmInF)
      if(SinkF) then
        call SinkIn(NMat,NumEl,NumNP,NumElD,KX,x,y,z,P0,POptm,P2H,P2L,
     !              P3,r2H,r2L,Beta)
        call SetSnk(NumNP,NMat,MatNum,hNew,rRoot,Sink,P0,POptm,P2H,P2L,
     !              P3,r2H,r2L,Beta,rLen)
      end if
      if(SeepF)
     !  call SeepIn(NSeepD,NumSPD,NSeep,NSP,NP)
      if(DrainF)
     !  call DrainIn(NDr,NDrD,NElDrD,NumEl,ND,NED,KElDr,KNoDr,EfDim,
     !               ConAxx,ConAyy,ConAzz,ConAxy,ConAxz,ConAyz)
      if(lChem) then
        call ChemIn(NMat,NumBP,cBound,ChPar,epsi,tPulse,KodCB,NLevel,
     !              lUpW,lArtD,PeCr)
        if(lWat)
     !  call ChInit(NumNP,NumEl,NumElD,NMat,x,y,z,KX,MatNum,NLevel,Con,
     !              hNew,Sink,cBound(12),Vx,Vy,Vz,ConAxx,ConAyy,ConAzz,
     !              ConAxy,ConAxz,ConAyz,Dispxx,Dispyy,Dispzz,Dispxy,
     !              Dispxz,Dispyz,ChPar,ThOld,thSat,Conc,Fc,Gc,ListNE,
     !              lUpW,WeTab,dt,dtMaxC,Peclet,Courant,lArtD,PeCr,
     !              PeCrMx,ConO)
        open(83,file='SWMS_3D.OUT\Conc.out'  ,status='unknown')
        open(74,file='SWMS_3D.OUT\Solute.out',status='unknown')
        call cOut  (NumNP,Conc,x,y,z,tInit,IJ)
      end if
      close(30)
      close(50)
      open(73,file='SWMS_3D.OUT\Q.out',       status='unknown')
      open(70,file='SWMS_3D.OUT\Run_Inf.out', status='unknown')

      call hOut  (hNew,x,y,z,NumNP,tInit,IJ)
      call thOut (ThOld,x,y,z,NumNP,tInit,IJ)
      call SubReg(NumEl,NumElD,NumNP,NMat,hNew,ThOld,ThOld,x,y,z,MatNum,
     !            LayNum,KX,tInit,dt,NLay,0,lWat,lChem,Conc,ChPar,
     !            wCumA,wCumT,cCumA,cCumT,wVolI,cVolI,WatIn,SolIn)
      if(NObs.gt.0) call ObsNod(tInit,NumNP,NObs,Node,hNew,ThOld,Conc)

      write(*,*)'beginning of numerical solution'
      call getdat(i,i,iday)
      call gettim(ihours,mins,isecs,i)
      Rtime1=iday*24.*60.*60.+ihours*60.*60.+mins*60.+isecs

* --- Beginning of time loop -------------------------------------------

11    continue

*     Calculate water flow
      if(lWat.or.TLevel.eq.1)
     !call WatFlow(NumNP,NumEl,NumElD,NTab,NTabD,MBand,MBandD,NMat,
     !             NSeep,NSeepD,NumSPD,NSP,NP,NumBP,ItCum,MaxIt,Iter,
     !             Kode,t,dt,dtMin,dtOpt,dtOld,tOld,hCritA,hCritS,
     !             TolTh,TolH,rLen,Width,rTop,vMeanR,hMeanR,AtmInf,
     !             SinkF,SeepF,qGWLF,FreeD,Par,hTab,ConTab,CapTab,
     !             TheTab,hNew,hOld,hTemp,thR,thSat,ThNew,ThOld,MatNum,
     !             Con,Cap,ConSat,Axz,Bxz,Dxz,hSat,A,B,Q,F,x,y,z,KX,
     !             Sink,DS,Beta,ConAxx,ConAyy,ConAzz,ConAxy,ConAxz,
     !             ConAyz,KXB,Explic,GWL0L,Aqh,Bqh,lWat,TLevel,lOrt,
     !             DrainF,ND,NDr,NDrD,NElDrD,KNoDr,rRoot,P0,POptm,P2H,
     !             P2L,P3,r2H,r2L,ConO,
     !             A1,B1,NumNPD,IAD,IADN,IADD,VRV,RES,RQI,RQ,QQ,QI,
     !             RQIDOT,ECNVRG,RCNVRG,ACNVRG,MNorth,MaxItO)
      if(.not.lWat.and.TLevel.eq.1) then
        if(lChem) then
          call ChInit(NumNP,NumEl,NumElD,NMat,x,y,z,KX,MatNum,NLevel,
     !                Con,hNew,Sink,cBound(12),Vx,Vy,Vz,ConAxx,ConAyy,
     !                ConAzz,ConAxy,ConAxz,ConAyz,Dispxx,Dispyy,Dispzz,
     !                Dispxy,Dispxz,Dispyz,ChPar,ThNew,thSat,Conc,Fc,Gc,
     !                ListNE,lUpW,WeTab,dt,dtMaxC,Peclet,Courant,lArtD,
     !                PeCr,PeCrMx,ConO)
        else
          call Veloc(NumNP,NumEl,NumElD,hNew,x,y,z,KX,ListNE,Con,
     !               ConAxx,ConAyy,ConAzz,ConAxy,ConAxz,ConAyz,Vx,Vy,Vz)
        end if
        Iter=1
      end if

*     Calculate solute transport
      if(lChem)
     !  call Solute(NumNP,NumEl,NumElD,MBand,MBandD,NMat,t,Kode,A,B,Q,
     !              hNew,hOld,F,x,y,z,KX,dt,DS,Sink,MatNum,Con,ConO,
     !              ConAxx,ConAyy,ConAzz,ConAxy,ConAxz,ConAyz,Vx,Vy,Vz,
     !              Dispxx,Dispyy,Dispzz,Dispxy,Dispxz,Dispyz,ChPar,
     !              ThNew,ThOld,thSat,Ac,Fc,Gc,Qc,Conc,ListNE,cBound,
     !              tPulse,NumBP,KodCB,KXB,NLevel,cPrec,crt,cht,lWat,
     !              lUpW,WeTab,epsi,CumCh0,CumCh1,CumChR,dtMaxC,Peclet,
     !              Courant,lArtD,PeCr,PeCrMx,lOrt,
     !              A1,B1,NumNPD,IAD,IADN,IADD,VRV,RES,RQI,RQ,QQ,QI,
     !              RQIDOT,ECNVRG,RCNVRG,ACNVRG,MNorth,MaxItO)

*     T-Level information
      call TLInf   (NumNP,NumBP,Kode,Q,hNew,CumQ,Width,SWidth,KXB,t,dt,
     !              TLevel,ShortF,TPrint(PLevel),Iter,ItCum,rTop,rRoot,
     !              vMeanR,hMeanT,hMeanR,hMeanG,AtmInF,SinkF,CumQrT,
     !              CumQrR,CumQvR,NumKD,hMean,vMean,lWat,lChem,rLen,
     !              Peclet,Courant,wCumT,wCumA,PeCrMx)
      if(lChem)
     !  call SolInf(1,NumNP,Kode,Qc,t,dt,TLevel,ShortF,TPrint(PLevel),
     !              AtmInF,NumKD,SMean,ChemS,CumCh0,CumCh1,CumChR,cCumA,
     !              cCumT,lWat,PeCrMx)
      if(NObs.gt.0) call ObsNod(t,NumNP,NObs,Node,hNew,ThNew,Conc)

*     P-Level information
      if(abs(TPrint(PLevel)-t).lt.0.001*dt) then
        if(lWat.or.(.not.lWat.and.PLevel.eq.1)) then
          call hOut (hNew,x,y,z,NumNP,t,IJ)
          call thOut(ThNew,x,y,z,NumNP,t,IJ)
          if(FluxF) then
            if(.not.lChem)
     !        call Veloc(NumNP,NumEl,NumElD,hNew,x,y,z,KX,ListNE,Con,
     !                   ConAxx,ConAyy,ConAzz,ConAxy,ConAxz,ConAyz,Vx,
     !                   Vy,Vz)
            call FlxOut(Vx,Vy,Vz,x,y,z,NumNP,t,IJ)
            call QOut  (Q,x,y,z,NumNP,t,IJ)
          end if
        end if
        call SubReg(NumEl,NumElD,NumNP,NMat,hNew,ThOld,ThNew,x,y,z,
     !              MatNum,LayNum,KX,t,dt,NLay,PLevel,lWat,lChem,Conc,
     !              ChPar,wCumA,wCumT,cCumA,cCumT,wVolI,cVolI,WatIn,
     !              SolIn)
        call BouOut(NumNP,NumBP,t,hNew,ThNew,Q,Width,KXB,Kode,x,y,z,
     !              Qc,Conc)
        if(lChem)
     !    call cOut(NumNP,Conc,x,y,z,t,IJ)
        PLevel=PLevel+1
      end if

*     A-level information
      if(abs(t-tAtm).le.0.001*dt.and.AtmInF) then
        if(lWat)
     !    call ALInf (t,CumQ,hMeanT,hMeanR,hMeanG,ALevel,CumQrT,CumQrR,
     !                CumQvR,NumKD)
        if(lChem)
     !    call SolInf(0,NumNP,Kode,Qc,t,dt,ALevel,ShortF,0.,AtmInF,
     !                NumKD,SMean,ChemS,CumCh0,CumCh1,CumChR,cCumA,
     !                cCumT,lWat,PeCrMx)
        if(ALevel.lt.MaxAL) then
          call SetAtm(tAtm,rTop,rRoot,hCritA,Width,KXB,NumBP,Kode,hNew,
     !                Q,NumNP,GWL0L,qGWLF,FreeD,cPrec,cht,crt)
          ALevel=ALevel+1
        end if
      end if

*     Root extraction
      if(SinkF)
     !  call SetSnk(NumNP,NMat,MatNum,hNew,rRoot,Sink,P0,POptm,P2H,P2L,
     !              P3,r2H,r2L,Beta,rLen)

*     Time governing
      if (abs(t-tMax).le.0.001*dt) then
        call getdat(i,i,iday)
        call gettim(ihours,mins,isecs,i)
        Rtime2=iday*24.*60.*60.+ihours*60.*60.+mins*60.+isecs
        write(70,*)
        write(70,*) 'Real time [sec]',Rtime2-RTime1
        write( *,*) 'Real time [sec]',Rtime2-RTime1
        stop
      end if
      tOld=t
      dtOld=dt
      call TmCont(dt,dtMax,dtOpt,dMul,dMul2,dtMin,Iter,TPrint(PLevel),
     !            tAtm,t,tMax,dtMaxC)
      TLevel=TLevel+1
      t=t+dt
      goto 11

* --- end of time loop -------------------------------------------------

      end

*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||