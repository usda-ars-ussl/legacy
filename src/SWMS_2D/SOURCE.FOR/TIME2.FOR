* Source file TIME2.FOR ||||||||||||||||||||||||||||||||||||||||||||||||

      subroutine TmCont(dt,dtMaxW,dtOpt,dMul,dMul2,dtMin,Iter,tPrint,
     !                  tAtm,t,tMax,dtMaxC,lMinStep,dtInit)
      logical lMinStep

      if(lMinStep) then
        dtMax=amin1(dtMaxW,dtMaxC,dtInit)
        lMinStep=.false.
      else
        dtMax=amin1(dtMaxW,dtMaxC)
      end if
      tFix=amin1(tPrint,tAtm,tMax)
      if(Iter.le.3.and.(tFix-t).ge.dMul*dtOpt) 
     !  dtOpt=amin1(dtMax,dMul*dtOpt)
      if(Iter.ge.7)
     !  dtOpt=amax1(dtMin,dMul2*dtOpt)
      dt=amin1(dtOpt,tFix-t)
      dt=amin1((tFix-t)/anint((tFix-t)/dt),dtMax)
      if(tFix-t.ne.dt.and.dt.gt.(tFix-t)/2.) dt=(tFix-t)/2.
      return
      end

************************************************************************

      subroutine SetAtm(tAtm,rTop,rRoot,hCritA,Width,KXB,NumBP,Kode,
     !                  hNew,Q,NumNP,GWL0L,qGWLF,FreeD,cPrec,cht,crt,
     !                  lMinStep)

      logical qGWLF,FreeD,lMinStep
      dimension Width(NumBP),KXB(NumBP),Kode(NumNP),hNew(NumNP),Q(NumNP)

      rTopOld=rTop
      read (31,*) tAtm,Prec,cPrec,rSoil,rRoot,hCritA,rGWL,GWL,crt,cht
      Prec=abs(Prec)
      rSoil=abs(rSoil)
      rRoot=abs(rRoot)
      hCritA=-abs(hCritA)
      hGWL=GWL+GWL0L
      rTop=rSoil-Prec
      do 11 i=1,NumBP
        n=KXB(i)
        K=Kode(n)
        if(K.eq.4.or.K.eq.-4) then
          Kode(n)=-4
          Q(n)=-Width(i)*rTop
          goto 11
        end if
        if(K.eq. 3) then
          if(abs(hNew(n)-hGWL).gt.1.e-8) lMinStep=.true.
          hNew(n)=hGWL
        end if
        if(K.eq.-3.and..not.qGWLF.and..not.FreeD) then
          if(Width(i).gt.0.) rGWLOld=-Q(n)/Width(i)
          if(abs(rGWLOld-rGWL).gt.1.e-8) lMinStep=.true.
          Q(n)=-Width(i)*rGWL
        end if
11    continue
      if((Prec-rSoil).gt.0.) then
        cPrec=Prec/(Prec-rSoil)*cPrec
      else
        cPrec=0.
      end if
      if(abs(rTop-rTopOld).gt.1.e-8.and.abs(rTop).gt.0.) lMinStep=.true.
      return
      end
      
************************************************************************

      real function Fqh(GWL,Aqh,Bqh)
      Fqh=-Aqh*exp(Bqh*abs(GWL))
      return
      end

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||