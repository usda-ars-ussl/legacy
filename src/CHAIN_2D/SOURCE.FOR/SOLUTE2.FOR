* Source file SOLUTE2.FOR ||||||||||||||||||||||||||||||||||||||||||||||

      subroutine Solute(NumNP,NumEl,NumElD,MBand,MBandD,NMat,NS,NSD,t,
     !                  Kode,A,B,Q,F,x,y,KX,KAT,dt,DS,Sink,MatNum,VxO,
     !                  VzO,VxN,VzN,Dispxx,Dispzz,Dispxz,ChPar,ThNew,
     !                  ThOld,thSat,Ac,Fc,Gc,Qc,Conc,cBound,tPulse,
     !                  NumBP,KodCB,KXB,NLevel,cPrec,crt,cht,lUpW,WeTab,
     !                  epsi,CumCh0,CumCh1,CumChR,dtMaxC,Peclet,Courant,
     !                  cPrevO,TempO,TempN,TDep,Width,cTolA,cTolR,
     !                  MaxItC,IterC,cTemp,cNew,Sorb,SorbN,Fc1,Gc1,
     !                  lLinear,lEquil,lArtD,PeCr,lOrt,
     !                  A1,B1,NumNPD,IAD,IADN,IADD,VRV,RES,RQI,RQ,QQ,QI,
     !                  RQIDOT,ECNVRG,RCNVRG,ACNVRG,MNorth,MaxItO)

      double precision A,B,A1,B1,VRV,RES,RQI,RQ,QQ,QI,RQIDOT,ECNVRG,
     !                 RCNVRG,ACNVRG
      logical lUpW,lLinear,lConv,lEquil,lOrt,lArtD
      dimension A(MBandD,NumNP),B(NumNP),Q(NumNP),F(NumNP),KX(NumElD,4),
     !          MatNum(NumNP),Sink(NumNP),DS(NumNP),x(NumNP),y(NumNP),
     !          Kode(NumNP),VxO(NumNP),VzO(NumNP),VxN(NumNP),VzN(NumNP),
     !          ThNew(NumNP),ThOld(NumNP),Conc(NSD,NumNP),Qc(NSD,NumNP),
     !          KXB(NumBP),Ac(NumNP),Fc(NumNP),Gc(NumNP),KodCB(NumBP),
     !          Dispxx(NumNP),Dispzz(NumNP),Dispxz(NumNP),thSat(NMat),
     !          ChPar(NSD*16+4,NMat),CumCh0(NS),CumCh1(NS),CumChR(NS),
     !          cPrevO(NumNP),Width(NumBP),cBound(NSD,9),TDep(NSD*16+4),
     !          TempO(NumNP),TempN(NumNP),WeTab(3,2*NumEl),cPrec(NS),
     !          cht(NS),crt(NS),cNew(NumNP),cTemp(NumNP),SorbN(NumNP),
     !          Sorb(NSD,NumNP),Fc1(NumNP),Gc1(NumNP),
     !          A1(MBandD,NumNP),B1(NumNP),RES(NumNP),IAD(MBandD,NumNP),
     !          IADN(NumNP),IADD(NumNP),VRV(NumNP),RQI(NumNPD,MNorth),
     !          RQ(NumNP),QQ(NumNP),RQIDOT(MNorth),QI(NumNPD,MNorth)

*     Initialisation

      IterC=1.
      Peclet=0.
      Courant=0.
      dtMaxC=1.e+30

*     Loop on species in the chain

      do 21 jj=1,NS
        Iter=0
        if(t.gt.tPulse) then
          do 12 i=1,4
            cBound(jj,i)=0.
12        continue
        end if

        if(.not.lLinear) then
          do 13 i=1,NumNP
            cNew(i)=Conc(jj,i)
            if(.not.lEquil) SorbN(i)=Sorb(jj,i)
13        continue
        end if

*       Iterative loop for a nonlinear adsorption isotherm

14      Iter=Iter+1
        Cum0=0.
        Cum1=0.
        CumR=0.
        do 16 i=1,NumNP
          if(.not.lLinear) then
            cTemp(i)=cNew(i)
            if(lOrt) B1(i)=cNew(i)
          else
            if(lOrt) B1(i)=Conc(jj,i)
          end if
          B(i)=0.d0
          Qc(jj,i)=0.
          if(epsi.lt.0.001) then
            if(lOrt) jjj=IADD(i)
            A(jjj,i)=0.d0
          else
            do 15 j=1,MBandD
              A(j,i)=0.d0
15          continue
          end if
16      continue

        do 17 Level=1,NLevel
          if(Level.eq.NLevel) then
            call Disper(jj,NSD,NumNP,NMat,Dispxx,Dispzz,Dispxz,VxN,VzN,
     !                  ThNew,thSat,ChPar,MatNum,TempN,TDep,lArtD,PeCr,
     !                  dt,lLinear,Conc)
            call PeCour(jj,NSD,NumNP,NumEl,NumElD,NMat,x,y,VxN,VzN,KX,
     !                  MatNum,Dispxx,Dispzz,ChPar,ThNew,dt,dtMxC,Pecl,
     !                  Cour,thSat,TempN,TDep,Conc,lLinear,lUpW,lArtD,
     !                  PeCr)
            Peclet=amax1(Peclet,Pecl)
            Courant=amax1(Courant,Cour)
            dtMaxC=amin1(dtMaxC,dtMxC)
          else
            call Disper(jj,NSD,NumNP,NMat,Dispxx,Dispzz,Dispxz,VxO,VzO,
     !                  ThOld,thSat,ChPar,MatNum,TempO,TDep,lArtD,PeCr,
     !                  dt,lLinear,Conc)
            if(lUpW.and.Iter.eq.1)
     !        call WeFact(NumNP,NumEl,NumElD,x,y,KX,WeTab,VxO,VzO,
     !                    Dispxx,Dispzz,Dispxz)
          end if

*         Calculate the coefficients
          call Coeff(jj,Level,NLevel,NumNP,NMat,NSD,lLinear,lEquil,epsi,
     !               dt,MatNum,ChPar,TDep,Conc,cNew,cPrevO,TempN,TempO,
     !               ThNew,ThOld,ThSat,Ac,Fc,Gc,Dispxx,Dispzz,Dispxz,
     !               VxN,VzN,VxO,VzO,SorbN,Sorb,cBound,Sink,Fc1,Gc1,
     !               lUpW,lArtD)

*          Calculate the coefficient matrix
           call SolMat(jj,NumNP,NumEl,NumElD,NLevel,Level,MBand,MBandD,
     !                 NSD,NumBP,KAT,epsi,dt,lUpW,lOrt,Cum0,Cum1,CumR,
     !                 A,B,F,DS,KX,x,y,Sink,Kode,VxO,VzO,VxN,VzN,Conc,
     !                 cBound,Gc,Gc1,Fc,Fc1,Qc,Ac,Dispxx,Dispzz,Dispxz,
     !                 WeTab,KodCB,KXB,IAD,IADN,IADD)
17      continue

*       Boundary condition
        call c_Bound(jj,NumNP,NMat,MBand,MBandD,NumBP,NS,NSD,A,B,Q,Qc,
     !               Conc,ChPar,Kode,KXB,KodCB,MatNum,Width,cBound,
     !               cPrec,crt,cht,TempN,TDep,epsi,dt,DS,lOrt,IADD)

*       Solve the global matrix equation for transport
        if(epsi.lt.0.001) then
          do 18 i=1,NumNP
            if(lOrt) jjj=IADD(i)
            B(i)=B(i)/A(jjj,i)
18        continue
        else
          if(lOrt) then
            call ILU     (A,NumNP,MBandD,IAD,IADN,IADD,A1)
            North=4
            call OrthoMin(A,B1,B,NumNP,MBandD,NumNPD,IAD,IADN,IADD,A1,
     !                    VRV,RES,RQI,RQ,QQ,QI,RQIDOT,ECNVRG,RCNVRG,
     !                    ACNVRG,North,MNorth,MaxItO)
          else
            call SolveT(A,B,MBand,MBandD,NumNP)
          end if
        end if

*       Test for convergence for nonlinear problem
        lConv=.true.
        do 19 i=1,NumNP
          if(NS.gt.1.and.Iter.eq.1) cPrevO(i)=Conc(jj,i)
          if(lOrt) B(i)=B1(i)
          if(lLinear) then
            Conc(jj,i)=sngl(B(i))
          else
            cNew(i)=sngl(B(i))
            if(abs(cNew(i)-cTemp(i)).gt.cTolA+cTolR*Conc(jj,i))
     !        lConv=.false.
          end if
19      continue
        if(.not.lLinear) then
          if(.not.lConv.and.iter.lt.MaxItC) goto 14
          do 20 i=1,NumNP
            Conc(jj,i)=cNew(i)
            if(.not.lEquil) Sorb(jj,i)=SorbN(i)
20        continue
        else if(.not.lEquil) then
*         Calculate sorbed concentration for linear noneq. adsorption
          Tr=293.15
          R=8.314
          jjj=(jj-1)*16
          do 22 i=1,NumNP
            M=MatNum(i)
            TT=(TempN(i)+273.15-Tr)/R/(TempN(i)+273.15)/Tr
            Frac=ChPar(4,M)*exp(TDep(4)*TT)
            xKs=ChPar(jjj+7,M)*exp(TDep(jjj+7)*TT)
            GamS=ChPar(jjj+12,M)*exp(TDep(jjj+12)*TT)
            GamS1=ChPar(jjj+15,M)*exp(TDep(jjj+15)*TT)
            Omega=ChPar(jjj+20,M)*exp(TDep(jjj+20)*TT)
            Sorb(jj,i)=Sorb(jj,i)+dt*Omega*(1.-Frac)*xKs*Conc(jj,i)/
     !                             (2.+dt*(Omega+GamS+GamS1))
22        continue
        end if
        CumCh0(jj)=CumCh0(jj)+Cum0
        CumCh1(jj)=CumCh1(jj)+Cum1
        CumChR(jj)=CumChR(jj)+CumR
        IterC=max0(IterC,Iter)
21    continue
      return
      end

************************************************************************

      subroutine Coeff(jj,Level,NLevel,NumNP,NMat,NSD,lLinear,lEquil,
     !                 epsi,dt,MatNum,ChPar,TDep,Conc,cNew,cPrevO,TempN,
     !                 TempO,ThNew,ThOld,ThSat,Ac,Fc,Gc,Dispxx,Dispzz,
     !                 Dispxz,VxN,VzN,VxO,VzO,SorbN,Sorb,cBound,Sink,
     !                 Fc1,Gc1,lUpW,lArtD)

      logical lLinear,lEquil,lUpW,lArtD
      dimension MatNum(NumNP),Sink(NumNP),VxO(NumNP),VzO(NumNP),
     !          VxN(NumNP),VzN(NumNP),ThNew(NumNP),ThOld(NumNP),
     !          Conc(NSD,NumNP),Ac(NumNP),Fc(NumNP),Gc(NumNP),
     !          Dispxx(NumNP),Dispzz(NumNP),Dispxz(NumNP),thSat(NMat),
     !          ChPar(NSD*16+4,NMat),cPrevO(NumNP),cBound(NSD,9),
     !          TDep(NSD*16+4),TempO(NumNP),TempN(NumNP),cNew(NumNP),
     !          SorbN(NumNP),Sorb(NSD,NumNP),Fc1(NumNP),Gc1(NumNP)

      alf=1.-epsi
      Tr=293.15
      R=8.314
      jjj=(jj-1)*16
      if(jj.gt.1) jj1=jjj-16

      do 11 i=1,NumNP
        M=MatNum(i)

*       Temperature dependent transport parameters
        if(Level.eq.NLevel) then
          TT=(TempN(i)+273.15-Tr)/R/(TempN(i)+273.15)/Tr
          if(jj.gt.1) cPrev=Conc(jj-1,i)
        else
          TT=(TempO(i)+273.15-Tr)/R/(TempO(i)+273.15)/Tr
          if(jj.gt.1) cPrev=cPrevO(i)
        end if
        ro=ChPar(1,M)*exp(TDep(1)*TT)
        Frac=ChPar(4,M)*exp(TDep(4)*TT)
        xKs=ChPar(jjj+7,M)*exp(TDep(jjj+7)*TT)
        xNu=ChPar(jjj+8,M)*exp(TDep(jjj+8)*TT)
        fExp=ChPar(jjj+9,M)*exp(TDep(jjj+9)*TT)
        Henry=ChPar(jjj+10,M)*exp(TDep(jjj+10)*TT)
        GamL=ChPar(jjj+11,M)*exp(TDep(jjj+11)*TT)
        GamS=ChPar(jjj+12,M)*exp(TDep(jjj+12)*TT)
        GamG=ChPar(jjj+13,M)*exp(TDep(jjj+13)*TT)
        GamL1=ChPar(jjj+14,M)*exp(TDep(jjj+14)*TT)
        GamS1=ChPar(jjj+15,M)*exp(TDep(jjj+15)*TT)
        GamG1=ChPar(jjj+16,M)*exp(TDep(jjj+16)*TT)
        xMuL=ChPar(jjj+17,M)*exp(TDep(jjj+17)*TT)
        xMuS=ChPar(jjj+18,M)*exp(TDep(jjj+18)*TT)
        xMuG=ChPar(jjj+19,M)*exp(TDep(jjj+19)*TT)
        Omega=ChPar(jjj+20,M)*exp(TDep(jjj+20)*TT)
        if(jj.gt.1) then
          xKsP=ChPar(jj1+7,M)*exp(TDep(jj1+7)*TT)
          xNuP=ChPar(jj1+8,M)*exp(TDep(jj1+8)*TT)
          fExpP=ChPar(jj1+9,M)*exp(TDep(jj1+9)*TT)
          HenryP=ChPar(jj1+10,M)*exp(TDep(jj1+10)*TT)
          GamL1P=ChPar(jj1+14,M)*exp(TDep(jj1+14)*TT)
          GamS1P=ChPar(jj1+15,M)*exp(TDep(jj1+15)*TT)
          GamG1P=ChPar(jj1+16,M)*exp(TDep(jj1+16)*TT)
        end if
        if(Level.eq.NLevel) then
          TTO=(TempO(i)+273.15-Tr)/R/(TempO(i)+273.15)/Tr
          xKsO=ChPar(jjj+7,M)*exp(TDep(jjj+7)*TTO)
          xNuO=ChPar(jjj+8,M)*exp(TDep(jjj+8)*TTO)
          fExpO=ChPar(jjj+9,M)*exp(TDep(jjj+9)*TTO)
          HenryO=ChPar(jjj+10,M)*exp(TDep(jjj+10)*TTO)
          GamSO=ChPar(jjj+12,M)*exp(TDep(jjj+12)*TTO)
          GamS1O=ChPar(jjj+15,M)*exp(TDep(jjj+15)*TTO)
          xMuSO=ChPar(jjj+18,M)*exp(TDep(jjj+18)*TTO)
          OmegaO=ChPar(jjj+20,M)*exp(TDep(jjj+20)*TTO)
          dKs=(xKs-xKsO)/dt
          dNu=(xNu-xNuO)/dt
          ddExp=(fExp-fExpO)/dt
          dHenry=(Henry-HenryO)/dt
        else
          TTN=(TempN(i)+273.15-Tr)/R/(TempN(i)+273.15)/Tr
          xKsN=ChPar(jjj+7,M)*exp(TDep(jjj+7)*TTN)
          xNuN=ChPar(jjj+8,M)*exp(TDep(jjj+8)*TTN)
          fExpN=ChPar(jjj+9,M)*exp(TDep(jjj+9)*TTN)
          HenryN=ChPar(jjj+10,M)*exp(TDep(jjj+10)*TTN)
          dKs=(xKsN-xKs)/dt
          dNu=(xNuN-xNu)/dt
          ddExp=(fExpN-fExp)/dt
          dHenry=(HenryN-Henry)/dt
        end if

        dSConc=1.
        dConc=1.
        SConcP=1.
        SConc=1.
        if(.not.lEquil) then
          SEqui=0.
          SEquiO=0.
        end if

*       Effects of nonlinear adsorption
        if(.not.lLinear) then
          cc=Conc(jj,i)
          if(Level.eq.NLevel) then
            cc=cNew(i)
            cMid=(Conc(jj,i)+cNew(i))/2.
          end if
          if(cc.gt.0.) then
            dSConc=fExp*cc**(fExp-1.)/(1+xNu*cc**fExp)**2
            SConc=cc**(fExp-1.)/(1.+xNu*cc**fExp)
            if(.not.lEquil) SEqui=xKs*cc**(fExp-1.)/(1.+xNu*cc**fExp)
          end if
          if(Level.eq.NLevel) then
            if(cMid.gt.0.) dConc=fExp*cMid**(fExp-1.)/
     !                           (1.+xNu*cMid**fExp)**2
            if(.not.lEquil.and.Conc(jj,i).gt.0.)
     !        SEquiO=xKsO*Conc(jj,i)**fExpO/(1.+xNuO*Conc(jj,i)**fExpO)
          end if
          if(jj.gt.1.and.cPrev.gt.0.) SConcP=cPrev**(fExpP-1.)/
     !                                       (1.+xNuP*cPrev**fExpP)
        end if

*       Correction of the dispersion coefficient
        if(Level.ne.NLevel) then
          ThW=ThOld(i)
          ThG=(thSat(M)-ThW)
          if(.not.lArtD.and..not.lUpW) then
            DPom=dt/6./(ThW+ro*Frac*xKs*dSConc+ThG*Henry)
            Dispxx(i)=Dispxx(i)+VxO(i)*VxO(i)*DPom
            Dispzz(i)=Dispzz(i)+VzO(i)*VzO(i)*DPom
            Dispxz(i)=Dispxz(i)+VxO(i)*VzO(i)*DPom
          end if
        else
          ThW=ThNew(i)
          ThG=(thSat(M)-ThW)
          Ac(i)=-(alf*ThOld(i)+epsi*ThW)-
     !           (alf*xKsO+epsi*xKs)*ro*Frac*dConc-
     !           (alf*(thSat(M)-ThOld(i))*HenryO+epsi*ThG*Henry)
          if(.not.lArtD.and..not.lUpW) then
            DPom=dt/6./(ThW+ro*Frac*xKs*dSConc+ThG*Henry)
            Dispxx(i)=Dispxx(i)-VxN(i)*VxN(i)*DPom
            Dispzz(i)=Dispzz(i)-VzN(i)*VzN(i)*DPom
            Dispxz(i)=Dispxz(i)-VxN(i)*VzN(i)*DPom
          end if
        end if

*       Calculate sorbed concentration on type-2 sorption sites
        if(.not.lEquil) then
          SSorb=Sorb(jj,i)
          if(Level.eq.NLevel) then
            if(.not.lLinear) then
              SorbN(i)=Sorb(jj,i)+dt*
     !             (epsi*    (Omega*((1.-Frac)*SEqui*cc-SorbN(i))-
     !                       (GamS+GamS1)*SorbN(i)+(1.-Frac)*xMuS)+
     !             (1.-epsi)*(OmegaO*((1.-Frac)*SEquiO-SSorb)-
     !                       (GamSO+GamS1O)*SSorb+(1.-Frac)*xMuSO))
              SSorb=SorbN(i)
            else
              Sorb(jj,i)=((2.-(OmegaO+GamSO+GamS1O)*dt)*Sorb(jj,i)+
     !                    dt*(1.-Frac)*OmegaO*xKsO*Conc(jj,i)+
     !                    dt*(1.-Frac)*(xMuSO+xMuS))/
     !                    (2.+dt*(Omega+GamS+GamS1))
              SSorb=Sorb(jj,i)
            end if
          end if
        end if

*       Calculate zero-order coefficient Gc
        cS=amin1(cBound(jj,5),Conc(jj,i))
        Gc(i)=xMuL*ThW+Frac*ro*xMuS+ThG*xMuG-Sink(i)*cS
        Gc1(i)=xMuL*ThW+ro*xMuS+ThG*xMuG
        if(.not.lEquil) Gc(i)=Gc(i)+Omega*ro*SSorb
        if(jj.gt.1) then
          cG=cPrev*(GamL1P*ThW+ro*Frac*xKsP*GamS1P*SConcP+
     !              ThG*HenryP*GamG1P)
          if(.not.lEquil) cG=cG+GamS1P*ro*Sorb(jj-1,i)
          Gc(i)=Gc(i)+cG
          Gc1(i)=Gc1(i)+cG
        end if

*       Calculate first-order coefficient Fc
        Fc(i)=-(GamL+GamL1)*ThW-(GamS+GamS1)*ro*Frac*xKs*SConc-
     !         (GamG+GamG1)*ThG*Henry+
     !         Sink(i)+Henry*(ThNew(i)-ThOld(i))/dt-ThG*dHenry
        if(.not.lLinear) then
          if(cMid.gt.0.) Fc(i)=Fc(i)-ro*Frac*
     !          (cMid**(fExp-1.)/(1.+xNu*cMid**fExp)*dKs-
     !           xKs*cMid**(2.*fExp-1.)/(1.+xNu*cMid**fExp)**2*dNu+
     !           xKs*alog(cMid)*cMid**(fExp-1.)/
     !           (1.+xNu*cMid**fExp)**2*ddExp)
        else
          Fc(i)=Fc(i)-ro*Frac*dKs
        end if
        if(.not.lEquil) then
          if(lLinear) then
            Fc(i)=Fc(i)-Omega*ro*(1.-Frac)*xKs
            if(Level.eq.NLevel) Fc(i)=Fc(i)+Omega*ro*
     !               (dt*Omega*(1.-Frac)*xKs/(2.+dt*(Omega+GamS+GamS1)))
          else
            Fc(i)=Fc(i)-Omega*ro*(1.-Frac)*SEqui
          end if
        end if
        Fc1(i)=(-(GamL+GamL1)*ThW-(GamS+GamS1)*ro*Frac*xKs*SConc-
     !           (GamG+GamG1)*ThG*Henry)*Conc(jj,i)
        if(.not.lEquil) Fc1(i)=Fc1(i)-(GamS+GamS1)*ro*Sorb(jj,i)
11    continue
      return
      end

************************************************************************

      subroutine SolMat(jj,NumNP,NumEl,NumElD,NLevel,Level,MBand,MBandD,
     !                  NSD,NumBP,KAT,epsi,dt,lUpW,lOrt,Cum0,Cum1,CumR,
     !                  A,B,F,DS,KX,x,y,Sink,Kode,VxO,VzO,VxN,VzN,Conc,
     !                  cBound,Gc,Gc1,Fc,Fc1,Qc,Ac,Dispxx,Dispzz,Dispxz,
     !                  WeTab,KodCB,KXB,IAD,IADN,IADD)

      logical lUpW,lOrt
      double precision A,B
      dimension A(MBandD,NumNP),B(NumNP),F(NumNP),KX(NumElD,4),y(NumNP),
     !          x(NumNP),Sink(NumNP),DS(NumNP),Kode(NumNP),VxO(NumNP),
     !          VzO(NumNP),VxN(NumNP),VzN(NumNP),Ac(NumNP),Fc(NumNP),
     !          Qc(NSD,NumNP),Conc(NSD,NumNP),Fc1(NumNP),Gc1(NumNP),
     !          Gc(NumNP),Dispxx(NumNP),Dispzz(NumNP),Dispxz(NumNP),
     !          cBound(NSD,9),WeTab(3,2*NumEl),KXB(NumBP),KodCB(NumBP),
     !          IAD(MBandD,NumNP),IADN(NumNP),IADD(NumNP),
     !          Bi(3),Ci(3),Wx(3),Wz(3),Vx(3),Vz(3),S(3,3),List(3)

      xMul=1.
      alf=1.-epsi
      Esp=epsi
      if(Level.ne.NLevel) Eps=alf
      jjj=MBand

      do 11 i=1,NumNP
        F(i)=0.
        if(Level.eq.NLevel) DS(i)=0.
11    continue

*     Loop on elements
      NumSEl=0
      do 15 n=1,NumEl
        NUS=4
        if(KX(n,3).eq.KX(n,4)) NUS=3

*       Loop on subelements
        do 14 k=1,NUS-2
          NumSEl=NumSEl+1
          i=KX(n,1)
          j=KX(n,k+1)
          l=KX(n,k+2)
          List(1)=i
          List(2)=j
          List(3)=l
          if(Level.ne.NLevel) then
            Vx(1)=VxO(i)
            Vx(2)=VxO(j)
            Vx(3)=VxO(l)
            Vz(1)=VzO(i)
            Vz(2)=VzO(j)
            Vz(3)=VzO(l)
          else
            Vx(1)=VxN(i)
            Vx(2)=VxN(j)
            Vx(3)=VxN(l)
            Vz(1)=VzN(i)
            Vz(2)=VzN(j)
            Vz(3)=VzN(l)
          end if
          Ci(1)=x(l)-x(j)
          Ci(2)=x(i)-x(l)
          Ci(3)=x(j)-x(i)
          Bi(1)=y(j)-y(l)
          Bi(2)=y(l)-y(i)
          Bi(3)=y(i)-y(j)
          AE=(Ci(3)*Bi(2)-Ci(2)*Bi(3))/2.
          if(KAT.eq.1) xMul=2.*3.1416*(x(i)+x(j)+x(l))/3.
          cS=cBound(jj,5)
          cS=(amin1(cS,Conc(jj,i))+amin1(cS,Conc(jj,j))+
     !        amin1(cS,Conc(jj,l)))/3.
          if(Level.eq.NLevel) Eps=epsi
          if(Level.ne.NLevel) Eps=alf
          RootCh=xMul*AE*dt*Eps*cS*(Sink(i)+Sink(j)+Sink(l))/3.
          Cum0=Cum0-xMul*AE*dt*Eps*(Gc1(i)+Gc1(j)+Gc1(l))/3.
          Cum1=Cum1-xMul*AE*dt*Eps*(Fc1(i)+Fc1(j)+Fc1(l))/3.
          CumR=CumR+RootCh
          Ec1=(Dispxx(i)+Dispxx(j)+Dispxx(l))/3.
          Ec2=(Dispxz(i)+Dispxz(j)+Dispxz(l))/3.
          Ec3=(Dispzz(i)+Dispzz(j)+Dispzz(l))/3.
          if(Level.eq.NLevel) AcE=(Ac(i)+Ac(j)+Ac(l))/3.
          FcE=(Fc(i)+Fc(j)+Fc(l))/3.
          GcE=(Gc(i)+Gc(j)+Gc(l))/3.
          Vxe=(Vx(1)+Vx(2)+Vx(3))/3.
          VzE=(Vz(1)+Vz(2)+Vz(3))/3.
          FMul=xMul*AE/4.
          SMul1=-1./AE/4.*xMul
          SMul2=AE/20.*xMul
          if(lUpW) then
            W1=WeTab(1,NumSEl)
            W2=WeTab(2,NumSEl)
            W3=WeTab(3,NumSEl)
            Wx(1)=2.*Vx(1)*(W2-W3)+Vx(2)*(W2-2.*W3)+Vx(3)*(2.*W2-W3)
            Wx(2)=Vx(1)*(2.*W3-W1)+2.*Vx(2)*(W3-W1)+Vx(3)*(W3-2.*W1)
            Wx(3)=Vx(1)*(W1-2.*W2)+Vx(2)*(2.*W1-W2)+2.*Vx(3)*(W1-W2)
            Wz(1)=2.*Vz(1)*(W2-W3)+Vz(2)*(W2-2.*W3)+Vz(3)*(2.*W2-W3)
            Wz(2)=Vz(1)*(2.*W3-W1)+2.*Vz(2)*(W3-W1)+Vz(3)*(W3-2.*W1)
            Wz(3)=Vz(1)*(W1-2.*W2)+Vz(2)*(2.*W1-W2)+2.*Vz(3)*(W1-W2)
          end if
          do 13 j1=1,3
            i1=List(j1)
            F(i1)=F(i1)+FMul*(GcE+Gc(i1)/3.)
            if(Level.eq.NLevel) DS(i1)=DS(i1)+FMul*(AcE+Ac(i1)/3.)
            iBound=0
            if(Kode(i).ne.0) then
              do 17 id=1,NumBP
                if(KXB(id).eq.i1.and.KodCB(id).gt.0) iBound=1
17            continue
            end if
            if(iBound.eq.1) Qc(jj,i1)=Qc(jj,i1)-Eps*FMul*(GcE+Gc(i1)/3.)
            do 12 j2=1,3
              i2=List(j2)
              S(j1,j2)=SMul1*(Ec1*Bi(j1)*Bi(j2)+Ec3*Ci(j1)*Ci(j2)+
     !                        Ec2*(Bi(j1)*Ci(j2)+Ci(j1)*Bi(j2)))
              S(j1,j2)=S(j1,j2)-(Bi(j2)/8.*(VxE+Vx(j1)/3.)+
     !                           Ci(j2)/8.*(VzE+Vz(j1)/3.))*xMul
              if(lUpW) S(j1,j2)=S(j1,j2)-xMul*
     !                         (Bi(j2)/40.*Wx(j1)+Ci(j2)/40.*Wz(j1))
              ic=1
              if(i1.eq.i2) ic=2
              S(j1,j2)=S(j1,j2)+SMul2*ic*(FcE+(Fc(i1)+Fc(i2))/3.)
              if(iBound.eq.1)
     !          Qc(jj,i1)=Qc(jj,i1)-Eps*S(j1,j2)*Conc(jj,i2)
              if(Level.ne.NLevel) then
                B(i1)=B(i1)-alf*S(j1,j2)*Conc(jj,i2)
              else
                if(lOrt) then
                  call Find(i1,i2,kk,NumNP,MBandD,IAD,IADN)
                  iB=kk
                else
                  iB=MBand+i2-i1
                end if
                A(iB,i1)=A(iB,i1)+epsi*S(j1,j2)
              end if
12          continue
13        continue
14      continue
15    continue
      do 16 i=1,NumNP
        if(Level.ne.NLevel) then
          B(i)=B(i)-alf*F(i)
        else
          if(lOrt) jjj=IADD(i)
          A(jjj,i)=A(jjj,i)+DS(i)/dt
          B(i)=B(i)+DS(i)/dt*Conc(jj,i)-epsi*F(i)
        end if
16    continue
      return
      end

************************************************************************

      subroutine c_Bound(jj,NumNP,NMat,MBand,MBandD,NumBP,NS,NSD,A,B,Q,
     !                   Qc,Conc,ChPar,Kode,KXB,KodCB,MatNum,Width,
     !                   cBound,cPrec,crt,cht,Temp,TDep,epsi,dt,DS,lOrt,
     !                   IADD)

      double precision A,B
      integer cKod
      logical lOrt
      dimension A(MBandD,NumNP),B(NumNP),MatNum(NumNP),Conc(NSD,NumNP),
     !          Qc(NSD,NumNP),Kode(NumNP),KXB(NumBP),KodCB(NumBP),
     !          Width(NumBP),Q(NumNP),Temp(NumNP),ChPar(NSD*16+4,NMat),
     !          TDep(NSD*16+4),cBound(NSD,9),cPrec(NS),cht(NS),crt(NS),
     !          DS(NumNP),IADD(NumNP)

      alf=1.-epsi
      iii=MBand
      do 14 i=1,NumNP
        if(Kode(i).ne.0) then
          M=MatNum(i)
          In=0
          do 11 j=1,NumBP
            if(KXB(j).eq.i) then
              if(KodCB(j).gt.0) then
                cKod=1
                if(abs(Kode(i)).le.2.or.abs(Kode(i)).ge.5)
     !                                    cBnd=cBound(jj,KodCB(j))
                if(abs(Kode(i)).eq.3) cBnd=cht(jj)
                if(abs(Kode(i)).eq.4) cBnd=cPrec(jj)
              else
                if(Q(i).gt.0.) then
                  cKod=3
                  if(abs(Kode(i)).eq.1.or.abs(Kode(i)).ge.5)
     !                                    cBnd=cBound(jj,-KodCB(j))
                  if(abs(Kode(i)).eq.3) cBnd=crt(jj)
                  if(abs(Kode(i)).eq.4) cBnd=cPrec(jj)
                else
                  cKod=2
                  if(abs(Kode(i)).eq.4) then
                    cKod=3
                    cBnd=0.
                  end if
                end if
              end if
              if(abs(Kode(i)).eq.2) cKod=2
              goto 12
            end if
11        continue

*         Point source or sink
          In=1
          if(Q(i).lt.0.) then
            cKod=2
          else
            cBnd=cBound(jj,6)
            cKod=3
          end if

12        continue

*         Dirichlet boundary condition
          if(cKod.eq.1) then
            Qc(jj,i)=Qc(jj,i)+Q(i)*(epsi*cBnd+alf*Conc(jj,i))-
     !                              DS(i)*(cBnd-Conc(jj,i))/dt
            if(lOrt) then
              A(IADD(i),i)=1.d100
              B(i)=1.d100*cBnd
            else
              do 13 k=1,2*MBand-1
                A(k,i)=0.
13            continue
              A(MBand,i)=1.
              B(i)=cBnd
            end if
          end if

*         Neumann boundary condition
          if(cKod.eq.2) then
            Qc(jj,i)=Q(i)*Conc(jj,i)
          end if

*         Cauchy boundary condition
          if(cKod.eq.3) then
            alf=1.-epsi
            B(i)=B(i)-Q(i)*(cBnd-alf*Conc(jj,i))
            if(lOrt) iii=IADD(i)
            A(iii,i)=A(iii,i)-epsi*Q(i)
            Qc(jj,i)=Q(i)*cBnd
          end if

*         Cauchy boundary condition for volatile solutes
          if(cKod.ge.2.and.KodCB(j).eq.-7.and.In.eq.0) then
            jjj=(jj-1)*16
            Tr=293.15
            R=8.314
            TT=(Temp(i)+273.15-Tr)/R/(Temp(i)+273.15)/Tr
            DispG=ChPar(jjj+6,M)*exp(TDep(jjj+6)*TT)
            Henry=ChPar(jjj+10,M)*exp(TDep(jjj+10)*TT)
            aa=DispG*Width(j)/cBound(jj,9)
            M=MatNum(i)
            B(i)=B(i)+alf*aa*Henry*Conc(jj,i)-
     !                aa*cBound(jj,8)
            if(lOrt) iii=IADD(i)
            A(iii,i)=A(iii,i)-epsi*aa*Henry
            Qc(jj,i)=Qc(jj,i)-aa*(Henry*Conc(jj,i)-CBound(jj,8))
          end if

        end if
14    continue
      return
      end

************************************************************************

*     Calculate the dispersion coefficient

      subroutine Disper(jj,NSD,NumNP,NMat,Dispxx,Dispzz,Dispxz,Vx,Vz,
     !                  theta,thSat,ChPar,MatNum,Temp,TDep,lArtD,PeCr,
     !                  dt,lLinear,Conc)

      logical lLinear,lArtD
      dimension Vx(NumNP),Vz(NumNP),theta(NumNP),ChPar(NSD*16+4,NMat),
     !          Dispxx(NumNP),Dispzz(NumNP),Dispxz(NumNP),MatNum(NumNP),
     !          thSat(NMat),Temp(NumNP),TDep(NSD*16+4),Conc(NSD,NumNP)

      jjj=(jj-1)*16
      Tr=293.15
      R=8.314
      do 11 i=1,NumNP
        M=MatNum(i)
        TT=(Temp(i)+273.15-Tr)/R/(Temp(i)+273.15)/Tr
        DispW=ChPar(jjj+5,M)*exp(TDep(jjj+5)*TT)
        DispG=ChPar(jjj+6,M)*exp(TDep(jjj+6)*TT)
        Henry=ChPar(jjj+10,M)*exp(TDep(jjj+10)*TT)
        TauW=theta(i)**(7./3.)/thSat(M)**2
        TauG=(thSat(M)-theta(i))**(7./3.)/thSat(M)**2
        Dif=theta(i)*DispW*TauW+(thSat(M)-theta(i))*DispG*Henry*TauG
        Dispxx(i)=Dif
        Dispzz(i)=Dif
        Dispxz(i)=0.
        Vabs=sqrt(Vx(i)*Vx(i)+Vz(i)*Vz(i))
        DispL=ChPar(2,M)
        DispT=ChPar(3,M)
        if(lArtD.and.Vabs.gt.0.) then
          ro=ChPar(1,M)*exp(TDep(1)*TT)
          Frac=ChPar(4,M)*exp(TDep(4)*TT)
          xKs=ChPar(jjj+7,M)*exp(TDep(jjj+7)*TT)
          xNu=ChPar(jjj+8,M)*exp(TDep(jjj+8)*TT)
          fExp=ChPar(jjj+9,M)*exp(TDep(jjj+9)*TT)
          aa=1.
          if(.not.lLinear.and.Conc(jj,i).gt.0.)
     !      aa=fExp*Conc(jj,i)**(fExp-1.)/(1.+xNu*Conc(jj,i)**fExp)**2
          ThRet=theta(i)+ro*xKs*aa*Frac+(thSat(M)-theta(i))*Henry
          DispL=amax1(DispL,Vabs*dt/ThRet/PeCr-Dif/Vabs)
        end if
        if(Vabs.gt.0.) then
          Dispxx(i)=DispL*Vx(i)*Vx(i)/Vabs+DispT*Vz(i)*Vz(i)/Vabs+Dif
          Dispzz(i)=DispL*Vz(i)*Vz(i)/Vabs+DispT*Vx(i)*Vx(i)/Vabs+Dif
          Dispxz(i)=(DispL-DispT)*Vx(i)*Vz(i)/Vabs
        end if
11    continue
      return
      end


************************************************************************

*     Calculate upstream weighing factors

      subroutine WeFact(NumNP,NumEl,NumElD,x,y,KX,WeTab,Vx,Vz,Dispxx,
     !                  Dispzz,Dispxz)

      integer e
      dimension x(NumNP),y(NumNP),KX(NumElD,4),Vx(NumNP),Vz(NumNP),
     !          Dispxx(NumNP),Dispzz(NumNP),Dispxz(NumNP),
     !          WeTab(3,2*NumEl),Beta(3),List(3)

      TanH(z)=(exp(z)-exp(-z))/(exp(z)+exp(-z))

      NumSEl=0
      do 13 e=1,NumEl
        NCorn=4
        if(KX(e,3).eq.KX(e,4)) NCorn=3
        do 12 n=1,NCorn-2
          NumSEl=NumSEl+1
          M1=KX(e,1)
          M2=KX(e,n+1)
          M3=KX(e,n+2)
          A=y(M2)-y(M1)
          B=x(M2)-x(M1)
          Beta(1)=atan2(A,B)
          A=y(M3)-y(M2)
          B=x(M3)-x(M2)
          Beta(2)=atan2(A,B)
          A=y(M1)-y(M3)
          B=x(M1)-x(M3)
          Beta(3)=atan2(A,B)
          List(1)=M1
          List(2)=M2
          List(3)=M3
          do 11 j=1,3
            k=j-1
            if(k.eq.0) k=3
            WeTab(k,NumSEl)=0.
            M1=List(j)
            jp1=j+1
            if(j.eq.3) jp1=1
            M2=List(jp1)
            Vxx=(Vx(M1)+Vx(M2))/2.
            Vzz=(Vz(M1)+Vz(M2))/2.
            if(abs(Vxx).lt.1.e-30.and.abs(Vzz).lt.1.e-30) goto 11
            BetaV=atan2(Vzz,Vxx)
            Delta=abs(BetaV-Beta(j))
            if(Delta.gt.0.314.and.abs(Delta-3.1416).gt.0.314) goto 11
            ALeng=sqrt((x(M2)-x(M1))**2+(y(M2)-y(M1))**2)
            CBeta=cos(Beta(j))
            SBeta=sin(Beta(j))
            Val=Vxx*CBeta+Vzz*SBeta
            VV=sqrt(Vxx*Vxx+Vzz*Vzz)
            DLL=(Dispxx(M1)+Dispxx(M2))/2.
            DLT=(Dispxz(M1)+Dispxz(M2))/2.
            DTT=(Dispzz(M1)+Dispzz(M2))/2.
            DAL=abs(DLL*CBeta*CBeta+2.0*CBeta*SBeta*DLT+DTT*SBeta*SBeta)
            Vel=VAL*ALeng
            Disp=2.0*DAL
            aa=11.
            if(Disp.gt.0.) aa=abs(Vel/Disp)
            if(Disp.lt.1.e-30.or.abs(Vel).lt.0.001*VV.or.aa.gt.10.) then
              if(abs(Vel).lt.0.001*VV) WeTab(k,NumSEl)=0.0
              if(Vel.gt.0.001*VV) WeTab(k,NumSEl)=1.0
              if(Vel.lt.-0.001*VV) WeTab(k,NumSEl)=-1.0
            else
              WeTab(k,NumSEl)=1.0/TanH(Vel/Disp)-Disp/Vel
            end if
11        continue
12      continue
13    continue
      return
      end

*************************************************************************

*     Calculate the maximum local Peclet and Courant numbers

      subroutine PeCour(jj,NSD,NumNP,NumEl,NumElD,NMat,x,y,Vx,Vz,KX,
     !                  MatNum,Dispxx,Dispzz,ChPar,theta,dt,dtMaxC,
     !                  Peclet,Courant,thSat,Temp,TDep,Conc,lLinear,
     !                  lUpW,lArtD,PeCr)

      logical lLinear,lUpW,lArtD
      dimension KX(NumElD,4),MatNum(NumNP),x(NumNP),y(NumNP),Vx(NumNP),
     !          Vz(NumNP),Dispxx(NumNP),Dispzz(NumNP),theta(NumNP),
     !          ChPar(NSD*16+4,NMat),thSat(NMat),TDep(NSD*16+4),
     !          Temp(NumNP),Conc(NSD,NumNP),Bi(3),Ci(3)

      Peclet=0.
      Courant=0.
      dtMaxC=1.e+30
      Tr=293.15
      R=8.314
      jjj=(jj-1)*16
      do 12 n=1,NumEl
        NUS=4
        if(KX(n,3).eq.KX(n,4)) NUS=3
        do 11 k=1,NUS-2
          PecX=99999.
          PecY=99999.
          dt1=1.e+30
          dt2=1.e+30
          i=KX(n,1)
          j=KX(n,k+1)
          l=KX(n,k+2)
          Ci(1)=x(l)-x(j)
          Ci(2)=x(i)-x(l)
          Ci(3)=x(j)-x(i)
          Bi(1)=y(j)-y(l)
          Bi(2)=y(l)-y(i)
          Bi(3)=y(i)-y(j)
          delX=amax1(abs(Ci(1)),abs(Ci(2)),abs(Ci(3)))
          delY=amax1(abs(Bi(1)),abs(Bi(2)),abs(Bi(3)))
          DxE=(Dispxx(i)+Dispxx(j)+Dispxx(l))/3.
          DzE=(Dispzz(i)+Dispzz(j)+Dispzz(l))/3.
          VxE=abs(Vx(i)+Vx(j)+Vx(l))/3.
          VzE=abs(Vz(i)+Vz(j)+Vz(l))/3.
          if(DxE.gt.0.) PecX=VxE*delX/DxE
          if(DzE.gt.0.) PecY=VzE*delY/DzE
          if(abs(VxE).lt.1.e-20) PecX=0.
          if(abs(VzE).lt.1.e-20) PecY=0.
          Peclet=amax1(Peclet,PecX,PecY)
          Peclet=amin1(Peclet,99999.)

          VxMax=amax1(abs(Vx(i))/theta(i),abs(Vx(j))/theta(j),abs(Vx(l))
     !                /theta(l))
          VzMax=amax1(abs(Vz(i))/theta(i),abs(Vz(j))/theta(j),abs(Vz(l))
     !                /theta(l))
          M1=MatNum(i)
          M2=MatNum(j)
          M3=MatNum(l)
          TT=(Temp(i)+Temp(j)+Temp(l))/3.+273.15
          f7=exp(TDep(jjj+7)*(TT-Tr)/R/TT/Tr)
          f8=exp(TDep(jjj+8)*(TT-Tr)/R/TT/Tr)
          f9=exp(TDep(jjj+9)*(TT-Tr)/R/TT/Tr)
          f10=exp(TDep(jjj+10)*(TT-Tr)/R/TT/Tr)
          a1=1.
          a2=1.
          a3=1.
          if(.not.lLinear) then
            if(Conc(jj,i).gt.0.) a1=ChPar(jjj+9,M1)*f9*Conc(jj,i)**
     !              (ChPar(jjj+9,M1)*f9-1.)/(1.+ChPar(jjj+8,M1)*f8*
     !              Conc(jj,i)**(ChPar(jjj+9,M1)*f9))**2
            if(Conc(jj,j).gt.0.) a2=ChPar(jjj+9,M2)*f9*Conc(jj,j)**
     !              (ChPar(jjj+9,M2)*f9-1.)/(1.+ChPar(jjj+8,M2)*f8*
     !              Conc(jj,j)**(ChPar(jjj+9,M2)*f9))**2
            if(Conc(jj,l).gt.0.) a3=ChPar(jjj+9,M3)*f9*Conc(jj,l)**
     !              (ChPar(jjj+9,M3)*f9-1.)/(1.+ChPar(jjj+8,M3)*f8*
     !              Conc(jj,l)**(ChPar(jjj+9,M3)*f9))**2
          end if
          R1=1.+(ChPar(1,M1)*ChPar(jjj+7,M1)*f7*a1*ChPar(4,M1)+
     !           (thSat(M1)-theta(i))*ChPar(jjj+10,M1)*f10)/theta(i)
          R2=1.+(ChPar(1,M2)*ChPar(jjj+7,M2)*f7*a2*ChPar(4,M2)+
     !           (thSat(M2)-theta(j))*ChPar(jjj+10,M2)*f10)/theta(j)
          R3=1.+(ChPar(1,M3)*ChPar(jjj+7,M3)*f7*a3*ChPar(4,M3)+
     !           (thSat(M3)-theta(l))*ChPar(jjj+10,M3)*f10)/theta(l)
          RMin=amin1(R1,R2,R3)
          if(RMin.gt.0.) CourX=VxMax*dt/delX/RMin
          if(RMin.gt.0.) CourY=VzMax*dt/delY/RMin
          Courant=amax1(Courant,CourX,CourY)

          CourMax=1.0
          if(.not.lUpW.and..not.lArtD)
     !      CourMax=amin1(1.,PeCr/amax1(0.5,PecX),PeCr/amax1(0.5,PecY))
          if(VxMax.gt.0.) dt1=CourMax*delX*RMin/VxMax
          if(VzMax.gt.0.) dt2=CourMax*delY*RMin/VzMax
          dtMaxC=amin1(dtMaxC,dt1,dt2)

11      continue
12    continue
      return
      end

*************************************************************************

*     Solve the global matrix equation for transport

      subroutine SolveT(A,B,MBand,MBandD,NumNP)

      double precision A,B,P,C,Sum
      dimension A(MBandD,NumNP),B(NumNP)

*     Reduction
      N1=NumNP-1
      do 12 k=1,N1
        P=1.d0/A(MBand,k)
        kk=k+1
        kc=MBand
        do 11 i=kk,NumNP
          kc=kc-1
          if(kc.le.0) goto 12
          C=-P*A(kc,i)
          A(kc,i)=C
          ii=kc+1
          L=kc+MBand-1
          do 11 j=ii,L
            jj=j+MBand-kc
            A(j,i)=A(j,i)+C*A(jj,k)
11      continue
12    continue

*     Back substitution
      do 14 i=2,NumNP
        jj=MBand+1-i
        ii=1
        if(jj.le.0) then
          jj=1
          ii=i-MBand+1
        end if
        Sum=0.
        do 13 j=jj,MBand-1
          Sum=Sum+A(j,i)*B(ii)
          ii=ii+1
13      continue
        B(i)=B(i)+Sum
14    continue
      B(NumNP)=B(NumNP)/A(MBand,NumNP)
      do 16 k=1,N1
        i=NumNP-k
        jj=i
        m=min0(2*MBand-1,MBand+k)
        Sum=0.
        do 15 j=MBand+1,m
          jj=jj+1
          Sum=Sum+A(j,i)*B(jj)
15      continue
        B(i)=(B(i)-Sum)/A(MBand,i)
16    continue
      return
      end

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
