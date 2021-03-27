* Source file OUTPUT2.FOR ||||||||||||||||||||||||||||||||||||||||||||||

      subroutine TLInf(NumNP,NumBP,Kode,Q,hNew,CumQ,Width,SWidth,KXB,t,
     !                 dt,TLevel,ShortF,TPrint,IterW,IterC,ItCum,rTop,
     !                 rRoot,vMeanR,hMeanT,hMeanR,hMeanG,AtmInF,SinkF,
     !                 CumQrT,CumQrR,CumQvR,NumKD,hMean,vMean,lWat,
     !                 lChem,rLen,Peclet,Courant,wCumT,wCumA)

      integer TLevel
      logical ShortF,SinkF,AtmInF,lWat,lChem
      dimension Q(NumNP),Kode(NumNP),KXB(NumBP),SWidth(NumKD),
     !          hNew(NumNP),Width(NumBP),CumQ(NumKD),hMean(NumKD),
     !          vMean(NumKD)

      if(TLevel.eq.1) then
        if(.not.lChem) then
          write(70,110)
        else
          if(lWat) then
            write(70,120)
          else
            write(70,130)
          end if
        end if
        write(71,140)
        write(77,150)
        if(lWat) write(78,160)
      end if
      if(lWat.or.TLevel.eq.1) then
        do 11 i=1,NumKD
          vMean(i)=0.
          hMean(i)=0.
11      continue
        do 12 i=1,NumBP
          n=KXB(i)
          j=iabs(Kode(n))
          if(j.eq.0) goto 12
          hMean(j)=hMean(j)+hNew(n)*Width(i)/SWidth(j)
          if(j.eq.4) vMean(j)=vMean(j)-Q(n)/SWidth(j)
12      continue
        hMeanT=hMean(4)
        hMeanG=hMean(3)
        do 13 i=1,NumNP
          j=iabs(Kode(i))
          wCumA=wCumA+abs(Q(i))*dt
          if(j.ne.0.and.j.ne.4) then
            vMean(j)=vMean(j)-Q(i)
          end if
13      continue
        if(.not.lWat.and.TLevel.eq.1) then
          write(71,170) t,rTop,rRoot,vMean(4),vMeanR,vMean(3),vMean(1),
     !                  vMean(2),(vMean(i),i=5,NumKD)
          write(77,180) t,hMean(4),hMeanR,hMean(3),hMean(1),hMean(2),
     !                  (hMean(i),i=5,NumKD)
        end if
      end if
      if(lWat) then
        wCumA=wCumA+abs(vMeanR*dt*rLen)
        wCumT=CumQvR
        do 14 j=1,NumKD
          if(j.eq.4) then
            CumQ(j)=CumQ(j)+vMean(j)*dt*SWidth(4)
          else
            CumQ(j)=CumQ(j)+vMean(j)*dt
          end if
          wCumT=wCumT+CumQ(j)
14      continue
        CumQrT=CumQrT+rTop  *dt*Swidth(4)
        CumQrR=CumQrR+rRoot *dt*rLen
        CumQvR=CumQvR+vMeanR*dt*rLen
      end if
      if(.not.ShortF.or.abs(TPrint-t).lt.0.001*dt) then
        if(lWat) then
          write(71,170) t,rTop,rRoot,vMean(4),vMeanR,vMean(3),vMean(1),
     !                  vMean(2),(vMean(i),i=5,NumKD)
          write(77,180) t,hMean(4),hMeanR,hMean(3),hMean(1),hMean(2),
     !                  (hMean(i),i=5,NumKD)
          write(78,190) t,CumQrT,CumQrR,CumQ(4),CumQvR,CumQ(3),
     !                  CumQ(1),CumQ(2),(CumQ(i),i=5,NumKD)
        end if
        if(lChem) then
          if(lWat) then
            write(70,200) TLevel,t,dt,IterW,IterC,ItCum,Peclet,Courant
          else
            write(70,210) TLevel,t,dt,IterC,Peclet,Courant
          end if
        else
          write(70,220) TLevel,t,dt,IterW,ItCum
        end if
      end if
      if(lWat) then
        if(SinkF) then
          write(*,230) t,IterW,ItCum,CumQ(4),CumQvR,CumQ(3),hMean(4),
     !                 hMeanR,hMean(3)
        else
          if(AtmInF) then
            write(*,230) t,IterW,ItCum,CumQ(4),CumQ(1),CumQ(3),hMean(4),
     !                   hMean(1),hMean(3)
          else
            write(*,230) t,IterW,ItCum,vMean(1),CumQ(1),CumQ(2),
     !                   hMean(1),hMean(2)
          end if
        end if
      end if

110   format(//' TLevel   Time         dt      Iter  ItCum'/)
120   format(//' TLevel   Time         dt     IterW IterC   ItCum  Pecle
     !t   Courant'/)
130   format(//' TLevel   Time         dt     IterC   Peclet  Courant'/)
140   format(' All fluxes (v) are positive out of the region'//
     !'     Time       rAtm       rRoot      vAtm       vRoot      ',
     !'vKode3     vKode1     vSeep      vKode5     vKode6 ...'/
     !'      [T]      [L/T]       [L/T]     [L/T]       [L/T]      ',
     !' [V/T]      [V/T]     [V/T]       [V/T]      [V/T]'/)
150   format(//
     !'     Time         hAtm       hRoot     hKode3     hKode1 ',
     !'     hSeep     hKode5     hKode6 ... '/
     !'      [T]          [L]        [L]        [L]        [L]  ',
     !'      [L]        [L]        [L]'/)
160   format(
     !' All cumulative fluxes (CumQ) are positive out of the region'//
     !'     Time     CumQAP      CumQRP     CumQA      CumQR     CumQ3',
     !'       CumQ1      CumQS      CumQ5       CumQ6 ....'/
     !'      [T]       [V]         [V]       [V]        [V]       [V] ',
     !'        [V]        [V]        [V]         [V]'/)
170   format(f12.4,9e11.3)
180   format(f12.4,9f11.1)
190   format(f12.4,9e11.3)
200   format(i5,2e12.3,2i5,i7,2f10.3)
210   format(i5,2e12.3,i5,2f10.3)
220   format(i5,2e12.3,i5,i7)
230   format(f14.4,i3,i6,1x,3e11.3,3f7.0)
      return
      end

************************************************************************

      subroutine ALInf(t,CumQ,hMeanT,hMeanR,hMeanG,ALevel,CumQrT,CumQrR,
     !                 CumQvR,NumKD)

      integer ALevel
      dimension CumQ(NumKD)

      if(ALevel.eq.1) write(72,110)
      write(72,120) t,CumQrT,CumQrR,CumQ(4),CumQvR,CumQ(3),hMeanT,
     !              hMeanR,hMeanG,ALevel

110   format(
     !' All cumulative fluxes (CumQ) are positive out of the region'//
     !'      Time      CumQAP     CumQRP     CumQA      CumQR ',
     !'     CumQ3        hAtm       hRoot     hKode3    A-level'/
     !'      [T]         [V]        [V]       [V]        [V]  ',
     !'      [V]          [L]        [L]        [L] '/)
120   format(f12.4,5e11.3,3f11.1,i8)
      return
      end

************************************************************************

      subroutine SubReg(NumEl,NumElD,NumNP,NMat,hNew,ThO,ThN,x,y,MatNum,
     !                  LayNum,KX,KAT,t,dt,NLay,PLevel,lWat,lChem,Conc,
     !                  ChPar,wCumA,wCumT,cCumA,cCumT,wVolI,cVolI,WatIn,
     !                  SolIn,thSat,lTemp,Temp,TPar,TDep,NS,NSD,Sorb,
     !                  lLinear,lEquil)

      logical lWat,lChem,lTemp,lLinear,lEquil
      integer PLevel,e
      dimension hNew(NumNP),x(NumNP),y(NumNP),MatNum(NumNP),Temp(NumNP),
     !          ThO(NumNP),ChPar(NSD*16+4,NMat),ThN(NumNP),WatIn(NumEl),
     !          SolIn(NumEl),KX(NumElD,4),LayNum(NumEl),Conc(NSD,NumNP),
     !          thSat(NMat),TDep(NSD*16+4),Area(10),hMean(10),ConVol(8),
     !          SubVol(10),SubCha(10),cMean(8,10),ConSub(8,10),cTot(8),
     !          cCumA(NS),cCumT(NS),cVolI(NS),TPar(10,NMat),SubT(10),
     !          TMean(10),Sorb(NSD,NumNP)

      xMul=1.
      ATot=0.
      Tr=293.15
      R=8.314
      if(lWat.or.PLevel.eq.0) then
        Volume=0.
        Change=0.
        hTot=0.
        DeltW=0.
      end if
      if(lTemp) then
        TTot=0.
        TVol=0.
      end if
      if(lChem) then
        do 11 jj=1,NS
          cTot(jj)=0.
          ConVol(jj)=0.
11      continue
        DeltC=0.
      end if
      do 13 i=1,NLay
        Area(i)=0.
        if(lWat.or.PLevel.eq.0) then
          SubVol(i)=0.
          SubCha(i)=0.
          hMean(i)=0.
        end if
        if(lTemp) then
          SubT(i)=0.
          TMean(i)=0.
        end if
        if(lChem) then
          do 12 jj=1,NS
            ConSub(jj,i)=0.
            cMean(jj,i)=0.
12        continue
        end if
13    continue
      do 16 e=1,NumEl
        Lay=LayNum(e)
        wEl=0.
        cEl=0.
        NUS=4
        if(KX(e,3).eq.KX(e,4)) NUS=3
        do 15 k=1,NUS-2
          i=KX(e,1)
          j=KX(e,k+1)
          l=KX(e,k+2)
          Mi=MatNum(i)
          Mj=MatNum(j)
          Mk=MatNum(l)
          Cj=x(i)-x(l)
          Ck=x(j)-x(i)
          Bj=y(l)-y(i)
          Bk=y(i)-y(j)
          TT=(Temp(i)+Temp(j)+Temp(l))/3.+273.15
          if(KAT.eq.1) xMul=2.*3.1416*(x(i)+x(j)+x(l))/3.
          AE=xMul*(Ck*Bj-Cj*Bk)/2.
          Area(Lay)=Area(Lay)+AE
          if(lWat.or.PLevel.eq.0) then
            hE=(hNew(i)+hNew(j)+hNew(l))/3.
            VNewE=AE*(thN(i)+thN(j)+thN(l))/3.
            VOldE=AE*(thO(i)+thO(j)+thO(l))/3.
            Volume=Volume+VNewE
            wEl=wEl+VNewE
            Change=Change+(VNewE-VOldE)/dt
            SubVol(Lay)=SubVol(Lay)+VNewE
            SubCha(Lay)=SubCha(Lay)+(VNewE-VOldE)/dt
            hTot=hTot+hE*AE
            hMean(Lay)=hMean(Lay)+hE*AE
          end if
          if(lTemp) then
            TE=(Temp(i)+Temp(j)+Temp(l))/3.
            TNewE=AE*((Temp(i)+273.15)*(TPar(1,Mi)*TPar(8,Mi)+
     !                     TPar(2,Mi)*TPar(9,Mi)+TPar(10,Mi)*ThN(i))+
     !                (Temp(j)+273.15)*(TPar(1,Mj)*TPar(8,Mj)+
     !                     TPar(2,Mj)*TPar(9,Mj)+TPar(10,Mj)*ThN(j))+
     !                (Temp(l)+273.15)*(TPar(1,Mk)*TPar(8,Mk)+
     !                     TPar(2,Mk)*TPar(9,Mk)+TPar(10,Mk)*ThN(l)))/3.
            TVol=TVol+TNewE
            SubT(Lay)=SubT(Lay)+TNewE
            TTot=TTot+TE*AE
            TMean(Lay)=TMean(Lay)+TE*AE
          end if
          if(lChem) then
            do 14 jj=1,NS
              jjj=(jj-1)*16
              cE=(Conc(jj,i)+Conc(jj,j)+Conc(jj,l))/3.
              f7=exp(TDep(jjj+7)*(TT-Tr)/R/TT/Tr)
              f8=exp(TDep(jjj+8)*(TT-Tr)/R/TT/Tr)
              f9=exp(TDep(jjj+9)*(TT-Tr)/R/TT/Tr)
              f10=exp(TDep(jjj+10)*(TT-Tr)/R/TT/Tr)
              a1=1.
              a2=1.
              a3=1.
              if(.not.lLinear) then
                if(Conc(jj,i).gt.0.) a1=Conc(jj,i)**
     !               (ChPar(jjj+9,Mi)*f9-1.)/(1.+ChPar(jjj+8,Mi)*f8*
     !                Conc(jj,i)**(ChPar(jjj+9,Mi)*f9))
                if(Conc(jj,j).gt.0.) a2=Conc(jj,j)**
     !               (ChPar(jjj+9,Mj)*f9-1.)/(1.+ChPar(jjj+8,Mj)*f8*
     !                Conc(jj,j)**(ChPar(jjj+9,Mj)*f9))
                if(Conc(jj,l).gt.0.) a3=Conc(jj,l)**
     !               (ChPar(jjj+9,Mk)*f9-1.)/(1.+ChPar(jjj+8,Mk)*f8*
     !                Conc(jj,l)**(ChPar(jjj+9,Mk)*f9))
              end if
              cNewE=AE*
     !            (Conc(jj,i)*(thN(i)+
     !                ChPar(4,Mi)*ChPar(1,Mi)*ChPar(jjj+7,Mi)*f7*a1+
     !                (thSat(Mi)-ThN(i))*ChPar(jjj+10,Mi)*f10)+
     !             Conc(jj,j)*(thN(j)+
     !                ChPar(4,Mj)*ChPar(1,Mj)*ChPar(jjj+7,Mj)*f7*a2+
     !                (thSat(Mj)-ThN(j))*ChPar(jjj+10,Mj)*f10)+
     !             Conc(jj,l)*(thN(l)+
     !                ChPar(4,Mk)*ChPar(1,Mk)*ChPar(jjj+7,Mk)*f7*a3+
     !                (thSat(Mk)-ThN(l))*ChPar(jjj+10,Mk)*f10))/3.
              if(.not.lEquil) cNewE=cNewE+AE*(ChPar(1,Mi)*Sorb(jj,i)+
     !                 ChPar(1,Mj)*Sorb(jj,j)+ChPar(1,Mk)*Sorb(jj,l))/3.
              ConVol(jj)=ConVol(jj)+cNewE
              if(jj.eq.1) cEl=cEl+cNewE
              ConSub(jj,Lay)=ConSub(jj,Lay)+cNewE
              cTot(jj)=cTot(jj)+cE*AE
              cMean(jj,Lay)=cMean(jj,Lay)+cE*AE
14          continue
          end if
          if(k.eq.NUS-2) then
            if(PLevel.eq.0) then
              if(lWat) WatIn(e)=wEl
              if(lChem) SolIn(e)=cEl
            else
              if(lWat) DeltW=DeltW+abs(WatIn(e)-wEl)
              if(lChem) DeltC=DeltC+abs(SolIn(e)-cEl)
            end if
          end if
15      continue
16    continue
      do 18 Lay=1,NLay
        if(lWat.or.PLevel.eq.0) hMean(Lay)=hMean(Lay)/Area(Lay)
        if(lTemp) TMean(Lay)=TMean(Lay)/Area(Lay)
        do 17 jj=1,NS
          if(lChem) cMean(jj,Lay)=cMean(jj,Lay)/Area(Lay)
17      continue
        ATot=ATot+Area(Lay)
18    continue
      if(lWat.or.PLevel.eq.0) hTot=hTot/ATot
      if(lTemp) TTot=TTot/ATot
      do 19 jj=1,NS
        if(lChem) cTot(jj)=cTot(jj)/ATot
19    continue
      if(PLevel.eq.0) write(80,110)
      write(80,120) t,     (        i,i=1,NLay)
      write(80,130) ATot,  (  Area(i),i=1,NLay)
      if(lWat.or.PLevel.eq.0) then
        write(80,140) Volume,(SubVol(i),i=1,NLay)
        write(80,150) Change,(SubCha(i),i=1,NLay)
        write(80,160) hTot,  ( hMean(i),i=1,NLay)
      end if
      if(lTemp) then
        write(80,170) TVol,  (  SubT(i),i=1,NLay)
        write(80,180) TTot,  ( TMean(i),i=1,NLay)
      end if
       if(lChem) then
        do 20 jj=1,NS
          write(80,190) jj,ConVol(jj),(ConSub(jj,i),i=1,NLay)
          write(80,200) jj,cTot(jj),  ( cMean(jj,i),i=1,NLay)
20      continue
      end if

*     Mass balance calculation
      if(PLevel.eq.0) then
        wVolI=Volume
        do 21 jj=1,NS
          cVolI(jj)=ConVol(jj)
21      continue
      else
        if(lWat) then
          wBalT=Volume-wVolI+wCumT
          write(80,210) wBalT
          ww=amax1(DeltW,wCumA)
          if(ww.ge.1.e-25) then
            wBalR=abs(wBalT)/ww*100.
            write(80,220) wBalR
          end if
        end if
        if(lChem) then
          do 22 jj=1,NS
            cBalT=ConVol(jj)-cVolI(jj)+cCumT(jj)
            write(80,230) jj,cBalT
            cc=amax1(DeltC,cCumA(jj))
            if(cc.ge.1.e-25) then
              cBalR=abs(cBalT)/cc*100.
              write(80,240) jj,cBalR
            end if
22        continue
        end if
      end if

110   format(/' Time [T]             Total     Sub-region number ...')
120   format(/f12.4,16x,10(i7,4x))
130   format(' Area    [V]       ',e11.3,10e11.3)
140   format(' Volume  [V]       ',e11.3,10e11.3)
150   format(' InFlow  [V/T]     ',e11.3,10e11.3)
160   format(' hMean   [L]       ',e11.3,10f11.1)
170   format(' HeatVol [VM/L/T2] ',e11.3,10e11.3)
180   format(' tMean   [K]       ',f11.3,10f11.3)
190   format(' ConcVol [VM/L3] ',i1,1x,e11.3,10e11.3)
200   format(' cMean   [M/L3]  ',i1,1x,e11.3,10e11.3)
210   format(' WatBalT [V]       ',e11.3)
220   format(' WatBalR [%]       ',f11.3)
230   format(' CncBalT [VM/L3] ',i1,1x,e11.3)
240   format(' CncBalR [%]     ',i1,1x,f11.3)
      return
      end

************************************************************************

      subroutine BouOut(NumNP,NumBP,NS,NSD,t,hNew,theta,Q,Width,KXB,
     !                  Kode,x,y,Temp,Conc)

      dimension hNew(NumNP),Q(NumNP),Width(NumBP),theta(NumNP),
     !          KXB(NumBP),Kode(NumNP),x(NumNP),y(NumNP),Temp(NumNP),
     !          Conc(NSD,NumNP)

      write(79,110) t
      ii=0
      do 12 i=1,NumNP
        if(Kode(i).ne.0) then
          do 11 j=1,NumBP
            n=KXB(j)
            if(n.eq.i) then
              ii=ii+1
              v=-Q(i)/Width(j)
              write(79,120) ii,i,x(i),y(i),Kode(i),Q(i),v,hNew(i),
     !                      theta(i),Temp(i),(Conc(jj,i),jj=1,NS)
              goto 12
            end if
11        continue
          ii=ii+1
          write(79,130) ii,i,x(i),y(i),Kode(i),Q(i),hNew(i),theta(i),
     !                  Temp(i),(Conc(jj,i),jj=1,NS)
        end if
12    continue

110   format(//' Time:',f12.4//
     !'    i    n    x      z    Code     Q          v       ',
     !'    h       th       Temp     Conc'/
     !'                                 [V/T]      [L/T]     ',
     !'   [L]     [-]        [C]    [M/L3]'/)
120   format(2i5,2f7.1,i5,2e11.3,f11.1,f7.3,f10.2,8e10.3)
130   format(2i5,2f7.1,i5,e11.3,11x,f11.1,f7.3,f10.2,8e10.3)
      return
      end

************************************************************************

      subroutine SolInf(K,NumNP,NS,NSD,Kode,Qc,t,dt,IterC,TLevel,ShortF,
     !                  TPrint,NumKD,SMean,ChemS,CumCh0,CumCh1,CumChR,
     !                  cCumA,cCumT,lWat)

      integer TLevel
      logical ShortF,lWat
      dimension Qc(NSD,NumNP),Kode(NumNP),CumCh0(NS),CumCh1(NS),
     !          CumChR(NS),ChemS(NSD,NumKD),SMean(NSD,NumKD),cCumA(NS),
     !          cCumT(NS)

      do 14 jj=1,NS
        if(K.eq.1) then
          do 11 i=1,NumKD
            SMean(jj,i)=0.
11        continue
          do 12 i=1,NumNP
            j=iabs(Kode(i))
            if(j.ne.0) then
              SMean(jj,j)=Smean(jj,j)-Qc(jj,i)
            end if
12        continue
          cCumA(jj)=abs(CumCh0(jj))+abs(CumCh1(jj))+abs(CumChR(jj))
          cCumT(jj)=CumCh0(jj)+CumCh1(jj)+CumChR(jj)
          do 13 j=1,NumKD
            ChemS(jj,j)=ChemS(jj,j)+SMean(jj,j)*dt
            cCumT(jj)=cCumT(jj)+ChemS(jj,j)
            cCumA(jj)=cCumA(jj)+abs(ChemS(jj,j))
13        continue
          if(TLevel.eq.1) write(83+jj,110)
          if(.not.ShortF.or.abs(TPrint-t).lt.0.001*dt)
     !      write(83+jj,120) t,CumCh0(jj),CumCh1(jj),CumChR(jj),
     !                   (ChemS(jj,j),j=1,NumKD),(SMean(jj,i),i=1,NumKD)
        else if(ShortF) then
          write(83+jj,120) t,CumCh0(jj),CumCh1(jj),CumChR(jj),
     !                   (ChemS(jj,j),j=1,NumKD),(SMean(jj,i),i=1,NumKD)
       end if
14    continue
      if(.not.lWat) write(*,130) t,IterC,TLevel,CumCh0(1),CumCh1(1),
     !                           ChemS(1,1)

110   format(' All solute fluxes (SMean) and cumulative solute fluxes',
     !' (ChemS) are positive out of the region'//
     !'     Time     CumCh0     CumCh1     CumChR   ',20('-'),'  ChemS(i
     !),i=1,NumKD  ',22('-'),'  ',21('-'),'  SMean(j),j=1,NumKD ',22('-'
     !)/'      [T]    [VM/L3]    [VM/L3]    [VM/L3]',31(' '),'[VM/L3]',5
     !9(' '),'[VM/T/L3]'/)
120   format(f10.2,15e11.3)
130   format(f14.4,2i6,1x,3e11.3,2x,2e11.3)
      return
      end

************************************************************************

      subroutine ObsNod(t,NumNP,NS,NSD,NObs,NObsD,Node,hNew,ThNew,Conc,
     !                  Sorb,Temp)

      dimension Node(NObsD),hNew(NumNP),ThNew(NumNP),Conc(NSD,NumNP),
     !          Sorb(NSD,NumNP),Temp(NumNP)

      if(NS.eq.0)
     ! write(92,100) t,(hNew(Node(i)),ThNew(Node(i)),Temp(Node(i)),
     !               i=1,NObs)
      if(NS.eq.1)
     ! write(92,110) t,(hNew(Node(i)),ThNew(Node(i)),Temp(Node(i)),
     !   (Conc(j,Node(i)),j=1,NS),(Sorb(j,Node(i)),j=1,NS),i=1,NObs)
      if(NS.eq.2)
     ! write(92,120) t,(hNew(Node(i)),ThNew(Node(i)),Temp(Node(i)),
     !   (Conc(j,Node(i)),j=1,NS),(Sorb(j,Node(i)),j=1,NS),i=1,NObs)
      if(NS.eq.3)
     ! write(92,130) t,(hNew(Node(i)),ThNew(Node(i)),Temp(Node(i)),
     !   (Conc(j,Node(i)),j=1,NS),(Sorb(j,Node(i)),j=1,NS),i=1,NObs)
      if(NS.eq.4)
     ! write(92,140) t,(hNew(Node(i)),ThNew(Node(i)),Temp(Node(i)),
     !   (Conc(j,Node(i)),j=1,NS),(Sorb(j,Node(i)),j=1,NS),i=1,NObs)
      if(NS.eq.5)
     ! write(92,150) t,(hNew(Node(i)),ThNew(Node(i)),Temp(Node(i)),
     !   (Conc(j,Node(i)),j=1,NS),(Sorb(j,Node(i)),j=1,NS),i=1,NObs)
      if(NS.eq.6)
     ! write(92,160) t,(hNew(Node(i)),ThNew(Node(i)),Temp(Node(i)),
     !   (Conc(j,Node(i)),j=1,NS),(Sorb(j,Node(i)),j=1,NS),i=1,NObs)
      if(NS.eq.7)
     ! write(92,170) t,(hNew(Node(i)),ThNew(Node(i)),Temp(Node(i)),
     !   (Conc(j,Node(i)),j=1,NS),(Sorb(j,Node(i)),j=1,NS),i=1,NObs)
      if(NS.eq.8)
     ! write(92,180) t,(hNew(Node(i)),ThNew(Node(i)),Temp(Node(i)),
     !   (Conc(j,Node(i)),j=1,NS),(Sorb(j,Node(i)),j=1,NS),i=1,NObs)

100   format(f11.3,5(f11.3,f8.4,f9.3,        2x))
110   format(f11.3,5(f11.3,f8.4,f9.3, 2e11.3,2x))
120   format(f11.3,5(f11.3,f8.4,f9.3, 4e11.3,2x))
130   format(f11.3,5(f11.3,f8.4,f9.3, 6e11.3,2x))
140   format(f11.3,5(f11.3,f8.4,f9.3, 8e11.3,2x))
150   format(f11.3,5(f11.3,f8.4,f9.3,10e11.3,2x))
160   format(f11.3,5(f11.3,f8.4,f9.3,12e11.3,2x))
170   format(f11.3,5(f11.3,f8.4,f9.3,14e11.3,2x))
180   format(f11.3,5(f11.3,f8.4,f9.3,16e11.3,2x))
      return
      end

************************************************************************

      subroutine hOut(hNew,x,y,NumNP,t,IJ)

      dimension hNew(NumNP),x(NumNP),y(NumNP)

      write(75,110) t
      L1=(IJ-1)/10+1
      do 12 n=1,NumNP,IJ
        do 11 L=1,L1
          m=n+(L-1)*10
          k=m+9
          if(L.eq.L1) k=n+IJ-1
          write(75,120) m,x(m),y(m),(hNew(j),j=m,k)
11      continue
12    continue

110   format(//' Time  ***',f12.4,' ***'//
     !         '    n    x(n)   z(n)       h(n)      h(n+1) ...'/)
120   format(i5,2f8.1,10f10.1)
      return
      end

************************************************************************

      subroutine thOut(theta,x,y,NumNP,t,IJ)

      dimension theta(NumNP),x(NumNP),y(NumNP)

      write(76,110) t
      L1=(IJ-1)/16+1
      do 12 n=1,NumNP,IJ
        do 11 L=1,L1
          m=n+(L-1)*16
          k=m+15
          if(L.eq.L1) k=n+IJ-1
          write(76,120) m,x(m),y(m),(theta(j),j=m,k)
11      continue
12    continue

110   format(//' Time  ***',f12.4,' ***'//
     !         '    n    x(n)   z(n)      th(n)     th(n+1) ...'/)
120   format(i5,2f8.1,16f6.3)
      return
      end

************************************************************************

      subroutine QOut(Q,x,y,NumNP,t,IJ)

      dimension Q(NumNP),x(NumNP),y(NumNP)

      write(73,110) t
      L1=(IJ-1)/10+1
      do 12 n=1,NumNP,IJ
        do 11 L=1,L1
          m=n+(L-1)*10
          k=m+9
          if(L.eq.L1) k=n+IJ-1
          write(73,120) m,x(m),y(m),(Q(j),j=m,k)
11      continue
12    continue

110   format(//' Time  ***',f12.4,' ***'//
     !         '    n    x(n)   z(n)       Q(n)      Q(n+1) ...'/)
120   format(i5,2f8.1,10e11.3)
      return
      end

************************************************************************

      subroutine FlxOut(Vx,Vz,x,y,NumNP,t,IJ)

      dimension x(NumNP),y(NumNP),Vx(NumNP),Vz(NumNP)

      write(81,110) t
      write(82,120) t
      L1=(IJ-1)/10+1
      do 12 n=1,NumNP,IJ
        do 11 L=1,L1
          m=n+(L-1)*10
          k=m+9
          if(L.eq.L1) k=n+IJ-1
          write(81,130) m,x(m),y(m),(Vz(j),j=m,k)
          write(82,130) m,x(m),y(m),(Vx(j),j=m,k)
11      continue
12    continue

110   format(//' Time  ***',f12.4,' ***'//
     !         '    n    x(n)  z(n)     vz(n)     vz(n+1) ...'/)
120   format(//' Time  ***',f12.4,' ***'//
     !         '    n    x(n)  z(n)     vx(n)     vx(n+1) ...'/)
130   format(i5,2f8.1,10e10.2)
      return
      end

************************************************************************

      subroutine cOut(jj,NumNP,Conc,x,y,t,IJ,NSD)

      dimension Conc(NSD,NumNP),x(NumNP),y(NumNP)

      write(74,110) t,jj
      L1=(IJ-1)/10+1
      do 12 n=1,NumNP,IJ
        do 11 L=1,L1
          m=n+(L-1)*10
          k=m+9
          if(L.eq.L1) k=n+IJ-1
          write(74,120) m,x(m),y(m),(Conc(jj,j),j=m,k)
11      continue
12    continue

110   format(//' Time  ***',f12.4,' *** Species ',i3//
     !         '    n    x(n)   z(n)      Conc(n)   Conc(n+1)  ...'/)
120   format(i5,2f8.1,15f7.3) !10e11.3)
      return
      end

************************************************************************

      subroutine TOut(NumNP,Temp,x,y,t,IJ)

      dimension Temp(NumNP),x(NumNP),y(NumNP)

      write(83,110) t
      L1=(IJ-1)/10+1
      do 12 n=1,NumNP,IJ
        do 11 L=1,L1
          m=n+(L-1)*10
         k=m+9
          if(L.eq.L1) k=n+IJ-1
          write(83,120) m,x(m),y(m),(Temp(j),j=m,k)
11      continue
12    continue

110   format(//' Time  ***',f12.4,' ***'//
     !         '    n    x(n)   z(n)      Temp(n)   Temp(n+1)  ...'/)
120   format(i5,2f8.1,10e10.3)
      return
      end

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||