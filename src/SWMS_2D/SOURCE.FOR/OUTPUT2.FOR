* Source file OUTPUT2.FOR ||||||||||||||||||||||||||||||||||||||||||||||

      subroutine TLInf(NumNP,NumBP,Kode,Q,hNew,CumQ,Width,SWidth,KXB,t,
     !                 dt,TLevel,ShortF,TPrint,Iter,ItCum,rTop,rRoot,
     !                 vMeanR,hMeanT,hMeanR,hMeanG,AtmInF,SinkF,CumQrT,
     !                 CumQrR,CumQvR,NumKD,hMean,vMean,lWat,lChem,rLen,
     !                 Peclet,Courant,wCumT,wCumA)

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
            write(70,200) TLevel,t,dt,Iter,ItCum,Peclet,Courant
          else
            write(70,210) TLevel,t,dt,Peclet,Courant
          end if
        else
          write(70,200) TLevel,t,dt,Iter,ItCum
        end if
      end if
      if(lWat) then
        if(SinkF) then
          write(*,220) t,Iter,ItCum,CumQ(4),CumQvR,CumQ(3),hMean(4),
     !                 hMeanR,hMean(3)
        else
          if(AtmInF) then
            write(*,220) t,Iter,ItCum,CumQ(4),CumQ(1),CumQ(3),hMean(4),
     !                   hMean(1),hMean(3)
          else
            write(*,220) t,Iter,ItCum,vMean(1),CumQ(1),CumQ(2),hMean(1),
     !                   hMean(2)
          end if
        end if
      end if

110   format(//' TLevel   Time         dt      Iter  ItCum'/)
120   format(//' TLevel   Time         dt      Iter  ItCum  Peclet   ',
     !'Courant'/)
130   format(//' TLevel   Time         dt        Peclet   Courant'/)
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
200   format(i5,2e12.3,i5,i6,2f10.3)
210   format(i5,2e12.3,2f10.3)
220   format(f14.4,i3,i6,1x,3e11.3,3f7.0)
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
     !                  SolIn)

      logical lWat,lChem
      integer PLevel,e
      dimension hNew(NumNP),x(NumNP),y(NumNP),MatNum(NumNP),Conc(NumNP),
     !          KX(NumElD,4),ChPar(10,NMat),LayNum(NumEl),ThO(NumNP),
     !          ThN(NumNP),WatIn(NumEl),SolIn(NumEl),Area(10),hMean(10),
     !          SubVol(10),SubCha(10),cMean(10),ConSub(10)

      xMul=1.
      ATot=0.
      if(lWat.or.PLevel.le.1) then
        Volume=0.
        Change=0.
        hTot=0.
        DeltW=0.
      end if
      if(lChem) then
        cTot=0.
        ConVol=0.
        DeltC=0.
      end if
      do 11 i=1,NLay
        Area(i)=0.
        if(lWat.or.PLevel.le.1) then
          SubVol(i)=0.
          SubCha(i)=0.
          hMean(i)=0.
        end if
        if(lChem) then
          ConSub(i)=0.
          cMean(i)=0.
        end if
11    continue
      do 13 e=1,NumEl
        Lay=LayNum(e)
        wEl=0.
        cEl=0.
        NUS=4
        if(KX(e,3).eq.KX(e,4)) NUS=3
        do 12 k=1,NUS-2
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
          if(KAT.eq.1) xMul=2.*3.1416*(x(i)+x(j)+x(l))/3.
          AE=xMul*(Ck*Bj-Cj*Bk)/2.
          Area(Lay)=Area(Lay)+AE
          if(lWat.or.PLevel.le.1) then
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
          if(lChem) then
            cE=(Conc(i)+Conc(j)+Conc(l))/3.
            cNewE=AE*((thN(i)+ChPar(1,Mi)*ChPar(5,Mi))*Conc(i)+
     !                (thN(j)+ChPar(1,Mj)*ChPar(5,Mj))*Conc(j)+
     !                (thN(l)+ChPar(1,Mk)*ChPar(5,Mk))*Conc(l))/3.
            ConVol=ConVol+cNewE
            cEl=cEl+cNewE
            ConSub(Lay)=ConSub(Lay)+cNewE
            cTot=cTot+cE*AE
            cMean(Lay)=cMean(Lay)+cE*AE
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
12      continue
13    continue
      do 14 Lay=1,NLay
        if(lWat.or.PLevel.le.1) hMean(Lay)=hMean(Lay)/Area(Lay)
        if(lChem) cMean(Lay)=cMean(Lay)/Area(Lay)
        ATot=ATot+Area(Lay)
14    continue
      if(lWat.or.PLevel.le.1) hTot=hTot/ATot
      if(lChem) cTot=cTot/ATot
      if(PLevel.eq.0) write(80,110)
      write(80,120) t,     (        i,i=1,NLay)
      write(80,130) ATot,  (  Area(i),i=1,NLay)
      if(lWat.or.PLevel.le.1) then
        write(80,140) Volume,(SubVol(i),i=1,NLay)
        write(80,150) Change,(SubCha(i),i=1,NLay)
        write(80,160) hTot,  ( hMean(i),i=1,NLay)
      end if
      if(lChem) then
        write(80,170) ConVol,(ConSub(i),i=1,NLay)
        write(80,180) cTot,  ( cMean(i),i=1,NLay)
      end if

*     Mass balance calculation
      if(PLevel.eq.0) then
        wVolI=Volume
        if(lChem) cVolI=ConVol
      else
        if(lWat) then
          wBalT=Volume-wVolI+wCumT
          write(80,190) wBalT
          ww=amax1(DeltW,wCumA)
          if(ww.ge.1.e-25) then
            wBalR=abs(wBalT)/ww*100.
            write(80,200) wBalR
          end if
        end if
        if(lChem) then
          cBalT=ConVol-cVolI+cCumT
          write(80,210) cBalT
          cc=amax1(DeltC,cCumA)
          if(cc.ge.1.e-25) then
            cBalR=abs(cBalT)/cc*100.
            write(80,220) cBalR
          end if
        end if
      end if

110   format(/' Time [T]             Total     Sub-region number ...')
120   format(/f12.4,16x,10(i7,4x))
130   format(' Area    [V]       ',e11.3,10e11.3)
140   format(' Volume  [V]       ',e11.3,10e11.3)
150   format(' InFlow  [V/T]     ',e11.3,10e11.3)
160   format(' hMean   [L]       ',e11.3,10f11.1)
170   format(' ConcVol [VM/L3]   ',e11.3,10e11.3)
180   format(' cMean   [M/L3]    ',e11.3,10e11.3)
190   format(' WatBalT [V]       ',e11.3)
200   format(' WatBalR [%]       ',f11.3)
210   format(' CncBalT [VM/L3]   ',e11.3)
220   format(' CncBalR [%]       ',f11.3)
      return
      end

************************************************************************

      subroutine BouOut(NumNP,NumBP,t,hNew,theta,Q,Width,KXB,Kode,x,y,
     !                  Conc)

      dimension hNew(NumNP),Q(NumNP),Width(NumBP),theta(NumNP),
     !          KXB(NumBP),Kode(NumNP),x(NumNP),y(NumNP),Conc(NumNP)

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
     !                      theta(i),Conc(i)
              goto 12
            end if
11        continue
          ii=ii+1
          write(79,130) ii,i,x(i),y(i),Kode(i),Q(i),hNew(i),theta(i),
     !                  Conc(i)
        end if
12    continue

110   format(//' Time:',f12.4//
     !'    i    n    x      z    Code     Q          v       ',
     !'    h       th      Conc'/
     !'                                 [V/T]      [L/T]     ',
     !'   [L]     [-]     [M/L3]'/)
120   format(2i5,2f7.1,i5,2e11.3,f11.1,f7.3,e10.3)
130   format(2i5,2f7.1,i5,e11.3,11x,f11.1,f7.3,e10.3)
      return
      end

************************************************************************

      subroutine SolInf(NumNP,Kode,Qc,t,dt,TLevel,ShortF,TPrint,NumKD,
     !                  SMean,ChemS,CumCh0,CumCh1,CumChR,cCumA,cCumT,
     !                  lWat)

      integer TLevel
      logical ShortF,lWat
      dimension Qc(NumNP),Kode(NumNP),ChemS(NumKD),SMean(NumKD)

      do 11 i=1,NumKD
        SMean(i)=0.
11    continue
      do 12 i=1,NumNP
        j=iabs(Kode(i))
        if(j.ne.0) then
          SMean(j)=Smean(j)-Qc(i)
        end if
12    continue
      cCumA=abs(CumCh0)+abs(CumCh1)+abs(CumChR)
      cCumT=CumCh0+CumCh1+CumChR
      do 13 j=1,NumKD
        ChemS(j)=ChemS(j)+SMean(j)*dt
        cCumT=cCumT+ChemS(j)
        cCumA=cCumA+abs(ChemS(j))
13    continue
      if(TLevel.eq.1) write(74,110)
      if(.not.ShortF.or.abs(TPrint-t).lt.0.001*dt)
     !  write(74,120) t,CumCh0,CumCh1,CumChR,(ChemS(j),j=1,NumKD),
     !                (SMean(i),i=1,NumKD)
      if(.not.lWat) write(*,130) t,TLevel,CumCh0,CumCh1,ChemS(1)

110   format(' All solute fluxes (SMean) and cumulative solute fluxes',
     !' (ChemS) are positive out of the region'//
     !'     Time     CumCh0     CumCh1     CumChR   ',20('-'),'  ChemS(i
     !),i=1,NumKD  ',22('-'),'  ',21('-'),'  SMean(j),j=1,NumKD ',22('-'
     !)/'      [T]    [VM/L3]    [VM/L3]    [VM/L3]',31(' '),'[VM/L3]',5
     !9(' '),'[VM/T/L3]'/)
120   format(f10.2,15e11.3)
130   format(f14.4,i6,1x,3e11.3,2x,2e11.3)
      return
      end

***********************************************************************

      subroutine ObsNod(t,NumNP,NObs,NObsD,Node,hNew,ThNew,Conc)

      dimension Node(NObsD),hNew(NumNP),ThNew(NumNP),Conc(NumNP)

      write(92,110) t,(hNew(Node(i)),ThNew(Node(i)),Conc(Node(i)),
     !              i=1,NObs)

110   format(f11.3,5(f11.3,f9.4,e11.3))
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

***********************************************************************

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

      subroutine cOut(NumNP,Conc,x,y,t,IJ)

      dimension Conc(NumNP),x(NumNP),y(NumNP)

      write(83,110) t
      L1=(IJ-1)/10+1
      do 12 n=1,NumNP,IJ
        do 11 L=1,L1
          m=n+(L-1)*10
          k=m+9
          if(L.eq.L1) k=n+IJ-1
          write(83,120) m,x(m),y(m),(Conc(j),j=m,k)
11      continue
12    continue

110   format(//' Time  ***',f12.4,' ***'//
     !         '    n    x(n)   z(n)      Conc(n)   Conc(n+1)  ...'/)
120   format(i5,2f8.1,10e11.3)
      return
      end

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
