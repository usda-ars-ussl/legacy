* Source file OUTPUT3.FOR ||||||||||||||||||||||||||||||||||||||||||||||

      subroutine TLInf(NumNP,NumBP,Kode,Q,hNew,CumQ,Width,SWidth,KXB,t,
     !                 dt,TLevel,ShortF,TPrint,Iter,ItCum,rTop,rRoot,
     !                 vMeanR,hMeanT,hMeanR,hMeanG,AtmInF,SinkF,CumQrT,
     !                 CumQrR,CumQvR,NumKD,hMean,vMean,lWat,lChem,rLen,
     !                 Peclet,Courant,wCumT,wCumA,PeCrMx)

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
            write(70,200) TLevel,t,dt,Iter,ItCum,Peclet,Courant,PeCrMx
          else
            write(70,210) TLevel,t,dt,Peclet,Courant,PeCrMx
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
     !'Courant   PeCrMax'/)
130   format(//' TLevel   Time         dt        Peclet   Courant   PeCr
     !Max'/)
140   format(' All fluxes (v) are positive out of the region'//
     !'     Time       rAtm       rRoot      vAtm       vRoot      ',
     !'vKode3     vKode1     vSeep      vKode5     vKode6 ...'/
     !'      [T]      [L/T]       [L/T]     [L/T]       [L/T]      ',
     !'[L3/T]     [L3/T]    [L3/T]      [L3/T]     [L3/T]'/)
150   format(//
     !'     Time         hAtm       hRoot     hKode3     hKode1 ',
     !'     hSeep     hKode5     hKode6 ... '/
     !'      [T]          [L]        [L]        [L]        [L]  ',
     !'      [L]        [L]        [L]'/)
160   format(
     !' All cumulative fluxes (CumQ) are positive out of the region'//
     !'     Time     CumQAP      CumQRP     CumQA      CumQR     CumQ3',
     !'       CumQ1      CumQS      CumQ5       CumQ6 ....'/
     !'      [T]      [L3]        [L3]       [L3]       [L3]      [L3]',
     !'        [L3]       [L3]       [L3]        [L3]'/)
170   format(f12.4,9e11.3)
180   format(f12.4,9f11.1)
190   format(f12.4,9e11.3)
200   format(i5,2e12.3,i5,i6,3f10.3)
210   format(i5,2e12.3,3f10.3)
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
     !'      [T]        [L3]       [L3]       [L3]       [L3]  ',
     !'      [L3]          [L]        [L]        [L] '/)
120   format(f12.4,5e11.3,3f11.1,i8)
      return
      end

************************************************************************

      subroutine SubReg(NumEl,NumElD,NumNP,NMat,hNew,ThO,ThN,x,y,z,
     !                  MatNum,LayNum,KX,t,dt,NLay,PLevel,lWat,lChem,
     !                  Conc,ChPar,wCumA,wCumT,cCumA,cCumT,wVolI,cVolI,
     !                  WatIn,SolIn)

      logical lWat,lChem
      integer PLevel,e
      dimension hNew(NumNP),x(NumNP),y(NumNP),MatNum(NumNP),Conc(NumNP),
     !          KX(NumElD,9),ChPar(10,NMat),LayNum(NumEl),ThO(NumNP),
     !          ThN(NumNP),WatIn(NumEl),SolIn(NumEl),Vol(10),hMean(10),
     !          SubVol(10),SubCha(10),cMean(10),ConSub(10),z(NumNP)

      VTot=0.
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
        Vol(i)=0.
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
        NUS=5
        if(KX(e,7).eq.0) NUS=3
        if(KX(e,5).eq.0) NUS=1
        do 12 n=1,NUS
          call Elem(NUS,NumElD,e,n,KX,i,j,k,l)
          Mi=MatNum(i)
          Mj=MatNum(j)
          Mk=MatNum(k)
          Ml=MatNum(l)
          Bl=+(y(j)-y(i))*(z(k)-z(i))-(y(k)-y(i))*(z(j)-z(i))
          Cl=-(x(j)-x(i))*(z(k)-z(i))+(x(k)-x(i))*(z(j)-z(i))
          Dl=+(x(j)-x(i))*(y(k)-y(i))-(x(k)-x(i))*(y(j)-y(i))
          Det=(x(l)-x(i))*Bl+(y(l)-y(i))*Cl+(z(l)-z(i))*Dl
          VE=Det/6.
          Vol(Lay)=Vol(Lay)+VE
          if(lWat.or.PLevel.le.1) then
            hE=(hNew(i)+hNew(j)+hNew(k)+hNew(l))/4.
            VNewE=VE*(thN(i)+thN(j)+thN(k)+thN(l))/4.
            VOldE=VE*(thO(i)+thO(j)+thO(k)+thO(l))/4.
            Volume=Volume+VNewE
            wEl=wEl+VNewE
            Change=Change+(VNewE-VOldE)/dt
            SubVol(Lay)=SubVol(Lay)+VNewE
            SubCha(Lay)=SubCha(Lay)+(VNewE-VOldE)/dt
            hTot=hTot+hE*VE
            hMean(Lay)=hMean(Lay)+hE*VE
          end if
          if(lChem) then
            cE=(Conc(i)+Conc(j)+Conc(k)+Conc(l))/4.
            cNewE=VE*((thN(i)+ChPar(1,Mi)*ChPar(5,Mi))*Conc(i)+
     !                (thN(j)+ChPar(1,Mj)*ChPar(5,Mj))*Conc(j)+
     !                (thN(k)+ChPar(1,Mk)*ChPar(5,Mk))*Conc(k)+
     !                (thN(l)+ChPar(1,Ml)*ChPar(5,Ml))*Conc(l))/4.
            ConVol=ConVol+cNewE
            cEl=cEl+cNewE
            ConSub(Lay)=ConSub(Lay)+cNewE
            cTot=cTot+cE*VE
            cMean(Lay)=cMean(Lay)+cE*VE
          end if
          if(n.eq.NUS) then
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
        if(lWat.or.PLevel.le.1) hMean(Lay)=hMean(Lay)/Vol(Lay)
        if(lChem) cMean(Lay)=cMean(Lay)/Vol(Lay)
        VTot=VTot+Vol(Lay)
14    continue
      if(lWat.or.PLevel.le.1) hTot=hTot/VTot
      if(lChem) cTot=cTot/VTot
      if(PLevel.eq.0) write(80,110)
      write(80,120) t,     (       i,i=1,NLay)
      write(80,130) VTot,  (  Vol(i),i=1,NLay)
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
130   format(' Volume  [L3]   ',e11.3,10e11.3)
140   format(' Water   [L3]   ',e11.3,10e11.3)
150   format(' InFlow  [L3/T] ',e11.3,10e11.3)
160   format(' hMean   [L]    ',e11.3,10f11.1)
170   format(' ConcVol [M]    ',e11.3,10e11.3)
180   format(' cMean   [M/L3] ',e11.3,10e11.3)
190   format(' WatBalT [L3]   ',e11.3)
200   format(' WatBalR [%]    ',f11.3)
210   format(' CncBalT [M]    ',e11.3)
220   format(' CncBalR [%]    ',f11.3)
      return
      end

************************************************************************

      subroutine BouOut(NumNP,NumBP,t,hNew,theta,Q,Width,KXB,Kode,x,y,z,
     !                  Qc,Conc)

      dimension hNew(NumNP),Q(NumNP),Width(NumBP),theta(NumNP),
     !          KXB(NumBP),Kode(NumNP),x(NumNP),y(NumNP),z(NumNP),
     !          Qc(NumNP),Conc(NumNP)

      write(79,110) t
      ii=0
      do 12 i=1,NumNP
        if(Kode(i).ne.0) then
          do 11 j=1,NumBP
            n=KXB(j)
            if(n.eq.i) then
              ii=ii+1
              v=-Q(i)/Width(j)
              vc=-Qc(i)/Width(j)
              write(79,120) ii,i,x(i),y(i),z(i),Kode(i),Q(i),v,hNew(i),
     !                      theta(i),Qc(i),vc,Conc(i)
              goto 12
            end if
11        continue
          ii=ii+1
          write(79,130) ii,i,x(i),y(i),z(i),Kode(i),Q(i),hNew(i),
     !                  theta(i),Qc(i),Conc(i)
        end if
12    continue

110   format(//' Time:',f12.4//
     !'    i    n    x      y      z    Code     Q          v       ',
     !'    h       th       Qc         vc       Conc'/
     !'                                       [L3/T]      [L/T]     ',
     !'   [L]     [-]     [M/T]     [M/L2/T]   [M/L3]'/)
120   format(2i5,3f7.1,i5,2e11.3,f11.1,f7.3,2e11.3,e10.3)
130   format(2i5,3f7.1,i5,e11.3,11x,f11.1,f7.3,e11.3,11x,e10.3)
      return
      end

************************************************************************

      subroutine SolInf(K,NumNP,Kode,Qc,t,dt,TLevel,ShortF,TPrint,
     !                  AtmInF,NumKD,SMean,ChemS,CumCh0,CumCh1,CumChR,
     !                  cCumA,cCumT,lWat,PeCrMx)

      integer TLevel
      logical ShortF,AtmInF,lWat
      dimension Qc(NumNP),Kode(NumNP),ChemS(NumKD),SMean(NumKD)

      if(K.eq.1) then
        do 11 i=1,NumKD
          SMean(i)=0.
11      continue
        do 12 i=1,NumNP
          j=iabs(Kode(i))
          if(j.ne.0) then
            SMean(j)=Smean(j)-Qc(i)
          end if
12      continue
        cCumA=abs(CumCh0)+abs(CumCh1)+abs(CumChR)
        cCumT=CumCh0+CumCh1+CumChR
        do 13 j=1,NumKD
          ChemS(j)=ChemS(j)+SMean(j)*dt
          cCumT=cCumT+ChemS(j)
          cCumA=cCumA+abs(ChemS(j))
13      continue
        if(TLevel.eq.1) write(74,110)
        if(.not.AtmInf.and.(.not.ShortF.or.abs(TPrint-t).lt.0.001*dt))
     !    write(74,120) t,CumCh0,CumCh1,CumChR,(ChemS(j),j=1,NumKD),
     !                  (SMean(i),i=1,NumKD)
      else
        write(74,120) t,CumCh0,CumCh1,CumChR,(ChemS(j),j=1,NumKD),
     !                (SMean(i),i=1,NumKD)
      end if
      if(.not.lWat) write(*,130) t,TLevel,CumCh0,CumCh1,ChemS(1),PeCrMx

110   format(' All solute fluxes (SMean) and cumulative solute fluxes',
     !' (ChemS) are positive out of the region'//
     !'     Time     CumCh0     CumCh1     CumChR   ',20('-'),'  ChemS(i
     !),i=1,NumKD  ',22('-'),'  ',21('-'),'  SMean(j),j=1,NumKD ',22('-'
     !)/'      [T]      [M]        [M]        [M]',35(' '),'[M]',63(' ')
     !,' [M/T]'/)
120   format(f10.2,15e11.3)
130   format(f14.4,i6,1x,3e11.3,2x,2e11.3)
      return
      end

***********************************************************************

      subroutine ObsNod(t,NumNP,NObs,Node,hNew,ThNew,Conc)

      dimension Node(NObs),hNew(NumNP),ThNew(NumNP),Conc(NumNP)

      write(92,110) t,(hNew(Node(i)),ThNew(Node(i)),Conc(Node(i)),
     !              i=1,NObs)

110   format(f11.3,5(f11.3,f9.4,e11.3))
      return
      end

************************************************************************

      subroutine hOut(hNew,x,y,z,NumNP,t,IJ)

      dimension hNew(NumNP),x(NumNP),y(NumNP),z(NumNP)

      write(75,110) t
      L1=(IJ-1)/10+1
      do 12 n=1,NumNP,IJ
        do 11 L=1,L1
          m=n+(L-1)*10
          k=m+9
          if(L.eq.L1) k=n+IJ-1
          write(75,120) m,x(m),y(m),z(m),(hNew(j),j=m,k)
11      continue
12    continue

110   format(//' Time  ***',f12.4,' ***'//
     !       '    n    x(n)   y(n)   z(n)       h(n)      h(n+1) ...'/)
120   format(i5,3f8.1,10f10.1)
      return
      end

***********************************************************************

      subroutine QOut(Q,x,y,z,NumNP,t,IJ)

      dimension Q(NumNP),x(NumNP),y(NumNP),z(NumNP)

      write(73,110) t
      L1=(IJ-1)/10+1
      do 12 n=1,NumNP,IJ
        do 11 L=1,L1
          m=n+(L-1)*10
          k=m+9
          if(L.eq.L1) k=n+IJ-1
          write(73,120) m,x(m),y(m),z(m),(Q(j),j=m,k)
11      continue
12    continue

110   format(//' Time  ***',f12.4,' ***'//
     !        '    n    x(n)   y(n)   z(n)       Q(n)      Q(n+1) ...'/)
120   format(i5,3f8.1,10e11.3)
      return
      end

************************************************************************

      subroutine thOut(theta,x,y,z,NumNP,t,IJ)

      dimension theta(NumNP),x(NumNP),y(NumNP),z(NumNP)

      write(76,110) t
      L1=(IJ-1)/16+1
      do 12 n=1,NumNP,IJ
        do 11 L=1,L1
          m=n+(L-1)*16
          k=m+15
          if(L.eq.L1) k=n+IJ-1
          write(76,120) m,x(m),y(m),z(m),(theta(j),j=m,k)
11      continue
12    continue

110   format(//' Time  ***',f12.4,' ***'//
     !        '    n    x(n)   y(n)   z(n)      th(n)     th(n+1) ...'/)
120   format(i5,3f8.1,16f6.3)
      return
      end

************************************************************************

      subroutine FlxOut(Vx,Vy,Vz,x,y,z,NumNP,t,IJ)

      dimension x(NumNP),y(NumNP),z(NumNP),Vx(NumNP),Vy(NumNP),Vz(NumNP)

      write(81,110) t
      write(82,120) t
      write(84,130) t
      L1=(IJ-1)/10+1
      do 12 n=1,NumNP,IJ
        do 11 L=1,L1
          m=n+(L-1)*10
          k=m+9
          if(L.eq.L1) k=n+IJ-1
          write(81,140) m,x(m),y(m),z(m),(Vz(j),j=m,k)
          write(82,140) m,x(m),y(m),z(m),(Vx(j),j=m,k)
          write(84,140) m,x(m),y(m),z(m),(Vy(j),j=m,k)
11      continue
12    continue

110   format(//' Time  ***',f12.4,' ***'//
     !         '    n    x(n)  y(n)  z(n)     vz(n)     vz(n+1) ...'/)
120   format(//' Time  ***',f12.4,' ***'//
     !         '    n    x(n)  y(n)  z(n)     vx(n)     vx(n+1) ...'/)
130   format(//' Time  ***',f12.4,' ***'//
     !         '    n    x(n)  y(n)  z(n)     vy(n)     vy(n+1) ...'/)
140   format(i5,3f8.1,10e10.2)
      return
      end

************************************************************************

      subroutine cOut(NumNP,Conc,x,y,z,t,IJ)

      dimension Conc(NumNP),x(NumNP),y(NumNP),z(NumNP)

      write(83,110) t
      L1=(IJ-1)/10+1
      do 12 n=1,NumNP,IJ
        do 11 L=1,L1
          m=n+(L-1)*10
          k=m+9
          if(L.eq.L1) k=n+IJ-1
          write(83,120) m,x(m),y(m),z(m),(Conc(j),j=m,k)
11      continue
12    continue

110   format(//' Time  ***',f12.4,' ***'//
     !     '    n    x(n)   y(n)   z(n)      Conc(n)   Conc(n+1)  ...'/)
120   format(i5,3f8.1,10e11.3)
      return
      end

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||