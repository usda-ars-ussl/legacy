* Source file INPUT2.FOR |||||||||||||||||||||||||||||||||||||||||||||||

      subroutine BasInf(KAT,MaxIt,TolTh,TolH,lWat,lChem,AtmInF,ShortF,
     !                  SeepF,CheckF,FluxF,FreeD,DrainF,lTemp,lWDep,
     !                  lEquil)

      character*100 Hed
      character*5  LUnit,TUnit,MUnit
      logical lWat,lChem,CheckF,AtmInF,ShortF,SeepF,FluxF,FreeD,DrainF,
     !        lTemp,lWDep,lEquil
      dimension IU(13)
      data IU /50,70,71,72,73,75,76,77,78,79,80,81,82/

      read(30,*)
      read(30,*)
      read(30,*) Hed
      read(30,*)
      read(30,*) LUnit,TUnit,MUnit
      read(30,*)
      read(30,*) KAT
      read(30,*)
      read(30,*) MaxIt,TolTh,TolH
      read(30,*)
      read(30,*) lWat,lChem,CheckF,ShortF,FluxF,AtmInF,SeepF,DrainF,
     !           FreeD,lTemp,lWDep,lEquil
      do 11 i=1,11
        write(IU(i),*) Hed
        write(IU(i),*)
        write(IU(i),*)'Program CHAIN_2D'
        call getdat(ii,imonth,iday)
        call gettim(ihours,mins,isecs,ii)
        write(iU(i),100) iday,imonth,ihours,mins,isecs
        if(AtmInF) then
          write(IU(i),*)'Time dependent boundary conditions'
        else
          write(IU(i),*)'Time independent boundary conditions'
        end if
        if(KAT.eq.0) write(IU(i),110)
        if(KAT.eq.1) write(IU(i),120)
        if(KAT.eq.2) write(IU(i),130)
        write(IU(i),*) 'Units: L = ',LUnit,', T = ',TUnit,', M = ',MUnit
11    continue
      write(*,*) '-----------------------------------------------------'
      write(*,*) '|                                                   |'
      write(*,*) '|                    CHAIN_2D                       |'
      write(*,*) '|                                                   |'
      write(*,*) '|    Code for simulating two-dimensional variably   |'
      write(*,*) '|     saturated water flow, heat transport, and     |'
      write(*,*) '|    transport of solutes involved in sequential    |'
      write(*,*) '|          first-order decay reactions              |'
      write(*,*) '|                                                   |'
      write(*,*) '|                   version 1.1                     |'
      write(*,*) '|                                                   |'
      write(*,*) '|            Last modified: October, 1994           |'
      write(*,*) '|                                                   |'
      write(*,*) '-----------------------------------------------------'
      write(*,*)
      write(*,*) Hed
      if(KAT.eq.0) write(*,110)
      if(KAT.eq.1) write(*,120)
      if(KAT.eq.2) write(*,130)
      write(50,140) MaxIt,TolTh,TolH
      if(lWDep) write(50,150)

100   format(' Date: ',i3,'.',i2,'.','    Time: ',i3,':',i2,':',i2)
110   format(' Horizontal plane flow, V = L*L')
120   format(' Axisymmetric flow, V = L*L*L')
130   format(' Vertical plane flow, V = L*L')
140   format(/' Max. number of iterations           ',i4/
     !        ' Absolute water content tolerance [-]',f8.5/
     !        ' Absolute pressure head tolerance [L]',f8.5/)
150   format(/' Temperature dependence of hydraulic properties is consid
     !ered'/)
      return
      end

************************************************************************

      subroutine MatIn(NMatD,NMat,NLay,Par,hTab1,hTabN)

      real K
      dimension Par(10,NMatD),Qe(10)
      data Qe /1.,.99,.90,.85,.75,.65,.50,.35,.20,.10/

      write(*,*) 'reading material information'
      Imax=10
      read(30,*)
      read(30,*)
      read(30,*) NMat,NLay,hTab1,hTabN,NPar
      if(NMat.gt.NMatD) then
        write(*,*) 'Dimension in NMatD is exceeded'
        stop
      end if
      hTab1=-amin1(abs(hTab1),abs(hTabN))
      hTabN=-amax1(abs(hTab1),abs(hTabN))
      read (30,*)
      write(50,110)
      do 11 M=1,NMat
        read (30,*)     (Par(i,M),i=1,NPar)
        write(50,120) M,(Par(i,M),i=1,NPar)
11    continue
      write(50,130)
      do 13 M=1,NMat
        write(50,*)
        do 12 i=1,Imax
          h=FH(Qe(i),Par(1,M))
          K=FK(h,Par(1,M))
          C=FC(h,Par(1,M))
          Q=FQ(h,Par(1,M))
          write(50,140) M,Qe(i),Q,h,C,K
12      continue
13    continue

110   format(/' MatNum         thR    thS    tha    thm      alpha      
     !   n          Ks          Kk          thk'/)
120   format(i5,8x,4f7.3,16e12.3)
130   format(//' MatNum          Qe     Q        h         C         K') 
140   format(i5,8x,2f7.3,f10.3,e10.2,e12.3)
      return
      end

************************************************************************

      subroutine GenMat(NTab,NTabD,NMat,thR,hSat,Par,hTab,ConTab,CapTab,
     !                  ConSat,TheTab,thSat)

      dimension Par(10,NMat),ConTab(NTabD,NMat),CapTab(NTabD,NMat),
     !          TheTab(NTabD,NMat),hTab(NTab),ConSat(NMat),hSat(NMat),
     !          thR(NMat),thSat(NMat)

      write(*,*) 'generating materials'
      hTab1=hTab(1)
      hTabN=hTab(NTab)
      dlh=(alog10(-hTabN)-alog10(-hTab1))/(NTab-1)
      do 11 i=1,NTab
        alh=alog10(-hTab1)+(i-1)*dlh
        hTab(i)=-10**alh
11    continue
      do 13 M=1,NMat
        Hr       =FH(0.0,Par(1,M))
        hSat(M)  =FH(1.0,Par(1,M))
        ConSat(M)=FK(0.0,Par(1,M))
        thR(M)   =FQ(Hr ,Par(1,M))
        thSat(M) =FQ(0.0,Par(1,M))
        do 12 i=1,NTab
          ConTab(i,M)=FK(hTab(i),Par(1,M))
          CapTab(i,M)=FC(hTab(i),Par(1,M))
          TheTab(i,M)=FQ(hTab(i),Par(1,M))
12      continue
13    continue
      return
      end

************************************************************************

      subroutine TmIn(tInit,tMax,tAtm,tOld,dt,dtMax,dMul,dMul2,dtMin,
     !                TPrint,t,dtOpt,dtOld,AtmInF)

      logical AtmInF
      dimension TPrint(50)

      write(*,*) 'reading time information'
      read(30,*)
      read(30,*)
      read(30,*) dt,dtMin,dtMax,dMul,dMul2,MPL
      read(30,*)
      read(30,*) (TPrint(i),i=1,MPL)
      dtOpt=dt
      dtOld=dt
      if(.not.AtmInF) then
        tMax=TPrint(MPL)
        tAtm=tMax
      end if
      TPrint(MPL+1)=tMax
      tOld=tInit
      t=tInit+dt
      return
      end

************************************************************************

      subroutine SeepIn(NSeepD,NumSPD,NSeep,NSP,NP)

      dimension NSP(NSeepD),NP(NSeepD,NumSPD)

      write(*,*) 'reading seepage face information'
      read(30,*)
      read(30,*)
      read(30,*) NSeep
      if(NSeep.gt.NSeepD) then
        write(*,*) 'Dimension in NSeepD is exceeded'
        stop
      end if
      read(30,*)
      read(30,*) (NSP(i),i=1,NSeep)
      do 11 i=1,NSeep
        if(NSP(i).gt.NumSPD) then
          write(*,*) 'Dimension in NumSPD is exceeded'
          stop
        end if
11    continue
      read(30,*)
      do 12 i=1,NSeep
        read(30,*) (NP(i,j),j=1,NSP(i))
12    continue
      return
      end

************************************************************************

      subroutine DrainIn(NDr,NDrD,NElDrD,NumEl,ND,NED,KElDr,EfDim,
     !                   ConAxx,ConAxz,ConAzz)

      integer e
      dimension ND(NDrD),NED(NDrD),KElDr(NDrD,NElDrD),EfDim(2,NDrD),
     !          ConAxx(NumEl),ConAzz(NumEl),ConAxz(NumEl)

      write(*,*) 'reading drainage information'
      read(30,*)
      read(30,*)
      read(30,*) NDr,DrCorr
      if(NDr.gt.NDrD) then
        write(*,*) 'Dimension in NDrD is exceeded'
        stop
      end if
      read(30,*)
      read(30,*) (ND(i),i=1,NDr)
      read(30,*)
      read(30,*) (NED(i),i=1,NDr)
      do 11 i=1,NDr
        if(NED(i).gt.NElDrD) then
          write(*,*) 'Dimension in NElDrD is exceeded'
          stop
        end if
11    continue
      read(30,*)
      do 12 i=1,NDr
        read(30,*) (EfDim(i,j),j=1,2)
12    continue
      read(30,*)
      do 13 i=1,NDr
        read(30,*) (KElDr(i,j),j=1,NED(i))
13    continue
      do 15 i=1,Ndr
        rho=EfDim(i,2)/EfDim(i,1)
        A=(1.+0.405*rho**(-4))/(1.-0.405*rho**(-4))
        B=(1.+0.163*rho**(-8))/(1.-0.163*rho**(-8))
        C=(1.+0.067*rho**(-12))/(1.-0.067*rho**(-12))
        Red=376.7/(138.*alog10(rho)+6.48-2.34*A-0.48*B-0.12*C)/DrCorr
        do 14 j=1,NED(i)
          e=KElDr(i,j)
          ConAxx(e)=ConAxx(e)*Red
          ConAxz(e)=ConAxz(e)*Red
          ConAzz(e)=ConAzz(e)*Red
14      continue
15    continue
      return
      end

************************************************************************

      subroutine NodInf(NumNP,NumEl,IJ,NumBP,NumNPD,NumElD,NumBPD,NumKD,
     !                  NS,NSD,NObs,NObsD,Kode,Q,Conc,hNew,hOld,hTemp,x,
     !                  y,MatNum,Beta,Axz,Bxz,Dxz,CheckF,TempN,TempO,
     !                  Sorb,lEquil)

      logical CheckF,lEquil
      dimension Kode(NumNPD),Q(NumNPD),hOld(NumNPD),x(NumNPD),y(NumNPD),
     !          hNew(NumNPD),hTemp(NumNPD),MatNum(NumNPD),Beta(NumNPD),
     !          Axz(NumNPD),Bxz(NumNPD),Dxz(NumNPD),Conc(NSD,NumNPD),
     !          TempN(NumNPD),TempO(NumNPD),Sorb(NSD,NumNPD),DConc(8),
     !          DSorb(8)

      write(*,*) 'reading nodal information'
      read(32,*)
      read(32,*)
      read(32,*) NumNP,NumEl,IJ,NumBP,NS,NObs
      if(NumNP.gt.NumNPD) then
        write(*,*) 'Dimension in NumNPD is exceeded'
        stop
      else if(NumEl.gt.NumElD) then
        write(*,*) 'Dimension in NumElD is exceeded'
        stop
      else if(NumBP.gt.NumBPD) then
        write(*,*) 'Dimension in NumBPD is exceeded'
        stop
      else if(NObs.gt.NObsD) then
        write(*,*) 'Dimension in NObsD is exceeded'
        stop
      end if
      read(32,*)
      NPR=0
      k=0

11    k=k+1
      if(lEquil) then
        read(32,*) n,Kode(n),x(n),y(n),hOld(n),Q(n),MatNum(n),Beta(n),
     !             Axz(n),Bxz(n),Dxz(n),TempO(n),(Conc(j,n),j=1,NS)
      else
        read(32,*) n,Kode(n),x(n),y(n),hOld(n),Q(n),MatNum(n),Beta(n),
     !             Axz(n),Bxz(n),Dxz(n),TempO(n),(Conc(j,n),j=1,NS),
     !             (Sorb(j,n),j=1,NS)
      end if
      if(Kode(n).gt.NumKD) then
        write(*,*) 'Dimension in NumKD is exceeded'
        stop
      end if
      if(n-k) 12,17,13
12    write(*,130) n
      stop
13    Deno=n-k+1
      DX=(x(n)-x(NPR))/Deno
      DY=(y(n)-y(NPR))/Deno
      DP=(hOld(n)-hOld(NPR))/Deno
      do 14 j=1,NS
        DConc(j)=(Conc(j,n)-Conc(j,NPR))/Deno
        if(.not.lEquil) DSorb(j)=(Sorb(j,n)-Sorb(j,NPR))/dx
14    continue
      DBeta=(Beta(n)-Beta(NPR))/Deno
      DA=(Axz(n)-Axz(NPR))/Deno
      DB=(Bxz(n)-Bxz(NPR))/Deno
      DD=(Dxz(n)-Dxz(NPR))/Deno
      DTemp=(TempO(n)-TempO(NPR))/Deno
15    x(k)=x(k-1)+DX
      y(k)=y(k-1)+DY
      hOld(k)=hOld(k-1)+DP
      do 16 j=1,NS
        Conc(j,k)=Conc(j,k-1)+DConc(j)
        if(.not.lEquil) Sorb(j,k)=Sorb(j,k-1)-DSorb(j)
16    continue
      Beta(k)=Beta(k-1)+DBeta
      Axz(k)=Axz(k-1)+DA
      Bxz(k)=Bxz(k-1)+DB
      Dxz(k)=Dxz(k-1)+DD
      TempO(k)=TempO(k-1)+DTemp
      MatNum(k)=MatNum(k-1)
      Kode(k)=Kode(k-1)
      Q(k)=Q(k-1)
      k=k+1
      if(k.lt.N) goto 15
17    NPR=N
      if(k.lt.NumNP) goto 11

      do 18 n=1,NumNP
        hNew(n)=hOld(n)
        hTemp(n)=hOld(n)
        TempN(n)=TempO(n)
18    continue
      if(CheckF) then
        write(50,110)
        do 19 n=1,NumNP
          if(lEquil) then
            write(50,120) n,Kode(n),x(n),y(n),hOld(n),Q(n),Conc(1,n),
     !                    TempO(n),MatNum(n),Beta(n),Axz(n),Bxz(n),
     !                    Dxz(n),(Conc(j,n),j=2,NS)
          else
            write(50,120) n,Kode(n),x(n),y(n),hOld(n),Q(n),Conc(1,n),
     !                    TempO(n),MatNum(n),Beta(n),Axz(n),Bxz(n),
     !                    Dxz(n),(Conc(j,n),j=2,NS),(Sorb(j,n),j=1,NS)
          end if
19      continue
      end if
      if(lEquil) then
        do 21 j=1,NS
          do 20 i=1,NumNP
            Sorb(j,i)=0.
20        continue
21      continue
      end if

110   format(////' Nodal Point Information'//'Node no   Kode',
     !       5x,'x,r',7x,'y,z',7x,'h',9x,'Q',8x,'Conc',6x,'Temp'/)
120   format(2i7,3f10.3,2e11.3,f8.2,i5,4f7.3,6e11.3)
130   format(' ERROR in NodInf at Node =',i5)
      return
      end

************************************************************************

      subroutine ElemIn(NumEl,NumElD,NumNP,KX,LayNum,ConAxx,ConAzz,
     !                  ConAxz,CheckF,ListNE,IJ,MBand,MBandD,lChem,lOrt)

      logical CheckF,lChem,lConst,lOrt
      integer e
      dimension KX(NumElD,4),ConAxx(NumEl),ConAzz(NumEl),ConAxz(NumEl),
     !          LayNum(NumEl),ListNE(NumNP)

      write(*,*) 'reading element information'
      Num=0
      read(32,*)
      read(32,*)
      do 14 e=1,NumEl
        IF (Num-e) 11,14,12
11      read(32,*) Num,(KX(Num,i),i=1,4),ConAxz(Num),ConAxx(Num),
     !             ConAzz(Num),LayNum(Num)
        if(KX(Num,4).eq.0) KX(Num,4)=KX(Num,3)
        if(Num.eq.e) goto 14
12      do 13 i=1,4
          KX(e,i)=KX(e-1,i)+1
13      continue
        ConAxx(e)=ConAxx(e-1)
        ConAzz(e)=ConAzz(e-1)
        ConAxz(e)=ConAxz(e-1)
        LayNum(e)=LayNum(e-1)
14    continue
      AA=3.141592654/180.
      do 15 e=1,NumEl
        Ang=AA*ConAxz(e)
        CAxx=ConAxx(e)
        CAzz=ConAzz(e)
        ConAxx(e)=CAxx*cos(Ang)*cos(Ang) + CAzz*sin(Ang)*sin(Ang)
        ConAzz(e)=CAxx*sin(Ang)*sin(Ang) + CAzz*cos(Ang)*cos(Ang)
        ConAxz(e)=(CAxx-CAzz)*sin(Ang)*cos(Ang)
15    continue
      if(CheckF) then
        write(50,110)
        do 16 e=1,NumEl
          write(50,120) e,(KX(e,i),i=1,4),ConAxz(e),ConAxx(e),ConAzz(e),
     !                  LayNum(e)
16      continue
      end if

      do 17 i=1,NumNP
        ListNE(i)=0
17    continue
      do 19 e=1,NumEl
        NCorn=4
        if(KX(e,3).eq.KX(e,4)) NCorn=3
        do 18 n=1,NCorn-2
          i=KX(e,1)
          j=KX(e,n+1)
          k=KX(e,n+2)
          ListNE(i)=ListNE(i)+1
          ListNE(j)=ListNE(j)+1
          ListNE(k)=ListNE(k)+1
18      continue
19    continue

      lOrt=.false.
      lConst=.true.
      MBand=1
      do 21 e=1,NumEl
        NUS=4
        if(KX(e,3).eq.KX(e,4)) NUS=3
        do 20 kk=1,NUS-2
          MB=1
          i=KX(e,1)
          j=KX(e,kk+1)
          k=KX(e,kk+2)
          if(abs(i-j).gt.MB) MB=abs(i-j)
          if(abs(i-k).gt.MB) MB=abs(i-k)
          if(abs(j-k).gt.MB) MB=abs(j-k)
          if(MB.gt.MBand) MBand=MB
          if(e.eq.1.and.kk.eq.1) then
            MB1=MB
          else
            if(MB1.ne.MB) lConst=.false.
          end if
20      continue
21    continue
      MBand=MBand+1
      if(MBand.gt.MBandD.or.(lChem.and.2*MBand-1.gt.MBandD)) lOrt=.true.
      if(.not.lConst) IJ=NumNP
      if(MBand.gt.10.or.NumNP.gt.200) lOrt=.true.

110   format (////' Element Information'//' Element    C O R N E R    N 
     !O D E S   ConAxz    ConAxx    ConAzz  LayNum'/)
120   format (i6,i9,3i6,e14.3,2f8.3,i5)
      return
      end

************************************************************************

      subroutine GeomIn(NumKD,NumNP,NumBP,NObs,NObsD,NS,SWidth,Width,
     !                  Kode,KXB,rLen,Node)

      character*30 Text1,Text2,Text3
      dimension KXB(NumBP),Width(NumBP),SWidth(NumKD),Kode(NumNP),
     !          Node(NObsD)

      write(*,*) 'reading geometric information'
      read(32,*)
      read(32,*)
      read(32,*) (KXB(i),i=1,NumBP)
      read(32,*)
      read(32,*) (Width(i),i=1,NumBP)
      read(32,*)
      read(32,*) rLen
      
      do 11 i=1,NumKD
        SWidth(i)=0.
11    continue
      do 12 i=1,NumBP
        n=KXB(i)
        j=iabs(Kode(n))
        if(j.eq.0) goto 12
        SWidth(j)=SWidth(j)+Width(i)
12    continue

      if(NObs.gt.0) then
        read(32,*)
        read(32,*) (Node(i),i=1,NObs)
        Text1='   hNew    theta    Temp    '
        Text2='   Conc       Sorb    '
        Text3='Node('
        if(NS.eq.0) write(92,110) (Text3,Node(j),j=1,NObs)
        if(NS.eq.1) write(92,120) (Text3,Node(j),j=1,NObs)
        if(NS.eq.2) write(92,130) (Text3,Node(j),j=1,NObs)
        if(NS.eq.3) write(92,140) (Text3,Node(j),j=1,NObs)
        if(NS.eq.4) write(92,150) (Text3,Node(j),j=1,NObs)
        if(NS.eq.5) write(92,160) (Text3,Node(j),j=1,NObs)
        if(NS.eq.6) write(92,170) (Text3,Node(j),j=1,NObs)
        if(NS.eq.0) write(92,180) (Text1,               i=1,NObs)
        if(NS.eq.1) write(92,190) (Text1,(Text2,j=1,NS),i=1,NObs)
        if(NS.eq.2) write(92,200) (Text1,(Text2,j=1,NS),i=1,NObs)
        if(NS.eq.3) write(92,210) (Text1,(Text2,j=1,NS),i=1,NObs)
        if(NS.eq.4) write(92,220) (Text1,(Text2,j=1,NS),i=1,NObs)
        if(NS.eq.5) write(92,230) (Text1,(Text2,j=1,NS),i=1,NObs)
        if(NS.eq.6) write(92,240) (Text1,(Text2,j=1,NS),i=1,NObs)
      end if

110   format (/// 6x,5(17x,a5,i3,')',  4x))
120   format (///14x,5(17x,a5,i3,')', 26x))
130   format (///14x,5(17x,a5,i3,')', 48x))
140   format (///14x,5(17x,a5,i3,')', 70x))
150   format (///14x,5(17x,a5,i3,')', 92x))
160   format (///14x,5(17x,a5,i3,')',114x))
170   format (///14x,5(17x,a5,i3,')',136x))
180   format (/'     time     ',5(a28,     2x))
190   format (/'     time     ',5(a28, a22,2x))
200   format (/'     time     ',5(a28,2a22,2x))
210   format (/'     time     ',5(a28,3a22,2x))
220   format (/'     time     ',5(a28,4a22,2x))
230   format (/'     time     ',5(a28,5a22,2x))
240   format (/'     time     ',5(a28,6a22,2x))
      return
      end

************************************************************************

      subroutine AtmIn(GWL0L,SinkF,qGWLF,tInit,tMax,Aqh,Bqh,hCritS,
     !                 MaxAL)

      logical SinkF,qGWLF

      write(*,*) 'reading atmospheric information'
      read(31,*)
      read(31,*)
      read(31,*)
      read(31,*)
      read(31,*) SinkF,qGWLF
      read(31,*)
      read(31,*) GWL0L,Aqh,Bqh
      read(31,*)
      read(31,*) tInit,MaxAL
      read(31,*)
      read(31,*) hCritS
      read(31,*)
      do 11 i=1,MaxAL-1
        read(31,*)
11    continue
      read(31,*) tMax
      rewind 31
      do 12 i=1,12
        read(31,*)
12    continue
      return
      end

************************************************************************

      subroutine SinkIn(NMat,NumEl,NumNP,NumElD,KAT,KX,x,y,P0,POptm,P2H,
     !                  P2L,P3,r2H,r2L,Beta)

      integer e
      dimension POptm(NMat),Beta(NumNP),KX(NumElD,4),x(NumNP),y(NumNP)

      write(*,*) 'reading sink information'
      read(30,*)
      read(30,*)
      read(30,*) P0,P2H,P2L,P3,r2H,r2L
      read(30,*)
      read(30,*) (POptm(i),i=1,NMat)
      P0 =-abs(P0)
      P2L=-abs(P2L)
      P2H=-abs(P2H)
      P3 =-abs(P3)
      xMul=1.
      SBeta=0.
      do 12 e=1,NumEl
        NUS=4
        IF(KX(e,3).eq.KX(e,4)) NUS=3
        do 11 k=1,NUS-2
          i=KX(e,1)
          j=KX(e,k+1)
          l=KX(e,k+2)
          CJ=x(i)-x(l)
          CK=x(j)-x(i)
          BJ=y(l)-y(i)
          BK=y(i)-y(j)
          AE=(CK*BJ-CJ*BK)/2.
          if(KAT.eq.1) xMul=2.*3.1416*(x(i)+x(j)+x(l))/3.
          BetaE=(Beta(i)+Beta(j)+Beta(l))/3.
          SBeta=SBeta+xMul*AE*BetaE
11      continue
12    continue
      do 13 i=1,NumNP
        Beta(i)=Beta(i)/SBeta
13    continue
      return
      end

************************************************************************

      subroutine ChemIn(NMat,NumBP,NS,NSD,NumKD,NLevel,cBound,ChPar,
     !                  epsi,tPulse,KodCB,ChemS,TDep,lUpW,cTolA,cTolR,
     !                  MaxItC,lLinear,lEquil,lArtD,PeCr)

      logical lUpW,lTDep,lLinear,lEquil,lArtD
      dimension ChPar(NSD*16+4,NMat),KodCB(NumBP),cBound(NSD,9),
     !          ChemS(NSD,NumKD),TDep(NSD*16+4)

      write(*,*) 'reading transport information'
      NLevel=1
      write(50,110)
      read(30,*)
      read(30,*)
      read(30,*) epsi,lUpW,lArtD,lTDep,cTolA,cTolR,MaxItC,PeCr
      PeCr=amax1(PeCr,0.1)
      if(epsi.lt.0.999) NLevel=2      
      if(lUpW) then
        write(50,120)
      else
        write(50,130)
        if(lArtD) write(50,140) PeCr
      end if
      write(50,150) NS
      read(30,*)
      do 11 M=1,NMat
        read(30,*) (ChPar(j,M),j=1,4)
        write(50,160) M,(ChPar(j,M),j=1,4)
        if(ChPar(4,M).lt.1.) lEquil=.false.
11    continue
      do 13 jj=1,NS
        jjj=(jj-1)*16
        write(50,170) jj
        read(30,*)
        read(30,*) (ChPar(jjj+j,1),j=5,6)
        write(50,180) (ChPar(jjj+j,1),j=5,6)
        read(30,*)
        lLinear=.true.
        do 12 M=1,NMat
          ChPar(jjj+5,M)=ChPar(jjj+5,1)
          ChPar(jjj+6,M)=ChPar(jjj+6,1)
          read(30,*) (ChPar(jjj+j,M),j=7,20)
          write(50,190) M,(ChPar(jjj+j,M),j=7,20)
          if(abs(ChPar(jjj+8,M)-0.0).gt.0.001) then
            write(50,200) M
          else if(abs(ChPar(jjj+9,M)-1.).gt.0.001) then
            write(50,210) M
          else 
            write(50,220) M
          end if
          if(abs(ChPar(jjj+8,M)-0.0).gt.0.001.or.
     !       abs(ChPar(jjj+9,M)-1.0).gt.0.001) lLinear=.false.
12      continue
13    continue
      if(lTDep) then
        read(30,*)
        do 14 jj=1,NS
          jjj=(jj-1)*16
          read(30,*)
          read(30,*) (TDep(jjj+j),j=5,6)
          read(30,*)
          read(30,*) (TDep(jjj+j),j=7,20)
14      continue
      else
       do 15 jj=1,NS*16+4
         TDep(jj)=0.
15     continue
      end if
      read(30,*)
      read(30,*) (KodCB(i),i=1,NumBP)
      read(30,*)
      write(50,230)
      do 16 j=1,NS
        read(30,*) (cBound(j,i),i=1,9)
        write(50,240) (cBound(j,i),i=1,9)
16    continue
      read(30,*)
      read(30,*) tPulse
      write(50,250) tPulse
      do 17 i=1,NS
        do 17 j=1,NumKD
          ChemS(i,j)=0.
17    continue
  
110   format (/' Solute transport information'/1x,28('='))
120   format (/' Upstream weighting finite-element method')
130   format (/' Galerkin finite-element method')
140   format (/' Artificial dispersion is added when Peclet number is',
     !         ' higher than',f10.3)
150   format (/' Number of species in the chain : ',i3//
     !' Mat.     Bulk.D.    DispL     DispT    Fraction')
160   format(i4,f13.4,3f10.4)
170   format(/'    Dif.w.      Dif.g.   ',50('-'),' (',i2,'.solute)')
180   format(2e12.4/' Mat.     KS         Nu         Beta      Henry    
     !  SinkL1     SinkS1     SinkG1     SinkL1`    SinkS1`    SinkG1`  
     !  SinkL0     SinkS0     SinkG0      Alfa')
190   format(i4,14e11.4)
200   format(/' Langmuir nonlinear adsorption isotherm for material ',
     !       i2)
210   format(/' Freundlich nonlinear adsorption isotherm for material ',
     !       i2)
220   format(/' No adsorption or linear adsorp. isotherm for material ',
     !       i2)
230   format (/'   Conc1     Conc2     Conc3     Conc4     cSink     ',
     !'cWell     Conc7     gAtm          d')
240   format (8e10.3,f10.4)
250   format (/' tPulse =   ',f15.3)
      return
      end

************************************************************************

      subroutine TempIn(NMat,NumBP,TPar,KodTB,TBound,Ampl,tPeriod)

      dimension TPar(10,NMat),KodTB(NumBP),TBound(6)

      write(*,*) 'reading heat transport information'
      write(50,110)
      read(30,*)     
      read(30,*)
      do 11 M=1,NMat    
        read(30,*) (TPar(j,M),j=1,10)
        write(50,120) (TPar(j,M),j=1,10)
11    continue   
      read(30,*)
      read(30,*) (KodTB(i),i=1,NumBP)
      read(30,*)
      read(30,*) (TBound(i),i=1,6)
      read(30,*)
      read(30,*) Ampl,tPeriod
      write(50,130)
      write(50,140) (TBound(i),i=1,6)
  
110   format (/' Heat transport information'/1x,26('=')//
     !'    Qn     Qo      Thermal. Disper.      B1         B2         B3
     !        Cn         Co         Cw')
120   format(2f7.3,2f10.3,6e11.3)
130   format(/'    BValT1    BValT2    BValT3    BValT4    BValT5...  ')
140   format(6f10.3)
      return
      end

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||