* Source file INPUT3.FOR |||||||||||||||||||||||||||||||||||||||||||||||

      subroutine BasInf(MaxIt,TolTh,TolH,lWat,lChem,AtmInF,ShortF,SeepF,
     !                  CheckF,FluxF,FreeD,DrainF)

      character*72 Hed
      character*5  LUnit,TUnit,MUnit
      logical lWat,lChem,CheckF,AtmInF,ShortF,SeepF,FluxF,FreeD,DrainF
      dimension IU(12)
      data IU /50,71,72,75,76,77,78,79,80,81,82,84/

      read(30,*)
      read(30,*)
      read(30,*) Hed
      read(30,*)
      read(30,*) LUnit,TUnit,MUnit
      read(30,*)
      read(30,*) MaxIt,TolTh,TolH
      read(30,*)
      read(30,*) lWat,lChem,CheckF,ShortF,FluxF,AtmInF,SeepF,FreeD,
     !           DrainF
      do 11 i=1,12
        write(IU(i),*) Hed
        write(IU(i),*)
        write(IU(i),*)'Program SWMS_3D'
        call getdat(ii,imonth,iday)
        call gettim(ihours,mins,isecs,ii)
        write(IU(i),100) iday,imonth,ihours,mins,isecs
        if(AtmInF) then
          write(IU(i),*)'Time dependent boundary conditions'
        else
          write(IU(i),*)'Time independent boundary conditions'
        end if
        write(IU(i),*) 'Units: L = ',LUnit,', T = ',TUnit,', M = ',MUnit
11    continue
      write(*,*) '-----------------------------------------------------'
      write(*,*) '|                                                   |'
      write(*,*) '|                     SWMS_3D                       |'
      write(*,*) '|                                                   |'
      write(*,*) '|     Code for simulating water flow and solute     |'
      write(*,*) '|      transport in three-dimensional variably      |'
      write(*,*) '|             saturated porous media                |'
      write(*,*) '|                                                   |'
      write(*,*) '|                   version 1.0                     |'
      write(*,*) '|           Last modified: August, 1994             |'
      write(*,*) '|                                                   |'
      write(*,*) '-----------------------------------------------------'
      write(*,*)
      write(*,*) Hed
      write(50,110) MaxIt,TolTh,TolH

100   format(' Date: ',i3,'.',i2,'.','    Time: ',i3,':',i2,':',i2)
110   format(/' Max. number of iterations           ',i4/
     !        ' Absolute water content tolerance [-]',f8.5/ 
     !        ' Absolute pressure head tolerance [L]',f8.5/)
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

110   format(/' MatNum,        thR    thS    tha    thm      alpha      
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

      subroutine DrainIn(NDr,NDrD,NElDrD,NumEl,ND,NED,KElDr,KNoDr,EfDim,
     !                   ConAxx,ConAyy,ConAzz,ConAxy,ConAxz,ConAyz)

      integer e
      dimension ND(NDrD),NED(NDrD),KElDr(NDrD,NElDrD),EfDim(2,NDrD),
     !          ConAxx(NumEl),ConAyy(NumEl),ConAzz(NumEl),
     !          ConAxy(NumEl),ConAxz(NumEl),ConAyz(NumEl),
     !          KNoDr(NDrD,NElDrD)

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
        read(30,*) (KNoDr(i,j),j=1,ND(i))
13    continue
      read(30,*)
      do 14 i=1,NDr
        read(30,*) (KElDr(i,j),j=1,NED(i))
14    continue
      do 16 i=1,Ndr
        rho=EfDim(i,2)/EfDim(i,1)
        A=(1.+0.405*rho**(-4))/(1.-0.405*rho**(-4))
        B=(1.+0.163*rho**(-8))/(1.-0.163*rho**(-8))
        C=(1.+0.067*rho**(-12))/(1.-0.067*rho**(-12))
        Red=376.7/(138.*alog10(rho)+6.48-2.34*A-0.48*B-0.12*C)/DrCorr
        do 16 j=1,NED(i)
          e=KElDr(i,j)
          ConAxx(e)=ConAxx(e)*Red
          ConAyy(e)=ConAyy(e)*Red
          ConAzz(e)=ConAzz(e)*Red
          ConAxy(e)=ConAxy(e)*Red
          ConAxz(e)=ConAxz(e)*Red
          ConAyz(e)=ConAyz(e)*Red
15      continue
16    continue
      return
      end

************************************************************************

      subroutine NodInf(NumNP,NumEl,IJ,NumBP,NumNPD,NumElD,NumBPD,NumKD,
     !                  NObs,NObsD,Kode,Q,Conc,hNew,hOld,hTemp,x,y,z,
     !                  MatNum,Beta,Axz,Bxz,Dxz,CheckF)

      logical CheckF
      dimension Kode(NumNPD),Q(NumNPD),hOld(NumNPD),x(NumNPD),y(NumNPD),
     !          z(NumNPD),hNew(NumNPD),hTemp(NumNPD),MatNum(NumNPD),
     !          Beta(NumNPD),Axz(NumNPD),Bxz(NumNPD),Dxz(NumNPD),
     !          Conc(NumNPD)

      write(*,*) 'reading nodal information'
      read(32,*)
      read(32,*)
      read(32,*) NumNP,NumEl,IJ,NumBP,NObs
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
      read(32,*) n,Kode(n),x(n),y(n),z(n),hOld(n),Conc(n),Q(n),
     !           MatNum(n),Beta(n),Axz(n),Bxz(n),Dxz(n)
      if(Kode(n).gt.NumKD) then
        write(*,*) 'Dimension in NumKD is exceeded'
        stop
      end if
      if(n-k) 12,15,13
12    write(*,130) n
      stop
13    Deno=n-k+1
      DX=(x(n)-x(NPR))/Deno
      DY=(y(n)-y(NPR))/Deno
      DZ=(z(n)-z(NPR))/Deno
      DP=(hOld(n)-hOld(NPR))/Deno
      DConc=(Conc(n)-Conc(NPR))/Deno
      DBeta=(Beta(n)-Beta(NPR))/Deno
      DA=(Axz(n)-Axz(NPR))/Deno
      DB=(Bxz(n)-Bxz(NPR))/Deno
      DD=(Dxz(n)-Dxz(NPR))/Deno
14    x(k)=x(k-1)+DX
      y(k)=y(k-1)+DY
      z(k)=z(k-1)+DZ
      hOld(k)=hOld(k-1)+DP
      Conc(k)=Conc(k-1)+DConc
      Beta(k)=Beta(k-1)+DBeta
      Axz(k)=Axz(k-1)+DA
      Bxz(k)=Bxz(k-1)+DB
      Dxz(k)=Dxz(k-1)+DD
      MatNum(k)=MatNum(k-1)
      Kode(k)=Kode(k-1)
      Q(k)=Q(k-1)
      k=k+1
      if(k.lt.N) goto 14
15    NPR=N
      if(k.lt.NumNP) goto 11

      do 16 n=1,NumNP
        hNew(n)=hOld(n)
        hTemp(n)=hOld(n)
16    continue
      if(CheckF) then
        write(50,110)
        do 17 n=1,NumNP
          write(50,120) n,Kode(n),x(n),y(n),z(n),hOld(n),Q(n),Conc(n),
     !                  MatNum(n),Beta(n),Axz(n),Bxz(n),Dxz(n)
17      continue
      end if

110   format(////' NODAL POINT INFORMATION'//'  NODE NO.',6x,'KODE',
     !       7x,'X,R',13x,'Y',13x,'Z',11x,'.PSI.',12x,'Q',11x,'Conc'/)
120   format(2i10,6e15.6,i10,4f7.3)
130   format(' ERROR IN NodInf AT N=',i5)
      return
      end

************************************************************************

      subroutine ElemIn(NumEl,NumElD,NumNP,KX,LayNum,ConAxx,ConAyy,
     !                  ConAzz,ConAxy,ConAxz,ConAyz,CheckF,ListNE,MBand,
     !                  MBandD,x,y,z,lChem,lOrt)

      logical CheckF,lChem,lOrt
      integer e
      dimension KX(NumElD,9),ConAxx(NumEl),ConAyy(NumEl),ConAzz(NumEl),
     !          ConAxy(NumEl),ConAxz(NumEl),ConAyz(NumEl),LayNum(NumEl),
     !          ListNE(NumNP),x(NumNP),y(NumNP),z(NumNP)

      write(*,*) 'reading element information'
      Num=0
      read(32,*)
      read(32,*)
      AA=3.141592654/180.
      do 14 e=1,NumEL
        if (Num-e) 11,14,12
c11      read(32,*) Num,(KX(Num,i),i=1,9),ConAxz(Num),ConAxx(Num),
c     !             ConAyy(Num),ConAzz(Num),LayNum(Num)
c        Ang=AA*ConAxz(Num)
c        CAxx=ConAxx(Num)
c        CAzz=ConAzz(Num)
c        ConAxx(Num)=CAxx*cos(Ang)*cos(Ang) + CAzz*sin(Ang)*sin(Ang)
c        ConAzz(Num)=CAxx*sin(Ang)*sin(Ang) + CAzz*cos(Ang)*cos(Ang)
c        ConAxz(Num)=(CAxx-CAzz)*sin(Ang)*cos(Ang)
c        ConAxy(Num)=0.
c        ConAyz(Num)=0.
11      read(32,*) Num,(KX(Num,i),i=1,9),c11,c22,c33,a11,a22,a33,a12,
     !             a13,a23,LayNum(Num)
        ConAxx(Num)=c11*a11*a11+c22*a12*a12+c33*a13*a13
        ConAyy(Num)=c11*a12*a22+c22*a22*a22+c33*a23*a23
        ConAzz(Num)=c11*a13*a13+c22*a23*a13+c33*a33*a33
        ConAxy(Num)=c11*a11*a12+c22*a12*a22+c33*a13*a23
        ConAxz(Num)=c11*a11*a13+c22*a12*a23+c33*a13*a33
        ConAyz(Num)=c11*a12*a13+c22*a22*a23+c33*a23*a33
        if(Num.eq.e) goto 14
12      do 13 i=1,8
          KX(e,i)=KX(e-1,i)+1
13      continue
        if(KX(e-1,9).eq.1) KX(e,9)=2
        if(KX(e-1,9).eq.2) KX(e,9)=1
        ConAxz(e)=ConAxz(e-1)
        ConAxy(e)=ConAxy(e-1)
        ConAyz(e)=ConAyz(e-1)
        ConAxx(e)=ConAxx(e-1)
        ConAyy(e)=ConAyy(e-1)
        ConAzz(e)=ConAzz(e-1)
        LayNum(e)=LayNum(e-1)
14    continue
      if(CheckF) then
        write(50,110)
        do 15 e=1,NumEl
          write(50,120) e,(KX(e,i),i=1,8),ConAxx(e),ConAyy(e),ConAzz(e),
     !                  ConAxy(e),ConAxz(e),ConAyz(e),LayNum(e)
15      continue
      end if

      do 16 i=1,NumNP
        ListNE(i)=0
16    continue
      do 18 e=1,NumEl
        NUS=5
        if(KX(e,7).eq.0) NUS=3
        if(KX(e,5).eq.0) NUS=1
        do 17 n=1,NUS
          call Elem(NUS,NumElD,e,n,KX,i,j,k,l)
          ListNE(i)=ListNE(i)+1
          ListNE(j)=ListNE(j)+1
          ListNE(k)=ListNE(k)+1
          ListNE(l)=ListNE(l)+1
17      continue
18    continue

      lOrt=.false.
      MBand=1
      do 20 e=1,NumEl
        NUS=5
        if(KX(e,7).eq.0) NUS=3
        if(KX(e,5).eq.0) NUS=1
        do 19 n=1,NUS
          call Elem(NUS,NumElD,e,n,KX,i,j,k,l)
          Bl=+(y(j)-y(i))*(z(k)-z(i))-(y(k)-y(i))*(z(j)-z(i))
          Cl=-(x(j)-x(i))*(z(k)-z(i))+(x(k)-x(i))*(z(j)-z(i))
          Dl=+(x(j)-x(i))*(y(k)-y(i))-(x(k)-x(i))*(y(j)-y(i))
          Det=(x(l)-x(i))*Bl+(y(l)-y(i))*Cl+(z(l)-z(i))*Dl
          if(Det.le.0.) then
            write( *,130) e,n
            write(50,130) e,n
            stop
          end if
          if(abs(i-j).gt.MB) MB=abs(i-j)
          if(abs(i-k).gt.MB) MB=abs(i-k)
          if(abs(i-l).gt.MB) MB=abs(i-l)
          if(abs(j-k).gt.MB) MB=abs(j-k)
          if(abs(j-l).gt.MB) MB=abs(j-l)
          if(abs(k-l).gt.MB) MB=abs(k-l)
          if(MB.gt.MBand) MBand=MB
19      continue
20    continue
      MBand=MBand+1
      if(MBand.gt.MBandD.or.(lChem.and.2*MBand-1.gt.MBandD)) lOrt=.true.
      if(MBand.gt.10.or.NumNP.gt.300) lOrt=.true.

110   format (////' Element Information'//' Element    C O R N E R    N 
     !O D E S   ConAxx    ConAyy    ConAzz   ConAxy    ConAxz    ConAyz 
     ! LayNum'/)
120   format (i6,i8,7i5,6f8.3,i5)
130   format (//' Wrong orientation of nodes of element: ',i7,
     !          ', subelement: ',i3,' !!!'
     !        //' Determinant is smaller than zero !!!')
      return
      end

************************************************************************

      subroutine GeomIn(NumKD,NumNP,NumBP,NObs,SWidth,Width,Kode,KXB,
     !                  rLen,Node)

      character*50 Text1
      dimension KXB(NumBP),Width(NumBP),SWidth(NumKD),Kode(NumNP),
     !          Node(NObs)

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
        Text1='    hNew     theta     conc    '
        write(92,110) (Node(j),j=1,NObs)
        write(92,120) (Text1,i=1,NObs)
      end if
110   format (///14x,5(11x,'Node(',i3,')',11x))
120   format (/'     time     ',5(a31))
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

      subroutine SinkIn(NMat,NumEl,NumNP,NumElD,KX,x,y,z,P0,POptm,P2H,
     !                  P2L,P3,r2H,r2L,Beta)

      integer e
      dimension POptm(NMat),Beta(NumNP),KX(NumElD,8),x(NumNP),y(NumNP),
     !          z(NumNP) 

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
      SBeta=0.
      do 12 e=1,NumEl
        NUS=5
        if(KX(e,7).eq.0) NUS=3
        if(KX(e,5).eq.0) NUS=1
        do 11 n=1,NUS
          call Elem(NUS,NumElD,e,n,KX,i,j,k,l)
          Bl=+(y(j)-y(i))*(z(k)-z(i))-(y(k)-y(i))*(z(j)-z(i))
          Cl=-(x(j)-x(i))*(z(k)-z(i))+(x(k)-x(i))*(z(j)-z(i))
          Dl=+(x(j)-x(i))*(y(k)-y(i))-(x(k)-x(i))*(y(j)-y(i))
          Det=(x(l)-x(i))*Bl+(y(l)-y(i))*Cl+(z(l)-z(i))*Dl
          VE=Det/6.
          BetaE=(Beta(i)+Beta(j)+Beta(k)+Beta(l))/4.
          SBeta=SBeta+VE*BetaE
11      continue
12    continue
      do 13 i=1,NumNP
        Beta(i)=Beta(i)/SBeta
13    continue
      return
      end

************************************************************************

      subroutine ChemIn(NMat,NumBP,cBound,ChPar,epsi,tPulse,KodCB,
     !                  NLevel,lUpW,lArtD,PeCr)

      logical lUpW,lArtD
      dimension ChPar(10,NMat),KodCB(NumBP),cBound(12)

      write(*,*) 'reading transport information'
      NLevel=1
      write(50,110)
      read(30,*)
      read(30,*)
      read(30,*) epsi,lUpW,lArtD,PeCr
      PeCr=amax1(PeCr,0.1)
      if(epsi.lt.0.999) NLevel=2      
      read(30,*) 
      if(lUpW) then
        write(50,120)
      else
        write(50,130)
        if(lArtD) write(50,140) PeCr
      end if
      write(50,150)
      do 11 M=1,NMat    
        read(30,*) (ChPar(j,M),j=1,9)
        write(50,160) (ChPar(j,M),j=1,9)
11    continue   
      read(30,*)
      read(30,*) (KodCB(i),i=1,NumBP)
      read(30,*)
      read(30,*) (cBound(i),i=1,12)
      write(50,170) (cBound(i),i=1,12)
      read(30,*)
      read(30,*) tPulse
      write(50,180) tPulse
  
110   format (/' Solute transport information'/1x,28('='))
120   format (/' Upstream weighting finite-element method')
130   format (/' Galerkin finite-element method')
140   format (/' Artificial dispersion is added when Peclet number is',
     !         ' higher than',f10.3)
150   format (/'  Bulk.d.   Difus.         Disper.        Adsorp.    ',
     !         'SinkL1    SinkS1    SinkL0    SinkS0')
160   format (10e10.3) 
170   format (/'   Conc1     Conc2     Conc3     Conc4     Conc5     C',
     !         'onc6'/6e10.3//'   Conc7     Conc8     Conc9     Conc10',
     !         '    cWell     cSink'/6e10.3)
180   format (/' tPulse =   ',f15.3)
      return
      end

************************************************************************

      subroutine Elem(NUS,NumElD,e,n,KX,i,j,k,l)
      dimension KX(NumElD,9)
      integer e
      
      if(NUS.eq.1) then
        i=KX(e,1)
        j=KX(e,2)
        k=KX(e,3)
        l=KX(e,4)
      else if(NUS.eq.3) then
        if(KX(e,9).eq.1) then
          if(n.eq.1) then
            i=KX(e,4)
            j=KX(e,2)
            k=KX(e,3)
            l=KX(e,1)
          else if(n.eq.2) then
            i=KX(e,4)
            j=KX(e,5)
            k=KX(e,6)
            l=KX(e,2)
          else if(n.eq.3) then
            i=KX(e,4)
            j=KX(e,2)
            k=KX(e,6)
            l=KX(e,3)
          end if
        else if(KX(e,9).eq.2) then
          if(n.eq.1) then
            i=KX(e,4)
            j=KX(e,2)
            k=KX(e,3)
            l=KX(e,1)
          else if(n.eq.2) then
            i=KX(e,4)
            j=KX(e,5)
            k=KX(e,3)
            l=KX(e,2)
          else if(n.eq.3) then
            i=KX(e,4)
            j=KX(e,5)
            k=KX(e,6)
            l=KX(e,3)
          end if
        else if(KX(e,9).eq.3) then
          if(n.eq.1) then
            i=KX(e,4)
            j=KX(e,5)
            k=KX(e,6)
            l=KX(e,1)
          else if(n.eq.2) then
            i=KX(e,1)
            j=KX(e,5)
            k=KX(e,6)
            l=KX(e,3)
          else if(n.eq.3) then
            i=KX(e,1)
            j=KX(e,5)
            k=KX(e,3)
            l=KX(e,2)
          end if
        else if(KX(e,9).eq.4) then
          if(n.eq.1) then
            i=KX(e,4)
            j=KX(e,5)
            k=KX(e,6)
            l=KX(e,1)
          else if(n.eq.2) then
            i=KX(e,1)
            j=KX(e,5)
            k=KX(e,6)
            l=KX(e,2)
          else if(n.eq.3) then
            i=KX(e,1)
            j=KX(e,2)
            k=KX(e,6)
            l=KX(e,3)
          end if
        else if(KX(e,9).eq.5) then
          if(n.eq.1) then
            i=KX(e,4)
            j=KX(e,5)
            k=KX(e,6)
            l=KX(e,2)
          else if(n.eq.2) then
            i=KX(e,4)
            j=KX(e,2)
            k=KX(e,6)
            l=KX(e,1)
          else if(n.eq.3) then
            i=KX(e,1)
            j=KX(e,2)
            k=KX(e,6)
            l=KX(e,3)
          end if
        else if(KX(e,9).eq.6) then
          if(n.eq.1) then
            i=KX(e,4)
            j=KX(e,5)
            k=KX(e,6)
            l=KX(e,3)
          else if(n.eq.2) then
            i=KX(e,1)
            j=KX(e,5)
            k=KX(e,3)
            l=KX(e,2)
          else if(n.eq.3) then
            i=KX(e,1)
            j=KX(e,4)
            k=KX(e,3)
            l=KX(e,5)
          end if
        end if
      else if(NUS.eq.5) then
        if(KX(e,9).eq.1) then
          if(n.eq.1) then
            i=KX(e,1)
            j=KX(e,4)
            k=KX(e,7)
            l=KX(e,3)
          else if(n.eq.2) then
            i=KX(e,1)
            j=KX(e,6)
            k=KX(e,4)
            l=KX(e,2)
          else if(n.eq.3) then
            i=KX(e,5)
            j=KX(e,6)
            k=KX(e,7)
            l=KX(e,1)
          else if(n.eq.4) then
            i=KX(e,7)
            j=KX(e,6)
            k=KX(e,8)
            l=KX(e,4)
          else if(n.eq.5) then
            i=KX(e,1)
            j=KX(e,6)
            k=KX(e,7)
            l=KX(e,4)
          end if
        else
          if(n.eq.1) then
            i=KX(e,5)
            j=KX(e,2)
            k=KX(e,3)
            l=KX(e,1)
          else if(n.eq.2) then
            i=KX(e,5)
            j=KX(e,6)
            k=KX(e,8)
            l=KX(e,2)
          else if(n.eq.3) then
            i=KX(e,5)
            j=KX(e,8)
            k=KX(e,7)
            l=KX(e,3)
          else if(n.eq.4) then
            i=KX(e,5)
            j=KX(e,2)
            k=KX(e,8)
            l=KX(e,3)
          else if(n.eq.5) then
            i=KX(e,2)
            j=KX(e,8)
            k=KX(e,3)
            l=KX(e,4)
          end if
        end if
      end if
      return
      end

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||