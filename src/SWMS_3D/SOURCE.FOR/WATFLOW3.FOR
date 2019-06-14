* Source file WATFLOW3.FOR |||||||||||||||||||||||||||||||||||||||||||||

      subroutine WatFlow(NumNP,NumEl,NumElD,NTab,NTabD,MBand,MBandD,
     !                   NMat,NSeep,NSeepD,NumSPD,NSP,NP,NumBP,ItCum,
     !                   MaxIt,Iter,Kode,t,dt,dtMin,dtOpt,dtOld,
     !                   tOld,hCritA,hCritS,TolTh,TolH,rLen,Width,rTop,
     !                   vMeanR,hMeanR,AtmInf,SinkF,SeepF,qGWLF,FreeD,
     !                   Par,hTab,ConTab,CapTab,TheTab,hNew,hOld,hTemp,
     !                   thR,thSat,ThNew,ThOld,MatNum,Con,Cap,ConSat,
     !                   Axz,Bxz,Dxz,hSat,A,B,Q,F,x,y,z,KX,Sink,DS,Beta,
     !                   ConAxx,ConAyy,ConAzz,ConAxy,ConAxz,ConAyz,KXB,
     !                   Explic,GWL0L,Aqh,Bqh,lWat,TLevel,lOrt,DrainF,
     !                   ND,NDr,NDrD,NElDrD,KNoDr,rRoot,P0,POptm,P2H,
     !                   P2L,P3,r2H,r2L,ConO,
     !                   A1,B1,NumNPD,IAD,IADN,IADD,VRV,RES,RQI,RQ,QQ,
     !                   QI,RQIDOT,ECNVRG,RCNVRG,ACNVRG,MNorth,MaxItO)

      logical AtmInf,SinkF,SeepF,Explic,ItCrit,lWat,qGWLF,FreeD,lOrt,
     !        DrainF
      double precision A,B,A1,B1,VRV,RES,RQI,RQ,QQ,QI,RQIDOT,ECNVRG,
     !                 RCNVRG,ACNVRG
      integer TLevel
      dimension A(MBandD,NumNP),B(NumNP),Kode(NumNP),Q(NumNP),F(NumNP),
     !          hNew(NumNP),hTemp(NumNP),hOld(NumNP),ConSat(NMat),
     !          hTab(NTab),ConTab(NTabD,NMat),CapTab(NTabD,NMat),
     !          Con(NumNP),Cap(NumNP),x(NumNP),y(NumNP),z(NumNP),
     !          MatNum(NumNP),KX(NumElD,9),KXB(NumBP),Par(10,NMat),
     !          Width(NumBP),ConAxx(NumEl),ConAyy(NumEl),ConAzz(NumEl),
     !          ConAxz(NumEl),ConAxy(NumEl),ConAyz(NumEl),hSat(NMat),
     !          NP(NSeepD,NumSPD),NSP(NSeepD),DS(NumNP),Beta(NumNP),
     !          Axz(NumNP),Bxz(NumNP),Dxz(NumNP),Sink(NumNP),thR(NMat),
     !          thSat(NMat),TheTab(NTabD,NMat),ThNew(NumNP),ND(NDrD),
     !          ThOld(NumNP),POptm(NMat),ConO(NumNP),KNoDr(NDrD,NElDrD),
     !          A1(MBandD,NumNP),B1(NumNP),RES(NumNP),IAD(MBandD,NumNP),
     !          IADN(NumNP),IADD(NumNP),VRV(NumNP),RQI(NumNPD,MNorth),
     !          RQ(NumNP),QQ(NumNP),RQIDOT(MNorth),QI(NumNPD,MNorth)

      if(lWat.and.TLevel.ne.1) then
        do 10 i=1,NumNP
          hOld(i) =hNew(i)
          ThOld(i)=ThNew(i)
          ConO(i)=Con(i)
          if(Kode(i).lt.1) then
            hTemp(i)=hNew(i)+(hNew(i)-hOld(i))*dt/dtOld
            hNew(i) =hTemp(i)
          else
            hTemp(i)=hNew(i)
          end if
10      continue
      end if

11    continue

      Iter=0
      Explic=.false.

* --- Beginning of iteration loop --------------------------------------

12    continue
      if(SinkF.and..not.lWat.and.Iter.ne.0)
     !  call SetSnk(NumNP,NMat,MatNum,hNew,rRoot,Sink,P0,POptm,P2H,P2L,
     !              P3,r2H,r2L,Beta,rLen)
      call SetMat(NumNP,NTab,NTabD,NMat,hTab,ConTab,CapTab,hNew,hOld,
     !            MatNum,Par,Con,Cap,ConSat,Axz,Bxz,Dxz,hSat,hTemp,
     !            Explic,TheTab,thSat,thR,ThNew)
      call Reset (A,B,Kode,Q,hNew,F,Con,Cap,x,y,z,KX,NumNP,NumEl,NumElD,
     !            MBand,MBandD,dt,Iter,SinkF,Sink,DS,Beta,rLen,
     !            vMeanR,hMeanR,ConAxx,ConAyy,ConAzz,ConAxy,ConAxz,
     !            ConAyz,ThNew,ThOld,lWat,lOrt,IAD,IADN,IADD,B1)
      call Shift (NumNP,NumBP,NSeepD,NumSPD,NSeep,NSP,NP,hNew,hOld,Q,
     !            Kode,rTop,Width,KXB,hCritA,hCritS,SeepF,AtmInF,Explic,
     !            qGWLF,FreeD,GWL0L,Aqh,Bqh,Con,ConO,DrainF,ND,NDr,NDrD,
     !            NElDrD,KNoDr,lWat)
      call Dirich(A,B,Kode,hNew,NumNP,MBand,MBandD,lOrt,IADD)
      if(lOrt) then
        call ILU (A,NumNP,MBandD,IAD,IADN,IADD,A1)
        North=0
        call OrthoMin(A,B1,B,NumNP,MBandD,NumNPD,IAD,IADN,IADD,A1,VRV,
     !                RES,RQI,RQ,QQ,QI,RQIDOT,ECNVRG,RCNVRG,ACNVRG,
     !                North,MNorth,MaxItO)
      else
        call Solve (A,B,NumNP,MBand,MBandD)
      end if

      do 13 i=1,NumNP
        hTemp(i)=hNew(i)
        if(lOrt) B(i)=B1(i)
        hNew(i)=sngl(B(i))
13    continue
      Iter =Iter+1
      ItCum=ItCum+1
      if(Explic) goto 18

*     Test for convergence
      ItCrit=.true.
      do 14 i=1,NumNP
        m=MatNum(i)
        EpsTh=0.
        EpsH=0.
        if(hTemp(i).lt.hSat(m).and.hNew(i).lt.hSat(m)) then
          Th=ThNew(i)+Cap(i)*(hNew(i)-hTemp(i))/(ThSat(m)-ThR(m))/Axz(i)
          EpsTh=abs(ThNew(i)-Th)
        else
          EpsH=abs(hNew(i)-hTemp(i))
        end if
        if(EpsTh.gt.TolTh.or.EpsH.gt.TolH) then
          ItCrit=.false.
          goto 15
        end if
14    continue
15    continue

      if(.not.ItCrit) then
        if(Iter.lt.MaxIt.or.(.not.lWat.and.Iter.lt.5*MaxIt)) then
          goto 12
        else if(dt.le.dtMin) then
          Explic=.true.
          do 16 i=1,NumNP
            hNew(i) =hOld(i)
            hTemp(i)=hOld(i)
16        continue
          goto 12
        else if(.not.lWat) then
          write(*,*) ' No steady state solution found'
          stop
        else
          do 17 i=1,NumNP
            hNew(i) =hOld(i)
            hTemp(i)=hOld(i)
17        continue
          dt=amax1(dt/3.,dtMin)
          dtOpt=dt
          t=tOld+dt
          goto 11
        end if
      end if

* --- End of iteration loop --------------------------------------------

18    continue

      if(.not.lWat.and.TLevel.eq.1) then
        write(*,101) Iter
        do 19 i=1,NumNP
          hOld(i) =hNew(i)
          ThOld(i)=ThNew(i)
19      continue
      end if

101   format(' Steady state was reached after',i6,' iterations.')
      return
      end

************************************************************************

      subroutine Reset(A,B,Kode,Q,hNew,F,Con,Cap,x,y,z,KX,NumNP,NumEl,
     !                 NumElD,MBand,MBandD,dt,Iter,SinkF,Sink,DS,
     !                 Beta,rLen,vMeanR,hMeanR,ConAxx,ConAyy,ConAzz,
     !                 ConAxy,ConAxz,ConAyz,ThNew,ThOld,lWat,lOrt,
     !                 IAD,IADN,IADD,B1)

      logical SinkF,lWat,lOrt
      double precision A,B,B1
      integer e
      dimension A(MBandD,NumNP),B(NumNP),Q(NumNP),hNew(NumNP),F(NumNP),
     !          Con(NumNP),Cap(NumNP),x(NumNP),y(NumNP),z(NumNP),
     !          KX(NumElD,9),ConAxx(NumEl),ConAyy(NumEl),ConAzz(NumEl),
     !          ConAxz(NumEl),ConAxy(NumEl),ConAyz(NumEl),Kode(NumNP),
     !          Sink(NumNP),DS(NumNP),Beta(NumNP),ThNew(NumNP),
     !          ThOld(NumNP),EE(4,4),iLoc(4),Bi(4),Ci(4),Di(4),
     !          B1(NumNP),IAD(MBandD,NumNP),IADN(NumNP),IADD(NumNP)

*     Initialisation
      if(Iter.eq.0) then
        vMeanR=0.
        hMeanR=0.
        VolR =0.
      end if
      do 12 i=1,NumNP
        B(i)=0.
        if(lOrt) B1(i)=hNew(i)
        F(i)=0.
        if(Iter.eq.0) DS(i)=0.
        do 11 j=1,MBandD
          A(j,i)=0.
11      continue
12    continue

*     Loop on elements
      do 16 e=1,NumEl
        CAxx=ConAxx(e)
        CAyy=ConAyy(e)
        CAzz=ConAzz(e)
        CAxy=ConAxy(e)
        CAxz=ConAxz(e)
        CAyz=ConAyz(e)
        NUS=5
        if(KX(e,7).eq.0) NUS=3
        if(KX(e,5).eq.0) NUS=1

*       Loop on subelements
        do 15 n=1,NUS
          call Elem(NUS,NumElD,e,n,KX,i,j,k,l)
          iLoc(1)=i
          iLoc(2)=j
          iLoc(3)=k
          iLoc(4)=l
          Bi(1)=-(y(k)-y(j))*(z(l)-z(j))+(y(l)-y(j))*(z(k)-z(j))
          Bi(2)=+(y(l)-y(k))*(z(i)-z(k))-(y(i)-y(k))*(z(l)-z(k))
          Bi(3)=-(y(i)-y(l))*(z(j)-z(l))+(y(j)-y(l))*(z(i)-z(l))
          Bi(4)=+(y(j)-y(i))*(z(k)-z(i))-(y(k)-y(i))*(z(j)-z(i))
          Ci(1)=+(x(k)-x(j))*(z(l)-z(j))-(x(l)-x(j))*(z(k)-z(j))
          Ci(2)=-(x(l)-x(k))*(z(i)-z(k))+(x(i)-x(k))*(z(l)-z(k))
          Ci(3)=+(x(i)-x(l))*(z(j)-z(l))-(x(j)-x(l))*(z(i)-z(l))
          Ci(4)=-(x(j)-x(i))*(z(k)-z(i))+(x(k)-x(i))*(z(j)-z(i))
          Di(1)=-(x(k)-x(j))*(y(l)-y(j))+(x(l)-x(j))*(y(k)-y(j))
          Di(2)=+(x(l)-x(k))*(y(i)-y(k))-(x(i)-x(k))*(y(l)-y(k))
          Di(3)=-(x(i)-x(l))*(y(j)-y(l))+(x(j)-x(l))*(y(i)-y(l))
          Di(4)=+(x(j)-x(i))*(y(k)-y(i))-(x(k)-x(i))*(y(j)-y(i))
          Det=(x(l)-x(i))*Bi(4)+(y(l)-y(i))*Ci(4)+(z(l)-z(i))*Di(4)
          VE=Det/6.
          CapE=(Cap(i)+Cap(j)+Cap(k)+Cap(l))/4.
          ConE=(Con(i)+Con(j)+Con(k)+Con(l))/4.
          AMul=ConE/36./VE
          BMul=ConE/6.
          FMul=VE/20.
          BetaE=(Beta(i)+Beta(j)+Beta(k)+Beta(l))/4.
          if(SinkF.and.BetaE.gt.0..and.Iter.eq.0) then
            SinkE=(Sink(i)+Sink(j)+Sink(k)+Sink(l))/4.
            DS(i)=DS(i)+FMul*(4.*SinkE+Sink(i))
            DS(j)=DS(j)+FMul*(4.*SinkE+Sink(j))
            DS(k)=DS(k)+FMul*(4.*SinkE+Sink(k))
            DS(l)=DS(l)+FMul*(4.*SinkE+Sink(l))
            hNewE=(hNew(i)+hNew(j)+hNew(k)+hNew(l))/4.
            vMeanR=vMeanR+VE*SinkE/rLen
            hMeanR=hMeanR+hNewE*VE
            VolR=VolR+VE
          end if
          do 14 ii=1,4
            iG=iLoc(ii)
            F(iG)=F(iG)+FMul*5.
            B(iG)=B(iG)+BMul*(CAxz*Bi(ii)+CAyz*Ci(ii)+CAzz*Di(ii))
            do 13 jj=1,4
              jG=iLoc(jj)
              EE(ii,jj)=CAxx*Bi(ii)*Bi(jj)+CAyy*Ci(ii)*Ci(jj)+
     !                  CAzz*Di(ii)*Di(jj)+
     !                  CAxy*(Bi(ii)*Ci(jj)+Bi(jj)*Ci(ii))+
     !                  CAxz*(Bi(ii)*Di(jj)+Bi(jj)*Di(ii))+
     !                  CAyz*(Ci(ii)*Di(jj)+Ci(jj)*Di(ii))
              if(lOrt) then
                call Find(iG,jG,kk,NumNP,MBandD,IAD,IADN)
                A(kk,iG)=A(kk,iG)+AMul*EE(ii,jj)
              else
                iB=iG-jG+1
                if(iB.ge.1) A(iB,jG)=A(iB,jG)+AMul*EE(ii,jj)
              end if
13          continue
14        continue
15      continue

16    continue
      if(VolR.gt.0..and.Iter.eq.0) hMeanR=hMeanR/VolR

*     Determine boundary fluxes
      do 19 n=1,NumNP
        if(Kode(n).lt.1) goto 19
        QN=sngl(B(n))+DS(n)+F(n)*(ThNew(n)-ThOld(n))/dt
        if(lOrt) then
          do 17 j=1,IADN(n)
            QN=QN+sngl(A(j,n))*hNew(IAD(j,n))
17        continue
        else
          QN=QN+sngl(A(1,n))*hNew(n)
          do 18 j=2,MBand
            k=n-j+1
            if(k.ge.1) then
              QN=QN+sngl(A(j,k))*hNew(k)
            end if
            k=n+j-1
            if(K.le.NumNP) then
              QN=QN+sngl(A(j,n))*hNew(k)
            end if
18        continue
        end if
        Q(n)=QN
19    continue

*     Complete construction of RHS vector and form effective matrix
      do 20 i=1,NumNP
        if(.not.lWat) F(i)=0.
        j=1
        if(lOrt) j=IADD(i)
        A(j,i)=A(j,i)+F(i)*Cap(i)/dt
        B(i)=F(i)*Cap(i)*hNew(i)/dt-F(i)*(ThNew(i)-ThOld(i))/dt+
     !       Q(i)-B(i)-DS(i)
20    continue
      return
      end

************************************************************************

      subroutine Dirich(A,B,Kode,hNew,NumNP,MBand,MBandD,lOrt,IADD)

      logical lOrt
      double precision A,B
      dimension A(MBandD,NumNP),B(NumNP),Kode(NumNP),hNew(NumNP),
     !          IADD(NumNP)

      do 12 n=1,NumNP
        if(Kode(n).lt.1) goto 12
        if(lOrt) then
          A(IADD(n),n)=10.d30
          B(n)=10.d30*hNew(n)
        else
          do 11 m=2,MBand
            k=n-m+1
            if(k.gt.0) then
              B(k)=B(k)-A(m,k)*hNew(n)
              A(m,k)=0.
            end if
            l=n+m-1
            if(NumNP-l.ge.0) then
              B(l)=B(l)-A(m,n)*hNew(n)
            end if
            A(m,n)=0.
11        continue
          A(1,n)=1.
          B(n)=hNew(n)
        end if
12    continue
      return
      end

************************************************************************

      subroutine Shift(NumNP,NumBP,NSeepD,NumSPD,NSeep,NSP,NP,hNew,hOld,
     !                 Q,Kode,EI,Width,KXB,hCritA,hCritS,SeepF,AtmInF,
     !                 Explic,qGWLF,FreeD,GWL0L,Aqh,Bqh,Con,ConO,DrainF,
     !                 ND,NDr,NDrD,NElDrD,KNoDr,lWat)

      logical SeepF,AtmInF,Explic,qGWLF,FreeD,DrainF,lWat
      dimension Kode(NumNP),Q(NumNP),hNew(NumNP),hOld(NumNP),
     !          Width(NumBP),KXB(NumBP),NP(NSeepD,NumSPD),NSP(NSeepD),
     !          Con(NumNP),ConO(NumNP),ND(NDrD),KNoDr(NDrD,NElDrD)

*     Modify conditions on seepage faces
      if(SeepF) then
        do 12 i=1,NSeep
          NS=NSP(i)
          do 11 j=1,NS
            n=NP(i,j)
            if(Kode(n).eq.-2) then
              if(hNew(n).ge.0.) then
                Kode(n)=2
                hNew(n)=0.
              end if
            else
              if(Q(n).ge.0.) then
                Kode(n)=-2
                Q(n)=0.
              end if
            end if
11        continue
12      continue
      end if

*     Modify conditions in drainage node
      if(DrainF) then
        do 13 i=1,NDr
          do 13 j=1,ND(i)
            n=KNoDr(i,ND(i))
            if(Kode(n).eq.-5) then
              if(hNew(n).ge.0.) then
                Kode(n)=5
                hNew(n)=0.
              end if
            else
              if(Q(n).ge.0.) then
                Kode(n)=-5
                Q(n)=0.
              end if
            end if
13      continue
      end if

*     Modify potential surface flux boundaries
      if(AtmInF) then
        do 14 i=1,NumBP
          n=KXB(i)
          k=Kode(n)
          if(Explic.and.iabs(k).eq.4) then
            Kode(n)=-iabs(k)
            goto 14
          end if

*         Critical surface pressure on ...
          if(K.eq.4) then
            if(abs(Q(n)).gt.abs(-EI*Width(i)).or.Q(n)*(-EI).le.0) then
              Kode(n)=-4
              Q(n)=-EI*Width(i)
            end if
            goto 14
          end if

*         Surface flux on ...
          if(K.eq.-4) then
            if(hNew(n).le.hCritA) then
              Kode(n)=4
              hNew(n)=hCritA
              goto 14
            end if
            if(hNew(n).ge.hCritS) then
              Kode(n)=4
              hNew(n)=hCritS
            end if
          end if

*         Time variable flux
          if(K.eq.-3) then
            if(lWat) then
              if(qGWLF) Q(n)=-Width(i)*Fqh(hOld(n)-GWL0L,Aqh,Bqh)
            else
              if(qGWLF) Q(n)=-Width(i)*Fqh(hNew(n)-GWL0L,Aqh,Bqh)
            end if
          end if
14      continue
      end if
      
*     Free Drainage      
      if(FreeD) then      
        do 15 i=1,NumBP
          n=KXB(i)
          k=Kode(n)
          if(lWat) then
            if(k.eq.-3) Q(n)=-Width(i)*ConO(n)
          else
            if(k.eq.-3) Q(n)=-Width(i)*Con(n)
          end if
15      continue
      end if
      return
      end

************************************************************************

      subroutine SetMat(NumNP,NTab,NTabD,NMat,hTab,ConTab,CapTab,hNew,
     !                  hOld,MatNum,Par,Con,Cap,ConSat,Axz,Bxz,Dxz,hSat,
     !                  hTemp,Explic,TheTab,thSat,thR,theta)

      logical Explic
      dimension hTab(NTab),ConTab(NTabD,NMat),CapTab(NTabD,NMat),
     !          hNew(NumNP),hOld(NumNP),MatNum(NumNP),Par(10,NMat),
     !          Con(NumNP),Cap(NumNP),ConSat(NMat),hSat(NMat),
     !          Axz(NumNP),Bxz(NumNP),Dxz(NumNP),hTemp(NumNP),
     !          TheTab(NTabD,NMat),thSat(NMat),thR(NMat),theta(NumNP)

      alh1=alog10(-hTab(1))
      dlh =(alog10(-hTab(NTab))-alh1)/(NTab-1)
      do 11 i=1,NumNP
        M=MatNum(i)
        hi1=amin1(hSat(M),hTemp(i)/Axz(i))
        hi2=amin1(hSat(M), hNew(i)/Axz(i))
        if(Explic) hi2=hi1
        hiM=0.1*hi1+0.9*hi2
        if(hi1.ge.hSat(M).and.hi2.ge.hSat(M)) then
          Ci=ConSat(M)
        else if(hiM.ge.hTab(NTab).and.hiM.le.hTab(1)) then
          iT=int((alog10(-hiM)-alh1)/dlh)+1
          S1=(ConTab(iT+1,M)-ConTab(iT,M))/(hTab(iT+1)-hTab(iT))
          Ci=ConTab(iT,M)+S1*(hiM-hTab(iT))
        else
          Ci=FK(hiM,Par(1,M))
        end if
        Con(i)=Bxz(i)*Ci
        hi1=hOld(i)/Axz(i)
        hi2=hNew(i)/Axz(i)
        if(Explic) hi2=hi1
        if(hi2.ge.hSat(M)) then
          Ci=0
          Ti=thSat(M)
        else if(hi2.ge.hTab(NTab).and.hi2.le.hTab(1)) then
          iT=int((alog10(-hi2)-alh1)/dlh)+1
          S1=(CapTab(iT+1,M)-CapTab(iT,M))/(hTab(iT+1)-hTab(iT))
          S2=(TheTab(iT+1,M)-TheTab(iT,M))/(hTab(iT+1)-hTab(iT))
          Ci=CapTab(iT,M)+S1*(hi2-hTab(iT))
          Ti=TheTab(iT,M)+S2*(hi2-hTab(iT))
        else
          Ci=FC(hi2,Par(1,M))
          Ti=FQ(hi2,Par(1,M))
        end if
        Cap(i)=Ci*Dxz(i)/Axz(i)
        theta(i)=thR(M)+(Ti-thR(M))*Dxz(i)
11    continue
      return
      end

************************************************************************

      subroutine Solve(A,B,NumNP,MBand,MBandD)

      double precision A,B,C
      dimension A(MBandD,NumNP),B(NumNP)

*     Reduction
      do 13 n=1,NumNP
        do 12 m=2,MBand
          if(abs(A(m,n)).lt.1.e-30) goto 12
          C=A(m,n)/A(1,n)
          i=n+m-1
          if(i.gt.NumNP) goto 12
          j=0
          do 11 k=m,MBand
            j=j+1
            A(j,i)=A(j,i)-C*A(k,n)
11        continue
          A(m,n)=C
          B(i)=B(i)-A(m,n)*B(n)
12      continue
        B(n)=B(n)/A(1,n)
13    continue

*     Back substitution
      n=NumNP
14    do 15 k=2,MBand
        l=n+k-1
        if(l.gt.NumNP) goto 16
        B(n)=B(n)-A(k,n)*B(l)
15    continue
16    n=n-1
      if(n.gt.0) goto 14
      return
      end

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||