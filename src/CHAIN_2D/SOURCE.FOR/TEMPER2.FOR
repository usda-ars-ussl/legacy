* Source file TEMPER2.FOR ||||||||||||||||||||||||||||||||||||||||||||||

      subroutine Temper(NumNP,NumEl,NumElD,MBand,MBandD,NMat,NumBP,t,A,
     !                  B,Q,x,y,KX,KAT,dt,DS,MatNum,VxO,VzO,VxN,VzN,
     !                  Dispxx,Dispzz,Dispxz,TPar,ThNew,ThOld,Ac,TempO,
     !                  TempN,TBound,KodTB,KXB,Kode,Th3,Th4,Ampl,
     !                  tPeriod,lOrt,A1,B1,NumNPD,IAD,IADN,IADD,VRV,RES,
     !                  RQI,RQ,QQ,QI,RQIDOT,ECNVRG,RCNVRG,ACNVRG,MNorth,
     !                  MaxItO)

      logical lOrt
      double precision A,B,A1,B1,VRV,RES,RQI,RQ,QQ,QI,RQIDOT,ECNVRG,
     !                 RCNVRG,ACNVRG
      dimension A(MBandD,NumNP),B(NumNP),Q(NumNP),KX(NumElD,4),x(NumNP),
     !          MatNum(NumNP),ThNew(NumNP),ThOld(NumNP),TPar(10,NMat),
     !          VxO(NumNP),VzO(NumNP),VxN(NumNP),VzN(NumNP),Ac(NumNP),
     !          y(NumNP),KXB(NumBP),DS(NumNP),TempO(NumNP),TempN(NumNP),
     !          KodTB(NumBP),Dispxx(NumNP),Dispzz(NumNP),Dispxz(NumNP),
     !          Kode(NumNP),TBound(6),S(3,3),Bi(3),Ci(3),List(3),Vx(3),
     !          Vz(3),
     !          A1(MBandD,NumNP),B1(NumNP),RES(NumNP),IAD(MBandD,NumNP),
     !          IADN(NumNP),IADD(NumNP),VRV(NumNP),RQI(NumNPD,MNorth),
     !          RQ(NumNP),QQ(NumNP),RQIDOT(MNorth),QI(NumNPD,MNorth)

*     Initialisation
      epsi=0.5
      alf=1.-epsi
      jjj=MBand
      NLevel=2
      xMul=1.
      do 12 i=1,NumNP
        B(i)=0.
        DS(i)=0.
        if(lOrt) B1(i)=TempN(i)
        M=MatNum(i)
        TempO(i)=TempN(i)
        Ac(i)=-(TPar(1,M)*TPar(8,M)+TPar(2,M)*TPar(9,M)+
     !          TPar(10,M)*(alf*ThOld(i)+epsi*ThNew(i)))
        do 11 j=1,MBandD
          A(j,i)=0.
11      continue
12    continue

      do 19 Level=1,NLevel
        if(Level.eq.NLevel) then
          call DispT(NumNP,NMat,Dispxx,Dispzz,Dispxz,VxN,VzN,ThNew,TPar,
     !               MatNum)
        else
          call DispT(NumNP,NMat,Dispxx,Dispzz,Dispxz,VxO,VzO,ThOld,TPar,
     !               MatNum)
        end if

        do 13 i=1,NumNP
          M=MatNum(i)
          if(Level.ne.NLevel) then
            DPom=dt*TPar(10,M)**2/6./(-Ac(i))
            Dispxx(i)=Dispxx(i)+VxO(i)*VxO(i)*DPom
            Dispzz(i)=Dispzz(i)+VzO(i)*VzO(i)*DPom
            Dispxz(i)=Dispxz(i)+VxO(i)*VzO(i)*DPom
          else
            DPom=dt*TPar(10,M)**2/6./(-Ac(i))
            Dispxx(i)=Dispxx(i)-VxN(i)*VxN(i)*DPom
            Dispzz(i)=Dispzz(i)-VzN(i)*VzN(i)*DPom
            Dispxz(i)=Dispxz(i)-VxN(i)*VzN(i)*DPom
          end if
13      continue

*       Loop on elements
        do 17 n=1,NumEl
          NUS=4
          if(KX(n,3).eq.KX(n,4)) NUS=3

*         Loop on subelements
          do 16 k=1,NUS-2
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
            FMul=xMul*AE/4.
            Ec1=(Dispxx(i)+Dispxx(j)+Dispxx(l))/3.
            Ec2=(Dispxz(i)+Dispxz(j)+Dispxz(l))/3.
            Ec3=(Dispzz(i)+Dispzz(j)+Dispzz(l))/3.
            AcE=(Ac(i)+Ac(j)+Ac(l))/3.
            VxE=(Vx(1)*TPar(10,MatNum(i))+Vx(2)*TPar(10,MatNum(j))+
     !           Vx(3)*TPar(10,MatNum(l)))/3.
            VzE=(Vz(1)*TPar(10,MatNum(i))+Vz(2)*TPar(10,MatNum(j))+
     !           Vz(3)*TPar(10,MatNum(l)))/3.
            SMul1=-1./AE/4.*xMul
            SMul2=AE/20.*xMul
            do 15 j1=1,3
              i1=List(j1)
              if(Level.eq.NLevel) DS(i1)=DS(i1)+FMul*(AcE+Ac(i1)/3.)
              do 14 j2=1,3
                i2=List(j2)
                S(j1,j2)=SMul1*(Ec1*Bi(j1)*Bi(j2)+Ec3*Ci(j1)*Ci(j2)+
     !                         Ec2*(Bi(j1)*Ci(j2)+Ci(j1)*Bi(j2)))
                S(j1,j2)=S(j1,j2)-
     !              (Bi(j2)/8.*(VxE+Vx(j1)*TPar(10,MatNum(j1))/3.)+
     !               Ci(j2)/8.*(VzE+Vz(j1)*TPar(10,MatNum(j1))/3.))*xMul
                if(Level.ne.NLevel) then
                  B(i1)=B(i1)-alf*S(j1,j2)*TempO(i2)
                else
                  if(lOrt) then
                    call Find(i1,i2,kk,NumNP,MBandD,IAD,IADN)
                    iB=kk
                  else
                    iB=MBand+i2-i1
                  end if
                  A(iB,i1)=A(iB,i1)+epsi*S(j1,j2)
                end if
14            continue
15          continue
16        continue
17      continue

        do 18 i=1,NumNP
          if(Level.eq.NLevel) then
            if(lOrt) jjj=IADD(i)
            A(jjj,i)=A(jjj,i)+DS(i)/dt
            B(i)=B(i)+DS(i)/dt*TempO(i)
          end if
18      continue
19    continue

*     Boundary condition
      call T_Bound(NumNP,MBand,MBandD,NumBP,NMat,A,B,Q,TPar,TempO,KXB,
     !             Kode,KodTB,TBound,Th3,Th4,Ampl,t,tPeriod,epsi,lOrt,
     !             IADD)

*     Solve the global matrix equation for transport
      if(lOrt) then
        call ILU     (A,NumNP,MBandD,IAD,IADN,IADD,A1)
        call OrthoMin(A,B1,B,NumNP,MBandD,NumNPD,IAD,IADN,IADD,A1,VRV,
     !                RES,RQI,RQ,QQ,QI,RQIDOT,ECNVRG,RCNVRG,ACNVRG,4,
     !                MNorth,MaxItO)
      else
        call SolveT(A,B,MBand,MBandD,NumNP)
      end if
      do 20 i=1,NumNP
        if(lOrt) B(i)=B1(i)
        TempN(i)=sngl(B(i))
20    continue
      return
      end

************************************************************************

      subroutine T_Bound(NumNP,MBand,MBandD,NumBP,NMat,A,B,Q,TPar,Temp,
     !                   KXB,Kode,KodTB,TBound,Th3,Th4,Ampl,t,tPeriod,
     !                   epsi,lOrt,IADD)

      double precision A,B
      integer TKod
      logical lOrt
      dimension A(MBandD,NumNP),TPar(10,NMat),Temp(NumNP),KXB(NumBP),
     !          Kode(NumNP),KodTB(NumBP),B(NumNP),Q(NumNP),TBound(6),
     !          IADD(NumNP)

      jjj=MBand
      alf=1.-epsi
      do 14 i=1,NumNP
        if(Kode(i).ne.0) then
          In=0
          do 11 j=1,NumBP
            if(KXB(j).eq.i) then
              if(KodTB(j).gt.0) then
                TKod=1
                if(abs(Kode(i)).le.2.or.abs(Kode(i)).ge.5)
     !                                    TBnd=TBound(KodTB(j))
                if(abs(Kode(i)).eq.3) TBnd=Th3
                if(abs(Kode(i)).eq.4) then
                  PI=3.141592654
                  tBnd=Th4
                  if(tPeriod.gt.0.)
     !              tBnd=Th4+Ampl*sin(2.*PI*t/tPeriod-7.*PI/12.)
                end if
              else
                if(Q(i).gt.0.) then
                  TKod=3
                  if(abs(Kode(i)).le.2.or.abs(Kode(i)).ge.5)
     !                                    TBnd=TBound(-KodTB(j))
                  if(abs(Kode(i)).eq.3) TBnd=Th3
                  if(abs(Kode(i)).eq.4) TBnd=Th4
                else
                  TKod=2
                end if
              end if
              goto 12
            end if
11        continue

*         Point source or sink
          In=1
          if(Q(i).lt.0.) then
            TKod=2
          else
            TBnd=TBound(6)
            TKod=3
          end if

12        continue

*         Dirichlet boundary condition
          if(TKod.eq.1) then
            if(lOrt) then
              A(IADD(i),i)=1.d30
              B(i)=1.d30*TBnd
            else
              do 13 j=1,2*MBand-1
                A(j,i)=0.
13            continue
              A(MBand,i)=1.
              B(i)=TBnd
            end if
          end if

*         Cauchy boundary condition
          if(TKod.eq.3) then
            alf=1.-epsi
            B(i)=B(i)-Q(i)*TPar(10,1)*(TBnd-alf*Temp(i))
            if(lOrt) jjj=IADD(i)
            A(jjj,i)=A(jjj,i)-epsi*Q(i)*TPar(10,1)
          end if
        end if
14    continue
      return
      end

************************************************************************

*     Calculate the dispersion coefficient

      Subroutine DispT(NumNP,NMat,Dispxx,Dispzz,Dispxz,Vx,Vz,theta,TPar,
     !                 MatNum)

      dimension Vx(NumNP),Vz(NumNP),theta(NumNP),TPar(10,NMat),
     !          Dispxx(NumNP),Dispzz(NumNP),Dispxz(NumNP),MatNum(NumNP)

      do 11 i=1,NumNP
        M=MatNum(i)
        XLamb=TPar(5,M)+TPar(6,M)*theta(i)+TPar(7,M)*sqrt(theta(i))
        Vabs=sqrt(Vx(i)*Vx(i)+Vz(i)*Vz(i))
        if(Vabs.gt.0.) then
          Dispxx(i)=TPar(3,M)*TPar(10,M)*Vx(i)*Vx(i)/Vabs+
     !              TPar(4,M)*TPar(10,M)*Vz(i)*Vz(i)/Vabs+XLamb
          Dispzz(i)=TPar(3,M)*TPar(10,M)*Vz(i)*Vz(i)/Vabs+
     !              TPar(4,M)*TPar(10,M)*Vx(i)*Vx(i)/Vabs+XLamb
          Dispxz(i)=(TPar(3,M)-TPar(4,M))*TPar(10,M)*Vx(i)*Vz(i)/Vabs
        else
          Dispxx(i)=XLamb
          Dispzz(i)=XLamb
          Dispxz(i)=0.
        end if
11    continue
      return
      end

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
