* Source file SOLUTE3.FOR ||||||||||||||||||||||||||||||||||||||||||||||

      subroutine Solute(NumNP,NumEl,NumElD,MBand,MBandD,NMat,t,Kode,A,B,
     !                  Q,hNew,hOld,F,x,y,z,KX,dt,DS,Sink,MatNum,Con,
     !                  ConO,ConAxx,ConAyy,ConAzz,ConAxy,ConAxz,ConAyz,
     !                  Vx,Vy,Vz,Dispxx,Dispyy,Dispzz,Dispxy,Dispxz,
     !                  Dispyz,ChPar,ThNew,ThOld,thSat,Ac,Fc,Gc,Qc,Conc,
     !                  ListNE,cBound,tPulse,NumBP,KodCB,KXB,NLevel,
     !                  cPrec,crt,cht,lWat,lUpW,WeTab,epsi,CumCh0,
     !                  CumCh1,CumChR,dtMaxC,Peclet,Courant,lArtD,PeCr,
     !                  PeCrMx,
     !                  lOrt,A1,B1,NumNPD,IAD,IADN,IADD,VRV,RES,RQI,RQ,
     !                  QQ,QI,RQIDOT,ECNVRG,RCNVRG,ACNVRG,MNorth,MaxItO)

      double precision A,B,A1,B1,VRV,RES,RQI,RQ,QQ,QI,RQIDOT,ECNVRG,
     !                 RCNVRG,ACNVRG
      logical lWat,lUpW,lOrt,lArtD
      integer e
      dimension A(MBandD,NumNP),B(NumNP),Q(NumNP),hNew(NumNP),F(NumNP),
     !          KX(NumElD,9),MatNum(NumNP),Sink(NumNP),DS(NumNP),
     !          x(NumNP),y(NumNP),z(NumNP),Kode(NumNP),Vx(NumNP),
     !          Vy(NumNP),Vz(NumNP),ThNew(NumNP),ThOld(NumNP),
     !          KXB(NumBP),Conc(NumNP),Con(NumNP),ConAxx(NumEl),
     !          ConAyy(NumEl),ConAzz(NumEl),ConAxy(NumEl),ConAxz(NumEl),
     !          ConAyz(NumEl),Qc(NumNP),Ac(NumNP),Fc(NumNP),Gc(NumNP),
     !          KodCB(NumBP),Dispxx(NumNP),Dispyy(NumNP),Dispzz(NumNP),
     !          Dispxy(NumNP),Dispxz(NumNP),Dispyz(NumNP),ListNE(NumNP),
     !          ChPar(10,NMat),cBound(12),thSat(NMat),S(4,4),Bi(4),
     !          Ci(4),Di(4),List(4),Wx(4),Wy(4),Wz(4),WeTab(6,5*NumEl),
     !          VxE(4),VyE(4),VzE(4),ConO(NumNP),hOld(NumNP),
     !          A1(MBandD,NumNP),B1(NumNP),RES(NumNP),IAD(MBandD,NumNP),
     !          IADN(NumNP),IADD(NumNP),VRV(NumNP),RQI(NumNPD,MNorth),
     !          RQ(NumNP),QQ(NumNP),RQIDOT(MNorth),QI(NumNPD,MNorth)

*     Initialisation
      alf=1.-epsi
      jjj=MBand
      if(t.gt.tPulse) then
        do 11 i=1,10
          cBound(i)=0.
11      continue
      end if
      do 13 i=1,NumNP
        B(i) =0.
        if(lOrt) B1(i)=Conc(i)
        Qc(i)=0.
        if(epsi.lt.0.001) then
          if(lOrt) jjj=IADD(i)
          A(jjj,i)=0.
        else
          do 12 j=1,MBandD
            A(j,i)=0.
12        continue
        end if
13    continue

      do 21 Level=1,NLevel
        if(Level.eq.NLevel) then
          Eps=epsi
          if(lWat)
     !      call Veloc(NumNP,NumEl,NumElD,hNew,x,y,z,KX,ListNE,Con,
     !               ConAxx,ConAyy,ConAzz,ConAxy,ConAxz,ConAyz,Vx,Vy,Vz)
          call Disper(NumNP,NMat,Dispxx,Dispyy,Dispzz,Dispxy,Dispxz,
     !                Dispyz,Vx,Vy,Vz,ThNew,thSat,ChPar,MatNum,lArtD,
     !                PeCr,dt)
          call PeCour(NumNP,NumEl,NumElD,NMat,x,y,z,Vx,Vy,Vz,KX,MatNum,
     !                Dispxx,Dispyy,Dispzz,ChPar,ThNew,dt,dtMaxC,Peclet,
     !                Courant,lUpW,lArtD,PeCr,PeCrMx)
        if(lUpW.and.lWat)
     !     call WeFact(NumNP,NumEl,NumElD,x,y,z,KX,WeTab,Vx,Vy,Vz,
     !                 Dispxx,Dispyy,Dispzz,Dispxy,Dispxz,Dispyz)
        else
          Eps=1.-epsi
          call Disper(NumNP,NMat,Dispxx,Dispyy,Dispzz,Dispxy,Dispxz,
     !                Dispyz,Vx,Vy,Vz,ThNew,thSat,ChPar,MatNum,lArtD,
     !                PeCr,dt)
        end if
        do 14 i=1,NumNP
          M=MatNum(i)
          if(Level.ne.NLevel) then
            if(.not.lArtD.and..not.lUpW) then
              DPom=dt/6./(ThOld(i)+ChPar(1,M)*ChPar(5,M))
              Dispxx(i)=Dispxx(i)+Vx(i)*Vx(i)*DPom
              Dispyy(i)=Dispyy(i)+Vy(i)*Vy(i)*DPom
              Dispzz(i)=Dispzz(i)+Vz(i)*Vz(i)*DPom
              Dispxy(i)=Dispxy(i)+Vx(i)*Vy(i)*DPom
              Dispxz(i)=Dispxz(i)+Vx(i)*Vz(i)*DPom
              Dispyz(i)=Dispyz(i)+Vy(i)*Vz(i)*DPom
            end if
          else
            Ac(i)=-(ThOld(i)*alf+ThNew(i)*epsi)-ChPar(1,M)*ChPar(5,M)
            if(.not.lArtD.and..not.lUpW) then
              DPom=dt/6./(ThNew(i)+ChPar(1,M)*ChPar(5,M))
              Dispxx(i)=Dispxx(i)-Vx(i)*Vx(i)*DPom
              Dispyy(i)=Dispyy(i)-Vy(i)*Vy(i)*DPom
              Dispzz(i)=Dispzz(i)-Vz(i)*Vz(i)*DPom
              Dispxy(i)=Dispxy(i)-Vx(i)*Vy(i)*DPom
              Dispxz(i)=Dispxz(i)-Vx(i)*Vz(i)*DPom
              Dispyz(i)=Dispyz(i)-Vy(i)*Vz(i)*DPom
            end if
            cS=cBound(12)
            if(cS.gt.Conc(i)) cS=Conc(i)
            Gc(i)=ChPar(8,M)*ThNew(i)+ChPar(1,M)*ChPar(9,M)-
     !            Sink(i)*cS
            Fc(i)=ChPar(6,M)*ThNew(i)+ChPar(1,M)*ChPar(7,M)*ChPar(5,M)+
     !            Sink(i)
          end if
14      continue
        do 15 i=1,NumNP
          F(i)=0.
          if(Level.eq.NLevel) Ds(i)=0.
15      continue

*       Loop on elements
        NumSEl=0
        do 19 e=1,NumEl
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
          do 18 n=1,NUS
            NumSEl=NumSEl+1
            call Elem(NUS,NumElD,e,n,KX,i,j,k,l)
            List(1)=i
            List(2)=j
            List(3)=k
            List(4)=l
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

*           Calculate Velocities
            AA=1./VE/6.
            Ai=CAxx*Bi(1)+CAxy*Ci(1)+CAxz*Di(1)
            Aj=CAxx*Bi(2)+CAxy*Ci(2)+CAxz*Di(2)
            Ak=CAxx*Bi(3)+CAxy*Ci(3)+CAxz*Di(3)
            Al=CAxx*Bi(4)+CAxy*Ci(4)+CAxz*Di(4)
            if(Level.eq.NLevel) then
              Vxx=AA*(Ai*hNew(i)+Aj*hNew(j)+Ak*hNew(k)+Al*hNew(l))+CAxz
            else
              Vxx=AA*(Ai*hOld(i)+Aj*hOld(j)+Ak*hOld(k)+Al*hOld(l))+CAxz
            end if
            Ai=CAxy*Bi(1)+CAyy*Ci(1)+CAyz*Di(1)
            Aj=CAxy*Bi(2)+CAyy*Ci(2)+CAyz*Di(2)
            Ak=CAxy*Bi(3)+CAyy*Ci(3)+CAyz*Di(3)
            Al=CAxy*Bi(4)+CAyy*Ci(4)+CAyz*Di(4)
            if(Level.eq.NLevel) then
              Vyy=AA*(Ai*hNew(i)+Aj*hNew(j)+Ak*hNew(k)+Al*hNew(l))+CAyz
            else
              Vyy=AA*(Ai*hOld(i)+Aj*hOld(j)+Ak*hOld(k)+Al*hOld(l))+CAyz
            end if
            Ai=CAxz*Bi(1)+CAyz*Ci(1)+CAzz*Di(1)
            Aj=CAxz*Bi(2)+CAyz*Ci(2)+CAzz*Di(2)
            Ak=CAxz*Bi(3)+CAyz*Ci(3)+CAzz*Di(3)
            Al=CAxz*Bi(4)+CAyz*Ci(4)+CAzz*Di(4)
            if(Level.eq.NLevel) then
              Vzz=AA*(Ai*hNew(i)+Aj*hNew(j)+Ak*hNew(k)+Al*hNew(l))+CAzz
            else
              Vzz=AA*(Ai*hOld(i)+Aj*hOld(j)+Ak*hOld(k)+Al*hOld(l))+CAzz
            end if
            if(Level.ne.NLevel) then
              ConE=(ConO(i)+ConO(j)+ConO(k)+ConO(l))/4.
              VxE(1)=-ConO(i)*Vxx
              VxE(2)=-ConO(j)*Vxx
              VxE(3)=-ConO(k)*Vxx
              VxE(4)=-ConO(l)*Vxx
              VyE(1)=-ConO(i)*Vyy
              VyE(2)=-ConO(j)*Vyy
              VyE(3)=-ConO(k)*Vyy
              VyE(4)=-ConO(l)*Vyy
              VzE(1)=-ConO(i)*Vzz
              VzE(2)=-ConO(j)*Vzz
              VzE(3)=-ConO(k)*Vzz
              VzE(4)=-ConO(l)*Vzz
            else
              ConE=(Con(i)+Con(j)+Con(k)+Con(l))/4.
              VxE(1)=-Con(i)*Vxx
              VxE(2)=-Con(j)*Vxx
              VxE(3)=-Con(k)*Vxx
              VxE(4)=-Con(l)*Vxx
              VyE(1)=-Con(i)*Vyy
              VyE(2)=-Con(j)*Vyy
              VyE(3)=-Con(k)*Vyy
              VyE(4)=-Con(l)*Vyy
              VzE(1)=-Con(i)*Vzz
              VzE(2)=-Con(j)*Vzz
              VzE(3)=-Con(k)*Vzz
              VzE(4)=-Con(l)*Vzz
            end if
            VxEE=-ConE*Vxx
            VyEE=-ConE*Vyy
            VzEE=-ConE*Vzz

            if(Level.eq.1) then
              cS=cBound(12)
              if(cBound(12).gt.conc(i)) cS=cS+(conc(i)-cBound(12))/4.
              if(cBound(12).gt.conc(j)) cS=cS+(conc(j)-cBound(12))/4.
              if(cBound(12).gt.conc(k)) cS=cS+(conc(k)-cBound(12))/4.
              if(cBound(12).gt.conc(l)) cS=cS+(conc(l)-cBound(12))/4.
              RootCh=VE*dt*cS*(Sink(i)+Sink(j)+Sink(k)+Sink(l))/4.
              CumCh0=CumCh0-VE*dt*(Gc(i)+Gc(j)+Gc(k)+Gc(l))/4.+RootCh
              CumCh1=CumCh1-VE*dt*((Fc(i)-Sink(i))*conc(i)+
     !               (Fc(j)-Sink(j))*conc(j)+(Fc(k)-Sink(k))*conc(k)+
     !               (Fc(l)-Sink(l))*conc(l))/4.
              CumChR=CumChR+RootCh
            end if
            FMul=VE/5.
            GcE=(Gc(i)+Gc(j)+Gc(k)+Gc(l))/4.
            Ec1=(Dispxx(i)+Dispxx(j)+Dispxx(k)+Dispxx(l))/4.
            Ec2=(Dispyy(i)+Dispyy(j)+Dispyy(k)+Dispyy(l))/4.
            Ec3=(Dispzz(i)+Dispzz(j)+Dispzz(k)+Dispzz(l))/4.
            Ec4=(Dispxy(i)+Dispxy(j)+Dispxy(k)+Dispxy(l))/4.
            Ec5=(Dispxz(i)+Dispxz(j)+Dispxz(k)+Dispxz(l))/4.
            Ec6=(Dispyz(i)+Dispyz(j)+Dispyz(k)+Dispyz(l))/4.
            if(Level.eq.NLevel) AcE=(Ac(i)+Ac(j)+Ac(k)+Ac(l))/4.
            FcE=(Fc(i)+Fc(j)+Fc(k)+Fc(l))/4.
            SMul1=-1./VE/36.
            SMul2=VE/30.
            if(lUpW) then
              W12=WeTab(1,NumSEl)
              W13=WeTab(2,NumSEl)
              W14=WeTab(3,NumSEl)
              W23=WeTab(4,NumSEl)
              W24=WeTab(5,NumSEl)
              W34=WeTab(6,NumSEl)
              A11=-2.*W12+2.*W14+2.*W13
              A12=-2.*W12+W14+W13
              A13=-W12+W14+2.*W13
              A14=-W12+2.*W14+W13
              A21=-W23+2.*W12+W24
              A22=-2.*W23+2.*W12+2.*W24
              A23=-2.*W23+W12+W24
              A24=-W23+W12+2.*W24
              A31=-W34+W23-2.*W13
              A32=-W34+2.*W23-W13
              A33=-2.*W34+2.*W23-2.*W13
              A34=-2.*W34+W23-W13
              A41=-2.*W14+W34-W24
              A42=-W14+W34-2.*W24
              A43=-W14+2.*W34-W24
              A44=-2.*W14+2.*W34-2.*W24
              Wx(1)=VxE(1)*A11+VxE(2)*A12+VxE(3)*A13+VxE(4)*A14
              Wx(2)=VxE(1)*A21+VxE(2)*A22+VxE(3)*A23+VxE(4)*A24
              Wx(3)=VxE(1)*A31+VxE(2)*A32+VxE(3)*A33+VxE(4)*A34
              Wx(4)=VxE(1)*A41+VxE(2)*A42+VxE(3)*A43+VxE(4)*A44
              Wy(1)=VyE(1)*A11+VyE(2)*A12+VyE(3)*A13+VyE(4)*A14
              Wy(2)=VyE(1)*A21+VyE(2)*A22+VyE(3)*A23+VyE(4)*A24
              Wy(3)=VyE(1)*A31+VyE(2)*A32+VyE(3)*A33+VyE(4)*A34
              Wy(4)=VyE(1)*A41+VyE(2)*A42+VyE(3)*A43+VyE(4)*A44
              Wz(1)=VzE(1)*A11+VzE(2)*A12+VzE(3)*A13+VzE(4)*A14
              Wz(2)=VzE(1)*A21+VzE(2)*A22+VzE(3)*A23+VzE(4)*A24
              Wz(3)=VzE(1)*A31+VzE(2)*A32+VzE(3)*A33+VzE(4)*A34
              Wz(4)=VzE(1)*A41+VzE(2)*A42+VzE(3)*A43+VzE(4)*A44
            end if
            do 17 j1=1,4
              i1=List(j1)
              F(i1)=F(i1)+FMul*(GcE+Gc(i1)/4.)
              if(Level.eq.NLevel) DS(i1)=DS(i1)+FMul*(AcE+Ac(i1)/4.)
              iBound=0
              if(Kode(i).ne.0) then
                do 24 id=1,NumBP
                  if(KXB(id).eq.i1.and.KodCB(id).gt.0) iBound=1
24              continue
              end if
              if(iBound.eq.1) Qc(i1)=Qc(i1)-Eps*FMul*(GcE+Gc(i1)/4.)
              do 16 j2=1,4
                i2=List(j2)
                S(j1,j2)=SMul1*(Ec1*Bi(j1)*Bi(j2)+Ec2*Ci(j1)*Ci(j2)+
     !                          Ec3*Di(j1)*Di(j2)+
     !                          Ec4*(Bi(j1)*Ci(j2)+Ci(j1)*Bi(j2))+
     !                          Ec5*(Bi(j1)*Di(j2)+Di(j1)*Bi(j2))+
     !                          Ec6*(Ci(j1)*Di(j2)+Di(j1)*Ci(j2)))
                S(j1,j2)=S(j1,j2)-(Bi(j2)/30.*(VxEE+VxE(j1)/4.)+
     !                             Ci(j2)/30.*(VyEE+VyE(j1)/4.)+
     !                             Di(j2)/30.*(VzEE+VzE(j1)/4.))
                if(lUpW) S(j1,j2)=S(j1,j2)-(Bi(j2)/240.*Wx(j1)+
     !                            Ci(j2)/240.*Wy(j1)+Di(j2)/240.*Wz(j1))
                ic=1
                if(i1.eq.i2) ic=2
                S(j1,j2)=S(j1,j2)+SMul2*ic*(FcE+(Fc(i1)+Fc(i2))/4.)
                if(Level.ne.NLevel) then
                  B(i1)=B(i1)-alf*S(j1,j2)*Conc(i2)
                else
                  if(lOrt) then
                    call Find(i1,i2,kk,NumNP,MBandD,IAD,IADN)
                    iB=kk
                  else
                    iB=MBand+i2-i1
                  end if
                  A(iB,i1)=A(iB,i1)+epsi*S(j1,j2)
                end if
                if(iBound.eq.1) Qc(i1)=Qc(i1)-Eps*S(j1,j2)*Conc(i2)
16            continue
17          continue
18        continue
19      continue

        do 20 i=1,NumNP
          M=MatNum(i)
          if(Level.ne.NLevel) then
            B(i)=B(i)-alf*F(i)
          else
            if(lOrt) jjj=IADD(i)
            A(jjj,i)=A(jjj,i)+DS(i)/dt
            B(i)=B(i)+DS(i)/dt*Conc(i)-epsi*F(i)
          end if
20      continue
21    continue

*     Boundary condition
      call c_Bound(NumNP,MBand,MBandD,NumBP,A,B,Q,Qc,Conc,Kode,KXB,
     !             KodCB,cBound,cPrec,crt,cht,epsi,dt,DS,lOrt,IADD)

*     Solve the global matrix equation for transport
      if(epsi.lt.0.001) then
        do 22 i=1,NumNP
          if(lOrt) jjj=IADD(i)
          B(i)=B(i)/A(jjj,i)
22      continue
      else
        if(lOrt) then
          call ILU     (A,NumNP,MBandD,IAD,IADN,IADD,A1)
          North=4
          call OrthoMin(A,B1,B,NumNP,MBandD,NumNPD,IAD,IADN,IADD,A1,VRV,
     !                  RES,RQI,RQ,QQ,QI,RQIDOT,ECNVRG,RCNVRG,ACNVRG,
     !                  North,MNorth,MaxItO)
        else
          call SolveT(A,B,MBand,MBandD,NumNP)
        end if
      end if
      do 23 i=1,NumNP
        if(lOrt) B(i)=B1(i)
        Conc(i)=sngl(B(i))
23    continue
      return
      end

************************************************************************

      subroutine c_Bound(NumNP,MBand,MBandD,NumBP,A,B,Q,Qc,Conc,Kode,
     !                   KXB,KodCB,cBound,cPrec,crt,cht,epsi,dt,DS,lOrt,
     !                   IADD)

      double precision A,B
      integer cKod
      logical lOrt
      dimension A(MBandD,NumNP),B(NumNP),Q(NumNP),Conc(NumNP),
     !          Qc(NumNP),Kode(NumNP),KXB(NumBP),KodCB(NumBP),
     !          cBound(12),DS(NumNP),IADD(NumNP)

      alf=1.-epsi
      jjj=MBand
      do 14 i=1,NumNP
        if(Kode(i).ne.0) then
          do 11 j=1,NumBP
            if(KXB(j).eq.i) then
              if(KodCB(j).gt.0) then
                cKod=1
                if(abs(Kode(i)).le.2.or.abs(Kode(i)).ge.5)
     !                                    cBnd=cBound(KodCB(j))
                if(abs(Kode(i)).eq.3) cBnd=cht
                if(abs(Kode(i)).eq.4) cBnd=cPrec
              else
                if(Q(i).gt.0.) then
                  cKod=3
                  if(abs(Kode(i)).eq.1.or.abs(Kode(i)).ge.5)
     !                                    cBnd=cBound(-KodCB(j))
                  if(abs(Kode(i)).eq.3) cBnd=crt
                  if(abs(Kode(i)).eq.4) cBnd=cPrec
                else
                  cKod=2
                  if(Kode(i).eq.-4) then
                    cKod=3
                    cBnd=0.
                  end if
                end if
              end if
              if(abs(Kode(i)).eq.2) cKod=2
              goto 12
            end if
11        continue

*    point source or sink
          if(Q(i).lt.0.) then
            cKod=2
          else
            cBnd=cBound(11)
            cKod=3
          end if

12        continue

*     Dirichlet boundary condition
          if(cKod.eq.1) then
            Qc(i)=Qc(i)+Q(i)*(epsi*cBnd+alf*Conc(i))-
     !                        DS(i)*(cBnd-Conc(i))/dt
            if(lOrt) then
              A(IADD(i),i)=1.d30
              B(i)=1.d30*cBnd
            else
              do 13 j=1,2*MBand-1
                A(j,i)=0.
13            continue
              A(MBand,i)=1.
              B(i)=cBnd
            end if
          end if

*     Neumann boundary condition
          if(cKod.eq.2) then
            Qc(i)=Q(i)*Conc(i)
          end if

*     Cauchy boundary condition
          if(cKod.eq.3) then
            B(i)=B(i)-Q(i)*(cBnd-alf*Conc(i))
            if(lOrt) jjj=IADD(i)
            A(jjj,i)=A(jjj,i)-epsi*Q(i)
            Qc(i)=Q(i)*cBnd
          end if

        end if
14    continue
      return
      end

************************************************************************

*     Initial values for solute transport calculation

      subroutine ChInit(NumNP,NumEl,NumElD,NMat,x,y,z,KX,MatNum,NLevel,
     !                  Con,hNew,Sink,cSink,Vx,Vy,Vz,ConAxx,ConAyy,
     !                  ConAzz,ConAxy,ConAxz,ConAyz,Dispxx,Dispyy,
     !                  Dispzz,Dispxy,Dispxz,Dispyz,ChPar,theta,thSat,
     !                  Conc,Fc,Gc,ListNE,lUpW,WeTab,dt,dtMaxC,Peclet,
     !                  Courant,lArtD,PeCr,PeCrMx,ConO)

      logical lUpW,lArtD
      dimension hNew(NumNP),x(NumNP),y(NumNP),z(NumNP),KX(NumElD,9),
     !          theta(NumNP),Sink(NumNP),ChPar(10,NMat),Vx(NumNP),
     !          Vy(NumNP),Vz(NumNP),MatNum(NumNP),Con(NumNP),
     !          ConAxx(NumEl),ConAyy(NumEl),ConAzz(NumEl),ConAxy(NumEl),
     !          ConAxz(NumEl),ConAyz(NumEl),Dispxx(NumNP),Dispyy(NumNP),
     !          Dispzz(NumNP),Dispxy(NumNP),Dispxz(NumNP),Dispyz(NumNP),
     !          ListNE(NumNP),Fc(NumNP),Gc(NumNP),WeTab(6,5*NumEl),
     !          thSat(NMat),Conc(NumNP),ConO(NumNP)

      do 11 i=1,NumNP
        M=MatNum(i)
        if(NLevel.eq.2) then
          cS=cSink
          if(cS.gt.Conc(i)) cS=Conc(i)
          Gc(i)=ChPar(8,M)*theta(i)+ChPar(1,M)*ChPar(9,M)-Sink(i)*cS
          Fc(i)=ChPar(6,M)*theta(i)+ChPar(1,M)*ChPar(7,M)*ChPar(5,M)+
     !          Sink(i)
        end if
        ConO(i)=Con(i)
11    continue
      call Veloc(NumNP,NumEl,NumElD,hNew,x,y,z,KX,ListNE,Con,
     !           ConAxx,ConAyy,ConAzz,ConAxy,ConAxz,ConAyz,Vx,Vy,Vz)
      call Disper(NumNP,NMat,Dispxx,Dispyy,Dispzz,Dispxy,Dispxz,Dispyz,
     !            Vx,Vy,Vz,theta,thSat,ChPar,MatNum,lArtD,PeCr,dt)
      call PeCour(NumNP,NumEl,NumElD,NMat,x,y,z,Vx,Vy,Vz,KX,MatNum,
     !            Dispxx,Dispyy,Dispzz,ChPar,theta,dt,dtMaxC,Peclet,
     !            Courant,lUpW,lArtD,PeCr,PeCrMx)
      if(lUpW)
     !  call WeFact(NumNP,NumEl,NumElD,x,y,z,KX,WeTab,Vx,Vy,Vz,Dispxx,
     !              Dispyy,Dispzz,Dispxy,Dispxz,Dispyz)
      return
      end

************************************************************************

*     Calculate  velocities

      subroutine Veloc(NumNP,NumEl,NumElD,hNew,x,y,z,KX,ListNE,Con,
     !                 ConAxx,ConAyy,ConAzz,ConAxy,ConAxz,ConAyz,Vx,Vy,
     !                 Vz)

      integer e
      dimension hNew(NumNP),x(NumNP),y(NumNP),z(NumNP),ListNE(NumNP),
     !          Con(NumNP),KX(NumElD,9),Vx(NumNP),Vy(NumNP),Vz(NumNP),
     !          ConAxx(NumEl),ConAyy(NumEl),ConAzz(NumEl),ConAxy(NumEl),
     !          ConAxz(NumEl),ConAyz(NumEl),List(4)

      do 11 i=1,NumNP
        Vx(i)=0.
        Vy(i)=0.
        Vz(i)=0.
11    continue
      do 14 e=1,NumEl
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
        do 13 n=1,NUS
          call Elem(NUS,NumElD,e,n,KX,i,j,k,l)
          List(1)=i
          List(2)=j
          List(3)=k
          List(4)=l
          Bi=-(y(k)-y(j))*(z(l)-z(j))+(y(l)-y(j))*(z(k)-z(j))
          Bj=+(y(l)-y(k))*(z(i)-z(k))-(y(i)-y(k))*(z(l)-z(k))
          Bk=-(y(i)-y(l))*(z(j)-z(l))+(y(j)-y(l))*(z(i)-z(l))
          Bl=+(y(j)-y(i))*(z(k)-z(i))-(y(k)-y(i))*(z(j)-z(i))
          Ci=+(x(k)-x(j))*(z(l)-z(j))-(x(l)-x(j))*(z(k)-z(j))
          Cj=-(x(l)-x(k))*(z(i)-z(k))+(x(i)-x(k))*(z(l)-z(k))
          Ck=+(x(i)-x(l))*(z(j)-z(l))-(x(j)-x(l))*(z(i)-z(l))
          Cl=-(x(j)-x(i))*(z(k)-z(i))+(x(k)-x(i))*(z(j)-z(i))
          Di=-(x(k)-x(j))*(y(l)-y(j))+(x(l)-x(j))*(y(k)-y(j))
          Dj=+(x(l)-x(k))*(y(i)-y(k))-(x(i)-x(k))*(y(l)-y(k))
          Dk=-(x(i)-x(l))*(y(j)-y(l))+(x(j)-x(l))*(y(i)-y(l))
          Dl=+(x(j)-x(i))*(y(k)-y(i))-(x(k)-x(i))*(y(j)-y(i))
          Det=(x(l)-x(i))*Bl+(y(l)-y(i))*Cl+(z(l)-z(i))*Dl
          VE=Det/6.
          A=1./VE/6.
          Ai=CAxx*Bi+CAxy*Ci+CAxz*Di
          Aj=CAxx*Bj+CAxy*Cj+CAxz*Dj
          Ak=CAxx*Bk+CAxy*Ck+CAxz*Dk
          Al=CAxx*Bl+CAxy*Cl+CAxz*Dl
          Vxx=A*(Ai*hNew(i)+Aj*hNew(j)+Ak*hNew(k)+Al*hNew(l))+CAxz
          Ai=CAxy*Bi+CAyy*Ci+CAyz*Di
          Aj=CAxy*Bj+CAyy*Cj+CAyz*Dj
          Ak=CAxy*Bk+CAyy*Ck+CAyz*Dk
          Al=CAxy*Bl+CAyy*Cl+CAyz*Dl
          Vyy=A*(Ai*hNew(i)+Aj*hNew(j)+Ak*hNew(k)+Al*hNew(l))+CAyz
          Ai=CAxz*Bi+CAyz*Ci+CAzz*Di
          Aj=CAxz*Bj+CAyz*Cj+CAzz*Dj
          Ak=CAxz*Bk+CAyz*Ck+CAzz*Dk
          Al=CAxz*Bl+CAyz*Cl+CAzz*Dl
          Vzz=A*(Ai*hNew(i)+Aj*hNew(j)+Ak*hNew(k)+Al*hNew(l))+CAzz
          do 12 m=1,4
            ii=List(m)
            Vx(ii)=Vx(ii)-Con(ii)*Vxx
            Vy(ii)=Vy(ii)-Con(ii)*Vyy
            Vz(ii)=Vz(ii)-Con(ii)*Vzz
12        continue
13      continue
14    continue
      do 15 i=1,NumNP
        Vx(i)=Vx(i)/ListNE(i)
        Vy(i)=Vy(i)/ListNE(i)
        Vz(i)=Vz(i)/ListNE(i)
15    continue
      return
      end

************************************************************************

*     Calculate the dispersion coefficient

      subroutine Disper(NumNP,NMat,Dispxx,Dispyy,Dispzz,Dispxy,Dispxz,
     !                  Dispyz,Vx,Vy,Vz,theta,thSat,ChPar,MatNum,lArtD,
     !                  PeCr,dt)

      logical lArtD
      dimension Vx(NumNP),Vy(NumNP),Vz(NumNP),theta(NumNP),thSat(NMat),
     !          Dispxx(NumNP),Dispyy(NumNP),Dispzz(NumNP),MatNum(NumNP),
     !          Dispxy(NumNP),Dispxz(NumNP),Dispyz(NumNP),ChPar(10,NMat)

      do 11 i=1,NumNP
        M=MatNum(i)
        Tau=theta(i)**(7./3.)/thSat(M)**2
        Vabs=sqrt(Vx(i)*Vx(i)+Vy(i)*Vy(i)+Vz(i)*Vz(i))
        Dif=theta(i)*ChPar(2,M)*Tau
        DispL=ChPar(3,M)
        DispT=ChPar(4,M)
        if(lArtD.and.VAbs.gt.0.) DispL=amax1(DispL,Vabs*dt
     !                 /(theta(i)+ChPar(1,M)*ChPar(5,M))/PeCr-Dif/Vabs)
        Dispxx(i)=Dif
        Dispyy(i)=Dif
        Dispzz(i)=Dif
        Dispxy(i)=0.
        Dispxz(i)=0.
        Dispyz(i)=0.
        if(Vabs.gt.0.) then
          Dispxx(i)=DispL*Vx(i)*Vx(i)/Vabs+DispT*Vy(i)*Vy(i)/Vabs+
     !                                     DispT*Vz(i)*Vz(i)/Vabs+Dif
          Dispyy(i)=DispL*Vy(i)*Vy(i)/Vabs+DispT*Vx(i)*Vx(i)/Vabs+
     !                                     DispT*Vz(i)*Vz(i)/Vabs+Dif
          Dispzz(i)=DispL*Vz(i)*Vz(i)/Vabs+DispT*Vx(i)*Vx(i)/Vabs+
     !                                     DispT*Vy(i)*Vy(i)/Vabs+Dif
          Dispxy(i)=(DispL-DispT)*Vx(i)*Vy(i)/Vabs
          Dispxz(i)=(DispL-DispT)*Vx(i)*Vz(i)/Vabs
          Dispyz(i)=(DispL-DispT)*Vy(i)*Vz(i)/Vabs
        end if
11    continue
      return
      end


************************************************************************

*     Calculate upstream weighing factors

      subroutine WeFact(NumNP,NumEl,NumElD,x,y,z,KX,WeTab,Vx,Vy,Vz,
     !                  Dispxx,Dispyy,Dispzz,Dispxy,Dispxz,Dispyz)

      integer e
      dimension x(NumNP),y(NumNP),z(NumNP),KX(NumElD,9),Vx(NumNP),
     !          Vy(NumNP),Vz(NumNP),Dispxx(NumNP),Dispyy(NumNP),
     !          Dispzz(NumNP),Dispxy(NumNP),Dispxz(NumNP),Dispyz(NumNP),
     !          WeTab(6,5*NumEl),List(4)

c      TanH(zz)=(exp(zz)-exp(-zz))/(exp(zz)+exp(-zz))

      NumSEl=0
      do 13 e=1,NumEl
        NUS=5
        if(KX(e,7).eq.0) NUS=3
        if(KX(e,5).eq.0) NUS=1
        do 12 n=1,NUS
          call Elem(NUS,NumElD,e,n,KX,M1,M2,M3,M4)
          NumSEl=NumSEl+1
          List(1)=M1
          List(2)=M2
          List(3)=M3
          List(4)=M4
          do 11 j=1,6
            WeTab(j,NumSEl)=0.
            if(j.eq.1) then
              M1=List(1)
              M2=List(2)
            else if(j.eq.2) then
              M1=List(3)
              M2=List(1)
            else if(j.eq.3) then
              M1=List(4)
              M2=List(1)
            else if(j.le.4) then
              M1=List(2)
              M2=List(3)
            else if(j.le.5) then
              M1=List(4)
              M2=List(2)
            else
              M1=List(3)
              M2=List(4)
            end if
            A=x(M2)-x(M1)
            B=y(M2)-y(M1)
            C=z(M2)-z(M1)
            ALeng=sqrt(A**2+B**2+C**2)
            cosx=1./Aleng*A
            cosy=1./Aleng*B
            cosz=1./Aleng*C
            Vxx=(Vx(M1)+Vx(M2))/2.
            Vyy=(Vy(M1)+Vy(M2))/2.
            Vzz=(Vz(M1)+Vz(M2))/2.
            VV=sqrt(Vxx*Vxx+Vyy*Vyy+Vzz*Vzz)
            if(abs(VV).lt.1.e-30) goto 11
            cosDel=1/VV/Aleng*(A*Vxx+B*Vyy+C*Vzz)
            Delta=acos(abs(cosDel))
            if(NUS.ne.5.and.Delta.gt.0.314) goto 11
            Val=Vxx*cosx+Vyy*cosy+Vzz*cosz
            Dxx=(Dispxx(M1)+Dispxx(M2))/2.
            Dyy=(Dispyy(M1)+Dispyy(M2))/2.
            Dzz=(Dispzz(M1)+Dispzz(M2))/2.
            Dxy=(Dispxy(M1)+Dispxy(M2))/2.
            Dxz=(Dispxz(M1)+Dispxz(M2))/2.
            Dyz=(Dispyz(M1)+Dispyz(M2))/2.
            DAL=abs(Dxx*cosx**2+Dyy*cosy**2+Dzz*cosz**2+
     !              2.*Dxy*cosx*cosy+2.*Dxz*cosx*cosz+2.*Dyz*cosy*cosz)
            Vel=VAL*ALeng
            Disp=2.0*DAL
            aa=11.
            if(Disp.gt.0.) aa=abs(Vel/Disp)
            if(Disp.lt.1.e-30.or.abs(Vel).lt.0.001*VV.or.
     !                           abs(aa).gt.10.) then
              if(abs(Vel).lt.0.001*VV) WeTab(j,NumSEl)=0.0
              if(Vel.gt.0.001*VV) WeTab(j,NumSEl)=1.0
              if(Vel.lt.-0.001*VV) WeTab(j,NumSEl)=-1.0
            else
              WeTab(j,NumSEl)=1.0/TanH(Vel/Disp)-Disp/Vel
            end if
11        continue
12      continue
13    continue
      return
      end

*************************************************************************

*     Calculate the maximum local Peclet and Courant numbers

      subroutine PeCour(NumNP,NumEl,NumElD,NMat,x,y,z,Vx,Vy,Vz,KX,
     !                  MatNum,Dispxx,Dispyy,Dispzz,ChPar,theta,dt,
     !                  dtMaxC,Peclet,Courant,lUpW,lArtD,PeCr,PeCrMx)

      logical lUpW,lArtD
      integer e
      dimension KX(NumElD,9),x(NumNP),y(NumNP),z(NumNP),Vx(NumNP),
     !          Vy(NumNP),Vz(NumNP),MatNum(NumNP),Dispxx(NumNP),
     !          Dispyy(NumNP),Dispzz(NumNP),theta(NumNP),
     !          ChPar(10,NMat)

      Peclet=0.
      Courant=0.
      PeCrMx=0.
      dtMaxC=1.e+30
      do 12 e=1,NumEl
        NUS=5
        if(KX(e,7).eq.0) NUS=3
        if(KX(e,5).eq.0) NUS=1
        do 11 n=1,NUS
          call Elem(NUS,NumElD,e,n,KX,i,j,k,l)
          PecX=99999.
          PecY=99999.
          PecZ=99999.
          dt1=1.e+30
          dt2=1.e+30
          dt3=1.e+30
          xmax=amax1(x(i),x(j),x(k),x(l))
          xmin=amin1(x(i),x(j),x(k),x(l))
          ymax=amax1(y(i),y(j),y(k),y(l))
          ymin=amin1(y(i),y(j),y(k),y(l))
          zmax=amax1(z(i),z(j),z(k),z(l))
          zmin=amin1(z(i),z(j),z(k),z(l))
          delX=xmax-xmin
          delY=ymax-ymin
          delZ=zmax-zmin
          DxE=(Dispxx(i)+Dispxx(j)+Dispxx(k)+Dispxx(l))/4.
          DyE=(Dispyy(i)+Dispyy(j)+Dispyy(k)+Dispyy(l))/4.
          DzE=(Dispzz(i)+Dispzz(j)+Dispzz(k)+Dispzz(l))/4.
          VxE=abs(Vx(i)+Vx(j)+Vx(k)+Vx(l))/4.
          VyE=abs(Vy(i)+Vy(j)+Vy(k)+Vy(l))/4.
          VzE=abs(Vz(i)+Vz(j)+Vz(k)+Vz(l))/4.
          if(DxE.gt.0.) PecX=VxE*delX/DxE
          if(DyE.gt.0.) PecY=VyE*delY/DyE
          if(DzE.gt.0.) PecZ=VzE*delZ/DzE
          if(PecX.ne.99999.) Peclet=amax1(Peclet,PecX)
          if(PecY.ne.99999.) Peclet=amax1(Peclet,PecY)
          if(PecZ.ne.99999.) Peclet=amax1(Peclet,PecZ)
          Peclet=amin1(Peclet,99999.)

          VxMax=amax1(abs(Vx(i))/theta(i),abs(Vx(j))/theta(j),abs(Vx(k))
     !                /theta(k),abs(Vx(l))/theta(l))
          VyMax=amax1(abs(Vy(i))/theta(i),abs(Vy(j))/theta(j),abs(Vy(k))
     !                /theta(k),abs(Vy(l))/theta(l))
          VzMax=amax1(abs(Vz(i))/theta(i),abs(Vz(j))/theta(j),abs(Vz(k))
     !                /theta(k),abs(Vz(l))/theta(l))
          R1=1.+ChPar(1,MatNum(i))*ChPar(5,MatNum(i))/theta(i)
          R2=1.+ChPar(1,MatNum(j))*ChPar(5,MatNum(j))/theta(j)
          R3=1.+ChPar(1,MatNum(k))*ChPar(5,MatNum(k))/theta(k)
          R4=1.+ChPar(1,MatNum(l))*ChPar(5,MatNum(l))/theta(l)
          RMin=amin1(R1,R2,R3,R4)
          CourX=VxMax*dt/delX/RMin
          CourY=VyMax*dt/delY/RMin
          CourZ=VzMax*dt/delZ/RMin
          Courant=amax1(Courant,CourX,CourY,CourZ)

          Cour1=1.0
          Cour2=1.0
          Cour3=1.0
          if(.not.lUpW.and..not.lArtD) then
            if(PecX.ne.99999.) Cour1=amin1(1.,PeCr/amax1(0.5,PecX))
            if(PecY.ne.99999.) Cour2=amin1(1.,PeCr/amax1(0.5,PecY))
            if(PecZ.ne.99999.) Cour3=amin1(1.,PeCr/amax1(0.5,PecZ))
          end if
          if(VxMax.gt.0.) dt1=Cour1*delX*RMin/VxMax
          if(VyMax.gt.0.) dt2=Cour2*delY*RMin/VyMax
          if(VzMax.gt.0.) dt3=Cour3*delZ*RMin/VzMax
          dtMaxC=amin1(dtMaxC,dt1,dt2,dt3)

          PeCrMx=amax1(PeCrMx,amax1(PecX*CourX,PecY*CourY,PecZ*CourZ))

11      continue
12    continue
      return
      end

************************************************************************

*     Solve the global matrix equation for transport

      subroutine SolveT(A,B,MBand,MBandD,NumNP)

      double precision A,B,P,C,Sum
      dimension A(MBandD,NumNP),B(NumNP)

      N1=NumNP-1
      do 12 k=1,N1
        P=1./A(MBand,k)
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