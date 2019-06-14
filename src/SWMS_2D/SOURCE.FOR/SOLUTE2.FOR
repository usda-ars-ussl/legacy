* Source file SOLUTE2.FOR ||||||||||||||||||||||||||||||||||||||||||||||

      subroutine Solute(NumNP,NumEl,NumElD,MBand,MBandD,NMat,t,Kode,A,B,
     !                  Q,hNew,hOld,F,x,y,KX,KAT,dt,DS,Sink,MatNum,Con,
     !                  ConO,ConAxx,ConAzz,ConAxz,Vx,Vz,Dispxx,Dispzz,
     !                  Dispxz,ChPar,ThNew,ThOld,thSat,Ac,Fc,Gc,Qc,Conc,
     !                  ListNE,cBound,tPulse,NumBP,KodCB,KXB,NLevel,
     !                  cPrec,crt,cht,lWat,lUpW,WeTab,epsi,CumCh0,
     !                  CumCh1,CumChR,dtMaxC,Peclet,Courant,lArtD,PeCr,
     !                  lOrt,A1,B1,NumNPD,IAD,IADN,IADD,VRV,RES,RQI,RQ,
     !                  QQ,QI,RQIDOT,ECNVRG,RCNVRG,ACNVRG,MNorth,MaxItO)

      double precision A,B,A1,B1,VRV,RES,RQI,RQ,QQ,QI,RQIDOT,ECNVRG,
     !                 RCNVRG,ACNVRG
      logical lWat,lUpW,lOrt,lArtD
      dimension A(MBandD,NumNP),B(NumNP),Q(NumNP),hNew(NumNP),F(NumNP),
     !          KX(NumElD,4),MatNum(NumNP),Sink(NumNP),DS(NumNP),
     !          x(NumNP),y(NumNP),Kode(NumNP),Vx(NumNP),Vz(NumNP),
     !          ThNew(NumNP),ThOld(NumNP),KXB(NumBP),Conc(NumNP),
     !          Con(NumNP),ConAxx(NumEl),ConAzz(NumEl),ConAxz(NumEl),
     !          Qc(NumNP),Ac(NumNP),Fc(NumNP),Gc(NumNP),KodCB(NumBP),
     !          Dispxx(NumNP),Dispzz(NumNP),Dispxz(NumNP),ListNE(NumNP),
     !          ChPar(10,NMat),cBound(6),thSat(NMat),S(3,3),Bi(3),Ci(3),
     !          List(3),Wx(3),Wz(3),WeTab(3,2*NumEl),VxE(3),VzE(3),
     !          ConO(NumNP),hOld(NumNP),
     !          A1(MBandD,NumNP),B1(NumNP),RES(NumNP),IAD(MBandD,NumNP),
     !          IADN(NumNP),IADD(NumNP),VRV(NumNP),RQI(NumNPD,MNorth),
     !          RQ(NumNP),QQ(NumNP),RQIDOT(MNorth),QI(NumNPD,MNorth)

*     Initialisation
      xMul=1.
      alf=1.-epsi
      jjj=MBand
      if(t.gt.tPulse) then
        do 11 i=1,4
          cBound(i)=0.
11      continue
      end if
      do 13 i=1,NumNP
        B(i) =0.d0
        if(lOrt) B1(i)=Conc(i)
        Qc(i)=0.
        if(epsi.lt.0.001) then
          if(lOrt) jjj=IADD(i)
          A(jjj,i)=0.d0
        else
          do 12 j=1,MBandD
            A(j,i)=0.d0
12        continue
        end if
13    continue

      do 21 Level=1,NLevel
        if(Level.eq.NLevel) then
          Eps=epsi
          if(lWat)
     !      call Veloc(KAT,NumNP,NumEl,NumElD,hNew,x,y,KX,ListNE,Con,
     !                 ConAxx,ConAzz,ConAxz,Vx,Vz)
          call Disper(NumNP,NMat,Dispxx,Dispzz,Dispxz,Vx,Vz,ThNew,thSat,
     !                ChPar,MatNum,lArtD,PeCr,dt)
          call PeCour(NumNP,NumEl,NumElD,NMat,x,y,Vx,Vz,KX,MatNum,
     !                Dispxx,Dispzz,ChPar,ThNew,dt,dtMaxC,Peclet,
     !                Courant,lUpW,lArtD,PeCr)
          if(lUpW.and.lWat)
     !      call WeFact(NumNP,NumEl,NumElD,x,y,KX,WeTab,Vx,Vz,Dispxx,
     !                  Dispzz,Dispxz)
        else
          Eps=1.-epsi
          call Disper(NumNP,NMat,Dispxx,Dispzz,Dispxz,Vx,Vz,ThNew,thSat,
     !                ChPar,MatNum,lArtD,PeCr,dt)
        end if
        do 14 i=1,NumNP
          M=MatNum(i)
          if(Level.ne.NLevel) then
            if(.not.lArtD.and..not.lUpW) then
              DPom=dt/6./(ThOld(i)+ChPar(1,M)*ChPar(5,M))
              Dispxx(i)=Dispxx(i)+Vx(i)*Vx(i)*DPom
              Dispzz(i)=Dispzz(i)+Vz(i)*Vz(i)*DPom
              Dispxz(i)=Dispxz(i)+Vx(i)*Vz(i)*DPom
            end if
          else
            Ac(i)=-(ThOld(i)*alf+ThNew(i)*epsi)-ChPar(1,M)*ChPar(5,M)
            if(.not.lArtD.and..not.lUpW) then
              DPom=dt/6./(ThNew(i)+ChPar(1,M)*ChPar(5,M))
              Dispxx(i)=Dispxx(i)-Vx(i)*Vx(i)*DPom
              Dispzz(i)=Dispzz(i)-Vz(i)*Vz(i)*DPom
              Dispxz(i)=Dispxz(i)-Vx(i)*Vz(i)*DPom
            end if
            cS=cBound(5)
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
        do 19 n=1,NumEl
          CAxx=ConAxx(n)
          CAzz=ConAzz(n)
          CAxz=ConAxz(n)
          NUS=4
          if(KX(n,3).eq.KX(n,4)) NUS=3

*         Loop on subelements
          do 18 k=1,NUS-2
            NumSEl=NumSEl+1
            i=KX(n,1)
            j=KX(n,k+1)
            l=KX(n,k+2)
            List(1)=i
            List(2)=j
            List(3)=l
            Ci(1)=x(l)-x(j)
            Ci(2)=x(i)-x(l)
            Ci(3)=x(j)-x(i)
            Bi(1)=y(j)-y(l)
            Bi(2)=y(l)-y(i)
            Bi(3)=y(i)-y(j)
            AE=(Ci(3)*Bi(2)-Ci(2)*Bi(3))/2.
            
*           Calculate Velocities            
            AE1=1./AE/2.
            Ai=CAxx*Bi(1)+CAxz*Ci(1)
            Aj=CAxx*Bi(2)+CAxz*Ci(2)
            Ak=CAxx*Bi(3)+CAxz*Ci(3)
            if(Level.eq.NLevel) then
              Vxx=AE1*(Ai*hNew(i)+Aj*hNew(j)+Ak*hNew(l))
            else
              Vxx=AE1*(Ai*hOld(i)+Aj*hOld(j)+Ak*hOld(l))
            end if
            if(KAT.gt.0) Vxx=Vxx+CAxz
            Ai=CAxz*Bi(1)+CAzz*Ci(1)
            Aj=CAxz*Bi(2)+CAzz*Ci(2)
            Ak=CAxz*Bi(3)+CAzz*Ci(3)
            if(Level.eq.NLevel) then
              Vzz=AE1*(Ai*hNew(i)+Aj*hNew(j)+Ak*hNew(l))
            else
              Vzz=AE1*(Ai*hOld(i)+Aj*hOld(j)+Ak*hOld(l))
            end if
            if(KAT.gt.0) Vzz=Vzz+CAzz
            if(Level.ne.NLevel) then
              ConE=(ConO(i)+ConO(j)+ConO(l))/3.
              VxE(1)=-ConO(i)*Vxx
              VxE(2)=-ConO(j)*Vxx
              VxE(3)=-ConO(l)*Vxx
              VzE(1)=-ConO(i)*Vzz
              VzE(2)=-ConO(j)*Vzz
              VzE(3)=-ConO(l)*Vzz
            else
              ConE=(Con(i)+Con(j)+Con(l))/3.
              VxE(1)=-Con(i)*Vxx
              VxE(2)=-Con(j)*Vxx
              VxE(3)=-Con(l)*Vxx
              VzE(1)=-Con(i)*Vzz
              VzE(2)=-Con(j)*Vzz
              VzE(3)=-Con(l)*Vzz
            end if
            VxEE=-ConE*Vxx
            VzEE=-ConE*Vzz
            
            if(KAT.eq.1) xMul=2.*3.1416*(x(i)+x(j)+x(l))/3.
            if(Level.eq.NLevel) then
              cS=cBound(5)
              if(abs(cBound(5)).gt.1.e-30) then
                if(cBound(5).gt.conc(i)) cS=cS+(conc(i)-cBound(5))/3.
                if(cBound(5).gt.conc(j)) cS=cS+(conc(j)-cBound(5))/3.
                if(cBound(5).gt.conc(l)) cS=cS+(conc(l)-cBound(5))/3.
              end if
              RootCh=xMul*AE*dt*cS*(Sink(i)+Sink(j)+Sink(l))/3.
              CumCh0=CumCh0-xMul*AE*dt*(Gc(i)+Gc(j)+Gc(l))/3.+RootCh
              CumCh1=CumCh1-xMul*AE*dt*((Fc(i)-Sink(i))*conc(i)+
     !               (Fc(j)-Sink(j))*conc(j)+(Fc(l)-Sink(l))*conc(l))/3.
              CumChR=CumChR+RootCh
            end if
            FMul=xMul*AE/4.
            GcE=(Gc(i)+Gc(j)+Gc(l))/3.
            Ec1=(Dispxx(i)+Dispxx(j)+Dispxx(l))/3.
            Ec2=(Dispxz(i)+Dispxz(j)+Dispxz(l))/3.
            Ec3=(Dispzz(i)+Dispzz(j)+Dispzz(l))/3.
            if(Level.eq.NLevel) AcE=(Ac(i)+Ac(j)+Ac(l))/3.
            FcE=(Fc(i)+Fc(j)+Fc(l))/3.
            SMul1=-1./AE/4.*xMul
            SMul2=AE/20.*xMul
            if(lUpW) then
              NS=NumSEl
              W1=WeTab(1,NS)
              W2=WeTab(2,NS)
              W3=WeTab(3,NS)
              Wx(1)=2*VxE(1)*(W2-W3)+VxE(2)*(W2-2.*W3)+VxE(3)*(2.*W2-W3)
              Wx(2)=VxE(1)*(2.*W3-W1)+2*VxE(2)*(W3-W1)+VxE(3)*(W3-2.*W1)
              Wx(3)=VxE(1)*(W1-2.*W2)+VxE(2)*(2.*W1-W2)+2*VxE(3)*(W1-W2)
              Wz(1)=2*VzE(1)*(W2-W3)+VzE(2)*(W2-2.*W3)+VzE(3)*(2.*W2-W3)
              Wz(2)=VzE(1)*(2.*W3-W1)+2*VzE(2)*(W3-W1)+VzE(3)*(W3-2.*W1)
              Wz(3)=VzE(1)*(W1-2.*W2)+VzE(2)*(2.*W1-W2)+2*VzE(3)*(W1-W2)
            end if
            do 17 j1=1,3
              i1=List(j1)
              F(i1)=F(i1)+FMul*(GcE+Gc(i1)/3.)
              if(Level.eq.NLevel) DS(i1)=DS(i1)+FMul*(AcE+Ac(i1)/3.)
              iBound=0
              if(Kode(i).ne.0) then
                do 24 id=1,NumBP
                  if(KXB(id).eq.i1.and.KodCB(id).gt.0) iBound=1
24              continue
              end if
              if(iBound.eq.1) Qc(i1)=Qc(i1)-Eps*FMul*(GcE+Gc(i1)/3.)
              do 16 j2=1,3
                i2=List(j2)
                S(j1,j2)=SMul1*(Ec1*Bi(j1)*Bi(j2)+Ec3*Ci(j1)*Ci(j2)+
     !                         Ec2*(Bi(j1)*Ci(j2)+Ci(j1)*Bi(j2)))
                S(j1,j2)=S(j1,j2)-(Bi(j2)/8.*(VxEE+VxE(j1)/3.)+
     !                            Ci(j2)/8.*(VzEE+VzE(j1)/3.))*xMul
                if(lUpW) S(j1,j2)=S(j1,j2)-xMul*
     !                            (Bi(j2)/40.*Wx(j1)+Ci(j2)/40.*Wz(j1))
                ic=1
                if(i1.eq.i2) ic=2
                S(j1,j2)=S(j1,j2)+SMul2*ic*(FcE+(Fc(i1)+Fc(i2))/3.)
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
        if(abs(Conc(i)).lt.1.e-38) Conc(i)=0.
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
     !          Qc(NumNP),Kode(NumNP),KXB(NumBP),KodCB(NumBP),cBound(6),
     !          DS(NumNP),IADD(NumNP)

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

*     Point source or sink
          if(Q(i).lt.0.) then
            cKod=2
          else
            cBnd=cBound(6)
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
                A(j,i)=0.d0
13            continue
              A(MBand,i)=1.d0
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

      subroutine ChInit(NumNP,NumEl,NumElD,NMat,x,y,KX,MatNum,NLevel,
     !                  Con,hNew,Sink,cSink,Vx,Vz,ConAxx,ConAzz,ConAxz,
     !                  Dispxx,Dispzz,Dispxz,ChPar,theta,thSat,Conc,Fc,
     !                  Gc,ListNE,lUpW,WeTab,dt,dtMaxC,Peclet,Courant,
     !                  KAT,lArtD,PeCr,ConO)

      logical lUpW,lArtD
      dimension hNew(NumNP),x(NumNP),y(NumNP),KX(NumElD,4),theta(NumNP),
     !          Sink(NumNP),ChPar(10,NMat),Vx(NumNP),Vz(NumNP),
     !          MatNum(NumNP),Con(NumNP),ConAxx(NumEl),ConAzz(NumEl),
     !          ConAxz(NumEl),Dispxx(NumNP),Dispzz(NumNP),Dispxz(NumNP),
     !          ListNE(NumNP),Fc(NumNP),Gc(NumNP),WeTab(3,2*NumEl),
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
      call Veloc(KAT,NumNP,NumEl,NumElD,hNew,x,y,KX,ListNE,Con,ConAxx,
     !           ConAzz,ConAxz,Vx,Vz)
      call Disper(NumNP,NMat,Dispxx,Dispzz,Dispxz,Vx,Vz,theta,thSat,
     !            ChPar,MatNum,lArtD,PeCr,dt)
      call PeCour(NumNP,NumEl,NumElD,NMat,x,y,Vx,Vz,KX,MatNum,Dispxx,
     !            Dispzz,ChPar,theta,dt,dtMaxC,Peclet,Courant,lUpW,
     !            lArtD,PeCr)
      if(lUpW) 
     !  call WeFact(NumNP,NumEl,NumElD,x,y,KX,WeTab,Vx,Vz,Dispxx,
     !              Dispzz,Dispxz)
      return
      end

************************************************************************

*     Calculate  velocities

      subroutine Veloc(KAT,NumNP,NumEl,NumElD,hNew,x,y,KX,ListNE,Con,
     !                 ConAxx,ConAzz,ConAxz,Vx,Vz)

      integer e
      dimension hNew(NumNP),x(NumNP),y(NumNP),ListNE(NumNP),Con(NumNP),
     !          KX(NumElD,4),Vx(NumNP),Vz(NumNP),ConAxx(NumEl),
     !          ConAzz(NumEl),ConAxz(NumEl),List(3)

      do 11 i=1,NumNP
        Vx(i)=0.
        Vz(i)=0.
11    continue    
      do 14 e=1,NumEl
        CAxx=ConAxx(e)
        CAzz=ConAzz(e)
        CAxz=ConAxz(e)
        NCorn=4
        if(KX(e,3).eq.KX(e,4)) NCorn=3
        do 13 n=1,NCorn-2
          i=KX(e,1)
          j=KX(e,n+1)
          k=KX(e,n+2)
          List(1)=i
          List(2)=j
          List(3)=k
          vi=y(j)-y(k)
          vj=y(k)-y(i)
          vk=y(i)-y(j)
          wi=x(k)-x(j)
          wj=x(i)-x(k)
          wk=x(j)-x(i)
          Area=.5*(wk*vj-wj*vk)
          A=1./Area/2.
          Ai=CAxx*vi+CAxz*wi
          Aj=CAxx*vj+CAxz*wj
          Ak=CAxx*vk+CAxz*wk  
          Vxx=A*(Ai*hNew(i)+Aj*hNew(j)+Ak*hNew(k))
          if(KAT.gt.0) Vxx=Vxx+CAxz
          Ai=CAxz*vi+CAzz*wi
          Aj=CAxz*vj+CAzz*wj
          Ak=CAxz*vk+CAzz*wk
          Vzz=A*(Ai*hNew(i)+Aj*hNew(j)+Ak*hNew(k))
          if(KAT.gt.0) Vzz=Vzz+CAzz
          do 12 m=1,3
            l=List(m)
            Vx(l)=Vx(l)-Con(l)*Vxx
            Vz(l)=Vz(l)-Con(l)*Vzz
12        continue
13      continue
14    continue
      do 15 i=1,NumNP
        Vx(i)=Vx(i)/ListNE(i)
        Vz(i)=Vz(i)/ListNE(i)
15    continue
      return
      end

************************************************************************    

*     Calculate the dispersion coefficient

      subroutine Disper(NumNP,NMat,Dispxx,Dispzz,Dispxz,Vx,Vz,theta,
     !                  thSat,ChPar,MatNum,lArtD,PeCr,dt)

      logical lArtD
      dimension Vx(NumNP),Vz(NumNP),theta(NumNP),ChPar(10,NMat),
     !          Dispxx(NumNP),Dispzz(NumNP),Dispxz(NumNP),MatNum(NumNP),
     !          thSat(NMat)

      do 11 i=1,NumNP
        M=MatNum(i)
        Tau=theta(i)**(7./3.)/thSat(M)**2
        Vabs=sqrt(Vx(i)*Vx(i)+Vz(i)*Vz(i))
        Dif=theta(i)*ChPar(2,M)*Tau
        DispL=ChPar(3,M)
        DispT=ChPar(4,M)
        if(lArtD.and.VAbs.gt.1.e-20) DispL=amax1(DispL,VAbs*dt
     !                 /(theta(i)+ChPar(1,M)*ChPar(5,M))/PeCr-Dif/VAbs)
        Dispxx(i)=Dif
        Dispzz(i)=Dif
        Dispxz(i)=0.
        if(VAbs.gt.1.e-20) then
          Dispxx(i)=DispL*Vx(i)*Vx(i)/VAbs+DispT*Vz(i)*Vz(i)/VAbs+Dif
          Dispzz(i)=DispL*Vz(i)*Vz(i)/VAbs+DispT*Vx(i)*Vx(i)/VAbs+Dif
          Dispxz(i)=(DispL-DispT)*Vx(i)*Vz(i)/VAbs
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
            if(abs(Disp).gt.1.e-30) aa=abs(Vel/Disp)
            if(abs(Disp).lt.1.e-30.or.abs(Vel).lt.0.001*VV.or.
     !                                abs(aa).gt.10.) then
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

      subroutine PeCour(NumNP,NumEl,NumElD,NMat,x,y,Vx,Vz,KX,MatNum,
     !                  Dispxx,Dispzz,ChPar,theta,dt,dtMaxC,Peclet,
     !                  Courant,lUpW,lArtD,PeCr)

      logical lUpW,lArtD
      dimension KX(NumElD,4),x(NumNP),y(NumNP),Vx(NumNP),Vz(NumNP),
     !          MatNum(NumNP),Dispxx(NumNP),Dispzz(NumNP),theta(NumNP),
     !          ChPar(10,NMat),Bi(3),Ci(3)

      Peclet=0.
      Courant=0.
      dtMaxC=1.e+30
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
          if(DxE.gt.1.e-20) PecX=VxE*delX/DxE
          if(DzE.gt.1.e-20) PecY=VzE*delY/DzE
          if(PecX.ne.99999.) Peclet=amax1(Peclet,PecX)
          if(PecY.ne.99999.) Peclet=amax1(Peclet,PecY)
          Peclet=amin1(Peclet,99999.)

          VxMax=amax1(abs(Vx(i))/theta(i),abs(Vx(j))/theta(j),abs(Vx(l))
     !                /theta(l))
          VzMax=amax1(abs(Vz(i))/theta(i),abs(Vz(j))/theta(j),abs(Vz(l))
     !                /theta(l))
          R1=1.+ChPar(1,MatNum(i))*ChPar(5,MatNum(i))/theta(i)
          R2=1.+ChPar(1,MatNum(j))*ChPar(5,MatNum(j))/theta(j)
          R3=1.+ChPar(1,MatNum(l))*ChPar(5,MatNum(l))/theta(l)
          RMin=amin1(R1,R2,R3)
          CourX=VxMax*dt/delX/RMin
          CourY=VzMax*dt/delY/RMin
          Courant=amax1(Courant,CourX,CourY)

          Cour1=1.0
          Cour2=1.0
          if(.not.lUpW.and..not.lArtD) then
            if(PecX.ne.99999.) Cour1=amin1(1.,PeCr/amax1(0.5,PecX))
            if(PecY.ne.99999.) Cour2=amin1(1.,PeCr/amax1(0.5,PecY))
          end if
          if(VxMax.gt.1.e-20) dt1=Cour1*delX*RMin/VxMax
          if(VzMax.gt.1.e-20) dt2=Cour2*delY*RMin/VzMax
          dtMaxC=amin1(dtMaxC,dt1,dt2)

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