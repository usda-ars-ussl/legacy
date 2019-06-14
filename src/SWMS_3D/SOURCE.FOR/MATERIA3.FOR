* Source file MATERIAL.FOR |||||||||||||||||||||||||||||||||||||||||||||

      real function FK(h,Par)

      implicit double precision(A-H,O-Z)
      double precision n,m,Ks,Kr,Kk
      real h,Par(9)
      integer PPar

      BPar=.5
      PPar=2.
      Qr=Par(1)
      Qs=Par(2)
      Qa=Par(3)
      Qm=Par(4)
      Alfa=Par(5)
      n=Par(6)
      Ks=Par(7)
      Kk=Par(8)
      Qk=Par(9)
      m=1.-1./n
      HMin=-1.d300**(1./n)/max(Alfa,1.d0)
      HH=max(dble(h),HMin)
      Qees=dmin1((Qs-Qa)/(Qm-Qa),.999999999999999d0)
      Qeek=dmin1((Qk-Qa)/(Qm-Qa),Qees)
      Hs=-1./Alfa*(Qees**(-1./m)-1.)**(1./n)
      Hk=-1./Alfa*(Qeek**(-1./m)-1.)**(1./n)
      if(h.lt.Hk) then
        Qee=(1.+(-Alfa*HH)**n)**(-m)
        Qe =(Qm-Qa)/(Qs-Qa)*Qee
        Qek=(Qm-Qa)/(Qs-Qa)*Qeek
        FFQ =1.-(1.-Qee**(1./m))**m
        FFQk=1.-(1.-Qeek**(1./m))**m
        if(FFQ.le.0.d0) FFQ=m*Qee**(1./m)
        Kr=(Qe/Qek)**Bpar*(FFQ/FFQk)**PPar*Kk/Ks
        FK=max(Ks*Kr,1.d-37)
        return
      end if
      if(h.ge.Hk.and.h.lt.Hs) then
        Kr=(1.-Kk/Ks)/(Hs-Hk)*(h-Hs)+1.
        FK=Ks*Kr
      end if
      if(h.ge.Hs) FK=Ks
      return
      end

************************************************************************

      real function FC(h,Par)

      implicit double precision(A-H,O-Z)
      double precision n,m
      real h,Par(9)

      Qr=Par(1)
      Qs=Par(2)
      Qa=Par(3)
      Qm=Par(4)
      Alfa=Par(5)
      n=Par(6)
      m=1.-1./n
      HMin=-1.d300**(1./n)/max(Alfa,1.d0)
      HH=max(dble(h),HMin)
      Qees=dmin1((Qs-Qa)/(Qm-Qa),.999999999999999d0)
      Hs=-1./Alfa*(Qees**(-1./m)-1.)**(1./n)
      if(h.lt.Hs) then
        C1=(1.+(-Alfa*HH)**n)**(-m-1.)
        C2=(Qm-Qa)*m*n*(Alfa**n)*(-HH)**(n-1.)*C1
        FC=max(C2,1.d-37)
        return
      else
        FC=0.0
      end if
      return
      end

************************************************************************

      real function FQ(h,Par)

      implicit double precision(A-H,O-Z)
      double precision n,m
      real h,Par(9)

      Qr=Par(1)
      Qs=Par(2)
      Qa=Par(3)
      Qm=Par(4)
      Alfa=Par(5)
      n=Par(6)
      m=1.-1./n
      HMin=-1.d300**(1./n)/max(Alfa,1.d0)
      HH=max(dble(h),HMin)
      Qees=dmin1((Qs-Qa)/(Qm-Qa),.999999999999999d0)
      Hs=-1./Alfa*(Qees**(-1./m)-1.)**(1./n)
      if(h.lt.Hs) then
        Qee=(1.+(-Alfa*HH)**n)**(-m)
        FQ=max(Qa+(Qm-Qa)*Qee,1.d-37)
        return
      else
        FQ=Qs
      end if
      return
      end

************************************************************************

      real function FH(Qe,Par)

      implicit double precision(A-H,O-Z)
      double precision n,m
      real Qe,Par(9)

      Qr=Par(1)
      Qs=Par(2)
      Qa=Par(3)
      Qm=Par(4)
      Alfa=Par(5)
      n=Par(6)
      m=1.-1./n
      HMin=-1.d300**(1./n)/max(Alfa,1.d0)
      QeeM=(1.+(-Alfa*HMin)**n)**(-m)
      Qee=dmin1(dmax1(Qe*(Qs-Qa)/(Qm-Qa),QeeM),.999999999999999d0)
      FH=max(-1./Alfa*(Qee**(-1./m)-1.)**(1./n),-1.d37)
      return
      end

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||