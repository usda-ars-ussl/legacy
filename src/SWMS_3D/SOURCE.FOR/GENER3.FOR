*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
*                                                                      *
*     *** Purpose:  Finite element grid generator for SWMS_3D program  *
*     Designed by J. Simunek, August 94                                *
*     Finite element mesh is composed of quadrilateral paralleleppipeds*
*                                                                      *
*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*

      program Gener3

      parameter (nzD=200,
     !           nxD=200,
     !           nyD=200)
 
      dimension xCol(nxD),yCol(nyD),zLin(nzD),KXB(2*nxD*nyD),
     !          Width(2*nxD*nyD),Kode(nzD),hInit(nzD),Q(nzD),M(nzD),
     !          B(nzD),Axz(nzD),Bxz(nzD),Dxz(nzD),Conc(nzD)
      logical lSink

      open(10,file='SWMS_3D.IN\Gener3.in',status='old')
      open(50,file='SWMS_3D.IN\Grid.in'  ,status='unknown')

      read(10,*) 
      read(10,*)
      read(10,*) Aniz1,Aniz2,Aniz3,a11,a22,a33,a12,a13,a23
c      read(10,*) Angle,Aniz1,Aniz2,Aniz3
      read(10,*)
      read(10,*) nz,nx,ny
      read(10,*)
      read(10,*) xCol(1),yCol(1),zLin(nz)
      read(10,*)
      read(10,*) 
      read(10,*) (xCol(i),i=2,nx)
      read(10,*)
      read(10,*) (yCol(i),i=2,ny)
      read(10,*)
      read(10,*) (zLin(i),i=1,nz-1)
      read(10,*)
      read(10,*)

*     Transformation from increments to coordinates
      do 11 i=2,nx
        xCol(i)=xCol(i-1)+xCol(i)
11    continue 
      do 12 i=2,ny
        yCol(i)=yCol(i-1)+yCol(i)
12    continue 
      do 13 i=nz-1,1,-1
        zLin(i)=zLin(i+1)+zLin(i)
13    continue

      nex=nx-1
      ney=ny-1
      nez=nz-1
      NLay=nx*ny
      NumNP=nz*nx*ny
      NumEl=(nz-1)*(nx-1)*(ny-1)
      NObs=0
      NumBP=NLay*2
      write(50,110)
      write(50,120)
      write(50,130) NumNP,NumEl,NLay,NumBP,NObs
      write(50,140)

*     Read and interpolate line information
      np=0
      i=0
      lSink=.false.
14    continue
        i=i+1
        read(10,*) n,Kode(n),hInit(n),Conc(n),Q(n),M(n),B(n),Axz(n),
     !             Bxz(n),Dxz(n)
        if(B(n).gt.0.) lSink=.true.

        if (n-i) 15,18,16
15      write(*,*)'error in line',n
        stop
16      continue
        dz=zLin(n)-zLin(np)
        hInitS=(hInit(n)-hInit(np))/dz
        BetaS =(B(n)-B(np))/dz
        ConcS=(Conc(n)-Conc(np))/dz 
        AS=(Axz(n)-Axz(np))/dz
        BS=(Bxz(n)-Bxz(np))/dz
        DS=(Dxz(n)-Dxz(np))/dz
17      continue
          dz=zLin(i)-zLin(np)
          hInit(i)=hInit(np)+hInitS*dz
          B(i)=B(np)+BetaS*dz
          Conc(i)=Conc(np)+ConcS*dz
          Axz(i)=Axz(np)+AS*dz
          Bxz(i)=Bxz(np)+BS*dz
          Dxz(i)=Dxz(np)+DS*dz
          M(i)=M(i-1)
          Kode(i)=Kode(i-1)
          Q(i)=Q(i-1)
          i=i+1
        if(i.lt.n) goto 17
18      continue
        np=n
      if(i.lt.nz) goto 14

      n=0
      do 21 iz=1,nz
        do 20 iy=1,ny
          do 19 ix=1,nx
            n=n+1
            write(50,150) n,Kode(iz),xCol(ix),yCol(iy),zLin(iz),
     !                    hInit(iz),Conc(iz),Q(iz),M(iz),B(iz),Axz(iz),
     !                    Bxz(iz),Dxz(iz)
19        continue            
20      continue
21    continue

      write(50,160)
      write(50,170)
      ie=0
      do 24 iz=1,nez
        do 23 iy=1,ney
          do 22 ix=1,nex
            ie=ie+1
            i1=ix+nx*(iy-1)+nx*ny*(iz-1)
            i2=i1+1
            i3=i1+nx
            i4=i1+nx+1
            i5=i1+nx*ny
            i6=i1+nx*ny+1
            i7=i1+nx*ny+nx
            i8=i1+nx*ny+nx+1
            iCord=ix+iy+iz
            iRest=mod(iCord,2)
            if(iRest.eq.1) then
c              write(50,180) ie,i1,i2,i3,i4,i5,i6,i7,i8,1,Angle,Aniz1,
c     !                      Aniz2,Aniz3,M(iz)
              write(50,180) ie,i1,i2,i3,i4,i5,i6,i7,i8,1,Aniz1,Aniz2,
     !                      Aniz3,a11,a22,a33,a12,a13,a23,M(iz)
            else
c              write(50,180) ie,i1,i2,i3,i4,i5,i6,i7,i8,2,Angle,Aniz1,
c     !                      Aniz2,Aniz3,M(iz)
              write(50,180) ie,i1,i2,i3,i4,i5,i6,i7,i8,2,Aniz1,Aniz2,
     !                      Aniz3,a11,a22,a33,a12,a13,a23,M(iz)
            end if            
22        continue            
23      continue
24    continue

      write(50,190)

*     Boundary geom. information
      iB=0
      do 26 iy=1,ny
        do 25 ix=1,nx
          iB=iB+1
          KXB(iB)=iB
          KXB(iB+ny*nx)=nx*ny*(nz-1)+iB
          if(ix.ne.1.and.ix.ne.nx) then
            WidthX=(xCol(ix+1)-xCol(ix-1))
          else if(ix.eq.1) then
            WidthX=(xCol(ix+1)-xCol(1))
          else if(ix.eq.nx) then
            WidthX=(xCol(nx)-xCol(nx-1))
          end if
          if(iy.ne.1.and.iy.ne.ny) then
            WidthY=(yCol(iy+1)-yCol(iy-1))
          else if(iy.eq.1) then
            WidthY=(yCol(iy+1)-yCol(1))
          else if(iy.eq.ny) then
            WidthY=(yCol(ny)-yCol(ny-1))
          end if
          iCord=ix+iy
          iRest=mod(iCord,2)
          iDiv=3
          if(iRest.eq.1) iDiv=6
          Width(iB)=WidthX*WidthY/iDiv
          iCord=ix+iy+nz
          iRest=mod(iCord,2)
          iDiv=3
          if(iRest.eq.0) iDiv=6
          Width(iB+ny*nx)=WidthX*WidthY/iDiv
25      continue
26    continue
      rLen=0.
      if(lSink) rLen=(xCol(nx)-xCol(1))*(yCol(ny)-yCol(1))

      write(50,200) (KXB(i),  i=1,2*nx*ny)
      write(50,210) (Width(i),i=1,2*nx*ny)
      write(50,220)  rLen
      write(50,230)

110   format('*** BLOCK H: NODAL INFORMATION ',78('*'))
120   format('   NumNP   NumEl      IJ   NumBP    NObs')
130   format(5(3x,i5))
140   format('    n Code       x         y         z          h        C
     !onc         Q        M    B    Axz   Bxz   Dxz')
150   format(i5,3x,i2,3(1x,f9.2),3(1x,e11.4),2x,i2,4(1x,f5.2))
160   format('*** BLOCK I: ELEMENT INFORMATION ',76('*'))
170   format('    e     i     j     k     l     m     n     o     p Sub 
     !Aniz1 Aniz2 Aniz3 ------- cosin of angle ------ Lay')
180   format(i5,1x,8(i5,1x),i2,1x,3(2x,f4.2),6(1x,f4.2),i3)
190   format('*** BLOCK J : BOUNDARY INFORMATION ',74('*'))
200   format(' Node number array:'/(8i13) )
210   format(' Width array:'/(8f13.5) )
220   format(' Length:'/f10.2)
230   format('***** End of file Grid.In ',83('*'))
      close (50)
      stop
      end
