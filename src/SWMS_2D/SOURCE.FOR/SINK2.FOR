* Source file SINK2.FOR ||||||||||||||||||||||||||||||||||||||||||||||||

      subroutine SetSnk(NumNP,NMat,MatNum,hNew,TPot,Sink,P0,POptm,P2H,
     !                  P2L,P3,r2H,r2L,Beta,Length)

      real Length
      dimension MatNum(NumNP),hNew(NumNP),POptm(NMat),Beta(NumNP),
     !          Sink(NumNP)

      do 11 i=1,NumNP
        if(Beta(i).gt.0.) then
          M=MatNum(i)
          Alfa=FAlfa(TPot,hNew(i),P0,POptm(M),P2H,P2L,P3,r2H,r2L)
          Sink(i)=Alfa*Beta(i)*Length*TPot
        end if
11    continue
      return
      end 
   
************************************************************************
 
      real function FAlfa(TPot,h,P0,P1,P2H,P2L,P3,r2H,r2L)

      if(TPot.lt.r2L) P2=P2L
      if(TPot.gt.r2H) P2=P2H
      if((TPot.ge.r2L).and.(TPot.le.r2H))
     !  P2=P2H+(r2H-TPot)/(r2H-r2L)*(P2L-P2H)
      FAlfa=0.0
      if((h.gt.P3).and.(h.lt.P2)) FAlfa=(h-P3)/(P2-P3)
      if((h.ge.P2).and.(h.le.P1)) FAlfa=1.0
      if((h.gt.P1).and.(h.lt.P0)) FAlfa=(h-P0)/(P1-P0)
      return
      end

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||