c  Part of R package KDE
c  Copyright (C) 2010 Bin Wang
c
c  Unlimited use and distribution (see LICENCE).

cccccccccc FORTRAN subroutine discdata.f cccccccccc

c Discretized the data to find the weight sequence xi.  A vector of
c frequency "F" is used for groupded data.


c Last changed: 20 MAR 2009

      subroutine discdata0(X,F,n,a,b,M,xi)
      double precision X(*),F(*),a,b,xi(*)
      double precision d1,d2,d3,delta,tsum,t,p
      integer n,M,i,j,k,l,ia,ib,ic

c     Initialize grid counts to zero

      do 10 i=1,M+1
         xi(i) = dble(0)
 10   continue
      delta = (b-a)/M

      t = delta
      k=1
      DO 80 i = 1,n
         IF(X(i) .LT. t) THEN 
            xi(k) = xi(k) + F(i)
            GOTO 70
         ELSEIF(X(i) .LT. t+delta) THEN
            p = (X(i)-t)/delta
            xi(k) = xi(k) + F(i)*(1-p)
            xi(k+1) = xi(k+1) + F(i)*p
         ENDIF 
 70      t = t+delta  
         k = k + 1
 80   ENDDO
 
      
      do 100 i=1,M+1
         xi(i) = xi(i)/n/delta
 100  continue
      
      return
      end

cccccccccc End of discdata.f cccccccccc
