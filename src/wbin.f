c  Modified by Bin Wang (2012/07/03)

c  Part of R package KernSmooth
c  Copyright (C) 1995  M. P. Wand
c
c  Unlimited use and distribution (see LICENCE).

cccccccccc FORTRAN subroutine linbin.f cccccccccc

c Obtains bin counts for univariate data
c via the linear binning strategy. If "trun=0" then
c weight from end observations is given to corresponding
c end grid points. If "trun=1" then end observations
c are truncated.

c Last changed: 20 MAR 2009

      subroutine wbin(X,W,n,a,b,M,trun,gcnts)
      double precision X(*),W(*),a,b,gcnts(*),lxi,delta,rem
      integer n,M,i,li,trun

c     Initialize grid counts to zero

      do 10 i=1,M
         gcnts(i) = dble(0)
10    continue

      delta = (b-a)/(M-1)
      do 20 i=1,n
         lxi = ((X(i)-a)/delta) + 1

c        Find integer part of "lxi"

         li = int(lxi) 

         if (li.ge.1.and.li.lt.M) then
            gcnts(li) = gcnts(li) + W(i)
         endif

         if (li.lt.1.and.trun.eq.0) then
            gcnts(1) = gcnts(1) +  W(i)
         endif

         if (li.ge.M.and.trun.eq.0) then
            gcnts(M-1) = gcnts(M-1) +  W(i)
         endif

20    continue

      return
      end

cccccccccc End of wlinbin.f cccccccccc
