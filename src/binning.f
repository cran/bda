C     This code is modified based on linbin subroutine by M.P.Wand in
C     KernSmooth.  (March 27, 2013 by Bin Wang)

c     Part of R package KernSmooth
c     Copyright (C) 1995  M. P. Wand
c     
c     Unlimited use and distribution (see LICENCE).
      
ccccccccccFORTRAN subroutine binning.f cccccccccc
      
c     Obtains bin counts for univariate data via different binning
c     strategy. 

c     Last changed: 19 APRIL 2013

 
      subroutine binning(X,F,W,a,b,n,xa,xb,M,gcnts,itype)
      double precision X(*),F(*),W(*),a, b,xa,xb, gcnts(*),
     *     lxi,delta,rem
      integer n,M,i,li,itype
      
c     Initialize grid counts to zero
      
      do i=1,M
         gcnts(i) = dble(0)
      end do
      
      delta = (xb - xa)/(M - 1)
c     KDE using FFT
      IF(itype .eq. 0) THEN
         do 20 i = 1,n
            lxi = ((X(i)-xa)/delta) + 1
            li = int(lxi) 
            rem = lxi - li
            if (li.ge.1.and.li.lt.M) then
               gcnts(li) = gcnts(li) + (1.0-rem)*F(i)*W(i)
               gcnts(li+1) = gcnts(li+1) + rem*F(i)*W(i)
            endif

            if (li.lt.1) then
               gcnts(1) = gcnts(1) + F(i)*W(i)
            endif
            
            if (li.ge.M) then
               gcnts(M) = gcnts(M) + F(i)*W(i)
            endif
 20      continue
c     Binning to compute E.D.F.
      ELSE IF (itype .eq. 1) THEN
         do 40 i = 1,n
            lxi = ((X(i)+b-xa)/delta) + 1
            li = int(lxi) 
            rem = lxi - li
            if (li.ge.1.and.li.lt.M) then
               if(rem .gt. 0) then
                  gcnts(li+1) = gcnts(li+1) + F(i)*W(i)
               else
                  gcnts(li) = gcnts(li) + F(i)*W(i)
               endif
            endif

            if (li.lt.1) then
               gcnts(1) = gcnts(1) + F(i)*W(i)
            endif

            if (li.eq.M) then
               gcnts(M) = gcnts(M) + F(i)*W(i)
            endif
 40      continue
c     binning to construct histogram
      ELSE
         do 60 i = 1,n
            lxi = ((X(i)+0.5*(b+a)-xa)/delta) + 1
            li = int(lxi) 
            rem = lxi - li
            if (li.ge.1.and.li.lt.M) then
               gcnts(li+1) = gcnts(li+1) + F(i)*W(i)
            endif
            
            if (li.eq.M) then
               gcnts(M) = gcnts(M) + F(i)*W(i)
            endif
 60      continue
         
      ENDIF
      
      return
      end


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

      subroutine linbin(X,n,a,b,M,trun,gcnts)
      double precision X(*),a,b,gcnts(*),lxi,delta,rem
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

         rem = lxi - li
         if (li.ge.1.and.li.lt.M) then
            gcnts(li) = gcnts(li) + (1-rem)
            gcnts(li+1) = gcnts(li+1) + rem
         endif

         if (li.lt.1.and.trun.eq.0) then
            gcnts(1) = gcnts(1) + 1
         endif

         if (li.ge.M.and.trun.eq.0) then
            gcnts(M) = gcnts(M) + 1
         endif

20    continue

      return
      end

cccccccccc End of linbin.f cccccccccc

c  Part of R package KernSmooth
c  Copyright (C) 1995  M. P. Wand
c
c  Unlimited use and distribution (see LICENCE).

cccccccccc FORTRAN subroutine linbin2D.f cccccccccc

c Obtains bin counts for bivariate data
c via the linear binning strategy. In this version
c observations outside the mesh are ignored.

      subroutine lbtwod(X,n,a1,a2,b1,b2,M1,M2,gcnts)
      integer n,M1,M2,i,li1,li2,ind1,ind2,ind3,ind4
      double precision X(*),a1,a2,b1,b2,gcnts(*)
      double precision lxi1,lxi2,delta1,delta2,rem1,rem2

c     Initialize grid cnts to zero

      do 10 i = 1,(M1*M2)
         gcnts(i) = dble(0)
10    continue

      delta1 = (b1 - a1)/(M1 - 1)
      delta2 = (b2 - a2)/(M2 - 1)
      do 20 i = 1,n
         lxi1 = ((X(i) - a1)/delta1) + 1
         lxi2 = ((X(n+i) - a2)/delta2) + 1

c        Find the integer part of "lxi1" and "lxi2"

         li1 = int(lxi1)
         li2 = int(lxi2)
         rem1 = lxi1 - li1
         rem2 = lxi2 - li2
         if (li1.ge.1) then
            if (li2.ge.1) then
               if (li1.lt.M1) then
                  if (li2.lt.M2) then
                     ind1 = M1*(li2-1) + li1
                     ind2 = M1*(li2-1) + li1 + 1
                     ind3 = M1*li2 + li1
                     ind4 = M1*li2 + li1 + 1
                     gcnts(ind1) = gcnts(ind1)+(1-rem1)*(1-rem2)
                     gcnts(ind2) = gcnts(ind2)+rem1*(1-rem2)
                     gcnts(ind3) = gcnts(ind3)+(1-rem1)*rem2
                     gcnts(ind4) = gcnts(ind4)+rem1*rem2
                  endif
               endif
            endif
         endif
20    continue

      return
      end

cccccccccc End of linbin2D.f cccccccccc

cccccccccc FORTRAN subroutine bintwod cccccccccc

c     Obtains bin counts for bivariate data over two grids the two grids
c     are both in increasing orders.
      
      subroutine bintwod(X,n,g1,g2,M1,M2,gcnts)
      integer n,M1,M2,i,j,li1,li2
      double precision X(*),g1(*),g2(*),gcnts(*)

c     Initialize grid cnts to zero

      do 10 i = 1,(M1*M2)
         gcnts(i) = dble(0)
10    continue

      do 100 i = 1,n

c     Find the subscripts of the counter
         li1 = 0
         do 20 j = 1,M1
            if(X(i).lt.g1(j)) then
               li1 = j
               exit
            endif
 20      continue
         
         li2 = 0
         do 30 j = 1,M2
            if(X(n+i).lt.g2(j)) then
               li2 = j
               exit
            endif
 30      continue
         
         if (li1.ge.1) then
            if (li2.ge.1) then
               j = M1*(li2-1) + li1
               gcnts(j) = gcnts(j) + 1
            endif
         endif
 100  continue

      return
      end

cccccccccc End of bintod cccccccccc

cccccccccc FORTRAN subroutine wlinbin.f cccccccccc

c Obtains bin counts for univariate data
c via the linear binning strategy. If "trun=0" then
c weight from end observations is given to corresponding
c end grid points. If "trun=1" then end observations
c are truncated.

c Last changed: 20 MAR 2009

      subroutine wlinbin(X,W,n,a,b,M,trun,gcnts)
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

         rem = lxi - li
         if (li.ge.1.and.li.lt.M) then
            gcnts(li) = gcnts(li) + (1-rem) * W(i)
            gcnts(li+1) = gcnts(li+1) + rem * W(i)
         endif

         if (li.lt.1.and.trun.eq.0) then
            gcnts(1) = gcnts(1) +  W(i)
         endif

         if (li.ge.M.and.trun.eq.0) then
            gcnts(M) = gcnts(M) +  W(i)
         endif

20    continue

      return
      end

cccccccccc End of wlinbin.f cccccccccc
