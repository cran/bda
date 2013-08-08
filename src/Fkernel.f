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
      
      do 10 i=1,M
         gcnts(i) = dble(0)
 10   continue
      
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


ccccccccccFORTRAN subroutine smoothkde cccccccccc

c     2013/04/01.

      subroutine smoothkde(fx,x0,n,x,f,m,w,h,iter)
      integer m,n,i,j,iter, nsum,k,i1,i2,loop
      double precision fx(n),f0(n),x0(n),x(m),f(m),h,w
      double precision wd, dx, fk(n,n),gk(n),df
      double precision t1,t2, xl,xu

      DOUBLE PRECISION PI,TWOPI, twoh, twopih
      PARAMETER(PI=3.141592653589793D0,TWOPI=2D0*PI)
      
      wd = 2.0 * w
      dx = x0(2) - x0(1)
      twoh = -0.5/h/h
      twopih = 1.0/sqrt(twopi)/h

      nsum = 0
      do 50 i=1, m
         nsum = nsum + f(i)
 50   enddo
      do 100 i=1, n
         f0(i) = fx(i)
 100  enddo
      do 150 i=1, n
         t1 = dx*(i-1.0)/h
         gk(i) = twopih * exp(-0.5*t1*t1)
 150  enddo


      do 250 i=1, n-1
         fk(i,i) = gk(1)
         do 200 j=i+1, n
            fk(i,j) = gk(j-i)
            fk(j,i) = fk(i,j)
 200      enddo
 250   enddo
       fk(n,n) = gk(1)

       loop = 0
       df = 1.0
       do while((loop.LT.iter) .AND. (df.GT.0.0001))
          df = 0.0
          do 500 k=1,n
             fx(k) = 0.0
             do 450 i=1, m
                xl = x(i) - w
                xu = x(i) + w
                i1 = max0(1, ceiling((xl-x0(1))/dx))
                i2 = min0(n, floor((xu-x0(1))/dx))
                t1 = 0.0
                t2 = 0.0
                do 400 j=i1,i2
                   t1 = t1 + fk(k,j)*f0(j)
                   t2 = t2 + f0(j)
 400            enddo
                fx(k) = fx(k) + f(i)*t1/t2/nsum
 450         enddo
             df = df + (fx(k) - f0(k))**2.0
             f0(k) = fx(k)
 500      enddo
       enddo
      iter = loop
      end
      
cccccccccc End of smoothkde cccccccccc

ccccccccccFORTRAN subroutine yldist.f cccccccccc
      
c     This subroutine is to compute the modules of y_ell based on the
c     fourier transformed data after binning.  It can be called to
c     compute the LSCV scores based on data FFT-transformed.
      
C     Last changed: 26 Oct 2012

      subroutine yldist(gcnts, M, Y)
      double precision Y(*),gcnts(*),theta, sumcos,sumsin
      integer n,M,l
      DATA PI2/6.2831853/

c     Initialize Y_ell  to zero
      do 10 l = 1,M/2
         Y(l) = dble(0)
 10   continue

      do 100 l = 1,M/2
         sumcos = 0.0
         sumsin = 0.0
         do 50 k = 1,M
            theta = PI2 * (k-1.0) * l / M
            sumcos = sumcos + gcnts(k) * COS(theta)
            sumsin = sumsin + gcnts(k) * SIN(theta)
 50      continue
         Y(l) = (sumcos*sumcos + sumsin * sumsin)/M/M
 100  continue

      return
      end

cccccccccc End of yldist.f cccccccccc

