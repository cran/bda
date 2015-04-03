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

c  Part of R package BDA
c  Copyright (C) 2009-2010 Bin Wang
c
c  Unlimited use and distribution (see LICENCE).

cccccccccc FORTRAN subroutine ofcpdf.f cccccccccc

c     To compute the pdf based on a Berkson measurement error with
c     uniform rounding errors

c     Last changed: 05 Oct 2010

      subroutine ofcpdf(y,f,a,b,ny,x,nx,h)
      integer ny,nx,i,j,nsum
      double precision y(ny),f(ny),w(ny),a(ny),b(ny),x(nx),h,
     *     tmp, t1,t2,sqrt2h,sqrt2
      sqrt2 = sqrt(2.0)
      sqrt2h = sqrt2 * h

      nsum=0
      do 10 i=1,ny
 10      nsum = nsum+f(i)

      do 100 i=1,nx
         tmp = 0.0
         do 90 j=1,ny
            t1 = (b(j) + y(j) - x(i)) / sqrt2h
            t2 = (x(i) - a(j) - y(j)) / sqrt2h
            tmp = tmp + 0.5 *f(j)*(erf(t1/sqrt2)+erf(t2/sqrt2))
     *           / (b(j)-a(j))
 90      continue
         x(i) = tmp/nsum
 100  continue
      end
      
cccccccccc End of ofcpdf.f cccccccccc

cccccccccc FORTRAN subroutine ofccdf.f cccccccccc

c     To compute the cdf based on a Berkson measurement error with
c     uniform rounding errors

c     Last changed: 06 Oct 2010

      subroutine ofccdf(y,f,a,b,ny,x,nx,h)
      integer ny,nx,i,j,nsum
      double precision y(ny),f(ny),w(ny),a(ny),b(ny),x(nx),h,
     *     tmp, t1,t2,sqrt2h,sqrt2,PI,PI2,h2
      PI = 3.1415926536D0
      PI2 = 1./DSQRT(2.*PI)
      sqrt2 = sqrt(2.0)
      sqrt2h = sqrt2 * h
      h2 = h * h

      nsum=0
      do 10 i=1,ny
 10      nsum = nsum+f(i)

      do 100 i=1,nx
         tmp = 0.0
         do 90 j=1,ny
            t1 = b(j) + y(j) - x(i)
            t2 = x(i) - a(j) - y(j)
            tmp = tmp + f(j) / (b(j) - a(j)) *
     *           (h * PI2 * (dexp(-0.5 * t2 * t2/h2)
     *           -dexp(-0.5 * t1 * t1/h2))
     *           - dabs(t1) * 0.5 * erf(dabs(t1)/h/sqrt2)
     *           + dabs(t2) * 0.5 * erf(dabs(t2)/h/sqrt2))
 90      continue
         x(i) = tmp /nsum + 0.5
 100  continue
      end
      
cccccccccc End of ofccdf.f cccccccccc

cccccccccc FORTRAN subroutine remp.f cccccccccc

c     To draw a random sample from a kernel estimate of a distribution
c     Last changed: 18 Oct 2010

      subroutine remp(ny,y,f,a,b,nx,Fx,x,u,n)
      integer n,ny,nx,i,j,k,l,icounter
      double precision y(ny),f(ny),a(ny),b(ny),x(nx),Fx(nx),u(n),
     *     F0,F1,y0,dy,t1

      icounter=0
      do 100 i=1,ny
c     search for F1 and y1
         y0 = y(i) + a(i)
         dy=x(nx)-x(1)
         do 10 j=1,nx
            t1 = abs(x(j)-y0)
            if(t1<dy) then
               dy = t1
               k=j
            endif
 10      continue
         F0 = Fx(k) 
c     search for F1 and y1
         y0 = y(i) + b(i)
         dy=x(nx)-x(1)
         do 30 j=1,nx
            t1 = abs(x(j)-y0)
            if(t1<dy) then
               dy = t1
               k=j
            endif
 30      continue
         F1 = Fx(k) 
c     temporary using y0 to keep information
         do 90 l=1,INT(f(i))
            icounter = icounter+1
            y0 = u(icounter)*F1+(1.0-u(icounter))*F0
            dy=1.0
            do 60 j=1,nx
               t1 = abs(Fx(j)-y0)
               if(t1<dy) then
                  dy = t1
                  k=j
               endif
 60         continue
            u(icounter) = x(k)
 90      continue
 100  continue
      end
      
cccccccccc End of remp.f cccccccccc


C ===========================================================================
      subroutine nrlogit(x1,beta,p,n)

      implicit integer (I-N) 
      double precision beta(3),p(n),fx,dfx,x1,logitp
      J=0
      fx = 1.0
      DO 100 I = 1, n
         logitp = LOG(p(I)/(1-p(I)))
 10      IF((J .LT. 100) .AND. (ABS(fx) .GT. 0.000001)) THEN
            fx = beta(1)-logitp+beta(2)*x1+beta(3)*LOG(x1)
            dfx = beta(2)+beta(3)/x1
            x1 = x1 - fx/dfx
            IF(x1<0) THEN x1 = 0.001
            J = J+1
            GOTO 10
         ENDIF

         p(I) = x1
 100  END DO
      END
