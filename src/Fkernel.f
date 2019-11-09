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
         nsum = nsum + INT(f(i))
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
      integer M,l
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
      double precision y(ny),f(ny),a(ny),b(ny),x(nx),h,
     *     tmp, t1,t2,sqrt2h,sqrt2
      sqrt2 = sqrt(2.0)
      sqrt2h = sqrt2 * h

      nsum=0
      do 10 i=1,ny
         nsum = nsum + INT(f(i))
 10   continue
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
      double precision y(ny),f(ny),a(ny),b(ny),x(nx),h,
     *     tmp, t1,t2,sqrt2h,sqrt2,PI,PI2,h2
      PI = 3.1415926536D0
      PI2 = 1./DSQRT(2.*PI)
      sqrt2 = sqrt(2.0)
      sqrt2h = sqrt2 * h
      h2 = h * h

      nsum=0
      do 10 i=1,ny
         nsum = nsum + INT(f(i))
 10   continue
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
      k = 1
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

c  Part of R package BDA
c  Copyright (C) 2009-2010 Bin Wang
c
c  Unlimited use and distribution (see LICENCE).

cccccccccc FORTRAN subroutine iterfx.f cccccccccc

c     2011/10/09.

      subroutine iterfx(fx,x0,n,x,f,m,w,h,iter)
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
         nsum = nsum + INT(f(i))
 50   enddo
      do 100 i=1, n
         f0(i) = fx(i)
 100   enddo
      do 150 i=1, n
         gk(i) = twopih * exp(dx*(i-1.0)**2.0)
 150   enddo


      do 250 i=1, n
         fk(i,i) = gk(1)
         do 200 j=i+1, n
            fk(i,j) = gk(j-i)
            fk(j,i) = fk(i,j)
 200      enddo
 250   enddo

       loop = 0
       df = 1.0
       do while(loop .LT. iter .AND. df .GT. 0.0001)
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
      
cccccccccc End of iterfx.f cccccccccc


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

