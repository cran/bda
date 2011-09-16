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

      subroutine remp(ny,y,f,a,b,nx,Fx,x,u)
      integer ny,nx,i,j,k,l,icounter
      double precision y(ny),f(ny),a(ny),b(ny),x(nx),Fx(nx),u(ny),
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
