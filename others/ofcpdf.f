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

