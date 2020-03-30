cccccccccc FORTRAN subroutine ofcpdf.f cccccccccc

c     To compute the pdf based on a Berkson measurement error with
c     uniform rounding errors

c     Last changed: 05 Oct 2010

      subroutine ofcpdf(y,f,a,b,ny,x,nx,h)
      integer ny,nx,i,j
      double precision nsum, y(ny),f(ny),a(ny),b(ny),x(nx),h,
     *     tmp, t1,t2,sqrt2h,sqrt2
      sqrt2 = sqrt(2.0)
      sqrt2h = sqrt2 * h

      nsum=0.
      do 10 i=1,ny
         nsum = nsum + f(i)
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
      integer ny,nx,i,j
      double precision nsum, y(ny),f(ny),a(ny),b(ny),x(nx),h,
     *     tmp, t1,t2,sqrt2h,sqrt2,PI,PI2,h2
      PI = 3.141592653589793238462643
      PI2 = 1./DSQRT(2.*PI)
      sqrt2 = sqrt(2.0)
      sqrt2h = sqrt2 * h
      h2 = h * h

      nsum=0.
      do 10 i=1,ny
         nsum = nsum + f(i)
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

cccccccccc FORTRAN subroutine mlensimp.f cccccccccc

c     To compute the pdf based on a Berkson measurement error with
c     uniform rounding errors.  In this version, both parameters will be
c     found numerically.  Check the simplified version.

c     Last changed: 05 Oct 2010

      subroutine mlensimp(w,f,a,b,n,theta)
      integer n,i,iter
      double precision w(n),f(n),a(n),b(n),theta(2),
     *     t2,t3,t4,t5,mu,sig,sig1,g0,g2,
     *     re2,za(n),zb(n),pdfa(n),pdfb(n),cdfa(n),cdfb(n)

c     Set initial values of theta
      mu = theta(1)
      sig = theta(2)

      iter=0
      re2 = 1.0
      g0=0.0
      g2=0.0

      do while(iter<10000 .and. re2 >0.000001)
         do 30 i=1,n
            za(i) = (w(i)+a(i)-mu)/sig
            zb(i) = (w(i)+b(i)-mu)/sig
            call dnorm(pdfa(i),za(i))
            call dnorm(pdfb(i),zb(i))
            call pnorm(cdfa(i),za(i))
            call pnorm(cdfb(i),zb(i))
            t2 = cdfb(i)-cdfa(i)
            t3 = zb(i)*pdfb(i)-za(i)*pdfa(i)
            t5 = zb(i)*zb(i)*zb(i)*pdfb(i)-za(i)*za(i)*za(i)*pdfa(i)
            g0 = g0 + f(i)*sig * t3/t2
            g2 = g2 + f(i)*(t5*t2+t3*t3)/t2/t2
 30      continue
         sig1 = sig
         sig = sig1 - g0/g2
         t2 = dabs((sig-sig1)/dmin1(sig,sig1))
         t4 = dabs((sig-sig1))
         re2 = dmax1(t2,t4)
         iter = iter + 1
      enddo
      theta(2) = sig
      n = iter
      end
      
cccccccccc End of mlensimp.f cccccccccc


cccccccccc FORTRAN subroutine mlenrom.f cccccccccc

c     To compute the pdf based on a Berkson measurement error with
c     uniform rounding errors.  In this version, both parameters will be
c     found numerically.  Check the simplified version.

c     Last changed: 05 Oct 2010


      subroutine mlenorm(w,f,a,b,n,theta)
      integer n,i,iter
      double precision w(n),f(n),a(n),b(n),theta(2),
     *     t1,t2,t3,t4,t5,mu,sig,mu1,sig1,f0,g0,f1,f2,g1,g2,
     *     re1,re2,za(n),zb(n),pdfa(n),pdfb(n),cdfa(n),cdfb(n)

c     Set initial values of theta
      mu = theta(1)
      sig = theta(2)

      iter=0
      re1 = 1.0
      re2 = 1.0
      f0= 0.
      f1=0.
      f2=0.0
      g0=0.0
      g1=0.0
      g2=0.0

      do while(iter<10000 .and. (re1>0.000001 .or. re2 >0.000001))
         do 30 i=1,n
            za(i) = (w(i)+a(i)-mu)/sig
            zb(i) = (w(i)+b(i)-mu)/sig
            call dnorm(pdfa(i),za(i))
            call dnorm(pdfb(i),zb(i))
            call pnorm(cdfa(i),za(i))
            call pnorm(cdfb(i),zb(i))
            t1 = pdfb(i)-pdfa(i)
            t2 = cdfb(i)-cdfa(i)
            t3 = zb(i)*pdfb(i)-za(i)*pdfa(i)
            t4 = zb(i)*zb(i)*pdfb(i)-za(i)*za(i)*pdfa(i)
            t5 = zb(i)*zb(i)*zb(i)*pdfb(i)-za(i)*za(i)*za(i)*pdfa(i)
            f0 = f0 + f(i)*t1/t2
            g0 = g0 + f(i)*sig * t3/t2
            f1 = f1 + f(i)*(t3*t2+t1*t1)/t2/t2
            f2 = f2 + f(i)*(t4*t2+t3*t1)/t2/t2
            g1 = g1 + f(i)*((t4-t1)*t2+t3*t1)/t2/t2
            g2 = g2 + f(i)*(t5*t2+t3*t3)/t2/t2
 30      continue
         f1 = f1/sig
         f2 = f2/sig
         mu1 = mu
         sig1 = sig
         mu = mu1 - (f0*g2-f2*g0)/(f1*g2-f2*g1)
         sig = sig1 - (f1*g0-f0*g1)/(f1*g2-f2*g1)
         t1 = dabs((mu-mu1)/dmin1(mu,mu1))
         t2 = dabs((sig-sig1)/dmin1(sig,sig1))
         t3 = dabs((mu-mu1))
         t4 = dabs((sig-sig1))
         re1 = dmax1(t1,t3)
         re2 = dmax1(t2,t4)
         iter = iter +1
      enddo

      theta(1) = mu
      theta(2) = sig
      n = iter
      end
      
cccccccccc End of mlenorm.f cccccccccc

CCCCCCCCCCCCCCCCCCCCC  Subroutines to be called 

      SUBROUTINE dnorm(fx,x)
      double precision fx, x,PI
      PI = 3.141592653589793238462643
      fx = dexp(-0.5*x*x)/sqrt(2.*PI)
      return
      end

      SUBROUTINE pnorm(Fx,x)
      double precision Fx, x, x2,PI
      PI = 3.141592653589793238462643
      x2 = x/sqrt(2.)
      Fx = 0.5 + 0.5*erf(x2)
      return
      end
