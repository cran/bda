c  Part of R package BDA
c  Copyright (C) 2009-2010 Bin Wang
c
c  Unlimited use and distribution (see LICENCE).

cccccccccc FORTRAN subroutine mlensimp.f cccccccccc

c     To compute the pdf based on a Berkson measurement error with
c     uniform rounding errors.  In this version, only the standard
c     deviation will be estimated.  We assume the sample mean is a good
c     one.

c     Last changed: 05 Oct 2010

      subroutine mlensimp(w,f,a,b,n,theta)
      integer n,i,j,iter,nsum
      double precision w(n),f(n),a(n),b(n),theta(2),
     *     t1,t2,t3,t4,t5,mu,sig,mu1,sig1,g0,g2
     *     re2,za(n),zb(n),pdfa(n),pdfb(n),cdfa(n),cdfb(n)

c     Set initial values of theta
      mu = theta(1)
      sig = theta(2)

      iter=0
      re2 = 1.0

      do while(iter<10000 .and. re2 >0.000001)
         g0=0.0
         g2=0.0
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
         t2 = abs((sig-sig1)/min(sig,sig1))
         t4 = abs((sig-sig1))
         re2 = max(t2,t4)
         iter = iter +1
      enddo
      theta(2) = sig
      n = iter
      end
      
cccccccccc End of mlensimp.f cccccccccc


cccccccccc FORTRAN subroutine mlenorm.f cccccccccc

c     Last changed: 17 Feb 2011

      subroutine remlenorm(w,f,b,n,theta)
      integer n,i,j,iter
      double precision w(n),f(n),b(n),theta(2),
     *     t1,t2,t3,t4,t5,mu,sig,mu1,sig1,f0,g0,f1,f2,g1,g2
     *     re1,re2,za(n),zb(n),pdfa(n),pdfb(n),cdfa(n),cdfb(n)

c     Set initial values of theta
      mu = theta(1)
      sig = theta(2)

      iter=0
      re1 = 1.0
      re2 = 1.0

      do while(iter<10000 .and. (re1>0.000001 .or. re2 >0.000001))
         f0 = 0.0
         f1 = 0.0
         f2 = 0.0
         g0 = 0.0
         g1 = 0.0
         g2 = 0.0
         do 30 i=1,n
            za(i) = (w(i)-b(i)-mu)/sig
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
            if(t2 .NE. 0.0) then
               f0 = f0 + f(i)*t1/t2
               g0 = g0 + f(i)*sig * t3/t2
               f1 = f1 + f(i)*(t3*t2+t1*t1)/t2/t2
               f2 = f2 + f(i)*(t4*t2+t3*t1)/t2/t2
               g1 = g1 + f(i)*((t4-t1)*t2+t3*t1)/t2/t2
               g2 = g2 + f(i)*(t5*t2+t3*t3)/t2/t2
            endif
 30      continue
         f1 = f1/sig
         f2 = f2/sig
         mu1 = mu
         sig1 = sig
         mu = mu1 - (f0*g2-f2*g0)/(f1*g2-f2*g1)
         sig = sig1 - (f1*g0-f0*g1)/(f1*g2-f2*g1)
         t1 = abs((mu-mu1)/min(mu,mu1))
         t2 = abs((sig-sig1)/min(sig,sig1))
         t3 = abs((mu-mu1))
         t4 = abs((sig-sig1))
         re1 = max(t1,t3)
         re2 = max(t2,t4)
         iter = iter +1
      enddo

      theta(1) = mu
      theta(2) = sig
      n = iter
      end
      
cccccccccc End of mlenorm.f cccccccccc


CCCCCCCCCCCCCCCCCCCCC  Subroutines to be called 

      SUBROUTINE dnorm(fx,x)
      double precision fx, x
      DATA PI/3.141592653589793238462643d0/
      fx = dexp(-0.5*x*x)/sqrt(2.*PI)
      return
      end

      SUBROUTINE pnorm(Fx,x)
      double precision Fx, x, x2
      DATA PI/3.141592653589793238462643d0/
      x2 = x/sqrt(2.)
      Fx = 0.5 + 0.5*erf(x2)
      return
      end
