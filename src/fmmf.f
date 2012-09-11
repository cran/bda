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
         nsum = nsum + f(i)
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
      
cccccccccc End of iterfx.f cccccccccc

c  Part of R package KernSmooth
c  Copyright (C) 1995  M. P. Wand
c
c  Unlimited use and distribution (see LICENCE).

cccccccccc FORTRAN subroutine rlbin.f cccccccccc

c Obtains bin counts for univariate regression data
c via the linear binning strategy. If "trun=0" then
c weight from end observations is given to corresponding
c end grid points. If "trun=1" then end observations
c are truncated.

c Last changed: 26 MAR 2009

      subroutine rlbin(X,Y,n,a,b,M,trun,xcnts,ycnts)
      double precision X(*),Y(*),a,b,xcnts(*),ycnts(*),lxi,delta,rem
      integer n,M,i,li,trun

c     Initialize grid counts to zero

      do 10 i=1,M
         xcnts(i) = dble(0)
         ycnts(i) = dble(0)
10    continue

      delta = (b-a)/(M-1)
      do 20 i=1,n
         lxi = ((X(i)-a)/delta) + 1

c        Find integer part of "lxi"

         li = int(lxi) 
         rem = lxi - li
         if (li.ge.1.and.li.lt.M) then
            xcnts(li) = xcnts(li) + (1-rem)
            xcnts(li+1) = xcnts(li+1) + rem
            ycnts(li) = ycnts(li) + (1-rem)*y(i)
            ycnts(li+1) = ycnts(li+1) + rem*y(i)
         endif

         if (li.lt.1.and.trun.eq.0) then
            xcnts(1) = xcnts(1) + 1
            ycnts(1) = ycnts(1) + y(i)
         endif      
  
         if (li.ge.M.and.trun.eq.0) then 
               xcnts(M) = xcnts(M) + 1
               ycnts(M) = ycnts(M) + y(i)
         endif

20    continue

      return
      end

cccccccccc End of rlbin.f cccccccccc


c  simplified to implement the kde with epnichkov kernel (BW, 2011/10/09)

cccccccccc FORTRAN subroutine linbin.f cccccccccc

c Obtains bin counts for univariate data
c via the linear binning strategy. If "trun=0" then
c weight from end observations is given to corresponding
c end grid points. If "trun=1" then end observations
c are truncated.

c Last changed: 20 MAR 2009

      subroutine linbin(X,n,a,b,M,gcnts)
      double precision X(*),a,b,gcnts(*),lxi,delta,rem
      integer n,M,i,li

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

         if (li.lt.1) then
            gcnts(1) = gcnts(1) + 1
         endif

         if (li.ge.M) then
            gcnts(M) = gcnts(M) + 1
         endif

20    continue

      return
      end

cccccccccc End of linbin.f cccccccccc



c  Modified by Bin Wang (2012/07/03)

c  Part of R package KernSmooth
c  Copyright (C) 1995  M. P. Wand
c
c  Unlimited use and distribution (see LICENCE).

C     ======Souroutines to be called =================================
C     To sampling from a KDE.  F(x) is known over x, where x and F(x)
C     are both non-decreasing. For another non-decreasing vector
C     y0~unif(0,1), find the corresponding vector x0.

      SUBROUTINE findxyz(y, x, n, y0, m)
      implicit integer (I-N) 
      double precision y(n), x(n), y0(m), z(m), t

      J = 1
      DO 200 I=1,m
 100     IF(J > n) THEN
            z(I) = x(n)
         ELSE IF(y0(I) .GT. y(J)) THEN
            J = J + 1
            GOTO 100
         ELSE
            IF(J .EQ. 1) THEN
               z(I) = x(1)
            ELSE
               z(I) = x(J) - (x(J)-x(J-1))*(y(J)-y0(I))/(y(J)-y(J-1))
            ENDIF
         ENDIF
 200  ENDDO

      DO 250 I=1,m
         y0(I) = z(I)
 250  ENDDO

      RETURN
      END
