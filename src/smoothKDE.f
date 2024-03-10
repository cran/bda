cccccccccc FORTRAN subroutine binPRO.f cccccccccc

c redistribute interval type PRO data via binning
c Last changed: 11 Oct 2023

      subroutine probin(X,Y,n,a,b,M,gcnts)
      double precision X(*),Y(*),a,b,gcnts(*),lxi,lyi,delta,rem,dx
      integer n,M,i,li,lj

c     Initialize grid counts to zero

      do 10 i=1,M
         gcnts(i) = dble(0)
10    continue

      delta = (b-a)/M
      do 30 i=1,n
         lxi = ((X(i)-a)/delta) + 1
         lyi = ((Y(i)-a)/delta) + 1
         dx = (Y(i) - X(i))/delta

c        Find integer part of "lxi"

         li = int(lxi) 
         lj = int(lyi)

         if(li.eq.lj) then
            gcnts(li) = gcnts(li) + 1
         else
            do 20 j=li,lj
               if (j.lt.lj.and.j.eq.li) then
                  rem = lxi - li
                  gcnts(j) = gcnts(j) + (1-rem)/dx
               elseif (j.lt.lj.and.j.gt.li) then
                  gcnts(j) = gcnts(j) + 1
               else
                  rem = lyi - lj
                  gcnts(j) = gcnts(j) + rem/dx
               endif
 20         continue

         endif

 30   continue

      return
      end
cccccccccc End of binPRO.f cccccccccc

ccccccccccFORTRAN subroutine smoothkde cccccccccc

c     2013/04/01.

      subroutine smoothkde(fx,x0,n,x,f,m,w,h,iter)
      integer m,n,i,j,iter, k,i1,i2,loop
      double precision nsum, fx(n),f0(n),x0(n),x(m),f(m),h,w
      double precision wd, dx, fk(n,n),gk(n),df
      double precision t1,t2, xl,xu

      DOUBLE PRECISION PI,TWOPI, twoh, twopih
      PARAMETER(PI=3.141592653589793D0,TWOPI=2D0*PI)
      
      wd = 2.0 * w
      dx = x0(2) - x0(1)
      twoh = -0.5/h/h
      twopih = 1.0/sqrt(twopi)/h

      nsum = 0.
      do i=1, m
         nsum = nsum + f(i)
      end do
      do i=1, n
         f0(i) = fx(i)
         gk(i) = 0.0
      end do
      do i=1, n
         t1 = dx*(i-1.0)/h
         gk(i) = twopih * exp(-0.5*t1*t1)
      end do

      do i=1, n
         do j=1, n
            fk(i,j) = 0.0
         enddo
      enddo

      do i=1, n-1
         fk(i,i) = gk(1)
         do j=i+1, n
            fk(i,j) = gk(j-i)
            fk(j,i) = gk(j-i)
         enddo
      enddo
*      fk(n,n) = gk(1)

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
 400           enddo
               fx(k) = fx(k) + f(i)*t1/t2/nsum
 450        enddo
            df = df + (fx(k) - f0(k))**2.0
            f0(k) = fx(k)
 500     enddo
      enddo
      iter = loop

      return
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


c  Part of R package BDA
c  Copyright (C) 2009-2010 Bin Wang
c
c  Unlimited use and distribution (see LICENCE).

cccccccccc FORTRAN subroutine iterfx.f cccccccccc

c     2011/10/09.

      subroutine iterfx(fx,x0,n,x,f,m,w,h,iter)
      integer m,n,i,j,iter,k,i1,i2,loop
      double precision nsum,fx(n),f0(n),x0(n),x(m),f(m),h,w
      double precision wd, dx, fk(n,n),gk(n),df
      double precision t1,t2, xl,xu

      DOUBLE PRECISION PI,TWOPI, twoh, twopih
      PARAMETER(PI=3.141592653589793D0,TWOPI=2D0*PI)

      wd = 2.0 * w
      dx = x0(2) - x0(1)
      twoh = -0.5/h/h
      twopih = 1.0/sqrt(twopi)/h

      nsum = 0.
      do 50 i=1, m
         nsum = nsum + f(i)
 50   enddo
      do 100 i=1, n
         f0(i) = fx(i)
 100   enddo
      do 150 i=1, n
         gk(i) = twopih * exp(dx*(i-1.0)**2.0)
 150   enddo

       do i=1, n
          do j=1, n
             fk(i,j) = 0.0
          enddo
       enddo


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

