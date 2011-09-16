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
