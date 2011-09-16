c  Part of R package BDA
c  Copyright (C) 2009-2010 Bin Wang
c
c  Unlimited use and distribution (see LICENCE).

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
         do 90 l=1,f(i)
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
