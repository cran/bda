c  Part of R package KDE
c  Copyright (C) 2010 Bin Wang
c
c  Unlimited use and distribution (see LICENCE).

cccccccccc FORTRAN subroutine discdata.f cccccccccc

c Discretized the data to find the weight sequence xi.  A vector of
c frequency "F" is used for groupded data.


c Last changed: 20 MAR 2009

      subroutine discdata(X,F,n,rb,a,b,M,xi)
      double precision X(*),F(*),rb(*),a,b,xi(*)
      double precision d1,d2,d3,delta,tsum,t,p
      integer n,M,i,j,k,l,ia,ib,ic

c     Initialize grid counts to zero

      do 10 i=1,M+1
         xi(i) = dble(0)
 10   continue
      delta = (b-a)/M

      DO 80 i = 1,n
         ia = floor((X(i)-rb(i)-a)/delta)
         ib = ceiling((X(i)+rb(i)-a)/delta)
         ic = ceiling((X(i)-a)/delta)
         IF(i .EQ. 1) THEN 
            ha = F(i)
            hb = F(i+1)
         ELSEIF(i .EQ. n) THEN
            ha = F(i-1)
            hb= F(n)
         ELSE
            ha = F(i-1)
            hb= F(i+1)
         ENDIF            
c     Compute the portion that will be redistributed uniformly
c     p = 0.5*(hb+ha+2.*F(i))/(max(hb,ha)+1.*F(i))
         p = (hb+ha+3.*F(i))/(hb+ha+3.*F(i)+abs(ha-hb))
         t = ib-ia+1.
         d2 = p*F(i)/t
         d1 = (1.-p)*F(i)*(ha)/(ha+hb)
         d3 = (1.-p)*F(i)*(hb)/(ha+hb)

         t = ic-ia+1.
         tsum = (1.+t)*t/2.
         if(ha .GE. F(i)) then
            DO 30 j=ia,ic
               xi(j) = xi(j) + (j-ia+1.)/tsum*d1 + d2
 30         ENDDO
         else
            DO 40 j=ia,ic
               xi(j) = xi(j) + (ic-j+1.)/tsum*d1 + d2
40         ENDDO
         endif

         t = ib-ic+1.
         tsum = (1.+t)*t/2.
         if(F(i) .GE. hb) then
            DO 50 j=ic,ib
               xi(j) = xi(j) + (j-ic+1.)/tsum*d3 + d2
 50         ENDDO
         else
            DO 60 j=ic,ib
               xi(j) = xi(j) + (ib-j+1.)/tsum*d3 + d2
 60         ENDDO
         endif

 80   ENDDO
 
      
      do 100 i=1,M+1
         xi(i) = xi(i)/n/delta
 100  continue
      
      return
      end

cccccccccc End of discdata.f cccccccccc

cccccccccc FORTRAN subroutine lscv.f cccccccccc

c     To select the optimal bandwidth using the FFT algorithm

c     void F77_SUB(lscv)(
c     double *bw, int *n, int *M, double *ab,
c     double *sl, double *dyl, int *ih, double *hs,double *lscv)

c     Last changed: 6 Jul 2011

      subroutine lscv(bw,n,m,ab,sl,dyl,ih,hs,cv)
      double precision sl(*),dyl(*),hs(*),cv(*)
      double precision bw,ab,t1,t2,h0,h,h2,hstep,fmin,adj
      integer n,m,ih,i,j,k,l
      real PI
      parameter (PI=3.1415927)


c     expand searching for up to 10 times.  If all failed, report
c     results as is.

      adj = 1.0
      h0 = bw
      do 100 k=1,10
         h =0.9*h0*adj
         hstep = 1.6*adj*h/ih
         do 80 j=1,ih
            fmin = -999.0
            t1 = 0.0
            h = h+hstep
            hs(j) = h
            h2 = h*h
            do 60 i=1,M/2
               t1=t1+(exp(-h2*sl(i))-2.*exp(-0.5*h2*sl(i)))*dyl(i)
 60         enddo
            t2 = n*h*sqrt(2.*PI)
            cv(j) = ab*t1 +1./t2
            IF(fmin .EQ. -999.0) THEN 
               fmin = cv(j)
               bw = h
               l=1
            ELSEIF(fmin .LT. cv(j)) THEN
               fmin = cv(j)
               bw = h
               l=j
            ENDIF
 80      enddo
         IF(l .EQ. 1) THEN
            h0 = 0.9 *bw
            adj = 1.
         ELSEIF(l .EQ. ih) THEN
            h0 = 1.4*bw
            adj = 1.
         ELSE
            h0 = bw
            adj = 0.95
         ENDIF
 100  continue
      
      return
      end

cccccccccc End of lscv.f cccccccccc

