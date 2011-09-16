cccccccccc FORTRAN subroutine edf.f cccccccccc

c     To compute the empirical distribution function

c     Last changed: 24 July 2012

      subroutine wedf(X,W,n,xgrid,fhat,m)
      double precision X(*),W(*),xgrid(*),fhat(*)
      integer n,m,i,j,k
      double precision tmp
      
c     Initialize grid counts to zero
      do 10 i=1,m
         fhat(i) = dble(0)
 10   enddo

      tmp = dble(0)
      k = 1
      do 100 i=1,m
         fhat(i) = tmp
         do 50 j=k,n
            if(X(j) .LE. xgrid(i)) then
               fhat(i) = fhat(i) + W(j)
            else
               goto 60
            endif
 50      enddo
 60      k = j
         tmp = fhat(i)
 100  enddo
      
      do 300 i=1,m
         fhat(i) = fhat(i)/fhat(m)
 300  enddo
      
      return
      end

cccccccccc End of edf.f cccccccccc
