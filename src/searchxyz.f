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
