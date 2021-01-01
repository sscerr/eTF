      SUBROUTINE alltoallv_real( a, ssize, sdisp, b, rsize, rdisp )
      IMPLICIT NONE
      REAL*8  a( * ), b( * )
      INTEGER ssize( * ), sdisp( * ), rsize( * ), rdisp( * )
      INTEGER ierr, i
      DO i = 1, ssize(1)
         b(i) = a(i)
      END DO
      RETURN
      END 
