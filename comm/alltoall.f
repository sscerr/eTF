      SUBROUTINE alltoallv_real( a, ssize, sdisp, b, rsize, rdisp )
      IMPLICIT NONE
      include 'mpif.h'
      REAL*8  a( * ), b( * )
      INTEGER ssize( * ), sdisp( * ), rsize( * ), rdisp( * )
      INTEGER ierr
      CALL MPI_ALLTOALLV( a, ssize, sdisp, MPI_DOUBLE_PRECISION, 
     &        b, rsize, rdisp, MPI_DOUBLE_PRECISION, 
     &        MPI_COMM_WORLD, ierr )

      RETURN
      END 
