      SUBROUTINE scatter_real(SEND,RECV,SIZE,root)
      IMPLICIT NONE
      include 'mpif.h'
      INTEGER SIZE, ROOT
      REAL*8 SEND(*), RECV(*)
      INTEGER IERR
      CALL MPI_SCATTER(send,size,MPI_DOUBLE_PRECISION,
     &  recv,size,MPI_DOUBLE_PRECISION,ROOT,MPI_COMM_WORLD,IERR)
      RETURN
      END

