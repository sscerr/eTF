      SUBROUTINE gather_real(SEND,RECV,SIZE,root)
      IMPLICIT NONE
      include 'mpif.h'
      INTEGER SIZE, ROOT
      REAL*8 SEND(*), RECV(*)
      INTEGER IERR
      CALL MPI_GATHER(send,size,MPI_DOUBLE_PRECISION,
     &  recv,size,MPI_DOUBLE_PRECISION,ROOT,MPI_COMM_WORLD,IERR)
      RETURN
      END SUBROUTINE

      SUBROUTINE gatherv_real( SEND, ssize, RECV, 
     &                         rsize, rdisp, root )
      IMPLICIT NONE
      include 'mpif.h'
      INTEGER ssize, ROOT
      INTEGER rsize(*), rdisp(*)
      REAL*8  SEND(*), RECV(*)
      INTEGER IERR
      CALL MPI_GATHERV( send, ssize, MPI_DOUBLE_PRECISION,
     &  recv, rsize, rdisp, MPI_DOUBLE_PRECISION, ROOT, 
     &  MPI_COMM_WORLD, IERR )
      RETURN
      END SUBROUTINE
