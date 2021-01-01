
       SUBROUTINE SYNCRONIZE 
       IMPLICIT NONE
       include 'mpif.h'
       INTEGER IERR
       CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
       RETURN 
       END 
