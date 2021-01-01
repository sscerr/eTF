! on T3e substitute MPI_DOUBLE_PRECISION with MPI_REAL
! on other architecture substitute MPI_REAL with MPI_DOUBLE_PRECISION

!------------------------------------------------------------------------------!

      SUBROUTINE BCAST_REAL(ARRAY,N,root)
        IMPLICIT NONE
        REAL*8 array(*)
        INTEGER n, ierr, root
        RETURN
      END SUBROUTINE BCAST_REAL 

!------------------------------------------------------------------------------!
!
      SUBROUTINE BCAST_INTEGER(ARRAY,N,root)
        IMPLICIT NONE
        INTEGER ARRAY(*), N
        INTEGER IERR, root
        RETURN
      END SUBROUTINE BCAST_INTEGER

!
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!

      SUBROUTINE BCAST_CHARACTER(ARRAY,N,root)
        IMPLICIT NONE
        CHARACTER*(*) ARRAY
        INTEGER N, root
        RETURN
      END  SUBROUTINE BCAST_CHARACTER  

!
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!


      SUBROUTINE BCAST_LOGICAL(ARRAY,N,root)
        IMPLICIT NONE
        LOGICAL ARRAY(*)
        INTEGER N
        INTEGER IERR,I,root
        RETURN
      END SUBROUTINE BCAST_LOGICAL 
