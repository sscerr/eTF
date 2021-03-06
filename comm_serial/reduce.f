
! on T3e substitute MPI_DOUBLE_PRECISION with MPI_REAL
! on other architecture substitute MPI_REAL with MPI_DOUBLE_PRECISION

      SUBROUTINE PARALLEL_SUM_REAL(ARRAY,N)
        IMPLICIT NONE
        INTEGER N, ERR, I
        REAL*8 ARRAY(N)
        RETURN
      END SUBROUTINE PARALLEL_SUM_REAL

      SUBROUTINE PARALLEL_SUM_COMPLEX(ARRAY,N)
        IMPLICIT NONE
        INTEGER N, ERR, I
        COMPLEX*16 ARRAY(N)
        RETURN
      END SUBROUTINE PARALLEL_SUM_COMPLEX

!
      SUBROUTINE PARALLEL_SUM_INTEGER(ARRAY,N)
        IMPLICIT NONE
        INTEGER N, ERR, I
        INTEGER ARRAY(N)
      RETURN
      END SUBROUTINE PARALLEL_SUM_INTEGER
!
!
      SUBROUTINE PARALLEL_SUM_REAL_TO(ARRAY_IN,ARRAY_OUT,N)
        IMPLICIT NONE
        INTEGER N,ERR
        REAL*8 ARRAY_IN(N)
        REAL*8 ARRAY_OUT(N)
        ARRAY_OUT = ARRAY_IN
        RETURN
      END SUBROUTINE PARALLEL_SUM_REAL_TO


      SUBROUTINE PARALLEL_MAX_INTEGER(ARRAY,N)
        IMPLICIT NONE
        INTEGER N, I, IERR
        INTEGER ARRAY(N)
        RETURN
      END SUBROUTINE PARALLEL_MAX_INTEGER

      SUBROUTINE PARALLEL_MIN_INTEGER(ARRAY,N)
        IMPLICIT NONE
        INTEGER N, I, IERR
        INTEGER ARRAY(N)
        RETURN
      END SUBROUTINE PARALLEL_MIN_INTEGER


      SUBROUTINE PARALLEL_MAX_REAL(ARRAY,N)
        IMPLICIT NONE
        INTEGER N, IERR
        REAL*8  ARRAY(N)
        RETURN
      END SUBROUTINE PARALLEL_MAX_REAL

