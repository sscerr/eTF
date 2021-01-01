! on T3e substitute MPI_DOUBLE_PRECISION with MPI_REAL
! on other architecture substitute MPI_REAL with MPI_DOUBLE_PRECISION

!------------------------------------------------------------------------------!

      SUBROUTINE BCAST_REAL(ARRAY,N,root)
        IMPLICIT NONE
        INCLUDE 'mpif.h'
        REAL*8 array(*)
        INTEGER n, ierr, root
        CALL MPI_BCAST(array,n,MPI_DOUBLE_PRECISION,root,
     &    MPI_COMM_WORLD,ierr)
        RETURN
      END SUBROUTINE BCAST_REAL 

!------------------------------------------------------------------------------!
!
      SUBROUTINE BCAST_INTEGER(ARRAY,N,root)
        IMPLICIT NONE
        INCLUDE 'mpif.h'
        INTEGER ARRAY(*), N
        INTEGER IERR, root
        CALL MPI_BCAST(ARRAY,N,MPI_INTEGER,root,MPI_COMM_WORLD,IERR)
        RETURN
      END SUBROUTINE BCAST_INTEGER

!
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!

      SUBROUTINE BCAST_CHARACTER(ARRAY,N,root)
        IMPLICIT NONE
        INCLUDE 'mpif.h'
        CHARACTER*(*) ARRAY
        INTEGER N, root
        INTEGER I, IERR, NTMP
        CHARACTER*256 TMP
        INTEGER*4 TMP_CLINK(256/4)
        EQUIVALENCE (TMP,TMP_CLINK) 
        IF(N.GT.256) THEN
          WRITE(*,*)'ROOT_BCAST_CHARACTER SIZE OUT OF RANGE ',N
          CALL hangup
          stop
        END IF
        DO I=1,N
          TMP(I:I) = ARRAY(I:I)
        END DO
        NTMP = N/4
        IF(MOD(N,4).GT.0) THEN
          NTMP = NTMP + 1
        END IF
        CALL MPI_BCAST(TMP_CLINK,NTMP,MPI_INTEGER4,root, 
     &    MPI_COMM_WORLD,IERR)
        DO I=1,N
          ARRAY(I:I) = TMP(I:I)
        END DO
        RETURN
      END  SUBROUTINE BCAST_CHARACTER  

!
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!


      SUBROUTINE BCAST_LOGICAL(ARRAY,N,root)
        IMPLICIT NONE
        INCLUDE 'mpif.h'
        LOGICAL ARRAY(*)
        INTEGER N
        INTEGER IERR,I,root
        INTEGER IARRAY(N)
        DO I=1,N
          IF(ARRAY(I)) THEN
            IARRAY(I) = 1
          ELSE
            IARRAY(I) = 0
          END IF
        END DO
        CALL MPI_BCAST(IARRAY,N,MPI_INTEGER,root,MPI_COMM_WORLD,IERR)
        DO I=1,N
          IF(IARRAY(I).EQ.1) THEN
            ARRAY(I) = .TRUE.
          ELSE
            ARRAY(I) = .FALSE.
          END IF
        END DO
        RETURN
      END SUBROUTINE BCAST_LOGICAL 
