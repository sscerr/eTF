       SUBROUTINE PARALLEL_STARTUP(NPROC,MPIME,ROOT,GROUP) 
       IMPLICIT NONE
       INTEGER NPROC
       INTEGER MPIME, I, ERR, ROOT, GROUP
!      ---------------------------
!      INITIALIZE MPI ENVIRONEMENT 
!      ---------------------------
       nproc = 1
       mpime = 0
       ROOT = 0
       GROUP = 0
       RETURN
       END

       SUBROUTINE PARALLEL_COMMLIB( parallal_build )
          LOGICAL :: parallal_build
          parallal_build = .false.
          return
       END SUBROUTINE

