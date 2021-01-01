       SUBROUTINE PARALLEL_STARTUP(NPROC,MPIME,ROOT,GROUP) 
       IMPLICIT NONE
       INCLUDE 'mpif.h'
       INTEGER NPROC
       INTEGER MPIME, I, ERR, ROOT, GROUP
!      ---------------------------
!      INITIALIZE MPI ENVIRONEMENT 
!      ---------------------------
       CALL MPI_INIT(ERR)  
       CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPROC,ERR)
       CALL MPI_COMM_RANK(MPI_COMM_WORLD,MPIME,ERR)
       ROOT = 0
       GROUP = MPI_COMM_WORLD
       RETURN
       END

       SUBROUTINE PARALLEL_COMMLIB( parallal_build )
          LOGICAL :: parallal_build
          parallal_build = .true.
          return
       END SUBROUTINE

