

! on T3e substitute MPI_DOUBLE_PRECISION with MPI_REAL
! on other architecture substitute MPI_REAL with MPI_DOUBLE_PRECISION

      subroutine sendrecv(array,size,dest,sour,ip)
        implicit none
        REAL*8 array(*)
        integer size, dest, sour, ip
        RETURN
      END SUBROUTINE

      subroutine fddsendrecv(fddin,fddout,dest,sour,ip)
        implicit none
        REAL*8 fddin(2), fddout(2)
        integer dest, sour, ip, ierr
        fddout(1) = fddin(1) 
        fddout(2) = fddin(2) 
        RETURN
      END SUBROUTINE
