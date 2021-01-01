SUBROUTINE traspdist( a, b, idir )

   use parameter_mod
   use dom_distr_mod
   use parallel_mod

   IMPLICIT NONE
   include 'mpif.h'

   REAL(dp) :: a( nxl, ny, nz ), b( nx, nyt, nz )
   INTEGER  :: idir

   INTEGER  :: ip, srsize, ib
   INTEGER  :: ix, iy, iz, ierr

   REAL(dp), ALLOCATABLE :: rbuf( : )
   REAL(dp), ALLOCATABLE :: sbuf( : )


   ALLOCATE( rbuf( nxlmax * nytmax * nz * nprow) )
   ALLOCATE( sbuf( nxlmax * nytmax * nz * nprow) )

srsize = nxlmax*nytmax*nz

!  do iy = 1, ny
!    do ix = 1, nxl
!      a(ix,iy,:) = (mpime + 1) * iy
!    enddo
!  enddo

!  do iy = 1, ny
!    do ix = 1, nxl
!      write(*,*) a(ix,iy,1)
!    enddo
!  enddo
!  write(*,*) 'puppa1'


   IF( idir > 0 ) THEN

do ip = 1, nprow
 ib = 0
  do iz = 1, nz
    do iy = 1, nytp(ip)
      do ix = 1, nxl
       ib = ib + 1
        sbuf(ib+(ip-1)*srsize) = a(ix,iy + iytgp(ip) - 1,iz)
      enddo
    enddo
  enddo
enddo


      CALL MPI_ALLTOALL( sbuf, srsize, MPI_DOUBLE_PRECISION, & 
             rbuf, srsize, MPI_DOUBLE_PRECISION, &
             MPI_COMM_WORLD, ierr )



     do ip = 1, nprow
      ib = 0  
       do iz = 1, nz
         do iy = 1, nyt
           do ix = 1, nxlp(ip)
            ib = ib + 1
             b(ix + ixlgp(ip) - 1, iy, iz) = rbuf(ib + (ip-1)*srsize)     
           enddo
         enddo 
       enddo
     enddo


!do iy = 1, nyt
!  do ix = 1, nx
!    write(*,*) b(ix,iy,1)
!  enddo
!enddo
!write(*,*) 'puppa2'


   ELSE

     do ip = 1, nprow
      ib = 0  
       do iz = 1, nz
         do iy = 1, nyt
           do ix = 1, nxlp(ip)
            ib = ib + 1
             rbuf(ib + (ip-1)*srsize) = b(ix + ixlgp(ip) - 1, iy, iz)
           enddo
         enddo 
       enddo
     enddo

      CALL MPI_ALLTOALL( rbuf, srsize, MPI_DOUBLE_PRECISION, & 
             sbuf, srsize, MPI_DOUBLE_PRECISION, &
             MPI_COMM_WORLD, ierr )


do ip = 1, nprow
 ib = 0
  do iz = 1, nz
    do iy = 1, nytp(ip)
      do ix = 1, nxl
       ib = ib + 1
        a(ix,iy + iytgp(ip) - 1,iz) = sbuf(ib+(ip-1)*srsize)
      enddo
    enddo
  enddo
enddo

!do iy = 1, ny
!  do ix = 1, nxl
!    write(*,*) a(ix,iy,1)
!  enddo
!enddo
!write(*,*) 'puppa3'
!call syncronize
!stop


   END IF

   DEALLOCATE( rbuf )
   DEALLOCATE( sbuf )

   RETURN
END SUBROUTINE
