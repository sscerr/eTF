subroutine  corrente

!****************************************************************!
! COMPUTE: current density from J = curl(B)                      !
!          [2D version]                                          !
! -------------------------------------------------------------- ! 
! original version:             F. Califano, 2006                !
! MPI parallel version:         M. Faganello/F. Valentini, 2008  !                  
! Anisotropic/FLR-MHD version:  S. S. Cerri, 2011                !
!****************************************************************!

use parameter_mod
use parallel_mod
use dom_distr_mod
use box_mod
use deriv_mod
use fields_UJ_mod
use fields_EB_mod

IMPLICIT NONE

integer :: ix, iy, iz

REAL(dp),ALLOCATABLE :: zx(:), zy(:), zz(:)
REAL(dp), ALLOCATABLE :: Bt(:,:,:)

allocate(zx(nx))
allocate(zy(ny))
allocate(zz(nz))

ALLOCATE( Bt( nx, nyt, nz ) )  ! nyt ~ ny/nproc


!--compute Jx
do iz = 1, nz 
  do ix = 1, nxl
    call  dery_1( Bz(ix,:,iz), Jx(ix,:,iz) )
  enddo
enddo

!--compute Jy
call traspdist( Bz, Bt, 1 ) 
do iz = 1, nz 
  do iy = 1, nyt
    call  derx_1(Bt(:,iy,iz), zx )
    Bt(:,iy,iz) = - zx
  enddo
enddo
call traspdist( Jy, Bt, -1 ) 

!--compute Jz
call traspdist( By, Bt, 1 ) 
do iz = 1, nz 
  do iy = 1, nyt
    call  derx_1( Bt(:,iy,iz), zx )
    Bt(:,iy,iz) = zx
  enddo
enddo
call traspdist( Jz, Bt, -1 ) 
!
do iz = 1, nz 
  do ix = 1, nxl
    call dery_1( Bx(ix,:,iz), zy )
    Jz(ix,:,iz) = Jz(ix,:,iz) - zy
  enddo
enddo

deallocate(zx)
deallocate(zy)
deallocate(zz)
deallocate(Bt)

end subroutine
