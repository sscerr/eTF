       subroutine  trac(RHS_x)

!----------------------------------- !
!     M. Faganello, 2007             !
!                                    !
!   Trasporto passivo: calcolo di    !
!                                    !
!     -(U * Nabla) Tracciante        !
!                                    !
!----------------------------------- !

!**************************************************
!  MPI PARALLEL VERSION: VALENTINI-FAGANELLO 2009
!   3D PARALLEL VERSION: FAGANELLO 2010
!**************************************************


use parameter_mod
use box_mod
use deriv_mod
use fields_DP_mod
use fields_UJ_mod
use dom_distr_mod, only: nxl, nyt, nzl

IMPLICIT NONE

integer :: ix, iy, iz

REAL(dp),DIMENSION (nxl,ny,nz) :: RHS_x 
REAL(dp), ALLOCATABLE :: Trt(:,:,:)
REAL(dp), ALLOCATABLE:: zy(:), zx(:), zz(:)

ALLOCATE( Trt( nx, nyt,nz ) )  ! nyt ~ ny/nproc 
allocate(zx(nx),zy(ny),zz(nz))


 CALL traspdist( Tracciante, Trt, 1 )

do iz = 1, nz 
 do iy = 1, nyt
  call  derx_1( Trt(:,iy,iz), zx )
  Trt(:,iy,iz) = zx
 enddo
enddo
 
 CALL traspdist( RHS_x, Trt, -1 )

RHS_x = - Ux * RHS_x


do iz = 1, nz 
 do ix = 1, nxl
  call  dery_1(Tracciante(ix,:,iz), zy)
  RHS_x(ix,:,iz) = RHS_x(ix,:,iz) - Uy(ix,:,iz) * zy
 enddo
enddo




deallocate(zx)
deallocate(zy)
deallocate(zz)
deallocate(Trt)



end subroutine

