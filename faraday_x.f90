       subroutine  faraday_x(RHS_x)

!------------------------------- !
!       F. Califano, 2006        !
!                                !
!   calcolo di - rot(E)_x        !
!                                !
!------------------------------- !

!**************************************************
!  MPI PARALLEL VERSION: VALENTINI-FAGANELLO 2009
!   3D PARALLEL VERSION: FAGANELLO 2010
!**************************************************

use parameter_mod
use box_mod
use deriv_mod
use fields_EB_mod
use dom_distr_mod

IMPLICIT NONE

integer :: ix, iy, iz

REAL(dp),DIMENSION (nxl,nyl,nzl) :: RHS_x

REAL(dp), ALLOCATABLE :: zz(:)

allocate(zz(nz))



! Calcolo componente x

do iz = 1, nz 
 do ix = 1, nxl
  call  dery_1(Ez(ix,:,iz), RHS_x(ix,:,iz))
 enddo
enddo

  RHS_x = - RHS_x


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! eventuale correz. su divB nel bulk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALL traspdist( Bx, At, 1 )
!do iz = 1, nz 
! do iy = 1, nyt
!  call  derx_1( At(:,iy,iz), zx )
!  At(:,iy,iz) = zx
! enddo
!enddo
! CALL traspdist( AA, At, -1 )
!do iz = 1, nz 
! do ix = 1, nxl
!  call  dery_1(By(ix,:,iz), zy)
!  AA(ix,:,iz) = AA(ix,:,iz) + zy
! enddo
!enddo
!RHS_x = RHS_x - Uex * AA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








! Iperviscosita' su Bx...........eventualmente aggiornare per versione 2D/3D parallela


deallocate(zz)

end subroutine
