       subroutine  faraday_y(RHS_y)

!------------------------------- !
!       F. Califano, 2006        !
!                                !
!   calcolo di - rot(E)_y        !
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
use fields_DP_mod
use fields_UJ_mod
use dom_distr_mod
use parallel_mod

IMPLICIT NONE

integer :: ix, iy, iz

REAL(dp),DIMENSION (nxl,nyl,nzl) :: RHS_y

REAL(dp), ALLOCATABLE :: Et(:,:,:)
REAL(dp),ALLOCATABLE :: zx(:), zy(:), zz(:)

allocate(zx(nx))
allocate(zy(ny))
allocate(zz(nz))
ALLOCATE( Et( nx, nyt, nz ) )  ! nyt ~ ny/nproc 

! Calcolo componente y 

 CALL traspdist( Ez, Et, 1 )

do iz = 1, nz 
 do iy = 1, nyt
  call  derx_1( Et(:,iy,iz), zx )
  Et(:,iy,iz) = zx
 enddo
enddo

 CALL traspdist( RHS_y, Et, -1 )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! eventuali correz. su divB   nel bulk
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
!RHS_y = RHS_y - Uey * AA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

deallocate(zx)
deallocate(Et)

IF (ibc .EQ. 4) then
 if ( mpime == 0 ) then 

 RHS_y(1,:,:) = 0.5 * dsqrt(Den(1,:,:)) * a_inv(1,:,:) * ( - segno(1,:,:) * az(1,:,:) *&
               (L_a_p(1,:,:) -  L_a_m(1,:,:)) &
           - ay(1,:,:) * (c_s(1,:,:) * alpha_1(1,:,:) * f_inv(1,:,:) *&
     (L_s_p(1,:,:) +  L_s_m(1,:,:)) + alpha_2(1,:,:) * (L_f_p(1,:,:) + L_f_m(1,:,:)))) 

do iz = 1, nz
 call dery_1(By(1,:,iz),zy)
 RHS_y(1,:,iz) = RHS_y(1,:,iz) - Uy(1,:,iz) * zy
enddo


 endif
 if ( mpime == nprow -1 ) then 

 RHS_y(nxl,:,:) = 0.5 * dsqrt(Den(nxl,:,:)) * a_inv(2,:,:) * (- segno(2,:,:) * az(2,:,:) *&
               (L_a_p(2,:,:) -  L_a_m(2,:,:)) &
           - ay(2,:,:) * (c_s(2,:,:) * alpha_1(2,:,:) * f_inv(2,:,:) *&
     (L_s_p(2,:,:) +  L_s_m(2,:,:)) + alpha_2(2,:,:) * (L_f_p(2,:,:) + L_f_m(2,:,:)))) 


do iz = 1, nz
 call dery_1(By(nxl,:,iz),zy)
 RHS_y(nxl,:,iz) = RHS_y(nxl,:,iz) - Uy(nxl,:,iz) * zy
enddo

 endif
ENDIF


deallocate(zy)
deallocate(zz)



end subroutine
