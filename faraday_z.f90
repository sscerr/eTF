       subroutine  faraday_z(RHS_x)

!------------------------------- !
!       F. Califano, 2006        !
!                                !
!   calcolo di - rot(E)_z        !
!                                !
!------------------------------- !

!**************************************************
!  MPI PARALLEL VERSION: VALENTINI-FAGANELLO 2009
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

REAL(dp),DIMENSION (nxl,nyl,nzl) :: RHS_x
REAL(dp), ALLOCATABLE :: Et(:,:,:), Ett(:,:,:)
REAL(dp), ALLOCATABLE :: zx(:), zy(:), zz(:)

allocate(zx(nx))
allocate(zy(ny))
allocate(zz(nz))
ALLOCATE( Et(nx,nyt,nz), Ett(nxl,ny,nz) )  ! nyt ~ ny/nproc 


! Calcolo componente z

do iz = 1, nz 
 do ix = 1, nxl
  call  dery_1(Ex(ix,:,iz), RHS_x(ix,:,iz))
 enddo
enddo


 CALL traspdist( Ey, Et, 1 )

do iz = 1, nz 
 do iy = 1, nyt
  call  derx_1( Et(:,iy,iz), zx )
  Et(:,iy,iz) = zx
 enddo
enddo

 CALL traspdist( Ett, Et, -1 )


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! eventuali correz. su divB   nel bulk  
!        N.B. Ett viene poi sommato con il meno                   !!!!!!!!!!
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
!Ett = Ett + Uez * AA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






IF (ibc .EQ. 4) then

! N.B. va sommato col meno

 if ( mpime == 0 ) then 
 Ett(1,:,:) = - 0.5 * dsqrt(Den(1,:,:)) * a_inv(1,:,:) * &
       ( segno(1,:,:) * ay(1,:,:) * (L_a_p(1,:,:) - L_a_m(1,:,:)) - az(1,:,:) * &
        (c_s(1,:,:) * alpha_1(1,:,:)* f_inv(1,:,:) * (L_s_p(1,:,:) +  L_s_m(1,:,:)) + &
               alpha_2(1,:,:) * (L_f_p(1,:,:) + L_f_m(1,:,:))))


do iz = 1, nz
 call dery_1(By(1,:,iz),zy)
 Ett(1,:,iz) = Ett(1,:,iz) + Uz(1,:,iz) * zy
enddo


 endif
 if ( mpime == nprow - 1 ) then
 Ett(nxl,:,:) = - 0.5 * dsqrt(Den(nxl,:,:)) * a_inv(2,:,:) * &
       ( segno(2,:,:) * ay(2,:,:) * (L_a_p(2,:,:) - L_a_m(2,:,:)) - az(2,:,:) * &
        (c_s(2,:,:) * alpha_1(2,:,:)* f_inv(2,:,:) * (L_s_p(2,:,:) +  L_s_m(2,:,:)) + &
               alpha_2(2,:,:) * (L_f_p(2,:,:) + L_f_m(2,:,:))))


do iz = 1, nz
 call dery_1(By(nxl,:,iz),zy)
 Ett(nxl,:,iz) = Ett(nxl,:,iz) + Uz(nxl,:,iz) * zy
enddo



endif

ENDIF
  




RHS_x = RHS_x - Ett 



deallocate(zx)
deallocate(zy)
deallocate(zz)
deallocate(Et)
deallocate(Ett)


end subroutine
