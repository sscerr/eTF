       subroutine  cont_n_U(RHS_x)

!----------------------------------- !
!       F. Califano, 2006            !
!                                    !
!   Eq. Continuita': calcolo di      !
!                                    !
!          - Div[ nU ]               !
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
use dom_distr_mod
use parallel_mod

IMPLICIT NONE

integer :: ix, iy, iz

REAL(dp),DIMENSION (nxl,nyl,nzl) :: RHS_x 
REAL(dp), ALLOCATABLE :: Ut(:,:,:)
REAL(dp), ALLOCATABLE :: zx(:), zy(:), zz(:)

ALLOCATE( Ut( nx, nyt, nz ) )  ! nyt ~ ny/nproc 

allocate(zy(ny))
allocate(zx(nx))
allocate(zz(nz))


! Il segno meno viene messo nel secondo pezzo 
!    RHS_x(ix,:) = - RHS_x(ix,:) - zy

 CALL traspdist( nU_x, Ut, 1 )

do iz = 1, nz 
 do iy = 1, nyt
  call  derx_1( Ut(:,iy, iz), zx )
  Ut(:,iy, iz) = zx
 enddo
enddo

 CALL traspdist( RHS_x, Ut, -1 )

IF (ibc .EQ. 4) then
 if ( mpime == 0 ) then 

! N.B. va sommato col meno

 RHS_x(1,:,:) = - Den(1,:,:) * (c_s_inv(1,:,:) * L_0(1,:,:) + 0.5d0 * &
                    (alpha_2(1,:,:) * c_s_inv(1,:,:) * &
          (L_s_p(1,:,:) + L_s_m(1,:,:)) - alpha_1(1,:,:) * f_inv(1,:,:) * &
                      (L_f_p(1,:,:) + L_f_m(1,:,:))))

 endif
 if ( mpime == nprow -1 ) then 

 RHS_x(nxl,:,:) = - Den(nxl,:,:) * (c_s_inv(2,:,:) * L_0(2,:,:) + 0.5d0 * &
                    (alpha_2(2,:,:) * c_s_inv(2,:,:) * &
          (L_s_p(2,:,:) + L_s_m(2,:,:)) - alpha_1(2,:,:) * f_inv(2,:,:) * &
                      (L_f_p(2,:,:) + L_f_m(2,:,:))))

 endif
ENDIF

do iz = 1, nz 
 do ix = 1, nxl
  call  dery_1(nU_y(ix,:,iz), zy)
  RHS_x(ix,:,iz) = - RHS_x(ix,:,iz) - zy
 enddo
enddo


deallocate(zz)
deallocate(zy)
deallocate(zx)
deallocate(Ut)

end subroutine
