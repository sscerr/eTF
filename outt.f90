       subroutine outt

!--------------------------------------------------------- !
!             F. Califano, 2006                            !
!                                                          !
!               Uscite campi integrati                     !
!                                                          !
!--------------------------------------------------------- !

!**************************************************
!  MPI PARALLEL VERSION: VALENTINI-FAGANELLO 2009
!**************************************************

use parameter_mod
use box_mod
use fields_UJ_mod
use fields_DP_mod
use fields_EB_mod
use dom_distr_mod
use parallel_mod

IMPLICIT NONE

INTEGER  :: ix, iy, iz
REAL(dp) :: E_ke, E_ki, E_cin, E_magn, E_elett, divB, divE

REAL(dp), ALLOCATABLE :: AA(:,:,:)
REAL(dp), ALLOCATABLE :: zx(:), zy(:), zz(:)
REAL(dp), ALLOCATABLE :: Bt(:,:,:), Et(:,:,:)

allocate(zx(nx))
allocate(zy(ny))
allocate(zz(nz))
allocate(AA(nxl, ny, nz))
ALLOCATE( Bt( nx, nyt, nz ) )  ! nyt ~ ny/nproc
ALLOCATE( Et( nx, nyt, nz ) )  ! nyt ~ ny/nproc

E_ke    = 0.0 
E_ki   = 0.0 
E_cin   = 0.0 
E_magn  = 0.0 
E_elett = 0.0


AA = Uex * Uex + Uey * Uey + Uez * Uez


do iz = 1, nz, nwz
 do iy = 1, ny, nwy 
  E_ke = E_ke + sum(AA(:,iy,iz))
 enddo
enddo

AA = Uix * Uix + Uiy * Uiy + Uiz * Uiz

do iz = 1, nz, nwz
 do iy = 1, ny, nwy
  E_ki = E_ki + sum(AA(:,iy,iz))
 enddo
enddo

AA = Ux * Ux + Uy * Uy + Uz * Uz

do iz = 1, nz, nwz
 do iy = 1, ny, nwy
 E_cin = E_cin +sum(AA(:,iy,iz))
 enddo
enddo

AA = Bx * Bx + By * By + Bz * Bz

do iz = 1, nz, nwz
 do iy = 1, ny, nwy
  E_magn = E_magn +sum(AA(:,iy,iz))
 enddo
enddo

AA = Ex * Ex + Ey * Ey + Ez * Ez

do iz = 1, nz, nwz
 do iy = 1, ny, nwy
  E_elett = E_elett + sum(AA(:,iy,iz))
 enddo
enddo


E_ke    = E_ke    * dxyz
E_ki    = E_ki    * dxyz
E_cin   = E_cin   * dxyz
E_magn  = E_magn  * dxyz
E_elett = E_elett * dxyz



! DIV(B)

 CALL traspdist( Bx, Bt, 1 )

do iz = 1, nz
 do iy = 1, nyt
  call  derx_1( Bt(:,iy,iz), zx )
  Bt(:,iy,iz) = zx
 enddo
enddo

 CALL traspdist( AA, Bt, -1 )

do iz = 1, nz
 do ix = 1, nxl
  call  dery_1(By(ix,:,iz), zy)
  AA(ix,:,iz) = AA(ix,:,iz) + zy
 enddo
enddo


do iz = 1, nz
 do iy = 1, ny
  zy(iy) = sum(abs(AA(:,iy,iz)))
 enddo
enddo

divB = sum(zy) * Tot_inv_nxnynz 
! DIV(E)

 CALL traspdist( Ex, Et, 1 )

do iz = 1, nz
 do iy = 1, nyt
  call  derx_1( Et(:,iy,iz), zx )
  Et(:,iy,iz) = zx
 enddo
enddo

 CALL traspdist( AA, Et, -1 )

do iz = 1, nz
 do ix = 1, nxl
  call  dery_1(Ey(ix,:,iz), zy)
  AA(ix,:,iz) = AA(ix,:,iz) + zy
 enddo
enddo


do iz = 1, nz
 do iy = 1, ny
  zy(iy) = sum(abs(AA(:,iy,iz)))
 enddo
enddo

divE = sum(zy) * Tot_inv_nxnynz 

deallocate(AA)
deallocate(zx)
deallocate(zy)
deallocate(zz)
deallocate(Et)
deallocate(Bt)


        call PARALLEL_SUM_REAL(E_ke,1)
        call PARALLEL_SUM_REAL(E_ki,1)
        call PARALLEL_SUM_REAL(E_cin,1)
        call PARALLEL_SUM_REAL(E_magn,1)
        call PARALLEL_SUM_REAL(E_elett,1)
        call PARALLEL_SUM_REAL(divB,1)
        call PARALLEL_SUM_REAL(divE,1)

if(mpime==0) then
write(*,  806) tempo, divB, divE, E_cin, E_magn, E_elett
write(l_div,  802) divB, divE
write(l_Eng,  805) E_ke, E_ki, E_cin, E_magn, E_elett
endif


802   format(1x, 2(1x, 1e11.4))
803   format(1x, 3(1x, 1e11.4))
804   format(1x, 4(1x, 1e11.4))
805   format(1x, 5(1x, 1e11.4))
806   format(1x, 6(1x, 1e11.4))
812   format(1x, 6(1x, 1e11.4))

end subroutine
