subroutine flr_i

!*****************************************************************!
! COMPUTE: ions' gyroviscous tensor G                             !
!          [2D version, assumes B "mainly" along z: Bperp << Bz]  !
! --------------------------------------------------------------- ! 
! S. S. Cerri, 2011                                               !
!                                                                 !
!                                                                 !
!                       Gi_ij  i,j=x,y,z                          !
!                                                                 !
! N.B.: Gi_ij is a SYMMETRIC tensor with ZERO TRACE, so one needs !  
!       only 5 linear independent components instead of the whole !
!       9 components of the tensor:                               ! 
!                                                                 !
!         i.e., Gi_xx, Gi_zz, Gi_xy, Gi_xz, and Gi_yz             !
!                                                                 !
!       Moreover, we consider a coordinate system such that the   !
!       z-axis is aligned with the mean magnetic field <B>, i.e.  !
!                                                                 !
!                  Gi_zz = 0  and  Gi_xx = -Gi_yy                 !
!                                                                 !
!    for further details, see:                                    !
!                                                                 !
!    *** Cerri et al., Phys. Plasmas 20, 112112 (2013) ***        !
!                                                                 !
!*****************************************************************!

use parameter_mod
use box_mod
use deriv_mod
use fields_UJ_mod
use fields_DP_mod
use fields_EB_mod
use dom_distr_mod
use parallel_mod

implicit none

integer :: ix, iy, iz

real(dp), allocatable :: dxui_j(:), dyui_j(:), dzui_j(:)
real(dp), allocatable :: Ut(:,:,:)


allocate(dxui_j(nx))    
allocate(dyui_j(ny))    
allocate(dzui_j(nz))    
allocate( Ut( nx, nyt, nz ) )  ! nyt ~ ny/nproc


!------------------------------------------------------!
!  Compute:  Gi_xx = - Gi_yy                           !
!                                                      !
!  Gi_xx = -0.5*(pi_perp/|B|)*[d(Uix)/dy + d(Uiy)/dx]  !
!                                                      !
!------------------------------------------------------!

call traspdist( Uiy, Ut, 1 )
do iz = 1, nz
  do iy = 1, nyt
    call derx_1( Ut(:,iy, iz), dxui_j )
    Ut(:,iy, iz) = dxui_j
  enddo
enddo
call traspdist( Gi_xx, Ut, -1 )
 
do iz = 1, nz
  do ix = 1, nxl
    call dery_1( Uix(ix,:,iz), dyui_j)
    Gi_xx(ix,:,iz) = - 0.5d0 * pi_perp(ix,:,iz) * ( Gi_xx(ix,:,iz) + dyui_j )
  enddo
enddo

Gi_xx = Gi_xx / dsqrt(B00)


!------------------------------------------------------!
!  Compute:  Gi_xy = Gi_yx                             !
!                                                      !
!  Gi_xy = -0.5*(pi_perp/|B|)*[d(Uiy)/dy - d(Uix)/dx]  !
!                                                      !
!------------------------------------------------------!

call traspdist( Uix, Ut, 1 )
do iz = 1, nz
  do iy = 1, nyt
    call derx_1( Ut(:,iy, iz), dxui_j )
    Ut(:,iy, iz) = dxui_j
  enddo
enddo
call traspdist( Gi_xy, Ut, -1 )

do iz = 1, nz
  do ix = 1, nxl
    call dery_1( Uiy(ix,:,iz), dyui_j)
    Gi_xy(ix,:,iz) = 0.5d0 * pi_perp(ix,:,iz) * (Gi_xy(ix,:,iz) - dyui_j)
  enddo
enddo

Gi_xy = Gi_xy / dsqrt(B00)


!--------------------------------------------------!
!  Compute:  Gi_xz = Gi_zx                         !
!                                                  !
!  Gi_xz = -2.0*(pi_para/|B|)*d(Uiy)/dz            !
!         - (pi_perp/|B|)*[d(Uiz)/dy - d(Uiy)/dz]  !
!                                                  !
! [ 2D version: Gi_xz = -(pi_perp/|B|)*d(Uiz)/dy ] !
!--------------------------------------------------!

do iz = 1, nz
  do ix = 1, nxl
    call dery_1( Uiz(ix,:,iz), dyui_j)
    Gi_xz(ix,:,iz) = -pi_perp(ix,:,iz)*dyui_j
  enddo
enddo

Gi_xz = Gi_xz / dsqrt(B00)


!--------------------------------------------------!
!  Compute:  Gi_yz = Gi_zy                         !
!                                                  !
!  Gi_yz = 2.0*(pi_para/|B|)*d(Uix)/dz             !
!         + (pi_perp/|B|)*[d(Uiz)/dx - d(Uix)/dz]  !
!                                                  ! 
! [ 2D version: Gi_yz = (pi_perp/|B|)*d(Uiz)/dx ]  ! 
!--------------------------------------------------!

call traspdist( Uiz, Ut, 1 )
do iz = 1, nz
  do iy = 1, nyt
    call derx_1( Ut(:,iy, iz), dxui_j )
    Ut(:,iy, iz) = dxui_j
  enddo
enddo
call traspdist( Gi_yz, Ut, -1 )

Gi_yz = pi_perp * Gi_yz
Gi_yz = Gi_yz / dsqrt(B00)


!----------------------------------------------------! 
! *** taking care of "omega*B" asymmetry ***         !
!  (arises from sign of Omega_ci)                    ! 
!                                                    !
! see: Cerri et al., Phys. Plasmas 20, 112112 (2013) !
!----------------------------------------------------! 
do iz = 1, nz
  do iy = 1, nyl
    do ix = 1, nxl
      if (Bz(ix,iy,iz) < 0.0) then
        Gi_xx(ix,iy,iz) = - Gi_xx(ix,iy,iz)
        Gi_xy(ix,iy,iz) = - Gi_xy(ix,iy,iz)
        Gi_xz(ix,iy,iz) = - Gi_xz(ix,iy,iz)
        Gi_yz(ix,iy,iz) = - Gi_yz(ix,iy,iz)
      endif
    enddo
  enddo
enddo


!-------------------------------------------------------!
! smoothing -> no gyroviscous terms at the x-boundaries !
!-------------------------------------------------------!
do ix = 1, nxl               
  Gi_xx(ix,:,:) = lambda(ixlg + ix - 1) * Gi_xx(ix,:,:)
  Gi_xy(ix,:,:) = lambda(ixlg + ix - 1) * Gi_xy(ix,:,:)
  Gi_xz(ix,:,:) = lambda(ixlg + ix - 1) * Gi_xz(ix,:,:)
  Gi_yz(ix,:,:) = lambda(ixlg + ix - 1) * Gi_yz(ix,:,:)
enddo  


deallocate(dxui_j)
deallocate(dyui_j)
deallocate(dzui_j)
deallocate(Ut)


end subroutine 
