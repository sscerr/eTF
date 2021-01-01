       subroutine u_ei

! --------------------------------------------------- !
!           F. Califano, February 2006                !
!                                                     !
!     Calcolo delle velocita' degli elettroni e ioni  !
!                                                     !
!     u_e = U - j / n;    u_i = U + de^2 j / n        !
!          ( de^2 = me / mi )             !           !
!                                                     !
! --------------------------------------------------- !

!**************************************************
!  MPI PARALLEL VERSION: VALENTINI-FAGANELLO 2009
!**************************************************

use box_mod
use parameter_mod
use parallel_mod
use dom_distr_mod
use fields_DP_mod
use fields_UJ_mod

IMPLICIT NONE

integer :: ix, iy, iz


! expressions without smoothing
!Uex = Ux - de12inv * Jx * Dinv
!Uey = Uy - de12inv * Jy * Dinv
!Uez = Uz - de12inv * Jz * Dinv
!Uix = Ux + de12inv * de2 * Jx * Dinv
!Uiy = Uy + de12inv * de2 * Jy * Dinv
!Uiz = Uz + de12inv * de2 * Jz * Dinv

! smoothing to ideal-MHD towards boundaries
do iz = 1, nzl
 do iy = 1, nyl
  do ix = 1, nxl
    Uex(ix,iy,iz) = Ux(ix,iy,iz) - de12inv * &
     lambda(ixlg + ix - 1) * Jx(ix,iy,iz) * Dinv(ix,iy,iz)
    Uey(ix,iy,iz) = Uy(ix,iy,iz) - de12inv * &
     lambda(ixlg + ix - 1) * Jy(ix,iy,iz) * Dinv(ix,iy,iz)
    Uez(ix,iy,iz) = Uz(ix,iy,iz) - de12inv * &
     lambda(ixlg + ix - 1) * Jz(ix,iy,iz) * Dinv(ix,iy,iz)
    Uix(ix,iy,iz) = Ux(ix,iy,iz) + de12inv * de2 * &
     lambda(ixlg + ix - 1) * Jx(ix,iy,iz) * Dinv(ix,iy,iz)
    Uiy(ix,iy,iz) = Uy(ix,iy,iz) + de12inv * de2 * &
     lambda(ixlg + ix - 1) * Jy(ix,iy,iz) * Dinv(ix,iy,iz)
    Uiz(ix,iy,iz) = Uz(ix,iy,iz) + de12inv * de2 * &
     lambda(ixlg + ix - 1) * Jz(ix,iy,iz) * Dinv(ix,iy,iz)
  enddo
 enddo
enddo





end subroutine
