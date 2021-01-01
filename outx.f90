subroutine outx

!*** OUTPUT: grid-based fields ***!


!----------------------------------------------------!
!  ORIGINAL MODEL:    F. Califano, 2006              !
!                     M. Faganello, 2008             !
!                                                    !
!  Faganello et al., New J. Phys. 11, 063008 (2009)  !
!                                                    !
!  ================================================  !
!  EXTENDED MODEL:    S. S. Cerri, 2011              !
!                                                    !
!  Cerri et al., Phys. Plasmas 20, 112112 (2013)     !
!                                                    !
!--------------------------------------------------- !

!**************************************************!
!  MPI PARALLEL VERSION: VALENTINI-FAGANELLO 2009  !
!  3D PARALLEL VERSION: FAGANELLO 2010             !
!  3D ANISOTROPIC/FLR H-MHD/TF: CERRI 2011         !
!**************************************************!


use parameter_mod
use box_mod
use fields_UJ_mod
use fields_DP_mod
use fields_EB_mod
use dom_distr_mod
use parallel_mod


IMPLICIT NONE

integer :: ix, iy, iz

write(l_U_me, 821) ixlg, ixlg+nxl-1, tempo
write(l_EB_me, 821) ixlg, ixlg+nxl-1, tempo
write(l_DPJ_me, 821) ixlg, ixlg+nxl-1, tempo
write(l_Gi_me, 821) ixlg, ixlg+nxl-1, tempo

do iz=1,nzl, nwz
 do iy=1,nyl, nwy
  do ix=1,nxl, nwx

   write(l_U_me, 806) Uex(ix,iy,iz), Uey(ix,iy,iz), Uez(ix,iy,iz) &
                    , Uix(ix,iy,iz), Uiy(ix,iy,iz), Uiz(ix,iy,iz)
  enddo
 enddo
enddo

do iz=1,nzl, nwz
 do iy=1,nyl, nwy
  do ix=1,nxl, nwx
   write(l_EB_me, 806) Ex(ix,iy,iz), Ey(ix,iy,iz), Ez(ix,iy,iz) &
                     , Bx(ix,iy,iz), By(ix,iy,iz), Bz(ix,iy,iz) 
  enddo
 enddo
enddo

do iz=1,nzl, nwz
 do iy=1,nyl, nwy
  do ix=1,nxl, nwx
   write(l_DPJ_me, 806) Den(ix,iy,iz), Tracciante(ix,iy,iz)  &
                      , pe_para(ix,iy,iz), pe_perp(ix,iy,iz) &
                      , pi_para(ix,iy,iz), pi_perp(ix,iy,iz)
  enddo
 enddo
enddo

do iz=1,nzl, nwz
 do iy=1,nyl, nwy
  do ix=1,nxl, nwx

   write(l_Gi_me, 804) Gi_xx(ix,iy,iz), Gi_xy(ix,iy,iz) &
                     , Gi_xz(ix,iy,iz), Gi_yz(ix,iy,iz)
  enddo
 enddo
enddo


802   format(1x, 2(1x, 1e11.4))
803   format(1x, 3(1x, 1e11.4))
804   format(1x, 4(1x, 1e11.4))
805   format(1x, 5(1x, 1e11.4))
806   format(1x, 6(1x, 1e11.4))
808   format(1x, 8(1x, 1e11.4))
809   format(1x, 9(1x, 1e11.4))
812   format(1x, 6(1x, 1e11.4))
821   format(1x, 2(1x, 1i5), 1x, 1e11.4)

end subroutine
