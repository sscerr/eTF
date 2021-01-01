                  subroutine moto_z(RHS_x)

!-------------------------------------------------------------!
!             F. Califano, 2006                               !
!                                                             !
!            Eq. MOTO componente z:                           !
!                                                             !
!          Div[ BB - n(UU + de^2 jj) - P_tot * Ident ]        !
!                                                             !
! Axz = - Den * (Uix * Uiz + de2 * Uex * Uez) + Bx Bz         !
! Ayz = - Den * (Uiy * Uiz + de2 * Uey * Uez) + By Bz         !
! Azz = - Den * (Uiz * Uiz + de2 * Uez * Uez) + Bz Bz - P_tot !
!-------------------------------------------------------------!


!**************************************************
!  MPI PARALLEL VERSION: VALENTINI-FAGANELLO 2009
!   3D PARALLE VERSION: FAGANELLO 2010
!  LANDAU-FLUID VERSION  FAGANELLO 2011
!
! FLR-Landau-Fluid VERSION: Cerri 2011
!
!**************************************************


use parameter_mod
use box_mod
use deriv_mod
use fields_UJ_mod
use fields_DP_mod
use fields_EB_mod
use dom_distr_mod
use parallel_mod

IMPLICIT NONE

integer :: ix, iy, iz

REAL(dp),DIMENSION (nxl,nyl,nzl) :: RHS_x

REAL(dp), ALLOCATABLE :: zy(:),zx(:), zz(:)
REAL(dp), ALLOCATABLE :: AA(:,:,:)
REAL(dp), ALLOCATABLE :: At(:,:,:), B_inv(:,:,:)
REAL(dp), ALLOCATABLE :: Ptot(:,:,:)

ALLOCATE( At( nx, nyt, nz ) )  ! nyt ~ ny/nproc 
allocate(zy(ny))
allocate(zx(nx))
allocate(zz(nz))
allocate(AA( nxl, ny, nz )) 
allocate(Ptot( nxl, ny, nz ))
allocate(B_inv(nxl, ny, nz))

Ptot = Bx**2 + By**2 + Bz**2
B_inv = 1.0d0 / Ptot 

 !*************************************************************!
 !  Calcolo del primo termine MOTO_z da derivare rispetto a x  !
 !*************************************************************!

AA = Bx * Bz * ( 1.0d0 - B_inv * ( pe_para + pi_para - pe_perp - pi_perp ) ) &
       - Den * (Uix * Uiz + de2 * Uex * Uez) - Gi_xz ! = Axz 

 CALL traspdist( AA, At, 1 )

do iz = 1, nz 
 do iy = 1, nyt
  call  derx_1( At(:,iy,iz), zx )
  At(:,iy,iz) = zx
 enddo
enddo

 CALL traspdist( RHS_x, At, -1 )



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! eventuali correz. su divB nel bulk
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
!RHS_x = RHS_x - Bz * AA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




! N.B per ora tralascio ogni termime de2....corretto solo MHD
! N.B. in questo caso ho messo, nei punti 1,nx U e non Ui/Ue

IF (ibc .EQ. 4) then
 if ( mpime == 0 ) then 

 RHS_x(1,:,:) = 0.5 * Den(1,:,:) * ( Uz(1,:,:) * (2. * L_0(1,:,:) *c_s_inv(1,:,:) + & 
                alpha_2(1,:,:) * c_s_inv(1,:,:) *&
               (L_s_p(1,:,:) +  L_s_m(1,:,:)) - alpha_1(1,:,:) * f_inv(1,:,:) *&
                 (L_f_p(1,:,:) +  L_f_m(1,:,:))) &
     + a_inv(1,:,:) * ( - ay(1,:,:) * (L_a_p(1,:,:) + L_a_m(1,:,:)) + az(1,:,:) * &
   (segno(1,:,:) * alpha_1(1,:,:)* (L_s_p(1,:,:) - L_s_m(1,:,:)) &
    + ax(1,:,:) * alpha_2(1,:,:) * f_inv(1,:,:) * (L_f_p(1,:,:) - L_f_m(1,:,:)))))  

!call derx_BC(Bz,grad)
!RHS_x(1,:,:) =  Bx(1,:,:) * grad(1,:,:)
do iz = 1, nz
 call dery_1(By(1,:,iz),zy)
 RHS_x(1,:,iz) = RHS_x(1,:,iz) - Bz(1,:,iz) * zy
enddo

 endif
 if ( mpime == nprow -1 ) then 

 RHS_x(nxl,:,:) = 0.5 * Den(nxl,:,:) * ( Uz(nxl,:,:) * (2. * L_0(2,:,:) *c_s_inv(2,:,:) + & 
                alpha_2(2,:,:) * c_s_inv(2,:,:) *&
               (L_s_p(2,:,:) +  L_s_m(2,:,:)) - alpha_1(2,:,:) * f_inv(2,:,:) *&
                 (L_f_p(2,:,:) +  L_f_m(2,:,:))) &
     + a_inv(2,:,:) * ( - ay(2,:,:) * (L_a_p(2,:,:) + L_a_m(2,:,:)) + az(2,:,:) * &
   (segno(2,:,:) * alpha_1(2,:,:)* (L_s_p(2,:,:) - L_s_m(2,:,:)) &
    + ax(2,:,:) * alpha_2(2,:,:) * f_inv(2,:,:) * (L_f_p(2,:,:) - L_f_m(2,:,:)))))  

!call derx_BC(Bz,grad)
!RHS_x(nxl,:,:) = Bx(nxl,:,:) * grad(nxl,:,:)
do iz = 1, nz
 call dery_1(By(nxl,:,iz),zy)
 RHS_x(nxl,:,iz) = RHS_x(nxl,:,iz) - Bz(nxl,:,iz) * zy
enddo

 endif
ENDIF

 !*************************************************************!
 ! Calcolo del secondo termine MOTO_z da derivare rispetto a y !
 !*************************************************************!

AA = By * Bz * ( 1.0d0 - B_inv * (pe_para + pi_para - pe_perp - pi_perp ) ) &
        - Den * (Uiy * Uiz + de2 * Uey * Uez) - Gi_yz ! = Ayz

do iz = 1, nz 
 do ix = 1, nxl
  call  dery_1(AA(ix,:,iz), zy)
  RHS_x(ix,:,iz) = RHS_x(ix,:,iz) + zy
 enddo
enddo

 !*************************************************************!
 ! Calcolo del terzo termine MOTO_z da derivare rispetto a z   !
 !                                                             !
 ! Ptot = pe_perp + pi_perp + B^2 / 2                          !
 !                                                             !
 ! PAI_zz = 0 per ipotesi (c.magnetico principale lungo z)     !
 ! (vedi commenti in flr_i.f90 per maggiori dettagli)          !
 !*************************************************************!


RHS_x = RHS_x * de12inv

deallocate(zx)
deallocate(zy)
deallocate(zz)
deallocate(AA)
deallocate(At)
deallocate(Ptot)
deallocate(B_inv)

end subroutine
