       subroutine moto_xy(RHS_x, RHS_y)

!--------------------------------------------------------------!
!             F. Califano, 2006                                !
!                                                              !
!            Eq. MOTO componente x:                            !
!                                                              !
!          Div[ BB - n(UU + de^2 jj) - P_tot * Ident ]         !
!                                                              !
! Axx = - Den * (Uix * Uix + de2 * Uex * Uex) + Bx Bx - P_tot  !
! Ayy = - Den * (Uiy * Uiy + de2 * Uey * Uey) + Bx By - P_tot  !
! Axy = - Den * (Uix * Uiy + de2 * Uex * Uey) + Bx By ! = Ayx  !
! Axz = - Den * (Uix * Uiz + de2 * Uex * Uez) + Bz Bz          !
! Ayz = - Den * (Uiy * Uiz + de2 * Uey * Uez) + By Bz          !
!--------------------------------------------------------------!


!**************************************************
!  MPI PARALLEL VERSION: VALENTINI-FAGANELLO 2009
!   3D PARALLEL VERSION: FAGANELLO 2010
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

REAL(dp),DIMENSION (nxl,nyl,nzl) :: RHS_x, RHS_y

REAL(dp),ALLOCATABLE :: zx(:)
REAL(dp),ALLOCATABLE :: zy(:)
REAL(dp),ALLOCATABLE :: zz(:)
REAL(dp),ALLOCATABLE :: AA(:,:,:)
REAL(dp),ALLOCATABLE :: Ptot(:,:,:), B_inv(:,:,:)
REAL(dp), ALLOCATABLE :: At(:,:,:)

ALLOCATE( At( nx, nyt, nz ) )  ! nyt ~ ny/nproc 

allocate(zx(nx))
allocate(zy(ny))
allocate(zz(nz))
allocate(AA(nxl, ny, nz))
allocate(Ptot(nxl, ny, nz))
allocate(B_inv(nxl, ny, nz))

! Ptot = P_e + p_perp + B^2 / 2 

Ptot = Bx**2 + By**2 + Bz**2
B_inv = 1.0d0 / Ptot 
Ptot = 0.5 * Ptot + pe_perp + pi_perp

 !************************************************************!
 ! Calcolo del primo termine MOTO_x da derivare rispetto a x  !
 !************************************************************!

AA = Bx * Bx * ( 1.0d0 - B_inv * ( pe_para + pi_para - pe_perp - pi_perp ) ) & 
      -  Ptot - Den * (Uix * Uix + de2 * Uex * Uex) - Gi_xx

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
!RHS_x = RHS_x - Bx * AA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





! N.B per ora tralascio ogni termime de2....corretto solo MHD
! N.B. in questo caso ho messo, nei punti 1,nx U e non Ui/Ue


IF (ibc .EQ. 4) then
 if ( mpime == 0 ) then 


 RHS_x(1,:,:) = 0.5d0 * Den(1,:,:) * (Ux(1,:,:)* (2.0d0 * L_0(1,:,:) *c_s_inv(1,:,:) + & 
          alpha_2(1,:,:) * c_s_inv(1,:,:) * &
            (L_s_p(1,:,:) + L_s_m(1,:,:)) - alpha_1(1,:,:) * f_inv(1,:,:) *&
                 (L_f_p(1,:,:) + L_f_m(1,:,:))) &
    + s(1,:,:) * alpha_2(1,:,:) * c_s_inv(1,:,:) * (L_s_p(1,:,:) - L_s_m(1,:,:)) &
     - alpha_1(1,:,:) * (L_f_p(1,:,:) - L_f_m(1,:,:)) )     


!call derx_BC(By,grad)
!RHS_x(1,:,:) = - By(1,:,:) * grad(1,:,:)

do iz = 1, nz
 call dery_1(By(1,:,iz),zy)
 RHS_x(1,:,iz) = RHS_x(1,:,iz) - Bx(1,:,iz) * zy
enddo

 endif
 if ( mpime == nprow -1 ) then 

 RHS_x(nxl,:,:) = 0.5 * Den(nxl,:,:) * (Ux(nxl,:,:)* ( 2. * L_0(2,:,:) *c_s_inv(2,:,:) + & 
         alpha_2(2,:,:) * c_s_inv(2,:,:) * &
               (L_s_p(2,:,:) + L_s_m(2,:,:)) - alpha_1(2,:,:) * f_inv(2,:,:) *&
                 (L_f_p(2,:,:) + L_f_m(2,:,:))) &
    + s(2,:,:) * alpha_2(2,:,:) * c_s_inv(2,:,:) * (L_s_p(2,:,:) - L_s_m(2,:,:)) &
     - alpha_1(2,:,:) * (L_f_p(2,:,:) - L_f_m(2,:,:)) )     

!call derx_BC(By,grad)
!RHS_x(nxl,:,:) = - By(nxl,:,:) * grad(nxl,:,:)
!IF (nz.GT.1) THEN
 !call derx_BC(Bz,grad)
 !RHS_x(nxl,:,:) = - Bz(nxl,:,:) * grad(nxl,:,:)
!ENDIF
do iz = 1, nz
 call dery_1(By(nxl,:,iz),zy)
 RHS_x(nxl,:,iz) = RHS_x(nxl,:,iz) - Bx(nxl,:,iz) * zy
enddo


 endif              
ENDIF


 !*************************************************************!
 ! Calcolo del primo termine MOTO_y  da derivare rispetto a x  !
 !*************************************************************!

AA = Bx * By * ( 1.0d0 - B_inv * ( pe_para + pi_para - pe_perp + pi_perp ) ) &
       - Den * (Uix * Uiy + de2 * Uex * Uey) - Gi_xy

 CALL traspdist( AA, At, 1 )

do iz = 1, nz 
 do iy = 1, nyt
  call  derx_1( At(:,iy,iz), zx )
  At(:,iy,iz) = zx
 enddo
enddo

 CALL traspdist( RHS_y, At, -1 )


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! eventuali correz. su divB nel bulk
! CALL traspdist( Bx, At, 1 )
!do iz = 1, nz 
! do iy = 1, nyt
!  call  derx_1( At(:,iy,iz), zx )
!  At(:,iy,iz) = zx
! enddo
!enddo
! CALL traspdist( AB, At, -1 )
!do iz = 1, nz 
! do ix = 1, nxl
!  call  dery_1(By(ix,:,iz), zy)
!  AB(ix,:,iz) = AB(ix,:,iz) + zy
! enddo
!enddo
!RHS_y = RHS_y - By * AB
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




IF (ibc .EQ. 4) then
 if ( mpime == 0 ) then

 RHS_y(1,:,:) = 0.5 * Den(1,:,:) * (Uy(1,:,:)* ( 2. * L_0(1,:,:) *c_s_inv(1,:,:) + & 
               alpha_2(1,:,:) * c_s_inv(1,:,:) * &
               (L_s_p(1,:,:) + L_s_m(1,:,:)) - alpha_1(1,:,:) * f_inv(1,:,:) * &
                 (L_f_p(1,:,:) + L_f_m(1,:,:))) &
     + a_inv(1,:,:) * (az(1,:,:) * (L_a_p(1,:,:) + L_a_m(1,:,:)) + ay(1,:,:) * &
   (segno(1,:,:) * alpha_1(1,:,:) * (L_s_p(1,:,:) - L_s_m(1,:,:)) &
      + ax(1,:,:) * alpha_2(1,:,:) * f_inv(1,:,:) * (L_f_p(1,:,:) - L_f_m(1,:,:)))))  

!call derx_BC(By,grad)
!RHS_y(1,:,:) =  Bx(1,:,:) * grad(1,:,:)
do iz = 1, nz
 call dery_1(By(1,:,iz),zy)
 RHS_y(1,:,iz) = RHS_y(1,:,iz) - By(1,:,iz) * zy
enddo

 endif
 if ( mpime == nprow -1 ) then 

 RHS_y(nxl,:,:) = 0.5 * Den(nxl,:,:) * (Uy(nxl,:,:)* (2.* L_0(2,:,:) *c_s_inv(2,:,:) + & 
                alpha_2(2,:,:) * c_s_inv(2,:,:) * &
               (L_s_p(2,:,:) + L_s_m(2,:,:)) - alpha_1(2,:,:) * f_inv(2,:,:) * &
                 (L_f_p(2,:,:) + L_f_m(2,:,:))) &
     + a_inv(2,:,:) * (az(2,:,:) * (L_a_p(2,:,:) + L_a_m(2,:,:)) + ay(2,:,:) * &
   (segno(2,:,:) * alpha_1(2,:,:) * (L_s_p(2,:,:) - L_s_m(2,:,:)) &
      + ax(2,:,:) * alpha_2(2,:,:) * f_inv(2,:,:) * (L_f_p(2,:,:) - L_f_m(2,:,:)))))  

!call derx_BC(By,grad)
!RHS_y(nxl,:,:) =  Bx(nxl,:,:) * grad(nxl,:,:)
do iz = 1, nz
 call dery_1(By(nxl,:,iz),zy)
 RHS_y(nxl,:,iz) = RHS_y(nxl,:,iz) - By(nxl,:,iz) * zy
enddo

 endif
ENDIF


 !**************************************************************!
 ! Calcolo del secondo termine MOTO_x da derivare rispetto a y  !
 !                                                              !
 ! Qui riutilizzo A_xy = A_yx calcolato prima                   !
 !                                                              !
 ! Ricorda che Gi_yx = Gi_xy                                    !
 ! (vedi commenti in flr_i.f90 per maggiori dettagli)           !
 !**************************************************************!

do iz = 1, nz 
 do ix = 1, nxl
  call  dery_1(AA(ix,:,iz), zy)
  RHS_x(ix,:,iz) = RHS_x(ix,:,iz) + zy
 enddo
enddo

 !*************************************************************!
 ! Calcolo del secondo termine MOTO_y da derivare rispetto a y !
 !                                                             !
 ! Ricorda che Gi_yy = -Gi_xx                                  !
 ! (vedi commenti in flr_i.f90 per maggiori dettagli)          !
 !*************************************************************!

AA = By * By * ( 1.0d0 - B_inv * ( pe_para + pi_para - pe_perp - pi_perp ) ) &
        -  Ptot - Den * (Uiy * Uiy + de2 * Uey * Uey) + Gi_xx

do iz = 1, nz 
 do ix = 1, nxl
  call  dery_1(AA(ix,:,iz), zy)
  RHS_y(ix,:,iz) = RHS_y(ix,:,iz) + zy
 enddo
enddo

 !***********************************************************!
 ! Calcolo del terzo termine MOTO_x da derivare rispetto a z !
 !                                                           !
 ! Ricorda che Gi_zx = Gi_xz                                 !
 ! (vedi commenti in flr_i.f90 per maggiori dettagli)        !
 !***********************************************************!



RHS_x = RHS_x * de12inv
RHS_y = RHS_y * de12inv

! Iperviscosita' su n * Ux............eventualmente aggiornare per versione parallela 2D/3D


deallocate(zx)
deallocate(zy)
deallocate(zz)
deallocate(AA)
deallocate(Ptot)
deallocate(At)
deallocate(B_inv)

end subroutine
